"""plantvelo.commands._run
~~~~~~~~~~~~~~~~~~~~~~~~~
Core pipeline function for plantvelo run.

Forked from velocyto/commands/_run.py and extended with:
  - ``ir_flanking`` parameter passed into the Logic constructor
  - Logic resolution also searches plantvelo.logic (PLANT_LOGICS registry)
  - loom file_attrs record plantvelo version and logic name
  - ``sample_name`` parameter: when provided, output filename becomes
    <sample_name>.loom and CellIDs are formatted as <sample_name>_<barcode>-1
    for Seurat compatibility
"""
import sys
import os
import glob
import re
import gzip
import array
import loompy
import numpy as np
import subprocess
import multiprocessing
import csv
import itertools
from collections import defaultdict
import logging
import h5py
from typing import *

import velocyto as vcy

from plantvelo.logic import PLANT_LOGICS


def _run(
    *,
    bamfile: Tuple[str, ...],
    gtffile: str,
    bcfile: str,
    outputfolder: str,
    sampleid: str,
    sample_name: Optional[str] = None,
    metadatatable: str,
    repmask: str,
    onefilepercell: bool,
    logic: str,
    ir_flanking: int,
    without_umi: str,
    umi_extension: str,
    multimap: bool,
    test: bool,
    samtools_threads: int,
    samtools_memory: int,
    loom_numeric_dtype: str,
    dump: bool,
    verbose: int,
    additional_ca: dict = {},
) -> None:
    """Run PlantVelocity analysis and output a loom file.

    Parameters
    ----------
    sample_name : str, optional
        When provided, used as the output loom filename and as prefix for
        CellIDs in Seurat-compatible format (<sample_name>_<barcode>-1).
        When None, CellIDs use the legacy format (<sampleid>:<barcode><gem_grp>)
        and the filename is derived from the BAM file.
    ir_flanking : int
        Minimum distance (bp) a read segment must lie from each intron edge
        to be classified as intron_retained rather than unspliced.
        Only used by plant-aware logic classes (PlantPermissive10X,
        PlantValidated10X).  Standard velocyto logics ignore this value.
    """

    ########################
    #    Resolve Inputs    #
    ########################

    logging.basicConfig(
        stream=sys.stdout,
        format='%(asctime)s - %(levelname)s - %(message)s',
        level=[logging.ERROR, logging.WARNING, logging.INFO, logging.DEBUG][verbose],
    )

    if isinstance(bamfile, tuple) and len(bamfile) > 1 and bamfile[-1][-4:] in [".bam", ".sam"]:
        multi = True
    elif isinstance(bamfile, tuple) and len(bamfile) == 1:
        multi = False
    else:
        raise IOError(f"Something went wrong in the argument parsing. You passed as bamfile: {bamfile}")

    if onefilepercell and multi:
        if bcfile is not None:
            raise ValueError("Inputs incompatibility. --bcfile/-b option was used together with --onefilepercell/-c option.")
        logging.warning("Each bam file will be interpreted as a DIFFERENT cell")
    elif not onefilepercell and multi:
        logging.warning("Several input files but --onefilepercell is False. Each bam file will be interpreted as containing a SET of cells!!!")

    if sampleid is None:
        assert metadatatable is None, "--metadatatable was specified but cannot fetch sample metadata without valid sampleid"
        if multi:
            logging.warning("When using multiple files you may want to use --sample-name option to specify the name of the output file")
        if multi and not onefilepercell:
            full_name = "_".join([os.path.basename(bamfile[i]).split(".")[0] for i in range(len(bamfile))])
            if len(full_name) > 50:
                sampleid = f'multi_input_{os.path.basename(bamfile[0]).split(".")[0]}'
            else:
                sampleid = f'multi_input_{full_name}'
        elif multi and onefilepercell:
            sampleid = f'onefilepercell_{os.path.basename(bamfile[0]).split(".")[0]}'
        else:
            sampleid = os.path.basename(bamfile[0]).split(".")[0]
        logging.info(f"No SAMPLEID specified, the sample will be called {sampleid}")

    if outputfolder is None:
        outputfolder = os.path.join(os.path.split(bamfile[0])[0], "plantvelo")
        logging.info(f"No OUTPUTFOLDER specified, find output files inside {outputfolder}")
    if not os.path.exists(outputfolder):
        os.mkdir(outputfolder)

    # ----------------------------------------------------------------
    # Resolve logic class
    # Search order: plantvelo.logic (PLANT_LOGICS) → velocyto namespace
    # ----------------------------------------------------------------
    if logic in PLANT_LOGICS:
        logic_class = PLANT_LOGICS[logic]
        logging.info(f"Using plantvelo logic: {logic}")
    else:
        logic_class = getattr(vcy, logic, None)
        if logic_class is None or not issubclass(logic_class, vcy.Logic):
            available = list(PLANT_LOGICS.keys()) + [
                k for k, v in vcy.logic.__dict__.items()
                if isinstance(v, type) and issubclass(v, vcy.Logic)
            ]
            raise ValueError(
                f"{logic} is not a valid logic. "
                f"Choose one among: {', '.join(available)}"
            )
        logging.info(f"Using velocyto logic: {logic}")

    # Instantiate logic: plant logics accept ir_flanking; velocyto logics do not
    try:
        logic_obj = logic_class(ir_flanking=ir_flanking)
        logging.debug(f"Logic {logic} instantiated with ir_flanking={ir_flanking}")
    except TypeError:
        logic_obj = logic_class()
        logging.debug(f"Logic {logic} instantiated (ir_flanking not supported, ignored)")

    ########################
    #     Barcodes         #
    ########################

    if bcfile is None:
        logging.debug("Cell barcodes will be determined while reading the .bam file")
        valid_bcset = None
    else:
        valid_bcs_list = (gzip.open(bcfile).read().decode() if bcfile.endswith(".gz") else open(bcfile).read()).rstrip().split()
        valid_cellid_list = np.array([f"{sampleid}:{v_bc}" for v_bc in valid_bcs_list])
        if len(set(bc.split('-')[0] for bc in valid_bcs_list)) == 1:
            gem_grp = f"-{valid_bcs_list[0].split('-')[-1]}"
        else:
            gem_grp = "x"
        valid_bcset = set(bc.split('-')[0] for bc in valid_bcs_list)
        logging.info(f"Read {len(valid_bcs_list)} cell barcodes from {bcfile}")
        logging.debug(f"Example barcode: {valid_bcs_list[0].split('-')[0]} cell_id: {valid_cellid_list[0]}")

    if metadatatable:
        try:
            sample_metadata = vcy.MetadataCollection(metadatatable)
            sample = sample_metadata.where("SampleID", sampleid)
            if len(sample) == 0:
                logging.error(f"Sample ID {sampleid} not found in sample sheet")
                sample = {}
            elif len(sample) > 1:
                logging.error(f"Sample ID {sampleid} has multiple lines in sample sheet")
                sys.exit(1)
            else:
                sample = sample[0].dict
        except (NameError, TypeError):
            logging.warning("SAMPLEFILE was not specified. add -s SAMPLEFILE to add metadata.")
            sample = {}
    else:
        sample = {}

    ########################
    #     Start Analysis   #
    ########################

    if without_umi:
        if umi_extension != "no":
            logging.warning("--umi-extension was specified but incompatible with --without-umi, it will be ignored!")
        umi_extension = "without_umi"

    exincounter = vcy.ExInCounter(
        sampleid=sampleid,
        logic=logic_class,
        valid_bcset=valid_bcset,
        umi_extension=umi_extension,
        onefilepercell=onefilepercell,
        dump_option=dump,
        outputfolder=outputfolder,
    )

    try:
        mb_available = int(subprocess.check_output('grep MemAvailable /proc/meminfo'.split()).split()[1]) / 1000
    except (subprocess.CalledProcessError, FileNotFoundError):
        logging.warning("Cannot determine available memory; assuming 32 GB")
        mb_available = 32000

    threads_to_use = min(samtools_threads, multiprocessing.cpu_count())
    mb_to_use = int(min(samtools_memory, mb_available / (len(bamfile) * threads_to_use)))
    compression = vcy.BAM_COMPRESSION

    if onefilepercell and without_umi:
        tagname = "NOTAG"
    elif onefilepercell:
        tagname = "NOTAG"
        exincounter.peek_umi_only(bamfile[0])
    else:
        exincounter.peek(bamfile[0])
        tagname = exincounter.cellbarcode_str

    if multi and onefilepercell:
        bamfile_cellsorted = list(bamfile)
    elif onefilepercell:
        bamfile_cellsorted = [bamfile[0]]
    else:
        bamfile_cellsorted = [
            f"{os.path.join(os.path.dirname(bmf), 'cellsorted_' + os.path.basename(bmf))}"
            for bmf in bamfile
        ]

    sorting_process: Dict[int, Any] = {}
    check_end_process = False
    for ni, bmf_cellsorted in enumerate(bamfile_cellsorted):
        command = f"samtools sort -l {compression} -m {mb_to_use}M -t {tagname} -O BAM -@ {threads_to_use} -o {bmf_cellsorted} {bamfile[ni]}"
        if os.path.exists(bmf_cellsorted):
            logging.warning(f"The file {bmf_cellsorted} already exists. Sorting step will be skipped.")
        else:
            sorting_process[ni] = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
            logging.info(f"Sorting {bamfile[ni]} → {bmf_cellsorted}")
            check_end_process = True

    logging.info(f"Load the annotation from {gtffile}")
    annotations_by_chrm_strand = exincounter.read_transcriptmodels(gtffile)
    chrs = list(v for k, v in annotations_by_chrm_strand.items())
    tms = list(itertools.chain.from_iterable((v.values() for v in chrs)))
    ivls = list(itertools.chain.from_iterable(tms))
    logging.debug(f"Generated {len(ivls)} features from {len(tms)} transcript models")
    del chrs, tms, ivls

    if repmask is not None:
        logging.info(f"Load repeat mask from {repmask}")
        exincounter.read_repeats(repmask)

    logging.info(f"Scan {' '.join(bamfile)} to validate intron intervals")
    if test:
        import pickle
        if os.path.exists("exincounter_dump.pickle"):
            exincounter = pickle.load(open("exincounter_dump.pickle", "rb"))
        else:
            pickle.dump(exincounter, open("exincounter_dump.pickle", "wb"))
            exincounter.mark_up_introns(bamfile=bamfile, multimap=multimap)
    else:
        exincounter.mark_up_introns(bamfile=bamfile, multimap=multimap)

    if check_end_process:
        logging.info("Waiting for bam sorting to finish…")
        for k in sorting_process.keys():
            returncode = sorting_process[k].wait()
            if returncode == 0:
                logging.info(f"bam file #{k} sorted")
            else:
                raise MemoryError(
                    f"bam file #{k} could not be sorted by cells. "
                    "Try upgrading samtools (>= 1.6) or increase --samtools-memory."
                )

    logging.debug("Start molecule counting!")
    results = exincounter.count(bamfile_cellsorted, multimap=multimap)
    dict_list_arrays, cell_bcs_order = results

    ########################
    #         Output       #
    ########################

    if not exincounter.filter_mode:
        valid_bcset = exincounter.valid_bcset
        valid_bcs_list = list(valid_bcset)
        gem_grp = ""
        valid_cellid_list = np.array([f"{sampleid}:{v_bc}" for v_bc in valid_bcs_list])

    if sample_name is not None:
        # Reformat barcodes for Seurat compatibility:
        # {sampleid}:{barcode}{gem_grp}  →  {sample_name}_{barcode}-1
        # gem_grp "x" (mixed gem groups) is normalised to "-1"; a numeric
        # suffix such as "-1" is kept as-is; empty suffix gets "-1".
        _gem_suffix = gem_grp if (gem_grp and gem_grp != "x") else "-1"
        cell_ids = np.array([f"{sample_name}_{v_bc}{_gem_suffix}" for v_bc in cell_bcs_order])
    else:
        cell_ids = np.array([f"{sampleid}:{v_bc}{gem_grp}" for v_bc in cell_bcs_order])
    ca = {"CellID": cell_ids}
    ca.update(additional_ca)
    for key, value in sample.items():
        ca[key] = np.full(len(cell_bcs_order), value)

    loom_name = sample_name if sample_name is not None else sampleid
    outfile = os.path.join(outputfolder, f"{loom_name}.loom")
    logging.debug(f"Generating output file {outfile}")

    atr_table = (
        ("Gene",       "genename", str),
        ("Accession",  "geneid",   str),
        ("Chromosome", "chrom",    str),
        ("Strand",     "strand",   str),
        ("Start",      "start",    int),
        ("End",        "end",      int),
    )

    ra: Dict[str, np.ndarray] = {}
    for name_col_attr, name_obj_attr, dtyp in atr_table:
        tmp_array = np.zeros((len(exincounter.genes),), dtype=object)
        for gene_id, gene_info in exincounter.genes.items():
            tmp_array[exincounter.geneid2ix[gene_id]] = getattr(gene_info, name_obj_attr)
        ra[name_col_attr] = tmp_array.astype(dtyp)

    layers: Dict[str, np.ndarray] = {}
    for layer_name in logic_obj.layers:
        layers[layer_name] = np.concatenate(dict_list_arrays[layer_name], axis=1)
        del dict_list_arrays[layer_name]

    total: np.ndarray
    for layer_name in logic_obj.layers:
        try:
            total += layers[layer_name]
        except NameError:
            total = np.array(layers[layer_name])

    tmp_layers = {"": total.astype("float32", order="C", copy=False)}
    tmp_layers.update({
        layer_name: layers[layer_name].astype(loom_numeric_dtype, order="C", copy=False)
        for layer_name in logic_obj.layers
    })

    import plantvelo
    loompy.create(
        outfile,
        tmp_layers,
        ra,
        ca,
        file_attrs={
            "plantvelo.__version__": plantvelo.__version__,
            "plantvelo.logic": logic,
            "plantvelo.ir_flanking": str(ir_flanking),
        },
    )
    logging.info(f"Loom file written: {outfile}")
    logging.debug("Terminated successfully!")
