"""plantvelo.commands.run
~~~~~~~~~~~~~~~~~~~~~~~~
``plantvelo run`` CLI command – a drop-in replacement for ``velocyto run``
with an additional ``--ir-flanking`` option for intron-retention detection.

Usage example:

    plantvelo run \\
        -b filtered_feature_bc_matrix/barcodes.tsv \\
        -o plantvelo_output \\
        -m /path/to/repeat_mask.gtf \\
        --logic PlantPermissive10X \\
        --ir-flanking 5 \\
        --sample-name samples \\
        /path/to/possorted_genome_bam.bam \\
        /path/to/genome_annotation.gtf
"""
import sys
import os
import logging
from typing import *

import click

from plantvelo.commands._run import _run


@click.command(short_help="Run PlantVelocity analysis outputting a loom file with intron_retained layer")
@click.argument(
    "bamfile",
    nargs=-1,
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
)
@click.argument(
    "gtffile",
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
)
@click.option(
    "--bcfile", "-b",
    help=(
        "Valid barcodes file, to filter the bam. "
        "If --bcfile is not specified all the cell barcodes will be included. "
        "Cell barcodes should be specified in the bcfile as the CB tag for each read."
    ),
    default=None,
    type=click.Path(resolve_path=True, file_okay=True, dir_okay=False, readable=True),
)
@click.option(
    "--outputfolder", "-o",
    help="Output folder, created if it does not exist.",
    default=None,
    type=click.Path(exists=False),
)
@click.option(
    "--sample-name",
    help=(
        "Used as the output loom filename (<sample_name>.loom) and to reformat "
        "cell barcodes for Seurat compatibility: <sample_name>_<barcode>-1. "
    ),
    default=None,
    type=str,
)
@click.option(
    "--metadatatable", "-s",
    help="CSV table containing metadata of the various samples (rows=samples, cols=entries).",
    default=None,
    type=click.Path(resolve_path=True, file_okay=True, dir_okay=False, readable=True),
)
@click.option(
    "--mask", "-m",
    help=".gtf file containing intervals to mask (repeat regions).",
    default=None,
    type=click.Path(resolve_path=True, file_okay=True, dir_okay=False, readable=True),
)
@click.option(
    "--onefilepercell", "-c",
    help=(
        "If set, every bamfile passed is interpreted as an independent cell; "
        "otherwise multiple files are a batch of cells (default: off)."
    ),
    default=False,
    is_flag=True,
)
@click.option(
    "--logic", "-l",
    help=(
        "The logic to use for read classification. "
        "Plant-aware logics: PlantPermissive10X (default), PlantValidated10X. "
    ),
    default="PlantPermissive10X",
)
@click.option(
    "--ir-flanking",
    help=(
        "Minimum distance (bp) a read segment must lie from each intron edge "
        "to be counted as intron_retained rather than unspliced. "
        "Only used by plant-aware logic classes. (default: 5)"
    ),
    default=5,
    type=int,
)
@click.option(
    "--without-umi", "-U",
    help="If set, data is assumed UMI-less and reads are counted instead of molecules (default: off).",
    default=False,
    is_flag=True,
)
@click.option(
    "--umi-extension", "-u",
    help=(
        "Extend UMI to guarantee uniqueness: `chr`, `Gene`, or `[N]bp`. "
        "See velocyto documentation for details. (default: no)"
    ),
    default="no",
)
@click.option(
    "--multimap", "-M",
    help="Consider non-unique mappings (not recommended).",
    default=False,
    is_flag=True,
)
@click.option(
    "--samtools-threads", "-@",
    help="Number of threads for samtools sort. (default: 16)",
    default=16,
)
@click.option(
    "--samtools-memory",
    help="MB of memory per thread for samtools sort. (default: 2048)",
    default=2048,
)
@click.option(
    "--dtype", "-t",
    help="dtype for loom layer matrices. Use uint32 if >6000 molecules/gene/cell are expected. (default: uint32)",
    default="uint32",
)
@click.option(
    "--dump", "-d",
    help="Debug: dump molecular mapping report to hdf5. --dump N saves every N cells. (default: 0)",
    default="0",
)
@click.option(
    "--verbose", "-v",
    help="Verbosity level: -v warnings, -vv info, -vvv debug.",
    count=True,
    default=1,
)
def run(
    bamfile: Tuple[str, ...],
    gtffile: str,
    bcfile: str,
    outputfolder: str,
    sample_name: str,
    metadatatable: str,
    mask: str,
    onefilepercell: bool,
    logic: str,
    ir_flanking: int,
    without_umi: str,
    umi_extension: str,
    multimap: bool,
    samtools_threads: int,
    samtools_memory: int,
    dtype: str,
    dump: str,
    verbose: int,
    additional_ca: dict = {},
) -> None:
    """Run PlantVelocity analysis outputting a loom file.

    BAMFILE  one or several position-sorted BAM files

    GTFFILE  genome annotation GTF file
    """
    return _run(
        bamfile=bamfile,
        gtffile=gtffile,
        bcfile=bcfile,
        outputfolder=outputfolder,
        sampleid=sample_name,
        sample_name=sample_name,
        metadatatable=metadatatable,
        repmask=mask,
        onefilepercell=onefilepercell,
        logic=logic,
        ir_flanking=ir_flanking,
        without_umi=without_umi,
        umi_extension=umi_extension,
        multimap=multimap,
        test=False,
        samtools_threads=samtools_threads,
        samtools_memory=samtools_memory,
        dump=dump,
        loom_numeric_dtype=dtype,
        verbose=verbose,
        additional_ca=additional_ca,
    )
