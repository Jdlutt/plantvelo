# PlantVelo

**PlantVelo** is a plant-specific RNA velocity analysis toolkit built on top of [velocyto.py](https://github.com/velocyto-team/velocyto.py). It extends the standard two-state splicing model (unspliced → spliced) with an explicit **intron retention (IR)** state, enabling a three-state kinetic framework tailored to the biology of plant cells.

---

## Background

### Why a plant-specific tool?

Standard RNA velocity methods model RNA dynamics as a two-state process:

```
Unspliced (U)  →  Spliced (S)
```

This model was designed primarily for animal cells. Plant cells, however, exhibit substantially higher rates of **intron retention (IR)** — estimated at 60% of transcripts under normal conditions. In plants, IR is not simply a splicing failure; it is an active regulatory mechanism that:

- Produces stable, cytoplasmic IR transcripts
- Is reversible: retained introns can be spliced post-transcriptionally

Lumping IR reads together with nascent pre-mRNA reads into a single "unspliced" category obscures this biology and leads to systematic misestimation of RNA velocity in plant data.

### The three-state model

PlantVelo separates the intron-only read population into a dedicated layer:

```
Unspliced (U)  →  Intron-Retained (IR)  →  Spliced (S)
```

| State | Biological meaning | Read signature |
|---|---|---|
| **Unspliced (U)** | Nascent pre-mRNA, actively being transcribed | Reads spanning exon–intron boundary |
| **Intron-Retained (IR)** | Stable IR isoform, intron fully retained | Reads lying entirely within intron, ≥ `ir_flanking` bp from each edge |
| **Spliced (S)** | Mature mRNA, intron removed | Reads mapping only to exon features |

This distinction allows downstream kinetic modelling to estimate separate rate constants for the U→IR and IR→S transitions, providing a more accurate description of gene expression dynamics in plant single-cell data.

---

## Features

- **Drop-in replacement for `velocyto run`** — same BAM/GTF inputs, same loom output format, fully backward-compatible
- **New `intron_retained` loom layer** — output loom files contain four layers (`spliced`, `unspliced`, `ambiguous`, `intron_retained`)
- **`--ir-flanking` parameter** — tunable boundary buffer (bp) to guard against mis-mapping at intron edges
- **Two plant-aware logic classes**:
  - `PlantPermissive10X` — counts both validated and non-validated intron-only reads as IR
  - `PlantValidated10X` — counts only validated intron-only reads as IR (more conservative)
- **Full compatibility with standard velocyto logics** — `Permissive10X`, `ValidatedIntrons10X`, available via `--logic`

---

## Installation

```bash
conda create -n plantvelo
conda activate plantvelo
conda install bioconda::velocyto.py
pip install git+https://github.com/Jdlutt/plantvelo.git
```

### Verify installation

```bash
plantvelo --version
plantvelo run --help
```

---

## Usage

### Basic usage (10x Genomics data)

```bash
plantvelo run \
    -b filtered_feature_bc_matrix/barcodes.tsv \
    -o plantvelo_output \
    -m /path/to/repeat_mask.gtf \
    --logic PlantPermissive10X \
    --ir-flanking 5 \
    --sample-name samples \
    /path/to/possorted_genome_bam.bam \
    /path/to/genome_annotation.gtf
```

This produces `plantvelo_output/<sample>.loom` with four layers:

| Layer | Content |
|-------|---------|
| `spliced` | Reads mapping only to exon features |
| `unspliced` | Reads spanning an exon–intron boundary (nascent pre-mRNA) |
| `ambiguous` | Reads compatible with both spliced and unspliced models |
| `intron_retained` | Reads lying fully within an intron (IR isoform) |

---

## Parameters

### `plantvelo run`

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `BAMFILE` | — | required | Position-sorted BAM file(s) |
| `GTFFILE` | — | required | Genome annotation GTF |
| `--bcfile` | `-b` | None | Valid cell barcodes file (e.g. `barcodes.tsv` from Cell Ranger). If omitted, all barcodes in the BAM are used. |
| `--outputfolder` | `-o` | `./plantvelo/` | Output directory (created if absent) |
| `--sample-name` | — | None | Sample name used for output filename and cell barcode |
| `--mask` | `-m` | None | Repeat masking GTF (strongly recommended) |
| `--logic` | `-l` | `PlantPermissive10X` | Read classification logic class |
| `--ir-flanking` | — | `5` | Minimum bp from intron edge required to classify a read as IR |
| `--dtype` | `-t` | `uint32` | Numeric dtype for loom layers. Use `uint32` if >6000 UMIs/gene/cell are expected. |
| `--samtools-threads` | `-@` | `16` | Threads for samtools sort |
| `--samtools-memory` | — | `2048` | MB per thread for samtools sort |
| `--without-umi` | `-U` | off | Count reads instead of UMI-collapsed molecules |
| `--multimap` | `-M` | off | Include non-uniquely mapped reads (not recommended) |
| `--verbose` | `-v` | 1 | Verbosity: `-v` warnings, `-vv` info, `-vvv` debug |

### Available logic classes

| Logic | IR layer | Non-validated introns |
|-------|----------|-----------------------|
| `PlantPermissive10X` | ✅ | Counted as IR |
| `PlantValidated10X` | ✅ strict | Discarded |

### `--ir-flanking` guidance

| Value | Recommended when |
|-------|-----------------|
| `5` | Default; standard 10x short reads (≥75 bp) |
| `10` | Read length < 75 bp, or alignment quality is uncertain |
| `20` | Long-read data |
| `0` | Maximum permissiveness; expect higher noise |

---

## Output format

The output is a standard [loom file](http://loompy.org/) readable by velocyto, scVelo, and other RNA velocity tools.
---

## How intron retention is detected

```
BAM reads
    │
    ▼
CIGAR parsing                    counter.py (velocyto)
    │  segments + ref_skipped flag
    ▼
Feature matching                 indexes.py (velocyto)
    │  each segment → exon / intron features
    ▼
PlantPermissive10X.count()       logic.py (plantvelo)
    │
    ├─ Only exon-compatible models
    │       └─► spliced
    │
    ├─ All models span exon–intron boundary
    │       └─► unspliced  (nascent pre-mRNA)
    │
    ├─ Only intron models
    │   ├─ all segments ≥ ir_flanking from intron edges
    │   │       └─► intron_retained  ◄── NEW
    │   └─ any segment close to edge (possible mis-map)
    │           └─► unspliced
    │
    └─ Ambiguous / multi-gene
            └─► ambiguous / discarded
```

---

## Project structure

```
plantvelo/
├── setup.py
├── README.md
└── plantvelo/
    ├── __init__.py
    ├── _version.py
    ├── logic.py              # PlantPermissive10X, PlantValidated10X
    └── commands/
        ├── __init__.py
        ├── plantvelo.py      # CLI entry point
        ├── run.py            # plantvelo run command
        └── _run.py           # pipeline core function
```

---

## Citation
PlantVelo builds on:

> La Manno G, Soldatov R, Zeisel A, et al. RNA velocity of single cells[J]. Nature, 2018, 560(7719): 494-498.

---

