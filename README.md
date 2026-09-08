# Code for *HOXA9 reinforces chromatin boundary insulation to maintain leukemia-associated genome architecture*

Wang, Yang, Feng et al. — analysis code for the **Hi-C / 3D genome** analyses in the
manuscript (Figs. 1–5, figs. S1–S4, tables S1, S2, S4, S8).

Everything here is short, plain Python / R / bash driving well-known published tools.
Nothing in this repository is a new algorithm; it is the exact parameterisation and glue
needed to reproduce the reported numbers from the deposited data.

---

## What is in scope

| Directory | Analysis | Manuscript |
|---|---|---|
| `00_hic_processing/` | Hi-C read alignment, valid-pair generation, library merging/downsampling, HiCRep similarity | Methods "Hi-C data processing"; table S8 |
| `01_compartments/` | A/B compartments (O/E → Pearson → PCA, sign fixed by H3K27ac) at 500 kb | fig. S1A; table S1 |
| `02_tads/` | TopDom TAD calling at 50 kb from `.hic` | Fig. 1A–C; figs. S1C, S1G, S4C; table S2 |
| `03_loops/` | FitHiChIP HOXA9-anchored loop calling; anchor annotation; random-loop background | Fig. 2A, C–E; Fig. 5C; fig. S3; table S4 |
| `04_apa/` | Aggregate peak analysis of the 5,824 HOXA9 loops across five contact maps | Fig. 2B |
| `05_insulation/` | Insulation scores; O/E contact pile-ups centred on loop anchors | Fig. 3B; Fig. 4B |

Each directory has its own README giving the exact parameters, the panels it produces, and
the numbers it should reproduce.

Note the split of labour on the TAD panels: this repository produces the TopDom TAD and
boundary calls; the downstream classification of each TAD as unchanged / fused / split /
shifted was done by a co-author (see below).

### Key reported numbers and where they come from

| Number | Source |
|---|---|
| 5,824 OCI-AML3 HOXA9 loops | `03_loops/configs/configfile_O3_WT_HOXA9_10000_0.05` → `FitHiChIP.interactions_FitHiC_Q0.05_MergeNearContacts.bed` |
| 13,041 AML4943 HOXA9 loops | `03_loops/configs/configfile_AML4943_HOXA9_10000_0.05` |
| APA centre enrichment 4.18 (OCI-AML3) / 2.22 (HSPC) | `04_apa/APA_summary.csv` |
| TopDom domain counts (OCI-AML3 3,180; HSPC 2,798; HOXA9-KO 3,064; dTAG-13 3,171) | `02_tads/` merged BEDs, rows where column 4 is `domain` |

## What is **not** in scope

These analyses were performed by other authors and are not in this repository. All of them
use standard published tools whose parameters are already stated in the Methods:

- CUT&RUN and ChIP-seq processing and peak calling (Trim Galore, Bowtie2, Picard, MACS2,
  HOMER, deepTools) — Figs. 2C, 3A, 5A, figs. S2, S4B
- RNA-seq processing and differential expression (Trim Galore, TopHat2, Cuffdiff) —
  Fig. 6A–B, Fig. 7F, tables S5, S6
- **Classification of TAD changes into fused / split / shifted** (BEDTools `windowBed` on the
  TopDom boundary calls produced by `02_tads/`) — Fig. 1B–G, fig. S1D, fig. S3C, fig. S4D
- E-P / E-E / E-S loop classification from histone-mark enrichment — Fig. 2C–D
- Motif enrichment and CTCF orientation calls — Fig. 2F–G, Fig. 5D, fig. S3B
- GO enrichment — Figs. 6C, 6E, 7F
- TCGA / BEAT AML expression correlation — Fig. 6F, fig. S6
- All wet-lab quantification — Fig. 7, fig. S5A, fig. S7

If Science Advances wants a single repository covering the whole paper, the co-authors
responsible for those sections need to contribute their scripts alongside these.

## Data

No data files are included here. Inputs come from:

- **GEO**: GSE301083, GSE301086, GSE306376, GSE302129, GSE302061, GSE302258
- **AML4943 high-resolution Hi-C**: https://doi.org/10.5281/zenodo.16930079
- **Previously published**: GSM5997594, GSM5997590, GSM3032904, GSM3032900, GSM6330185
  (OCI-AML3 histone marks / CTCF / ATAC) and GSM2861708 (CD34+ HSPC Hi-C)

All coordinates are **hg19 / GRCh37**.

## External software

| Tool | Version | Use |
|---|---|---|
| HiC-Pro | 3.1.0 | Hi-C alignment / valid pairs |
| bwa + pairtools + cooler + pairix | — | 4DN-style Hi-C alignment (see `NOTES`) |
| juicer_tools | 1.11.04 | `.hic` build and `dump` |
| hicstraw | — | `.hic` → locus-pair bed |
| TopDom | 0.0.2 | TAD calling (source inlined in `02_tads/02_run_TopDom.R`) |
| FitHiChIP | ay-lab/FitHiChIP | loop calling |
| hicrep (dejunlin) | — | contact-map similarity |
| scHiCTools | liu-bioinfo-lab/scHiCTools | 2D/3D embedding of similarity |
| BEDTools | 2.30.0 | interval overlaps |
| Python | 3 + numpy, pandas, scipy, scikit-learn, matplotlib | |
| R | 4.4.3 | TopDom driver |

## Paths

The scripts were run on a SLURM cluster and carry **absolute paths** to `.hic` files,
`juicer_tools`, chrom.sizes and conda environments. They are left as-run so that the exact
inputs of every reported number are unambiguous. To re-run elsewhere, edit the path
constants at the top of each script (they are all grouped there) and point them at the
GEO/Zenodo downloads.

`.hic` files are not redistributed here; `hg19.chrom.sizes` is available from UCSC.
