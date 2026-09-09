# HOXA9-anchored chromatin loops (FitHiChIP)

Fig. 2A/C–E, Fig. 5C, fig. S3, table S4. Methods §"Loop calling, aggregated peak analyses,
and protein enrichment at loop anchors".

FitHiChIP (ay-lab/FitHiChIP), unmodified upstream release. Every run uses the **same
parameter set**, which is exactly what the Methods describe:

| Parameter | Value | Meaning |
|---|---|---|
| `IntType` | 3 | peak-to-all |
| `BINSIZE` | 10000 | 10 kb resolution |
| `LowDistThr` / `UppDistThr` | 20000 / 3000000 | 20 kb – 3 Mb |
| `UseP2PBackgrnd` | 0 | peak-to-all background |
| `BiasType` | 1 | coverage bias correction |
| `MergeInt` | 1 | merge nearby interactions |
| `QVALUE` | 0.05 | FDR threshold |

Loops are read from `<OutDir>/.../Merge_Nearby_Interactions/FitHiChIP.interactions_FitHiC_Q0.05_MergeNearContacts.bed`
(tab-separated, one header line).

## Runs

| Config | Contact map | Peaks | Reported |
|---|---|---|---|
| `configs/configfile_O3_WT_HOXA9_10000_0.05` | OCI-AML3 WT (locus-pair bed dumped from the WT map) | 2,331 HOXA9 CUT&RUN peaks (HOXA9 ∩ HA, IgG-corrected) | **5,824 HOXA9 loops** — Fig. 2, table S4 |
| `configs/configfile_AML4943_HOXA9_10000_0.05` | AML4943 | AML4943 HOXA9 peaks (IgG-corrected) | **13,041 loops** — Fig. 5C |
| `configs/configfile_HSPC_HOXA9_10000_0.05` | CD34+ HSPC (GSM2861708) | same OCI-AML3 HOXA9 peaks | HSPC comparison |
| `configs/configfile_O3_HOXA9KO_HOXA9_10000_0.05` | HOXA9-KO, merged 4DN map | same | Fig. 2A |
| `configs/configfile_O3_24h_dTAG13_HOXA9_10000_0.05` | 24 h dTAG-13, merged 4DN map | same | Fig. 2A |
| `configs/configfile_O3_24h_dTAG13_dsWT_HOXA9_10000_0.05` | 24 h dTAG-13 downsampled to WT depth | same | depth control for Fig. 2A |
| `configs/configfile_O3_HOXA9KO_old_HOXA9_10000_0.05` | HOXA9-KO, earlier single-library map | same | superseded by the merged map above; kept for provenance |

## Scripts

| Script | What it does |
|---|---|
| `run_fithichip.sh` | `sbatch run_fithichip.sh configs/<configfile>`. Set `FITHICHIP_HOME` to your FitHiChIP checkout. |
| `dump_hic_to_bed.py` | `.hic` → FitHiChIP 7-column locus-pair bed (`chr s1 e1 chr s2 e2 count`), cis only, raw observed counts, via `hicstraw`. Needed when feeding a `.hic` through FitHiChIP's `Bed=` input. |
| `annotate_anchors.py` | Annotates each 10 kb anchor with four Yes/No features plus gene names: overlaps a HOXA9 peak; overlaps a TopDom `boundary` region of that condition's own TAD calls; overlaps a promoter (TSS ± 2 kb); overlaps a gene body. This is what defines the **enhancer–boundary (E-B)** class in Fig. 2E and the loop-to-gene assignments in Fig. 2D. Output: `<run>_merged.annotated.csv`. |
| `random_loop_background.py` | Enrichment background: builds 2 matched random loop sets per real set — (1) fully random anchors with the observed span distribution, (2) the real HOXA9 anchor kept, the partner anchor placed at a shuffled span. Fixed seed (0). Used for the overlap-enrichment significance described in the Methods. |

**Note on `Bed=` vs `HIC=`:** FitHiChIP overwrites the shared `ChrSizeFile` when the `HIC=`
input is used, which can clobber a chrom.sizes file shared between runs. Where that mattered,
the `.hic` was dumped to a 7-column bed with `dump_hic_to_bed.py` and passed via `Bed=`
instead. Both forms appear in the configs.

**Not included here:** E-P / E-E / E-S classification from H3K4me3 / H3K4me1 / H3K27ac /
H3K27me3 enrichment, motif enrichment, and CTCF orientation calls.
