# Aggregate peak analysis (APA)

Fig. 2B.

The **same** 5,824-loop OCI-AML3 WT HOXA9 loop set (from `../03_loops/`) is piled up in
every contact map, so the panels are directly comparable: they show how strongly each map
supports the WT loop calls, not how many loops each map calls on its own.

Method: `juicer_tools dump oe KR` at 10 kb → for each loop, the 11 × 11 window (±5 bins =
±50 kb) of O/E values centred on the loop apex → average over loops. The periphery sits at
~1 (expected background); the centre pixel is the fold enrichment at the apex.
P2LL is the centre against the mean of the lower-left 3 × 3 corner (Rao et al. 2014).

One of the 5,824 loops falls too close to a chromosome end for a full window, so 5,823 are
piled.

| Script | What it does |
|---|---|
| `01_compute_APA.py` | Steps above; writes `<map>_APA_matrix.{npy,csv}` and `APA_summary.csv`. Streams the juicer dump and keeps only the pixels inside a loop window — a dense chromosome matrix at 10 kb is ~7 GB and will OOM. |
| `02_plot_APA.py` | Renders the matrices as Reds heatmaps on a shared 0–5 colour scale, centre enrichment annotated. |
| `APA_summary.csv` | The published values, included as a checksum for a re-run. |

Results:

| Map | Centre O/E | P2LL |
|---|---|---|
| OCI-AML3 WT | 4.18 | 2.36 |
| OCI-AML3 HOXA9-KO | 3.19 | 1.72 |
| OCI-AML3 24 h dTAG-13 | 3.23 | 1.77 |
| 24 h dTAG-13, downsampled to WT depth | 3.22 | 1.76 |
| CD34+ HSPC | 2.22 | 1.26 |

Downsampling barely changes the dTAG value (3.23 → 3.22), so the reduction in APA signal
after HOXA9 removal is not a sequencing-depth artifact.
