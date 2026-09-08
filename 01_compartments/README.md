# A/B compartments

fig. S1A and table S1. Methods §"Compartment, TAD calling and analysis of TAD boundary".

Follows Lieberman-Aiden et al. (2009): O/E contact matrix → Pearson correlation matrix →
PCA. Because the sign of an eigenvector is arbitrary, A and B are assigned by H3K27ac
signal — the sign with the higher H3K27ac inner product is called A.

| Script | What it does |
|---|---|
| `00_process_bedGraph.py` | H3K27ac bedGraph → one length-normalised signal vector per chromosome at 500 kb. Writes `00_bg_vectors/<prefix>_<chr>.npy`. |
| `01_call_compartment.py` | Per chromosome: `juicer_tools dump oe NONE` at 500 kb → dense symmetric matrix → drop all-zero rows/columns → `np.corrcoef` → `PCA(n_components=2)` → re-insert the dropped bins as zeros → orient PC1/PC2 by the H3K27ac vector. Writes `01_AB_vec/<sample>_<chr>.npy` (PC1, PC2) and a per-chromosome summary CSV of the A/B H3K27ac totals. |

Run order: `00_process_bedGraph.py` then `01_call_compartment.py`.
Samples: HSPC (GSM2861708), OCI-AML3, OCI-AML3-HOXA9-KO, AML4943. Chromosomes 1–22, X, Y.
