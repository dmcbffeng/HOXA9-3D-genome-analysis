# Chromatin insulation

Fig. 3B and Fig. 4B. Methods §"Calculating chromatin insulation scores".

For a locus *k* on a 10 kb contact matrix *M* with window *w*:

    S1 = sum over {k-w <= i,j <= k+w, i != j} M_ij      all contacts in the window
    S2 = 2 * sum over {k-w <= i < k < j <= k+w} M_ij    contacts crossing locus k
    score_k = log2( S2 / (S1 - S2) )

Valleys mark candidate boundaries. Loci where `S2 == 0` or `S1 - S2 == 0` are undefined and
reported as NA (a `1e-8` term guards the 0/0 case).

| Script | What it does |
|---|---|
| `insulation.py` | `python insulation.py <sample> <observed\|oe>`. Dumps `dump <observed\|oe> KR` at 10 kb per chromosome and writes one CSV per window (w = 10 → 100 kb, w = 20 → 200 kb) with columns `chr,start,end,score,score_z`. |
| `run_insulation.sh` | SLURM wrapper: `sbatch run_insulation.sh <sample> <observed\|oe>`. |
| `pileup_at_anchors.py` | The pile-up heatmaps under Fig. 3B / 4B: `dump oe KR` at 10 kb, ±20 bins (41 × 41) around the centre of each site in an anchor list, averaged. Takes one peak/anchor BED per class and writes `APA_<class>.npy`. The classes in Fig. 3B are HOXA9 E-P, E-E, E-B anchors plus RAD21 and CTCF sites; Fig. 4B swaps in HOXA9 vs HOXA9/MEIS1 co-bound E-B anchors and MEIS1/PBX2/PBX3-only controls. |

## Read this before re-running

The matrix choice matters. On the **observed** matrix the "within" pairs counted in `S1 - S2` are on average closer in
genomic distance than the boundary-crossing pairs in `S2`. Contact frequency decays steeply
with distance, so `S2 << S1 - S2` at essentially every locus whether or not a boundary is
present, and >99.9% of raw scores come out negative. Z-scoring recentres them and makes the
signs look balanced, but it is an affine rescale — it hides the bias rather than removing it.
Dividing by the distance-expected value removes it properly, which is why the Methods
specify **O/E**; on O/E about 43% of bins are positive with no rescaling.

Also: each sample is z-scored against its own genome-wide mean and SD, so `score_z` is not
comparable across samples. Use the raw `score` column for anything differential.

Finally, a few hundred bins per sample sit in pericentromeric / segmental-duplication regions
where KR leaves the matrix near-degenerate, `S1 - S2` collapses and the `1e-8` floor dominates.
These are numerical artifacts, not boundaries (mostly chr9, chr1, chr6, chr10). Drop
`|score_z| > 5` — about 0.3% of valid bins — before plotting or boundary calling.
