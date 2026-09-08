"""Insulation scores (localNorm) on KR-normalized contact maps.

Method (matches ins_scores/Image20260803201952.png and the original
FitHiChIP/2025Nov_MEIS1/insulation/O3_ins2.py):

    S1 = sum_{k-w<=i,j<=k+w, i!=j} M_ij          (all contacts in the window)
    S2 = 2 * sum_{k-w<=i<k, k<j<=k+w} M_ij       (contacts crossing locus k)
    Score_k = log2( S2 / (S1 - S2 + 1e-8) )

MATRIX CHOICE -- this matters a lot:

  M must be the **O/E** (observed/expected) KR matrix, as the written method says.

  On the raw *observed* KR matrix the score is dominated by distance decay rather
  than by insulation: the "within" pairs counted in S1-S2 are on average much closer
  in genomic distance than the boundary-crossing pairs in S2, so S2 << S1-S2 at every
  locus and the log2 ratio is negative genome-wide (>99.9% of bins). Dividing by the
  distance-expected value removes that bias, and the score then varies around 0 with
  both signs, which is what an insulation score is supposed to do.

  The archived O3_ins2.py had a bug here: line 21 builds an `oe` dump command as a
  string, then lines 23-27 overwrite it with an `observed` argument list. So every
  archived CSV is observed-based. Z-scoring those (as was done for the file sent out
  as "OCI-AML3 insulation score.csv") recentres them to mean 0 / std 1 and makes the
  signs look balanced, but it is only an affine rescale -- it does not remove the
  distance-decay bias, it just hides it.

The 1e-8 in the denominator guards the 0/0 case; loci with S2 == 0 or S1 - S2 == 0
are undefined and reported as NA.

Usage:  python insulation.py <tag> <observed|oe>
Writes one CSV per window (w=10 -> 100Kb, w=20 -> 200Kb) with columns
chr,start,end,score,score_z  (score_z = genome-wide z-score of score).
"""

import os
import subprocess
import sys

import numpy as np
import pandas as pd

JUICER = '/nfs/turbo/umms-drjieliu/juicer_tools_1.11.04_jcuda.0.8.jar'
OUT_DIR = '/nfs/turbo/umms-drjieliu/usr/temp_Fan/AML_Analysis/ins_scores'
CHROM_SIZES = os.path.join(OUT_DIR, 'hg19.chrom.sizes')

WINDOWS = [10, 20]          # bins; at 10 kb resolution -> 100 Kb and 200 Kb
RESOLUTION = 10000

chroms = [f'chr{i}' for i in list(range(1, 23)) + ['X']]
chrom_sizes = {}
for line in open(CHROM_SIZES):
    [ch, length] = line.strip().split()
    if ch in chroms:
        chrom_sizes[ch] = int(length)

HIC_FILES = {
    'O3': '/nfs/turbo/umms-drjieliu/proj/4dn/data/bulkHiC/haley_HSPC_hic/O3.hic',
    'O3_HOXA9KO': '/nfs/turbo/umms-drjieliu/proj/4dn/data/bulkHiC/haley_HSPC_hic/O3_HOXA9KO.allValidPairs.hic',
    'O3_24h': '/nfs/turbo/umms-drjieliu/usr/temp_Fan/AML_Analysis/HiC/O3_24h/04.4DN_alignment/results/DW_HiC_24h.hic',
    # Merged 4DN maps (multiple libraries per biological sample, cross-library
    # deduplicated at the .pairs level). These supersede the two single-library
    # maps above: O3_HOXA9KO_merged >= O3_HOXA9KO, O3_24h_dTAG13_merged >= O3_24h.
    'O3_HOXA9KO_merged': '/nfs/turbo/umms-drjieliu/usr/temp_Fan/AML_Analysis/HiC/O3_4DN_merged/O3_HOXA9KO/O3_HOXA9KO.hic',
    'O3_24h_dTAG13_merged': '/nfs/turbo/umms-drjieliu/usr/temp_Fan/AML_Analysis/HiC/O3_4DN_merged/O3_24h_dTAG13/O3_24h_dTAG13.hic',
}


def window_scores(mat, max_size, window, ch, resolution):
    rows = []
    for j in range(max_size):
        if j < window or j + window >= max_size:
            rows.append([ch, j * resolution, (j + 1) * resolution, pd.NA])
            continue
        sub_mat = mat[j-window:j+window+1, j-window:j+window+1]
        _sum = np.sum(sub_mat)                              # S1
        _gap = 2 * np.sum(sub_mat[window+1:, :window])       # S2
        _non_gap = _sum - _gap
        if _gap == 0 or _non_gap == 0:                       # undefined
            score = pd.NA
        else:
            score = np.log2(_gap / (_non_gap + 1e-8))
        rows.append([ch, j * resolution, (j + 1) * resolution, score])
    return rows


def insulation(hic_file, tag, kind, chroms=chroms, resolution=RESOLUTION, windows=WINDOWS):
    tables = {w: [] for w in windows}
    for ch in chroms:
        tmp = os.path.join(OUT_DIR, f'temp_{tag}_{kind}_{ch}.txt')
        cmd = ["java", "-jar", JUICER, "dump", kind, "KR", hic_file,
               ch, ch, "BP", str(resolution), tmp]
        print(f'[{tag}/{kind}] dumping {ch}', flush=True)
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0 or not os.path.exists(tmp):
            print(f'[{tag}/{kind}] dump FAILED for {ch}: {result.stderr[-500:]}', flush=True)
            continue

        max_size = chrom_sizes[ch] // resolution + 1
        mat = np.zeros((max_size, max_size))
        for line in open(tmp):
            lst = line.strip().split()
            if lst[2].lower() == 'nan':
                continue
            p1, p2, v = int(lst[0]) // resolution, int(lst[1]) // resolution, float(lst[2])
            if p1 != p2:                      # diagonal excluded from S1
                mat[p1, p2] += v
                mat[p2, p1] += v
        os.remove(tmp)

        for w in windows:
            tables[w].extend(window_scores(mat, max_size, w, ch, resolution))

    for w in windows:
        df = pd.DataFrame(tables[w], columns=['chr', 'start', 'end', 'score'])
        s = pd.to_numeric(df['score'], errors='coerce')
        df['score_z'] = (s - s.mean()) / s.std(ddof=0)
        out = os.path.join(OUT_DIR, f'{tag}_insulation_localNorm_window{w*10}Kb_{kind}.csv')
        df.to_csv(out, index=False)
        print(f'[{tag}/{kind}] wrote {out}  '
              f'(valid={s.notna().sum()}, mean={s.mean():.3f}, %pos={(s > 0).mean():.2%})', flush=True)


if __name__ == '__main__':
    tag = sys.argv[1]
    kind = sys.argv[2] if len(sys.argv) > 2 else 'oe'
    insulation(HIC_FILES[tag], tag=tag, kind=kind)
