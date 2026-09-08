#!/usr/bin/env python
"""
Aggregate Peak Analysis (APA) for the 4 O3 contact maps.

Loops: the 5,824 O3_WT FitHiChIP interactions (the WT loop set) are used for ALL
four maps, so every panel piles up the SAME anchor pairs -> a fair comparison of
how strongly each contact map supports the WT loop calls.

Pileup value: observed/expected on KR-normalized matrices (`juicer dump oe KR`),
10 kb resolution, +/- 5 bins around each loop apex -> an 11 x 11 mean matrix.
Periphery ~1 (expected background); the center pixel is the fold-enrichment at the
loop apex (the number annotated top-right in the reference figure).

Contact maps (merged 4DN maps, same .hic used by the FitHiChIP runs; WT = the
OCI-AML3 bulk Hi-C that the precomputed WT locus-pair file was derived from):
  O3_WT                     -> haley_HSPC_hic/O3.hic
  O3_HOXA9KO                -> O3_4DN_merged/O3_HOXA9KO/O3_HOXA9KO.hic
  O3_24h_dTAG13             -> O3_4DN_merged/O3_24h_dTAG13/O3_24h_dTAG13.hic
  O3_24h_dTAG13_downsampled -> O3_4DN_downsampled_O3WT/O3_24h_dTAG13_dsWT/...hic
OE-KR normalization removes the 60x raw-count inflation of O3.hic (raw/KR-count
analyses would need /60, but OE is unaffected).
"""
import os
import subprocess
import numpy as np
import pandas as pd

RES = 10000
WINDOW = 5                       # +/- 5 bins -> 11 x 11
CHROMS = [f'chr{i}' for i in list(range(1, 23)) + ['X']]
JUICER = '/nfs/turbo/umms-drjieliu/juicer_tools_1.11.04_jcuda.0.8.jar'

B = '/nfs/turbo/umms-drjieliu/usr/temp_Fan/AML_Analysis'
ALLDIR = f'{B}/FitHiChIP_202608_newO3_dsWT/allFitHiChIP'
OUTDIR = f'{ALLDIR}/APA'
LOOP_FILE = f'{ALLDIR}/O3_WT_HOXA9_10000_0.05/FitHiChIP.interactions_FitHiC_Q0.05_MergeNearContacts.bed'
CHROM_SIZES = f'{B}/FitHiChIP/Nov21_all_summary/hg19.chrom.sizes'
TMP = os.environ.get('APA_TMP', '/tmp/apa_temp.txt')


def load_chrom_bins():
    d = {}
    for line in open(CHROM_SIZES):
        c = line.split()
        if len(c) >= 2:
            d[c[0]] = int(c[1]) // RES      # max bin index for the chromosome
    return d

MAPS = [
    ('O3_WT',                     f'/nfs/turbo/umms-drjieliu/proj/4dn/data/bulkHiC/haley_HSPC_hic/O3.hic'),
    ('O3_HOXA9KO',                f'{B}/HiC/O3_4DN_merged/O3_HOXA9KO/O3_HOXA9KO.hic'),
    ('O3_24h_dTAG13',             f'{B}/HiC/O3_4DN_merged/O3_24h_dTAG13/O3_24h_dTAG13.hic'),
    ('O3_24h_dTAG13_downsampled', f'{B}/HiC/O3_4DN_downsampled_O3WT/O3_24h_dTAG13_dsWT/O3_24h_dTAG13_dsWT.hic'),
    ('HSPC',                      f'/nfs/turbo/umms-drjieliu/proj/4dn/data/bulkHiC/haley_HSPC_hic/GSM2861708_HIC1258_aligned_inter_30.hic'),
]


def load_loops(loop_file):
    df = pd.read_csv(loop_file, sep='\t')
    loops = {}
    for _, r in df.iterrows():
        ch = r['chr1']
        p1, p2 = int(r['s1']) // RES, int(r['s2']) // RES
        if p1 > p2:
            p1, p2 = p2, p1
        loops.setdefault(ch, []).append((p1, p2))
    return loops


def apa(hic_file, loops, chrom_bins):
    """Memory-light pileup: keep only OE pixels that fall inside a loop window.

    The full-chromosome dense matrix is ~30k x 30k (7 GB) which OOMs, so instead
    we (1) list the pixels needed by all loop windows on the chromosome, then
    (2) stream `juicer dump oe KR` keeping only those pixels in a dict. The Hi-C
    matrix is symmetric; the dump emits i<=j, so pixels are keyed as (min,max).
    """
    n = 2 * WINDOW + 1
    pile_sum = np.zeros((n, n))
    used = 0
    for ch in CHROMS:
        if ch not in loops:
            continue
        max_bin = chrom_bins.get(ch, 0)

        # valid loops on this chromosome + the set of pixel keys we will need
        valid = []
        needed = set()
        for (p1, p2) in loops[ch]:
            if p1 < WINDOW or (max_bin and p2 + WINDOW > max_bin):
                continue
            valid.append((p1, p2))
            for di in range(-WINDOW, WINDOW + 1):
                for dj in range(-WINDOW, WINDOW + 1):
                    r, c = p1 + di, p2 + dj
                    needed.add((r, c) if r <= c else (c, r))
        if not valid:
            print(f'  {ch}: no valid loops, skipped')
            continue

        cmd = ["java", "-jar", JUICER, "dump", "oe", "KR", hic_file,
               ch, ch, "BP", str(RES), TMP]
        if os.path.exists(TMP):
            os.remove(TMP)
        res = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if res.returncode != 0 or not os.path.exists(TMP):
            print(f'  {ch}: dump failed, skipped')
            continue

        oe = {}
        for line in open(TMP):
            lst = line.split()
            if lst[2].lower() == 'nan':
                continue
            a, b = int(lst[0]) // RES, int(lst[1]) // RES
            key = (a, b) if a <= b else (b, a)
            if key in needed:
                oe[key] = oe.get(key, 0.0) + float(lst[2])
        os.remove(TMP)
        if not oe:
            print(f'  {ch}: empty KR matrix, skipped')
            continue

        cnt = 0
        for (p1, p2) in valid:
            sub = np.zeros((n, n))
            for di in range(-WINDOW, WINDOW + 1):
                for dj in range(-WINDOW, WINDOW + 1):
                    r, c = p1 + di, p2 + dj
                    sub[di + WINDOW, dj + WINDOW] = oe.get((r, c) if r <= c else (c, r), 0.0)
            pile_sum += sub
            cnt += 1
        used += cnt
        print(f'  {ch}: {cnt} loops piled')
    mean_mat = pile_sum / used if used else pile_sum
    return mean_mat, used


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    loops = load_loops(LOOP_FILE)
    chrom_bins = load_chrom_bins()
    n_total = sum(len(v) for v in loops.values())
    print(f'WT loop set: {n_total} loops on {len(loops)} chromosomes\n')

    summary = []
    c = WINDOW  # center index
    for stem, hic in MAPS:
        print(f'=== {stem} : {hic} ===')
        mean_mat, used = apa(hic, loops, chrom_bins)
        np.save(f'{OUTDIR}/{stem}_APA_matrix.npy', mean_mat)
        np.savetxt(f'{OUTDIR}/{stem}_APA_matrix.csv', mean_mat,
                   delimiter=',', fmt='%.6f')
        center = mean_mat[c, c]
        # P2LL: center vs mean of lower-left 3x3 corner (Rao et al. 2014)
        ll = mean_mat[-3:, :3].mean()
        p2ll = center / ll if ll else np.nan
        summary.append([stem, used, center, p2ll, mean_mat.min(), mean_mat.max()])
        print(f'  -> n_used={used}  center(OE)={center:.3f}  P2LL={p2ll:.3f}\n')

    sdf = pd.DataFrame(summary, columns=['map', 'n_loops', 'center_OE_enrichment',
                                         'P2LL', 'matrix_min', 'matrix_max'])
    sdf.to_csv(f'{OUTDIR}/APA_summary.csv', index=False)
    print(sdf.to_string(index=False))


if __name__ == '__main__':
    main()
