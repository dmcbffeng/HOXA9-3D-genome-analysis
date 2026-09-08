#!/usr/bin/env python
"""Dump a .hic map to FitHiChIP 7-col locus-pair bed (chr s1 e1 chr s2 e2 cc).
Intra-chromosomal (cis) only, raw observed counts ('NONE'), at the given resolution,
chr-prefixed. Matches the format of FitHiChIP/DS_HSPC_v2/O3_10000.bed."""
import sys, hicstraw
HIC = sys.argv[1]
OUT = sys.argv[2]
RES = int(sys.argv[3]) if len(sys.argv) > 3 else 10000
CHROMS = [str(i) for i in range(1, 23)] + ['X', 'Y']
hic = hicstraw.HiCFile(HIC)
assert RES in hic.getResolutions(), f'{RES} not in {hic.getResolutions()}'
n = 0
with open(OUT, 'w') as f:
    for ch in CHROMS:
        rec = hicstraw.straw('observed', 'NONE', HIC, ch, ch, 'BP', RES)
        for r in rec:
            c = int(r.counts)
            if c <= 0:
                continue
            f.write(f'chr{ch}\t{int(r.binX)}\t{int(r.binX)+RES}\t'
                    f'chr{ch}\t{int(r.binY)}\t{int(r.binY)+RES}\t{c}\n')
            n += 1
        print(f'  chr{ch}: done', flush=True)
print(f'total lines: {n}', flush=True)
