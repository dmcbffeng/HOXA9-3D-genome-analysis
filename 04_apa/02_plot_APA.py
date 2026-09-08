#!/usr/bin/env python
"""
Plot APA heatmaps from the 11 x 11 mean OE matrices produced by 01_compute_APA.py.

For each of the 4 maps, writes two PDFs:
  <stem>_APA_raw.pdf        : clean Reds heatmap, center enrichment annotated
                              top-right (matches the reference Figure 2B style).
  <stem>_APA_annotated.pdf  : same heatmap with every cell's OE value printed.

Colour scale is shared across all 4 maps (vmin=0, vmax=common) so the panels are
directly comparable. The center (apex) OE enrichment is shown top-right in red.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUTDIR = os.path.dirname(os.path.abspath(__file__))
WINDOW = 5
C = WINDOW

MAPS = ['O3_WT', 'O3_HOXA9KO', 'O3_24h_dTAG13', 'O3_24h_dTAG13_downsampled', 'HSPC']
TITLES = {
    'O3_WT': 'O3 WT',
    'O3_HOXA9KO': 'O3 HOXA9 KO',
    'O3_24h_dTAG13': 'O3 24h dTAG',
    'O3_24h_dTAG13_downsampled': 'O3 24h dTAG (downsampled)',
    'HSPC': 'HSPC',
}

mats = {m: np.load(f'{OUTDIR}/{m}_APA_matrix.npy') for m in MAPS}
vmax = max(4.0, np.ceil(max(m.max() for m in mats.values())))
vmin = 0.0
extent = [-WINDOW, WINDOW, -WINDOW, WINDOW]


def base_heatmap(ax, mat):
    im = ax.imshow(mat, cmap='Reds', vmin=vmin, vmax=vmax,
                   origin='lower', extent=extent, interpolation='nearest')
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_linewidth(0.8)
    return im


for m in MAPS:
    mat = mats[m]
    center = mat[C, C]

    # ---- raw heatmap ----
    fig, ax = plt.subplots(figsize=(3.2, 3.0))
    im = base_heatmap(ax, mat)
    ax.set_title(TITLES[m], fontsize=11)
    ax.text(0.97, 0.95, f'{center:.2f}', transform=ax.transAxes,
            ha='right', va='top', fontsize=12, color='black', fontweight='bold')
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, ticks=[vmin, vmax])
    cb.ax.tick_params(labelsize=9)
    fig.tight_layout()
    fig.savefig(f'{OUTDIR}/{m}_APA_raw.pdf', bbox_inches='tight')
    plt.close(fig)

    # ---- number-annotated heatmap ----
    fig, ax = plt.subplots(figsize=(5.2, 4.8))
    im = base_heatmap(ax, mat)
    ax.set_title(f'{TITLES[m]}   (center OE = {center:.2f})', fontsize=11)
    n = mat.shape[0]
    for i in range(n):        # row -> y (origin lower)
        for j in range(n):    # col -> x
            v = mat[i, j]
            txt_color = 'white' if v > 0.6 * vmax else 'black'
            ax.text(j - WINDOW, i - WINDOW, f'{v:.1f}',
                    ha='center', va='center', fontsize=6, color=txt_color)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, ticks=[vmin, vmax])
    cb.ax.tick_params(labelsize=9)
    fig.tight_layout()
    fig.savefig(f'{OUTDIR}/{m}_APA_annotated.pdf', bbox_inches='tight')
    plt.close(fig)

    print(f'{m}: center OE={center:.3f}  vmax={vmax}  -> raw + annotated PDF')

print('done')
