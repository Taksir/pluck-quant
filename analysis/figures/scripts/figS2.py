#!/usr/bin/env python3
"""Supplementary Fig. S2, the read simulation on the 2,062 extended-tier genes (reported separately). Output figures/FigS2.pdf/.svg.
Reads only figures/data/: fig2a_e2_per_sample_error.tsv (a), fig2a_e2_metrics_all.tsv (b, c)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from fig2_style import P, LAYOUT, GREY, MM, DATA, setup, mstyle, panel_letter, save

setup()
rd = lambda f: pd.read_csv(f'{DATA}/{f}', sep='\t', comment='#')
LAYS = ['se50', 'se101', 'pe']
ps = rd('fig2a_e2_per_sample_error.tsv'); ps = ps[(ps.cond == 'clean')]
mt = rd('fig2a_e2_metrics_all.tsv'); mt = mt[(mt.cond == 'clean') & (mt.truth == 'expected')]
code = {'core': 'hybrid', 'extended': 'v2ext'}
fig = plt.figure(figsize=(180 * MM, 68 * MM))
gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1, 1], wspace=0.38, left=0.07, right=0.99, top=0.84, bottom=0.2)
rng = np.random.default_rng(3)
M3 = ['hybrid', 'salmon_pan', 'salmon_ref']
hl = [Line2D([], [], ls='none', label=P[c][0], marker=P[c][2], ms=5.2 if c == 'hybrid' else 4.2, mfc=P[c][1] if P[c][3] else 'white',
             mec='white' if c == 'hybrid' else P[c][1], mew=0.5 if c == 'hybrid' else 0.8) for c in M3]
hl += [Line2D([], [], ls='none', marker='o', ms=4.2, mfc='white', mec='#444444', label='extended tier (b, c)'),
       Line2D([], [], ls='none', marker='o', ms=4.2, mfc='#444444', mec='#444444', label='core tier (b, c)')]
fig.legend(handles=hl, loc='upper center', ncol=5, frameon=False, bbox_to_anchor=(0.53, 1.0), columnspacing=1.4, handletextpad=0.3)

# a: per-sample error on extended-tier genes (data: fig2a_e2_per_sample_error.tsv, tier == extended)
ax = fig.add_subplot(gs[0])
for gi, lay in enumerate(LAYS):
    for k, m in enumerate(['v2ext', 'salmon_pan', 'salmon_ref']):
        sty = 'hybrid' if m == 'v2ext' else m
        v = ps[(ps.tier == 'extended') & (ps.layout == lay) & (ps.method == m)].err_B6.values
        x = gi * 4 + k
        ax.scatter(x + (rng.random(len(v)) - 0.5) * 0.4, v, **mstyle(sty, 9))
        ax.plot([x - 0.3, x + 0.3], [np.median(v)] * 2, color=P[sty][1], lw=1.0)
ax.set_xticks([gi * 4 + 1 for gi in range(3)]); ax.set_xticklabels([LAYOUT[l].replace(', ', ',\n') for l in LAYS])
ax.set_ylabel('Median |log$_2$ error| per sample\n(extended-tier B6 sample-genes)'); ax.set_ylim(0.05, 0.095)
ax.set_title('Extended tier, 2,062 genes')
panel_letter(ax, 'a', x=-0.25)

# b: within-gene r2, core -> extended tier; c: error ratio to Salmon-pan, core -> extended (data: fig2a_e2_metrics_all.tsv)
for k, (col, ylab, lim) in enumerate([('M1_within_r2', 'Within-gene r²', (0.9, 0.97)),
                                       ('M2_ratio_to_pan', 'Error ratio to Salmon-pan', (0.9, 1.65))]):
    ax = fig.add_subplot(gs[k + 1])
    for gi, lay in enumerate(LAYS):
        for j, c in enumerate(M3):
            if k and c == 'salmon_pan':
                continue
            x = gi + (j - 1) * 0.25
            a = mt[(mt.layout == lay) & (mt.tier == 'core') & (mt.method == c)][col].iloc[0]
            b = mt[(mt.layout == lay) & (mt.tier == 'extended') & (mt.method == ('v2ext' if c == 'hybrid' else c))][col].iloc[0]
            ax.plot([x, x], [a, b], color=P[c][1], lw=0.8)
            ax.scatter([x], [a], **mstyle(c, 11))
            kw = mstyle(c, 11); kw.update(facecolors='white', edgecolors=P[c][1], linewidths=0.8)
            ax.scatter([x], [b], **kw)
    ax.set_xticks(range(3)); ax.set_xticklabels([LAYOUT[l].replace(', ', ',\n') for l in LAYS])
    ax.set_ylim(*lim); ax.set_ylabel(ylab); ax.set_xlim(-0.5, 2.5)
    if k:
        ax.axhline(1, color='black', lw=0.5)
    ax.set_title('Core (filled) to extended (hollow)')
    panel_letter(ax, 'bc'[k], x=-0.28)
save(fig, 'FigS2')
print('saved FigS2')
