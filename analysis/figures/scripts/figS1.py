#!/usr/bin/env python3
"""Supplementary Fig. S1, the full read simulation (every estimator and condition). Output figures/FigS1.pdf/.svg.
Reads only figures/data/. Data file per panel in the comment above it.
Text at 7 pt (6.5 pt annotations); distinct markers for the two branches; ratios in e from the metrics table at three decimals."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from fig2_style import P, LAYOUT, GREY, MM, DATA, FS_ANN, setup, mstyle, panel_letter, save

setup()
rd = lambda f: pd.read_csv(f'{DATA}/{f}', sep='\t', comment='#')
LAYS = ['se50', 'se101', 'pe']
M6 = ['hybrid', 'vf', 'emg', 'salmon_pan', 'salmon_ref', 'ref']
ps = rd('fig2a_e2_per_sample_error.tsv'); ps = ps[ps.tier == 'core']
mt = rd('fig2a_e2_metrics_all.tsv'); mt = mt[(mt.truth == 'expected') & (mt.tier == 'core')]
fig = plt.figure(figsize=(180 * MM, 215 * MM))
gs = fig.add_gridspec(3, 1, height_ratios=[1, 0.95, 1.0], hspace=0.75, left=0.085, right=0.99, top=0.955, bottom=0.115)
rng = np.random.default_rng(2)
hl = [Line2D([], [], ls='none', label=P[c][0], marker=P[c][2], ms=5.2 if c == 'hybrid' else 4.2, mfc=P[c][1] if P[c][3] else 'white',
             mec='white' if c == 'hybrid' else P[c][1], mew=0.5 if c == 'hybrid' else 0.8) for c in M6]
fig.legend(handles=hl, loc='upper center', ncol=6, frameon=False, bbox_to_anchor=(0.53, 0.998), columnspacing=1.4, handletextpad=0.3)

# a: per-sample error of every estimator, clean condition (data: fig2a_e2_per_sample_error.tsv)
ga = gs[0].subgridspec(1, 3, wspace=0.1); axs = []
for j, lay in enumerate(LAYS):
    ax = fig.add_subplot(ga[j], sharey=axs[0] if axs else None); axs.append(ax)
    for i, c in enumerate(M6):
        v = ps[(ps.cond == 'clean') & (ps.layout == lay) & (ps.method == c)].err_B6.values
        ax.scatter(i + (rng.random(len(v)) - 0.5) * 0.4, v, **mstyle(c, 9))
        ax.plot([i - 0.3, i + 0.3], [np.median(v)] * 2, color=P[c][1], lw=1.0)
    ax.set_xticks(range(6)); ax.set_xticklabels([P[c][0] for c in M6], rotation=35, ha='right', rotation_mode='anchor')
    ax.set_title(LAYOUT[lay])
    if j:
        plt.setp(ax.get_yticklabels(), visible=False)
    else:
        ax.set_ylabel('Median |log$_2$ error| per sample\n(B6 sample-genes)')
axs[0].set_ylim(0.04, 0.095)
panel_letter(axs[0], 'a', x=-0.3)

# b-c: within-gene r2 and spurious strain-effect rate per estimator and layout (data: fig2a_e2_metrics_all.tsv)
gbc = gs[1].subgridspec(1, 2, wspace=0.25)
off = dict(zip(M6, np.linspace(-0.3, 0.3, 6)))
for k, (col, ylab) in enumerate([('M1_within_r2', 'Within-gene r² (median over genes)'), ('M4_spurious', 'Spurious strain-effect rate')]):
    ax = fig.add_subplot(gbc[k])
    for gi, lay in enumerate(LAYS):
        for c in M6[::-1]:
            r = mt[(mt.cond == 'clean') & (mt.layout == lay) & (mt.method == c)].iloc[0]
            v = r[col] * (100 if k else 1)
            ax.scatter([gi + off[c]], [v], **mstyle(c, 12))
    ax.set_xticks(range(3)); ax.set_xticklabels([LAYOUT[l].replace(', ', ',\n') for l in LAYS])
    ax.set_xlim(-0.5, 2.5); ax.set_ylabel(ylab)
    if k == 0:
        ax.set_ylim(0.9, 0.97)
        panel_letter(ax, 'b', x=-0.17)
    else:
        ax.axhline(5, color=GREY, lw=0.6, ls=(0, (3, 2)))
        ax.set_yscale('log'); ax.set_ylim(0.5, 40); ax.set_yticks([1, 5, 10, 30]); ax.set_yticklabels(['1%', '5%', '10%', '30%'])
        ax.text(-0.45, 0.62, 'PLUCK vs Salmon-pan at 50 bp: 11 vs 3 of 340 genes', fontsize=FS_ANN, ha='left', va='center')
        ax.text(1.0 + off['hybrid'], 3.3, '6 vs 3\nof 224', fontsize=FS_ANN, ha='center', va='bottom')
        ax.text(2.0 + off['hybrid'], 3.3, '6 vs 3\nof 224', fontsize=FS_ANN, ha='center', va='bottom')
        panel_letter(ax, 'c', x=-0.17)

# d: pre-mRNA (15% unspliced fragments, 101 bp PE): r2 clean -> pre-mRNA (data: fig2a_e2_metrics_all.tsv)
gde = gs[2].subgridspec(1, 4, width_ratios=[0.95, 0.75, 0.75, 0.75], wspace=0.3)
ax = fig.add_subplot(gde[0])
for i, c in enumerate(M6):
    a = mt[(mt.cond == 'clean') & (mt.layout == 'pe') & (mt.method == c)].M1_within_r2.iloc[0]
    b = mt[(mt.cond == 'premrna') & (mt.layout == 'pe') & (mt.method == c)].M1_within_r2.iloc[0]
    ax.plot([i, i], [a, b], color=P[c][1], lw=0.8)
    ax.scatter([i], [a], **mstyle(c, 10)); ax.scatter([i], [b], marker='_', s=30, c=P[c][1], lw=1.2, zorder=4)
ax.set_xticks(range(6)); ax.set_xticklabels([P[c][0] for c in M6], rotation=40, ha='right', rotation_mode='anchor')
ax.set_ylabel('Within-gene r²'); ax.set_ylim(0.9, 0.97)
ax.set_title('Pre-mRNA, 101 bp PE:\nclean (marker), 15% (bar)')
panel_letter(ax, 'd', x=-0.36)

# e: bias arm, error per sample against bias strength a_s = s/7 (data: fig2a_e2_per_sample_error.tsv)
M7 = ['hybrid', 'vf', 'emg', 'ref', 'salmon_pan', 'salmon_pan_flags', 'salmon_pan_flags_tpm']
axe = []
for k, (cond, ttl) in enumerate([('gc', 'GC bias'), ('pos3', '3′ positional bias'), ('hex', 'Hexamer bias')]):
    ax = fig.add_subplot(gde[k + 1], sharey=axe[0] if axe else None); axe.append(ax)
    for c in M7[::-1]:
        d = ps[(ps.cond == cond) & (ps.layout == 'pe') & (ps.method == c)].sort_values('sample')
        a = d['sample'].values / 7
        sty = 'salmon_pan' if c.startswith('salmon_pan') else c
        ls = {'salmon_pan': ':', 'salmon_pan_flags': '--', 'salmon_pan_flags_tpm': '-'}.get(c, P[sty][5])
        ax.plot(a, d.err_B6, color=P[sty][1], lw=1.3 if c == 'hybrid' else 0.8, ls=ls, zorder=4 if c == 'hybrid' else 2)
        if c not in ('salmon_pan', 'salmon_pan_flags'):   # markers on every estimator line except the two Salmon NumReads lines
            kw = mstyle(sty, 6)
            if c == 'salmon_pan_flags_tpm':
                kw.update(marker='D')
            ax.scatter(a, d.err_B6, **kw)
    m = mt[(mt.cond == cond) & (mt.layout == 'pe')].set_index('method').M2_err_B6   # ratios at 3 dp, as stored
    ax.text(0.03, 0.97, f'PLUCK / Salmon-pan flags:\nNumReads {m.hybrid / m.salmon_pan_flags:.3f}×\n'
            f'TPM† {m.hybrid / m.salmon_pan_flags_tpm:.3f}×', fontsize=FS_ANN, transform=ax.transAxes, va='top')
    ax.set_title(f'{ttl}, 101 bp PE'); ax.set_xlabel('Bias strength $a_s$')
    ax.set_xlim(-0.03, 1.03)
    if k:
        plt.setp(ax.get_yticklabels(), visible=False)
    else:
        ax.set_ylabel('Median |log$_2$ error|')
axe[0].set_ylim(0.04, 0.42)
panel_letter(axe[0], 'e', x=-0.3)
lh = [Line2D([], [], color=P['salmon_pan'][1], ls=':', lw=0.8, label='Salmon-pan, no flags'),
      Line2D([], [], color=P['salmon_pan'][1], ls='--', lw=0.8, label='Salmon-pan, flags, NumReads'),
      Line2D([], [], color=P['salmon_pan'][1], ls='-', lw=0.8, marker='D', ms=3, label='Salmon-pan, flags, TPM × eff. length (†)')]
fig.legend(handles=lh, loc='lower center', bbox_to_anchor=(0.55, 0.0), ncol=3, frameon=False, fontsize=FS_ANN, columnspacing=1.2)
fig.canvas.draw()
r = fig.canvas.get_renderer()
for a in fig.axes:
    bb = a.get_tightbbox(r)
    if bb.x0 < 0 or bb.x1 > fig.bbox.x1:
        print('outside', a.get_title(), round(bb.x0), round(bb.x1), round(fig.bbox.x1))
for l in fig.legends + [a.get_legend() for a in fig.axes if a.get_legend()]:
    bb = l.get_window_extent(r); print('legend', round(bb.x0), round(bb.x1))
save(fig, 'FigS1')
print('saved FigS1')
