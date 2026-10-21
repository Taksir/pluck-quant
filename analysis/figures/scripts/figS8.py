#!/usr/bin/env python3
"""Fig. S8: the reliability flag per cohort-tissue, against the secondary models and the plain baselines,
and how flag rate, recall and the unflagged tail move with the u_g threshold.
Reads only figures/data/. Writes figures/FigS8.pdf and FigS8.svg."""
import numpy as np, pandas as pd
from u5_style import *

setup()
rd = lambda f: pd.read_csv(f'{DATA}/{f}', sep='\t', comment='#')
fig = plt.figure(figsize=(180 * MM, 188 * MM))
ORANGE = PIPE['PLUCK'][0]

# a: data fig5a_roc_curves_u5.tsv (curves and AUCs), fig5a_operating_points_u5.tsv (u_g < 0.6 and < 0.7)
roc = rd('fig5a_roc_curves_u5.tsv'); ops = rd('fig5a_operating_points_u5.tsv')
MODELS = [('U', '$u_g$ (design flag, no fit)', ORANGE, '-', 1.1), ('compact', 'compact logistic (8 design features)', 'black', '--', 0.8),
          ('postcount', 'post-count logistic (+ expression, branch)', 'black', ':', 0.9), ('B2', 'baseline + multi-copy fraction', GREY_B2, '-', 0.7),
          ('B1', 'gene-feature baseline: length, transcripts, TSL1, GC, expression', GREY_B1, '-', 0.7)]
W, H, X0, Y0, GX, GY = 42, 42, 14, 8, 14, 15
for k, c in enumerate(COHORT_ORDER):
    col, row = k % 3, k // 3
    ax = ax_mm(fig, X0 + col * (W + GX), Y0 + row * (H + GY), W, H)
    ax.plot([0, 1], [0, 1], color='#CCCCCC', lw=0.5, ls=':')
    txt = []
    for m, lab, colr, ls, lw in MODELS[::-1]:
        r = roc[(roc.cohort == c) & (roc.model == m)]
        cur = r[r.fpr >= 0]; a = r[r.fpr < 0].tpr.iloc[0]
        ax.plot(cur.fpr, cur.tpr, color=colr, ls=ls, lw=lw, zorder=4 if m == 'U' else 3)
    for m, lab, colr, ls, lw in MODELS:
        a = roc[(roc.cohort == c) & (roc.model == m) & (roc.fpr < 0)].tpr.iloc[0]
        txt.append(({'U': '$u_g$', 'compact': 'compact', 'postcount': 'post-count', 'B1': 'gene-feature baseline', 'B2': 'baseline + multi-copy'}[m], a, colr))
    for thr, face in ((0.6, ORANGE), (0.7, 'white')):
        o = ops[(ops.cohort == c) & (ops.threshold == thr)].iloc[0]
        ax.plot(o.fpr, o.tpr_recall, 'o', ms=4, mfc=face, mec=ORANGE if face == 'white' else 'white', mew=0.8, zorder=6)
    for j, (lab, a, colr) in enumerate(txt):
        ax.text(0.97, 0.04 + j * 0.085, '%s %.3f' % (lab, a), ha='right', va='bottom', color=colr, bbox=dict(fc='white', ec='none', pad=0.3, alpha=0.85))
    ax.text(0.97, 0.04 + 5 * 0.085, 'AUC', ha='right', va='bottom', fontweight='bold')
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect('equal')
    ax.set_xticks([0, 0.5, 1]); ax.set_yticks([0, 0.5, 1]); ax.set_xticklabels(['0', '0.5', '1']); ax.set_yticklabels(['0', '0.5', '1'])
    ax.set_title(COHORT_LABEL[c] + ('  (in-sample)' if c == 'keele_liver' else ''), loc='left', pad=2)
    if col == 0: ax.set_ylabel('True-positive rate')
    if row == 1: ax.set_xlabel('False-positive rate')
letter(fig, 1, 4, 'a')
ky = Y0 + 2 * H + GY + 12.5
for j, (m, lab, colr, ls, lw) in enumerate(MODELS):
    xx = 14 + (j % 2) * 84; yy = ky + (j // 2) * 3.6
    lax = ax_mm(fig, xx, yy - 1.5, 7, 3); lax.axis('off'); lax.set_xlim(0, 1); lax.set_ylim(0, 1)
    lax.plot([0, 1], [0.5, 0.5], color=colr, ls=ls, lw=lw + 0.2)
    mtext(fig, xx + 8, yy, lab, va='center')
xx = 14; yy = ky + 3 * 3.6
mtext(fig, xx, yy, '● $u_g$ < 0.6 (primary)   ○ $u_g$ < 0.7 (sensitivity)', va='center', color=ORANGE)

# b: data figS8_threshold_sweep_u5.tsv
sw = rd('figS8_threshold_sweep_u5.tsv')
STY = {'keele_kidney': ('black', '-'), 'keele_liver': ('black', '--'), 'munger': ('black', ':'),
       'gse199702_liver': ('#8A8A8A', '-'), 'gse199702_heart': ('#8A8A8A', '--'), 'gse199702_kidney': ('#8A8A8A', ':')}
by = 140
for j, (col, lab, ylim, yt) in enumerate((('flag_rate', 'Genes flagged', (0, 0.6), [0, 0.2, 0.4, 0.6]),
                                         ('recall', 'Recall of the >2× tail', (0, 1), [0, 0.25, 0.5, 0.75, 1]),
                                         ('tail_unflagged', 'Tail rate among unflagged genes', (0, 0.1), [0, 0.025, 0.05, 0.075, 0.1]))):
    ax = ax_mm(fig, 14 + j * 56, by, 42, 30)
    for c in COHORT_ORDER:
        d = sw[sw.cohort == c]; colr, ls = STY[c]
        ax.plot(d.threshold, d[col], color=colr, ls=ls, lw=0.9)
    for t, ls in ((0.6, '-'), (0.7, '--')):
        ax.axvline(t, color=ORANGE, lw=0.7, ls=ls, zorder=0)
    ax.set_xlim(0.3, 0.95); ax.set_ylim(*ylim); ax.set_yticks(yt)
    ax.set_yticklabels([('%g%%' % (v * 100)) for v in yt])
    ax.set_xlabel('Threshold on $u_g$ (flag if below)'); ax.set_ylabel(lab)
    if j == 0:
        ax.text(0.6, ylim[1] * 0.97, ' 0.6', color=ORANGE, va='top', ha='left'); ax.text(0.7, ylim[1] * 0.97, ' 0.7', color=ORANGE, va='top', ha='left')
letter(fig, 1, by - 4, 'b')
for k, c in enumerate(COHORT_ORDER):
    xx = 14 + (k % 3) * 56; yy = by + 41 + (k // 3) * 3.6
    colr, ls = STY[c]
    lax = ax_mm(fig, xx, yy - 1.5, 7, 3); lax.axis('off'); lax.set_xlim(0, 1); lax.set_ylim(0, 1)
    lax.plot([0, 1], [0.5, 0.5], color=colr, ls=ls, lw=1.0)
    mtext(fig, xx + 8, yy, COHORT_LABEL[c], va='center')

save(fig, 'FigS8')
print('wrote FigS8')
