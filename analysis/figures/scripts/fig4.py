#!/usr/bin/env python3
"""Figure 4: founder-haplotype shift by founder, tool and read length, and the real-genotype DO check.

Writes figures/Fig4.pdf and Fig4.svg (vector, text as text). Reads only figures/data/.
Colour, marker, order and labels follow the shared figure style; PLUCK is the TSV tool 'hybrid'.
Run: python3 figures/scripts/fig4.py
"""
import json
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                      # noqa: E402
import numpy as np                                   # noqa: E402
import pandas as pd                                  # noqa: E402
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec   # noqa: E402
from matplotlib.lines import Line2D                  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
OUTD = os.path.join(HERE, '..')

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans'],
    'font.size': 7, 'axes.titlesize': 7, 'axes.labelsize': 7, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
    'legend.fontsize': 7, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
    'axes.spines.top': False, 'axes.spines.right': False,
})

# ---- pipeline style (label, colour, marker, filled, size factor); order is the dict order
ENC = {
    'PLUCK':              dict(c='#D55E00', m='o', filled=True, s=1.4),
    'Salmon-pan':         dict(c='#0072B2', m='s', filled=True, s=1.0),
    'Salmon-ref':         dict(c='#56B4E9', m='s', filled=False, s=1.0),
    'STAR+fC':            dict(c='#009E73', m='^', filled=True, s=1.0),
    'kallisto-pan':       dict(c='#882255', m='D', filled=True, s=0.85),
    'kallisto-ref':       dict(c='#CC79A7', m='D', filled=False, s=0.85),
    'Kmerator+countTags': dict(c='#7F7F7F', m='v', filled=True, s=1.0),
    'RNA-Skim':           dict(c='#333333', m='x', filled=True, s=1.0),
    'ref-unique probes':  dict(c='#332288', m='o', filled=False, s=1.0),
    'GEO reference RSEM': dict(c='#555555', m='^', filled=False, s=1.0),
}
TOOL2LABEL = {'hybrid': 'PLUCK'}
MS = 4.2          # base marker size (points)
FOUNDERS = ['A', 'C', 'D', 'E', 'F', 'G', 'H']                       # standard CC order as Fig. 2b: domestic A/J 129S1 NOD NZO, then wild-derived CAST PWK WSB
FNAME = dict(A='A/J', C='129S1', D='NOD', E='NZO', H='WSB', F='CAST', G='PWK')
CAST_PWK_X = (3.5, 5.5)                                               # x span of CAST and PWK


def read(name):
    return pd.read_csv(os.path.join(DATA, name), sep='\t', comment='#')


def lab(tool):
    return TOOL2LABEL.get(tool, tool)


def point(ax, x, y, label, lo=None, hi=None, horizontal=False, z=3, size=1.0, alpha=1.0):
    e = ENC[label]
    ms = MS * e['s'] * size
    face = e['c'] if e['filled'] else 'white'
    edge = 'white' if label == 'PLUCK' else e['c']
    z = 10 if label == 'PLUCK' else z                                  # PLUCK drawn on top
    if lo is not None:
        if horizontal:
            ax.plot([lo, hi], [y, y], color=e['c'], lw=0.9, zorder=z - 0.5, solid_capstyle='butt')
        else:
            ax.plot([x, x], [lo, hi], color=e['c'], lw=0.9, zorder=z - 0.5, solid_capstyle='butt')
    if e['m'] == 'x':
        ax.plot([x], [y], marker='x', ms=ms, mew=1.0, color=e['c'], ls='none', zorder=z, alpha=alpha)
    else:
        ax.plot([x], [y], marker=e['m'], ms=ms, mfc=face, mec=edge if e['filled'] else e['c'],
                mew=0.6 if label == 'PLUCK' else 0.8, ls='none', zorder=z, alpha=alpha)


def ci(s):
    v = json.loads(s) if isinstance(s, str) else s
    return float(v[0]), float(v[1])


def panel_letter(ax, letter, x=-0.12, y=1.06):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=10, fontweight='bold', va='bottom', ha='left')


def cast_pwk_shade(ax, ymin, ymax):
    """Shade CAST and PWK only: WSB is wild-derived but M. m. domesticus and behaves like the domestic founders."""
    ax.axvspan(*CAST_PWK_X, color='#EFEFEF', zorder=0, lw=0)
    ax.text(np.mean(CAST_PWK_X), ymax, 'CAST, PWK', ha='center', va='top', fontsize=7, color='#555555')


# ======================================================================= figure
fig = plt.figure(figsize=(7.09, 7.4))     # 180 x 188 mm; all text 7 pt at this size
gs = GridSpec(2, 2, figure=fig, width_ratios=[1.05, 1.0], height_ratios=[1, 1], left=0.075, right=0.985, top=0.885,
              bottom=0.10, wspace=0.30, hspace=0.42)

# ---------------------------------------------------------------- panel a: Munger founders, per founder
# data: fig4a_munger_founder_delta_frameP.tsv (per founder, scoring set P_kmerator, recomputed by fig4_founder_delta_recompute.py)
ax = fig.add_subplot(gs[0, 0])
pf = read('fig4a_munger_founder_delta_frameP.tsv')
pf = pf[pf.founder != 'all']
order_a = ['PLUCK', 'Salmon-pan', 'Salmon-ref', 'Kmerator+countTags', 'ref-unique probes']
off = dict(zip(order_a, np.linspace(-0.28, 0.28, len(order_a))))
YA = (-0.335, 0.03)
cast_pwk_shade(ax, *YA)
ax.axhline(0, color='black', lw=0.6, zorder=1)
for _, r in pf.iterrows():
    L = lab(r.tool)
    x = FOUNDERS.index(r.founder) + off[L]
    lo, hi = ci(r.delta_ci95)
    point(ax, x, r.median_delta, L, lo, hi)
ax.set_xticks(range(7))
ax.set_xticklabels([FNAME[f] for f in FOUNDERS])
ax.text(5.0, -0.105, 'wild-derived', transform=ax.get_xaxis_transform(), ha='center', va='top', fontsize=7, color='#555555')
ax.annotate('', xy=(3.6, -0.085), xytext=(6.4, -0.085), xycoords=ax.get_xaxis_transform(), textcoords=ax.get_xaxis_transform(),
            arrowprops=dict(arrowstyle='-', color='#999999', lw=0.6), annotation_clip=False)
ax.set_xlim(-0.6, 6.6)
ax.set_ylim(*YA)
ax.set_ylabel('Median founder-haplotype shift Δ (log2)')
ax.set_title('Munger founders, 100 bp (47 livers): Δ per founder', loc='left', pad=4)
panel_letter(ax, 'a', x=-0.13)

# ---------------------------------------------------------------- panel b: domestic vs wild-derived split, every pipeline, three cohorts
# data: fig4b_domestic_wild_split.tsv (scoring set A, all core genes; RNA-Skim from set P_rnaskim, its own gene set, as in Table 2)
sub = GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[0, 1], wspace=0.10)
fb = read('fig4b_domestic_wild_split.tsv')
cohorts_b = [('Keele kidney', 'Keele kidney, 50 bp'), ('Keele liver', 'Keele liver, 50 bp†'), ('Munger founders', 'Munger founders, 100 bp')]
order_b = ['PLUCK', 'Salmon-pan', 'Salmon-ref', 'STAR+fC', 'kallisto-pan', 'kallisto-ref', 'Kmerator+countTags', 'RNA-Skim', 'ref-unique probes']
XB = (-0.30, 0.03)
axes_b = []
for k, (coh, title) in enumerate(cohorts_b):
    axb = fig.add_subplot(sub[0, k])
    axes_b.append(axb)
    axb.axvline(0, color='black', lw=0.6, zorder=1)
    for i, L in enumerate(order_b):
        y = len(order_b) - 1 - i
        r = fb[(fb.cohort == coh) & (fb.tool.map(lab) == L)]
        if len(r) == 0:
            axb.text(XB[0] + 0.01, y, 'not run', va='center', ha='left', fontsize=7, color='#777777')
            continue
        r = r.iloc[0]
        e = ENC[L]
        # domestic (small, faint) -> wild-derived (full marker), joined by a thin line
        axb.plot([r.delta_domestic, r.delta_wild], [y, y], color=e['c'], lw=0.8, zorder=2, alpha=0.8)
        point(axb, r.delta_domestic, y, L, size=0.65, alpha=0.55)
        point(axb, r.delta_wild, y, L)
    axb.set_xlim(*XB)
    axb.set_ylim(-0.7, len(order_b) - 0.3)
    axb.set_yticks(range(len(order_b)))
    axb.set_yticklabels(([L + ('‡' if L == 'RNA-Skim' else '') for L in order_b])[::-1] if k == 0 else [])
    axb.set_xticks([-0.2, -0.1, 0])
    axb.set_xticklabels(['−0.2', '−0.1', '0'])
    axb.set_title(title.replace(', ', ',\n'), fontsize=7, pad=3, linespacing=1.1)
    axb.grid(axis='y', color='#E5E5E5', lw=0.5, zorder=0)
    axb.tick_params(axis='y', length=0)
axes_b[1].set_xlabel('Median Δ (log2): small faint marker, domestic founders;\nlarge marker, wild-derived founders')
panel_letter(axes_b[0], 'b', x=-1.05)

# ---------------------------------------------------------------- panel c: DO dose classes against GEO individualized RSEM
# data: fig4c_do_dose_classes.tsv (median Δ, class against the gene's D0 animals, with CI; placebo rows),
#       fig4c_do_per_gene_hybrid_minus_refunique.tsv (per-gene paired difference; medians checked against fig4c_do_paired_hybrid_minus_refunique.tsv)
subc = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1, 0], width_ratios=[1.75, 1.0], wspace=0.42)
axc = fig.add_subplot(subc[0, 0])
dc = read('fig4c_do_dose_classes.tsv')
order_c = ['PLUCK', 'Salmon-pan', 'Salmon-ref', 'ref-unique probes', 'GEO reference RSEM']
offc = dict(zip(order_c, np.linspace(-0.3, 0.3, len(order_c))))
classes = [('placebo', 'placebo\n(D0 split)'), ('D1', 'D1\nhet-like'), ('D2', 'D2\nhom.\ndivergent')]   # ordered by dose
axc.axhline(0, color='black', lw=0.6, zorder=1)
for j, (cl, _) in enumerate(classes):
    for _, r in dc[dc['class'] == cl].iterrows():
        L = lab(r.tool)
        point(axc, j + offc[L], r.median_delta, L, r.ci_lo, r.ci_hi)
axc.set_xticks(range(len(classes)))
axc.set_xticklabels([t for _, t in classes], fontsize=7)
axc.set_xlim(-0.55, 2.55)
axc.set_ylim(-0.03, 0.245)
axc.set_ylabel('Median Δ vs GEO individualized RSEM (log2)')
axc.set_title('Munger DO, 100 bp (48 animals)', loc='left', pad=4)
panel_letter(axc, 'c', x=-0.30)

axd2 = fig.add_subplot(subc[0, 1])
pg = read('fig4c_do_per_gene_hybrid_minus_refunique.tsv')
YC2 = (-0.25, 0.55)
cls2 = ['D1', 'D2']
vals = [pg.loc[pg['class'] == cl, 'diff'].to_numpy() for cl in cls2]
axd2.axhline(0, color='black', lw=0.6, zorder=1)
vp = axd2.violinplot(vals, positions=[0, 1], widths=0.75, showextrema=False, points=400)   # density on unclipped values; axis limits below
for b in vp['bodies']:
    b.set_facecolor('#D9D9D9'); b.set_edgecolor('#8C8C8C'); b.set_linewidth(0.5); b.set_alpha(1)
for x, v in enumerate(vals):
    q1, med, q3 = np.percentile(v, [25, 50, 75])
    axd2.plot([x, x], [q1, q3], color='black', lw=2.2, solid_capstyle='butt', zorder=3)
    axd2.plot([x - 0.17, x + 0.17], [med, med], color='white', lw=1.2, zorder=4)
    axd2.text(x, q3 + 0.015, f'{med:.3f}', va='bottom', ha='center', fontsize=7)      # median, above the IQR bar
axd2.set_xticks([0, 1])
axd2.set_xticklabels([f'{cl}\n{100 * (v > 0).mean():.0f}%\n> 0' for cl, v in zip(cls2, vals)], fontsize=7)
axd2.set_xlim(-0.5, 1.75)
axd2.set_ylim(*YC2)
axd2.set_ylabel('Paired difference per gene (log2)')
axd2.set_title('PLUCK − ref-unique', loc='left', pad=4)

# ---------------------------------------------------------------- panel d: DO, per founder
# data: fig4d_do_per_founder.tsv (carrier probability >= 0.35 with a hit, against the gene's D0 animals; no interval in the source)
axd = fig.add_subplot(gs[1, 1])
dd = read('fig4d_do_per_founder.tsv')
offd = dict(zip(order_c, np.linspace(-0.3, 0.3, len(order_c))))
YD = (-0.01, 0.235)
cast_pwk_shade(axd, *YD)
axd.axhline(0, color='black', lw=0.6, zorder=1)
for _, r in dd.iterrows():
    L = lab(r.tool)
    point(axd, FOUNDERS.index(r.founder) + offd[L], r.median_delta, L)
axd.set_xticks(range(7))
axd.set_xticklabels([FNAME[f] for f in FOUNDERS])
axd.set_xlim(-0.6, 6.6)
axd.set_ylim(*YD)
axd.set_ylabel('Median Δ vs GEO individualized RSEM (log2)')
axd.set_title('Munger DO, 100 bp (48 animals): Δ per founder allele', loc='left', pad=4)
axd.text(0.02, 0.97, 'founder carrier probability ≥ 0.35', transform=axd.transAxes, fontsize=7, color='#555555', va='top')
panel_letter(axd, 'd', x=-0.13)

# ---------------------------------------------------------------- shared legend (standard pipeline order)
handles = []
for L in ['PLUCK', 'Salmon-pan', 'Salmon-ref', 'STAR+fC', 'kallisto-pan', 'kallisto-ref', 'Kmerator+countTags', 'RNA-Skim',
          'ref-unique probes', 'GEO reference RSEM']:
    e = ENC[L]
    if e['m'] == 'x':
        h = Line2D([], [], marker='x', ms=MS * e['s'], mew=1.0, color=e['c'], ls='none', label=L)
    else:
        h = Line2D([], [], marker=e['m'], ms=MS * e['s'], mfc=e['c'] if e['filled'] else 'white',
                   mec='white' if L == 'PLUCK' else e['c'], mew=0.6 if L == 'PLUCK' else 0.8, ls='none', label=L)
        if L == 'PLUCK':
            h.set_markeredgecolor('white')
    handles.append(h)
fig.legend(handles=handles, loc='upper center', ncol=5, frameon=False, bbox_to_anchor=(0.53, 0.995), columnspacing=1.4,
           handletextpad=0.3, borderaxespad=0.2)
fig.text(0.075, 0.003, '† tuning cohort.  ‡ RNA-Skim scored on its own gene frame (protein-coding index only).\n'
         'Error bars: 95% bootstrap intervals over genes, mostly narrower than the markers (none in b and d).', fontsize=7, color='#444444')

for ext in ('pdf', 'svg'):
    fig.savefig(os.path.join(OUTD, f'Fig4.{ext}'))
print('wrote Fig4.pdf, Fig4.svg')
