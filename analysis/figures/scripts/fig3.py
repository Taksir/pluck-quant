#!/usr/bin/env python3
"""Fig. 3: held-out real-cohort accuracy of PLUCK against seven pipelines.
Reads only figures/data/ (file named above each panel). Writes figures/Fig3.pdf and Fig3.svg.
Run: python3 -I figures/scripts/fig3.py   (needs fig3_compute.py to have been run once)"""
import sys
# Folder constants below point to the original working folder; set them for your copy.
sys.path.insert(0, '/mnt/project-files/paper/figures/scripts')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, NullFormatter
from scipy.stats import gaussian_kde
from u3_encoding import (MM, ORDER, STYLE, COHORTS, COHORT_LABEL, NAME2KEY, setup, read, label, mk, ax_mm, save)

setup()
W, H = 180, 220
fig = plt.figure(figsize=(W * MM, H * MM))
LBL = 27            # left label column (mm)
NROW = len(ORDER)   # 9 pipelines
TSVNAME = {'PLUCK': 'hybrid'}


def letter(x, y, s):
    fig.text(x / W, 1 - y / H, s, fontsize=10, fontweight='bold', va='top', ha='left')


def rows_axis(ax, labels=True):
    ax.set_ylim(NROW - 0.5, -0.5)
    ax.set_yticks(range(NROW))
    ax.set_yticklabels(ORDER if labels else [])
    ax.tick_params(axis='y', length=0, pad=2)
    for i in range(NROW):
        if i % 2 == 0:
            ax.axhspan(i - 0.5, i + 0.5, color='#F3F3F3', zorder=0, lw=0)


def style_ylabels(ax):
    for t in ax.get_yticklabels():
        if t.get_text() == 'PLUCK':
            t.set_fontweight('bold')


# ---------------------------------------------------------------- panels a and b: per-gene distributions
# data: fig3ab_pergene_u3.tsv (per-gene within-gene Spearman rho_<tool> and noise SD sd_<tool>; from fig3_compute.py)
# data: table2_accuracy_rnaskim_frame.tsv (RNA-Skim median only, own gene frame)
# data: table2_accuracy_frameA.tsv (medians of the scoring set with all methods, checked against the per-gene medians below)
pg = read('fig3ab_pergene_u3.tsv')
rs = read('table2_accuracy_rnaskim_frame.tsv')
fa = read('table2_accuracy_frameA.tsv')


def violin_row(ax, vals, y, tool, lo, hi):
    s = STYLE[tool]
    v = vals[np.isfinite(vals)]
    grid = np.linspace(lo, hi, 300)
    kde = gaussian_kde(v[(v > lo - 0.5 * (hi - lo)) & (v < hi + 0.5 * (hi - lo))], bw_method=0.12)
    d = kde(grid)
    d = d / d.max() * 0.42
    ax.fill_between(grid, y + 0.05, y + 0.05 - d, color=s['color'], alpha=0.25 if tool != 'PLUCK' else 0.40, lw=0, zorder=2)
    q10, q25, q50, q75, q90 = np.percentile(v, [10, 25, 50, 75, 90])
    ax.plot([max(q10, lo), min(q90, hi)], [y + 0.18, y + 0.18], color=s['color'], lw=0.6, zorder=3, solid_capstyle='butt')
    ax.plot([max(q25, lo), min(q75, hi)], [y + 0.18, y + 0.18], color=s['color'], lw=2.0, zorder=3, solid_capstyle='butt')
    mk(ax, q50, y + 0.18, tool, ms=3.0)
    return q50


PH = 33   # height of each cohort axis in panels a and b (mm)


def dist_panel(top, metric, lo, hi, ticks, xlabel, fmt, letter_s, ref_rule=True, below_note=True):
    letter(1, top - 8, letter_s)
    for k, c in enumerate(COHORTS):
        x0 = LBL + k * 50
        ax = ax_mm(fig, x0, top, 38, PH, W, H)
        rows_axis(ax, labels=(k == 0))
        if k == 0:
            style_ylabels(ax)
        ax.set_xlim(lo, hi)
        ax.xaxis.set_major_locator(FixedLocator(ticks))
        if k == 1:
            ax.set_xlabel(xlabel, labelpad=1.5)
        ax.set_title(COHORT_LABEL[c], loc='left', fontsize=7)
        d = pg[pg.cohort == c]
        meds, below = {}, []
        for i, tool in enumerate(ORDER):
            tsv = TSVNAME.get(tool, tool)
            col = f'{metric}_{tsv}'
            if tool == 'RNA-Skim':
                r = rs[(rs.cohort_label.map(NAME2KEY) == c) & (rs.tool == 'RNA-Skim')].iloc[0]
                m = r['within_gene_spearman' if metric == 'rho' else 'noise_sd']
                mk(ax, m, i + 0.18, tool, ms=3.0)
                meds[tool] = m
                continue
            if col not in d or d[col].notna().sum() == 0:
                ax.text(lo + 0.02 * (hi - lo), i + 0.05, 'not run', va='center', ha='left', fontsize=7, color='#666666')
                continue
            v = d[col].to_numpy(dtype=float)
            meds[tool] = violin_row(ax, v, i, tool, lo, hi)
            vv = v[np.isfinite(v)]
            below.append(np.mean(vv < lo) if metric == 'rho' else np.mean(vv > hi))
            # the plotted median must equal the scorer's common-gene-set value (table2_accuracy_frameA.tsv)
            ref = fa[(fa.cohort.map(NAME2KEY) == c) & (fa.tool == tsv)]
            want = ref['within_gene_spearman' if metric == 'rho' else 'noise_sd'].iloc[0]
            assert round(float(np.median(vv)), 4) == round(float(want), 4), (c, tool, np.median(vv), want)
        if ref_rule and 'Salmon-pan' in meds:
            ax.axvline(meds['Salmon-pan'], color=STYLE['Salmon-pan']['color'], lw=0.5, ls=(0, (1, 1)), zorder=1)
        # median column to the right of the axes
        tx = ax_mm(fig, x0 + 38.5, top, 9, PH, W, H)
        tx.set_axis_off(); tx.set_ylim(NROW - 0.5, -0.5); tx.set_xlim(0, 1)
        tx.text(0.05, -0.75, 'median', fontsize=7, color='#444444', va='bottom')
        for i, tool in enumerate(ORDER):
            if tool in meds:
                tx.text(0.05, i + 0.18, fmt.format(meds[tool]), fontsize=7, va='center',
                        fontweight='bold' if tool == 'PLUCK' else 'normal')
        if below_note and below and max(below) > 0:
            pass
    return


# a: within-gene Spearman, per gene
dist_panel(12, 'rho', 0.90, 1.0, [0.90, 0.95, 1.0], 'Within-gene Spearman against the cohort truth (one value per gene)', '{:.4f}', 'a')
# b: within-gene noise SD, per gene
dist_panel(64, 'sd', 0.0, 0.30, [0, 0.1, 0.2, 0.3], 'Within-gene noise SD against the cohort truth (log2, one value per gene)', '{:.3f}', 'b')

# ---------------------------------------------------------------- panels c and d: spurious rate and median shift, all and unflagged genes
# data: fig3_metrics_all_frames.tsv (scoring set A and set P_rnaskim; subsets all and unflagged)
mt = read('fig3_metrics_all_frames.tsv')


def summary_row(c, tool, subset):
    tsv = TSVNAME.get(tool, tool)
    frame = 'P_rnaskim' if tool == 'RNA-Skim' else 'A'
    r = mt[(mt.cohort == c) & (mt.frame == frame) & (mt.subset == subset) & (mt.tool == tsv)]
    return None if r.empty else r.iloc[0]


TOPC = 114
letter(1, TOPC - 8, 'c')
letter(99, TOPC - 8, 'd')
for k, c in enumerate(COHORTS):
    for panel, x0 in (('c', LBL + k * 24.5), ('d', 104 + k * 24.5)):
        ax = ax_mm(fig, x0, TOPC, 21, 36, W, H)
        rows_axis(ax, labels=(panel == 'c' and k == 0))
        if panel == 'c' and k == 0:
            style_ylabels(ax)
        ax.set_title(COHORT_LABEL[c].replace(', ', '\n'), loc='left', fontsize=7, linespacing=1.0)
        if panel == 'c':
            ax.set_xscale('log'); ax.set_xlim(0.8, 125)
            ax.xaxis.set_major_locator(FixedLocator([1, 5, 20, 100]))
            ax.set_xticklabels(['1', '5', '20', '100'])
            ax.xaxis.set_minor_locator(FixedLocator([2, 3, 4, 6, 7, 8, 9, 20, 30, 40, 60, 70, 80, 90]))
            ax.xaxis.set_minor_formatter(NullFormatter())
            ax.axvline(5, color='#555555', lw=0.6, ls=(0, (3, 2)), zorder=1)
            if k == 1:
                ax.set_xlabel('Spurious strain-effect rate (%)', labelpad=1.5)
        else:
            ax.set_xlim(-0.30, 0.035)   # Munger unflagged Kmerator+countTags is -0.274; same limits in all three cohorts
            ax.xaxis.set_major_locator(FixedLocator([-0.2, -0.1, 0]))
            ax.set_xticklabels(['−0.2', '−0.1', '0'])
            ax.axvline(0, color='#555555', lw=0.5, zorder=1)
            if k == 1:
                ax.set_xlabel('Median founder shift Δ (log2)', labelpad=1.5)
        for i, tool in enumerate(ORDER):
            ra, ru = summary_row(c, tool, 'all'), summary_row(c, tool, 'unflagged')
            if ra is None:
                ax.text((1.2 if panel == 'c' else -0.29), i, 'not run', va='center', ha='left', fontsize=7, color='#666666')
                continue
            if panel == 'c':
                xa, xu = 100 * ra.spurious_rate, 100 * ru.spurious_rate
                ax.plot([xa, xu], [i - 0.17, i + 0.2], color=STYLE[tool]['color'], lw=0.5, zorder=2)
                mk(ax, xa, i - 0.17, tool, ms=3.0)
                mk(ax, xu, i + 0.2, tool, ms=3.0, scale=0.75, alpha=0.55)
            else:
                for rr, yy, sc, al in ((ra, i - 0.17, 1.0, 1.0), (ru, i + 0.2, 0.75, 0.55)):
                    lo_, hi_ = [float(v) for v in rr.delta_ci95.strip('[]').split(',')]
                    ax.plot([lo_, hi_], [yy, yy], color=STYLE[tool]['color'], lw=0.8, alpha=al, zorder=2)
                    mk(ax, rr.median_delta, yy, tool, ms=3.0, scale=sc, alpha=al)
# key for the two sub-rows
kx = ax_mm(fig, LBL, TOPC + 44.5, 120, 4, W, H); kx.set_axis_off(); kx.set_xlim(0, 120); kx.set_ylim(0, 1)
kx.text(0, 0.5, 'c, d: upper point, all core genes; lower point, unflagged genes (calibrated within the subset).\n'
        'Bars in d: bootstrap 95% CI over genes, mostly narrower than the markers.', va='center', fontsize=7, linespacing=1.15)

# ---------------------------------------------------------------- panel e: rank robustness across samples
# data: fig3e_persample_u3.tsv (per-sample Spearman across genes against the cohort truth, rank among the tools of the scoring set with all methods)
ps = read('fig3e_persample_u3.tsv')
ps['label'] = ps.tool.map(label)
TOPE = 176
ER = [t for t in ORDER if t != 'RNA-Skim']   # RNA-Skim is not in the scoring set with all methods
letter(1, TOPE - 11.3, 'e')
fig.text(LBL / W, 1 - (TOPE - 11.6) / H, 'Rank of each pipeline in each sample, by Spearman\ncorrelation across genes with the cohort truth',
         fontsize=7, va='top', linespacing=1.0)
for k, c in enumerate(COHORTS):
    d = ps[ps.cohort == c]
    n_s = d.srr.nunique()
    ax = ax_mm(fig, LBL + k * 21.5, TOPE, 19.5, 30, W, H)
    ax.set_ylim(len(ER) - 0.5, -0.5); ax.set_xlim(0.4, 8.6)
    ax.set_yticks(range(len(ER))); ax.set_yticklabels(ER if k == 0 else [])
    ax.tick_params(axis='y', length=0, pad=2)
    if k == 0:
        style_ylabels(ax)
    ax.set_xticks(range(1, 9)); ax.set_xticklabels([str(i) for i in range(1, 9)])
    ax.tick_params(axis='x', length=1.5, pad=1)
    for i in range(len(ER)):
        if i % 2 == 0:
            ax.axhspan(i - 0.5, i + 0.5, color='#F3F3F3', zorder=0, lw=0)
    ax.set_title(COHORT_LABEL[c].replace(', ', '\n'), loc='left', fontsize=7, linespacing=1.0)
    if k == 1:
        ax.set_xlabel('Rank (1 = highest)', labelpad=1.5)
    for i, tool in enumerate(ER):
        dt = d[d.label == tool]
        if dt.empty:
            ax.text(0.6, i, 'not run', va='center', ha='left', fontsize=7, color='#666666')
            continue
        for r_, n_ in dt['rank'].value_counts().items():
            ax.scatter([r_], [i], s=(n_ / n_s) * 22, marker='o', color=STYLE[tool]['color'], alpha=0.9, lw=0,
                       zorder=5 if tool == 'PLUCK' else 3)
# bubble size key (area proportional to the share of samples)
kx = ax_mm(fig, LBL, TOPE + 38.5, 62, 4, W, H); kx.set_axis_off(); kx.set_xlim(0, 62); kx.set_ylim(0, 1)
kx.text(0, 0.5, 'Share of samples at that rank:', va='center', fontsize=7)
for j, f in enumerate([0.1, 0.5, 1.0]):
    kx.scatter([37 + j * 9], [0.5], s=f * 22, color='#777777', lw=0)
    kx.text(38.3 + j * 9, 0.5, f'{int(f * 100)}%', va='center', fontsize=7)

# ---------------------------------------------------------------- panel f: GSE199702, three tissues, 101 bp
# data: fig3f_e4b_gates.tsv (values and gene-bootstrap 95% intervals)
g = read('fig3f_e4b_gates.tsv')
TIS = ['liver', 'heart', 'kidney']
TN = {'liver': 'Liver (53)', 'heart': 'Heart (30)', 'kidney': 'Kidney (30)'}
letter(93, TOPE - 11.3, 'f')
fig.text(97 / W, 1 - (TOPE - 10.6) / H, 'GSE199702 liver, heart and kidney, 101 bp', fontsize=7, va='top')
FX = 110
spec = [  # key in the data file, x0, y0, width, xlim, ticks, title, tool marker
    ('G2', FX, TOPE, 67, (1.0, 2.0), [1.0, 1.25, 1.5, 1.75, 2.0], 'Spurious-rate ratio, PLUCK / Salmon-pan', 'PLUCK'),
    ('G1', FX, TOPE + 21, 21, (0, 0.025), [0, 0.01, 0.02], 'Spearman gap,\nSalmon-pan − PLUCK', 'PLUCK'),
    ('G3', FX + 26, TOPE + 21, 14, (0, 0.035), [0, 0.03], '|Δ|, PLUCK', 'PLUCK'),
    ('G4', FX + 45, TOPE + 21, 19, (-0.13, -0.07), [-0.12, -0.08], 'Δ, ref-unique\nprobes', 'ref-unique probes'),
]
for gate, x0, y0, w, xl, tk, title, tool in spec:
    ax = ax_mm(fig, x0, y0, w, 10, W, H)
    ax.set_ylim(2.6, -0.6); ax.set_xlim(*xl)
    ax.set_yticks(range(3)); ax.set_yticklabels([TN[t] for t in TIS] if x0 == FX else [])
    ax.tick_params(axis='y', length=0, pad=2)
    ax.xaxis.set_major_locator(FixedLocator(tk))
    ax.set_xticklabels([('−' if v < 0 else '') + f'{abs(v):g}' for v in tk])
    ax.set_title(title, loc='left', fontsize=7)
    for i, t in enumerate(TIS):
        r = g[(g.tissue == t) & (g.gate == gate)].iloc[0]
        ax.plot([r.gene_bootstrap_lo, r.gene_bootstrap_hi], [i, i], color=STYLE[tool]['color'], lw=0.8, zorder=2)
        mk(ax, r.value, i, tool, ms=3.0)
        if gate == 'G2':
            ax.text(r.gene_bootstrap_hi + 0.03, i, f'{r.value:.3f}', va='center', ha='left', fontsize=7)
fig.text(FX / W, 1 - (TOPE + 35.5) / H, 'Bars: gene-bootstrap 95% interval.', fontsize=7, va='top')

save(fig, 'Fig3')
