#!/usr/bin/env python3
"""Supplementary Figs S4-S7. Reads only figures/data/ (file named above each panel).
Writes figures/FigS4..FigS7 (.pdf, .svg).  Run: python3 -I figures/scripts/figS4_S7.py"""
import sys
# Folder constants below point to the original working folder; set them for your copy.
sys.path.insert(0, '/mnt/project-files/paper/figures/scripts')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, NullFormatter
from matplotlib.colors import LinearSegmentedColormap
from u3_encoding import MM, ORDER, STYLE, COHORTS, COHORT_LABEL, setup, read, label, mk, ax_mm, save

setup()
MINUS = '−'


def letter(fig, W, H, x, y, s):
    fig.text(x / W, 1 - y / H, s, fontsize=10, fontweight='bold', va='top', ha='left')


def fmt_tick(v):
    return (MINUS if v < 0 else '') + f'{abs(v):g}'


# =============================================================== Fig. S4: pairwise agreement matrices
# data: figS4_pairwise_u3.tsv (per-sample Spearman across genes between two estimates, or an estimate and the truth)
def fig_s4():
    pw = read('figS4_pairwise_u3.tsv')
    pw['a'] = pw.tool_a.map(label); pw['b'] = pw.tool_b.map(label)
    names = ['cohort truth'] + [t for t in ORDER if t != 'RNA-Skim']
    n = len(names)
    W, H = 180, 89
    fig = plt.figure(figsize=(W * MM, H * MM))
    cm = LinearSegmentedColormap.from_list('agree', ['#FFFFFF', '#C6DBEF', '#6BAED6', '#2171B5', '#08306B'])
    vmin, vmax = 0.6, 1.0
    cell = 5.3
    for k, c in enumerate(COHORTS):
        x0 = 26 + k * (n * cell + 4.5)
        ax = ax_mm(fig, x0, 8, n * cell, n * cell, W, H)
        d = pw[pw.cohort == c]
        ax.set_xlim(-0.5, n - 0.5); ax.set_ylim(n - 0.5, -0.5)
        for i in range(n):
            for j in range(n):
                if i == j:
                    ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, color='#E6E6E6', lw=0))
                    continue
                a, b = names[min(i, j)], names[max(i, j)]
                v = d[((d.a == a) & (d.b == b)) | ((d.a == b) & (d.b == a))].spearman
                if v.empty:
                    ax.text(j, i, '–', ha='center', va='center', fontsize=7, color='#999999')
                    continue
                val = v.median() if i > j else v.min()
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, color=cm((val - vmin) / (vmax - vmin)), lw=0))
                ax.text(j, i, f'{val:.3f}'.lstrip('0'), ha='center', va='center', fontsize=5.5 if False else 7,
                        color='white' if val >= 0.83 else 'black', fontstretch='condensed')
        ax.set_xticks(range(n)); ax.set_yticks(range(n))
        ax.set_xticklabels(names, rotation=90, fontsize=7)
        ax.set_yticklabels(names if k == 0 else [], fontsize=7)
        ax.tick_params(length=0, pad=1.5)
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_title(COHORT_LABEL[c], loc='left', fontsize=7)
        for t in ax.get_yticklabels() + ax.get_xticklabels():
            if t.get_text() == 'PLUCK':
                t.set_fontweight('bold')
        letter(fig, W, H, x0 - (25 if k == 0 else 3.5), 2, 'abc'[k])
    fig.text(26 / W, 1 - 82 / H, 'Below the diagonal: median over samples. Above the diagonal: lowest sample. –, not run.', fontsize=7, va='top')
    # colour scale: linear from 0.6 (white) to 1.0 (dark blue)
    cb = ax_mm(fig, 142, 81.5, 32, 2.2, W, H)
    cb.imshow(np.linspace(vmin, vmax, 256)[None, :], aspect='auto', cmap=cm, vmin=vmin, vmax=vmax, extent=(vmin, vmax, 0, 1))
    cb.set_yticks([]); cb.set_xticks([0.6, 0.7, 0.8, 0.9, 1.0]); cb.set_xticklabels(['0.6', '0.7', '0.8', '0.9', '1'])
    cb.tick_params(axis='x', length=1.5, pad=1)
    for sp in cb.spines.values():
        sp.set_visible(False)
    fig.text(140 / W, 1 - 82.6 / H, 'Spearman', fontsize=7, va='center', ha='right')
    save(fig, 'FigS4')


# =============================================================== Fig. S5: GSE199702 per-tissue metrics
# data: figS5_e4b_tissue_metrics_u3.tsv (scorer compare block per tissue; GBRS/EMASE truth)
# data: fig3f_e4b_gates.tsv (G2 values and gene-bootstrap interval, for panel d)
def fig_s5():
    m = read('figS5_e4b_tissue_metrics_u3.tsv'); m['label'] = m.tool.map(label)
    g = read('fig3f_e4b_gates.tsv')
    tools = ['PLUCK', 'Salmon-pan', 'Salmon-ref', 'ref-unique probes']
    TIS = ['liver', 'heart', 'kidney']
    TN = {'liver': 'Liver, 101 bp (53)', 'heart': 'Heart, 101 bp (30)', 'kidney': 'Kidney, 101 bp (30)'}
    W, H = 180, 62
    fig = plt.figure(figsize=(W * MM, H * MM))
    panels = [('a', 'within_gene_spearman', 'Within-gene Spearman\n(GBRS truth)', (0.93, 1.0), [0.94, 0.96, 0.98, 1.0], False),
              ('b', 'spurious_rate', 'Spurious strain-effect\nrate (%)', (1, 100), [1, 5, 20, 100], True),
              ('c', 'median_delta', 'Median founder\nshift Δ (log2)', (-0.13, 0.02), [-0.1, -0.05, 0], False)]
    for pi, (L, col, xl, lim, tk, logx) in enumerate(panels):
        x0 = 30 + pi * 37
        ax = ax_mm(fig, x0, 9, 31, 36, W, H)
        rows = [(t, tl) for t in TIS for tl in tools]
        ax.set_ylim(len(rows) - 0.5, -0.5)
        ax.set_yticks([i * len(tools) + 1.5 for i in range(3)])
        ax.set_yticklabels([TN[t] for t in TIS] if pi == 0 else [])
        ax.tick_params(axis='y', length=0)
        for i in range(3):
            if i % 2 == 0:
                ax.axhspan(i * 4 - 0.5, i * 4 + 3.5, color='#F3F3F3', lw=0, zorder=0)
        if logx:
            ax.set_xscale('log'); ax.xaxis.set_minor_formatter(NullFormatter())
            ax.axvline(5, color='#555555', lw=0.6, ls=(0, (3, 2)))
        if col == 'median_delta':
            ax.axvline(0, color='#555555', lw=0.5)
        ax.set_xlim(*lim); ax.xaxis.set_major_locator(FixedLocator(tk)); ax.set_xticklabels([fmt_tick(v) for v in tk])
        ax.set_xlabel(xl, labelpad=1.5)
        for r, (t, tl) in enumerate(rows):
            x = m[(m.tissue == t) & (m.label == tl)].iloc[0]
            v = x[col] * (100 if col == 'spurious_rate' else 1)
            if col == 'median_delta':
                lo, hi = [float(z) for z in x.delta_ci95.strip('[]').split(',')]
                ax.plot([lo, hi], [r, r], color=STYLE[tl]['color'], lw=0.8)
            mk(ax, v, r, tl, ms=3.0)
        letter(fig, W, H, x0 - (29 if pi == 0 else 4), 2, L)
    # legend of pipelines
    lx = ax_mm(fig, 30, 55, 105, 4, W, H); lx.set_axis_off(); lx.set_xlim(0, 105); lx.set_ylim(0, 1)
    xs = 0
    for tl in tools:
        mk(lx, xs + 1, 0.5, tl, ms=3.0); lx.text(xs + 2.6, 0.5, tl, va='center', fontsize=7); xs += 26
    # d: spurious-rate ratio against both Salmon indexes
    ax = ax_mm(fig, 150, 9, 27, 36, W, H)
    ax.set_ylim(2.6, -0.6); ax.set_xlim(0.9, 1.8)
    ax.set_yticks(range(3)); ax.set_yticklabels(['Liver', 'Heart', 'Kidney']); ax.tick_params(axis='y', length=0)
    ax.axvline(1.0, color='#555555', lw=0.5)
    for i, t in enumerate(TIS):
        r = g[(g.tissue == t) & (g.gate == 'G2')].iloc[0]
        ax.plot([r.gene_bootstrap_lo, r.gene_bootstrap_hi], [i - 0.15, i - 0.15], color=STYLE['Salmon-pan']['color'], lw=0.8)
        mk(ax, r.value, i - 0.15, 'Salmon-pan', ms=3.0)
        h = m[(m.tissue == t) & (m.label == 'PLUCK')].spurious_rate.iloc[0]
        sr = m[(m.tissue == t) & (m.label == 'Salmon-ref')].spurious_rate.iloc[0]
        mk(ax, h / sr, i + 0.2, 'Salmon-ref', ms=3.0)
    ax.xaxis.set_major_locator(FixedLocator([1.0, 1.25, 1.5]))
    ax.set_xlabel('PLUCK spurious rate\n/ comparator rate', labelpad=1.5)
    ax.set_ylim(2.6, -0.6)
    letter(fig, W, H, 141, 2, 'd')
    save(fig, 'FigS5')


# =============================================================== Fig. S6: 13 low-mapping libraries
# data: figS6_hits_identity_u3.tsv (hits per read and identity depth per library, ratio to cohort median)
# data: figS6_sensitivity_u3.tsv (scoring set A metrics with all libraries and with the 13 excluded)
def fig_s6():
    hi = read('figS6_hits_identity_u3.tsv')
    se = read('figS6_sensitivity_u3.tsv'); se['label'] = se.tool.map(label)
    W, H = 180, 112
    fig = plt.figure(figsize=(W * MM, H * MM))
    CO = {'keele_kidney': 'Keele kidney\n50 bp', 'keele_liver': 'Keele liver\n50 bp†', 'munger': 'Munger\n100 bp'}
    for pi, (chk, title) in enumerate((('hits per read', 'Probe hits per read / cohort median'),
                                      ('identity depth', 'Identity-check depth / cohort median'))):
        x0 = 22 + pi * 86
        ax = ax_mm(fig, x0, 8, 66, 36, W, H)
        d = hi[hi.check == chk]
        ax.set_yscale('log'); ax.set_ylim(0.12, 5); ax.yaxis.set_major_locator(FixedLocator([0.25, 0.5, 1, 2, 4]))
        ax.set_yticklabels(['0.25', '0.5', '1', '2', '4']); ax.yaxis.set_minor_formatter(NullFormatter())
        ax.tick_params(axis='y', length=2)
        for yv in ((0.5, 2.0) if chk == 'hits per read' else (0.5,)):
            ax.axhline(yv, color='#555555', lw=0.6, ls=(0, (3, 2)))
        ax.set_xlim(-0.5, 2.5); ax.set_xticks(range(3)); ax.set_xticklabels([CO[c] for c in COHORTS])
        ax.tick_params(axis='x', length=0)
        rng = np.random.default_rng(1)
        for ci, c in enumerate(COHORTS):
            dc = d[d.cohort == c]
            if dc.empty:
                ax.text(ci, 1, 'no identity\ntable', ha='center', va='center', fontsize=7, color='#666666')
                continue
            jit = rng.uniform(-0.28, 0.28, len(dc))
            low = dc.in_low_mapping_13.to_numpy(bool)
            ax.scatter(ci + jit[~low], dc.ratio[~low], s=5, color='#9A9A9A', lw=0)
            ax.scatter(ci + jit[low], dc.ratio[low], s=9, color='black', marker='D', lw=0)
            fl = dc.flag.to_numpy(bool)
            ax.scatter(ci + jit[fl], dc.ratio[fl], s=24, facecolor='none', edgecolor='#D55E00', lw=0.6)
        ax.set_ylabel(title, labelpad=1.5)
        letter(fig, W, H, x0 - 20, 2, 'ab'[pi])
    kx = ax_mm(fig, 22, 50, 150, 4, W, H); kx.set_axis_off(); kx.set_xlim(0, 150); kx.set_ylim(0, 1)
    kx.scatter([1], [0.5], s=5, color='#9A9A9A', lw=0); kx.text(2.5, 0.5, 'library', va='center', fontsize=7)
    kx.scatter([16], [0.5], s=9, color='black', marker='D', lw=0); kx.text(17.5, 0.5, 'one of the 13 low-mapping libraries', va='center', fontsize=7)
    kx.scatter([64], [0.5], s=24, facecolor='none', edgecolor='#D55E00', lw=0.6); kx.text(65.8, 0.5, 'flagged by this check', va='center', fontsize=7)
    kx.text(96, 0.5, 'dashed: check bounds (0.5–2×; depth ≥ 0.5×)', va='center', fontsize=7)
    # c-f: primary against 13-excluded, Keele kidney and liver, scoring set A, all three gene subsets
    mets = [('within_gene_spearman', 'Within-gene Spearman', (0.92, 1.0)), ('noise_sd', 'Noise SD', (0.0, 0.2)),
            ('spurious_rate', 'Spurious rate', (0.0, 0.55)), ('median_delta', 'Median Δ', (-0.13, 0.02))]
    for mi, (col, tt, lim) in enumerate(mets):
        x0 = 14 + mi * 42
        ax = ax_mm(fig, x0, 64, 32, 32, W, H)
        ax.plot(lim, lim, color='#999999', lw=0.5, zorder=1)
        for c in ('keele_kidney', 'keele_liver'):
            a = se[(se.cohort == c) & (se.run.str.startswith('all'))]
            b = se[(se.cohort == c) & (se.run.str.startswith('13'))]
            j = a.merge(b, on=['subset', 'method', 'label'], suffixes=('_a', '_b'))
            for _, r in j.iterrows():
                mk(ax, r[col + '_a'], r[col + '_b'], r['label'], ms=2.6)
        ax.set_xlim(*lim); ax.set_ylim(*lim)
        ax.set_title(tt, loc='left', fontsize=7)
        ax.set_xlabel('All 47 libraries', labelpad=1.5)
        if mi == 0:
            ax.set_ylabel('Low-mapping libraries excluded', labelpad=1.5)
        ax.tick_params(axis='y', length=2)
        letter(fig, W, H, x0 - 13 if mi == 0 else x0 - 8, 58, 'cdef'[mi])
    lx = ax_mm(fig, 14, 106, 166, 4, W, H); lx.set_axis_off(); lx.set_xlim(0, 166); lx.set_ylim(0, 1)
    xs = 0
    for tl in [t for t in ORDER if t not in ('RNA-Skim', 'kallisto-pan')]:
        mk(lx, xs + 1, 0.5, tl, ms=2.6); lx.text(xs + 2.4, 0.5, tl, va='center', fontsize=7)
        xs += 6 + len(tl) * 1.25
    save(fig, 'FigS6')


# =============================================================== Fig. S7: per-sample switch
# data: figS7_switch_u3.tsv (set A, all genes: PLUCK under --switch cohort (primary), --switch sample (A) and --switch sample-scaled (B))
# data: figS7_switch_flagged_u3.tsv (flagged-gene median shift, both calibrations)
# data: figS7_switch_claims_u3.tsv (kallisto-pan Munger R1/R2 under each switch)
def fig_s7():
    h = read('figS7_switch_u3.tsv')
    f = read('figS7_switch_flagged_u3.tsv')
    cl = read('figS7_switch_claims_u3.tsv')
    V = ['primary', 'A', 'B']
    VL = {'primary': '', 'A': '', 'B': ''}   # variants are told apart by the marker key under the panels
    W, H = 180, 58
    fig = plt.figure(figsize=(W * MM, H * MM))
    CO = {'keele_kidney': 'Kidney\n50 bp', 'keele_liver': 'Liver\n50 bp†', 'munger': 'Munger\n100 bp'}

    def grouped(ax, getter, ylim, yt, ylab, ref=None, limit=None):
        xs, labs = [], []
        for ci, c in enumerate(COHORTS):
            for vi, v in enumerate(V):
                x = ci * 4 + vi
                y = getter(c, v)
                if y is None:
                    continue
                kw = dict(scale=1.0) if v == 'primary' else dict(scale=0.85, alpha=0.75 if v == 'A' else 0.45)
                mk(ax, x, y, 'PLUCK', ms=3.0, **kw)
                xs.append(x); labs.append(VL[v])
            ax.text(ci * 4 + 1, ylim[0] - (ylim[1] - ylim[0]) * 0.05, CO[c], ha='center', va='top', fontsize=7, linespacing=1.0)
        ax.set_xticks(xs); ax.set_xticklabels(labs, fontsize=7); ax.tick_params(axis='x', length=0, pad=1)
        ax.set_xlim(-0.8, 10.8); ax.set_ylim(*ylim); ax.yaxis.set_major_locator(FixedLocator(yt))
        ax.set_yticklabels([fmt_tick(v) for v in yt]); ax.tick_params(axis='y', length=2)
        ax.set_ylabel(ylab, labelpad=1.5)
        if ref is not None:
            ax.axhline(ref, color='#555555', lw=0.5)
        if limit is not None:
            ax.axhline(limit, color='#555555', lw=0.6, ls=(0, (3, 2)))

    def hv(col):
        return lambda c, v: float(h[(h.cohort == c) & (h.variant == v)][col].iloc[0])

    ax = ax_mm(fig, 14, 6, 34, 30, W, H)
    grouped(ax, hv('gap_to_salmon_pan'), (0, 0.022), [0, 0.01, 0.02], 'Spearman gap to Salmon-pan')
    letter(fig, W, H, 1, 2, 'a')
    ax = ax_mm(fig, 60, 6, 34, 30, W, H)
    grouped(ax, hv('noise_ratio_vs_salmon_pan'), (1.0, 1.8), [1.0, 1.2, 1.4, 1.6, 1.8], 'Noise SD ratio to Salmon-pan')
    letter(fig, W, H, 50, 2, 'b')

    def flag_get(cal):
        def g_(c, v):
            r = f[(f.cohort == c) & (f.variant == v) & (f.subset == 'flagged') & (f.calibration == cal) & (f.method == 'vf')]
            return None if r.empty else float(r.median_delta.iloc[0])
        return g_
    ax = ax_mm(fig, 106, 6, 30, 30, W, H)
    grouped(ax, flag_get('own'), (-0.01, 0.07), [0, 0.02, 0.04, 0.06], 'Flagged genes: median Δ', ref=0)
    ax.set_title('calibrated within subset', loc='left', fontsize=7)
    letter(fig, W, H, 96, 2, 'c')
    ax = ax_mm(fig, 145, 6, 30, 30, W, H)
    grouped(ax, flag_get('global'), (-0.01, 0.07), [0, 0.02, 0.04, 0.06], '', ref=0)
    ax.set_title('all-gene calibration', loc='left', fontsize=7)
    letter(fig, W, H, 139, 2, 'd')
    # marker key for the three switch variants
    kx = ax_mm(fig, 14, 45.5, 160, 4, W, H); kx.set_axis_off(); kx.set_xlim(0, 160); kx.set_ylim(0, 1)
    for xk, (kw, txt) in zip((1, 28, 66), ((dict(scale=1.0), 'cohort switch'), (dict(scale=0.85, alpha=0.75), 'per-sample switch'),
                                          (dict(scale=0.85, alpha=0.45), 'scaled per-sample switch'))):
        mk(kx, xk, 0.5, 'PLUCK', ms=3.0, **kw); kx.text(xk + 1.8, 0.5, txt, va='center', fontsize=7)
    # kallisto-pan note under the panels
    r2 = {v: cl[(cl.variant == v) & cl.claim.str.contains('ratio_hybrid_over_M')] for v in ('A', 'B')}
    txt = (f"Munger, PLUCK / kallisto-pan spurious-rate ratio: cohort switch {r2['A'].primary.iloc[0]}, "
           f"per-sample switch {r2['A'].per_sample.iloc[0]}, scaled per-sample switch {r2['B'].per_sample.iloc[0]}.")
    fig.text(14 / W, 1 - 52 / H, txt, fontsize=7, va='top')
    save(fig, 'FigS7')


if __name__ == '__main__':
    fig_s4(); fig_s5(); fig_s6(); fig_s7()
