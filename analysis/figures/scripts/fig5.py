#!/usr/bin/env python3
"""Fig. 5: a design-only reliability flag and PLUCK's error tail.
Reads only figures/data/ (each panel names its files below). Writes figures/Fig5.pdf and Fig5.svg."""
import numpy as np, pandas as pd
from scipy.stats import gaussian_kde
from u5_style import *

setup()
rd = lambda f: pd.read_csv(f'{DATA}/{f}', sep='\t', comment='#')
fig = plt.figure(figsize=(180 * MM, 217 * MM))
ORANGE = PIPE['PLUCK'][0]
SKY = PIPE['Salmon-ref'][0]
rng = np.random.default_rng(5)
from decimal import Decimal, ROUND_HALF_UP
half_up = lambda v: str(Decimal(repr(float(v))).quantize(Decimal('0.01'), rounding=ROUND_HALF_UP))   # 1.585 -> 1.59 (round half up)

# =============================================================== a: ROC and AUC
# data: fig5a_roc_curves_u5.tsv (curves), fig5a_operating_points_u5.tsv (u_g < 0.6 points),
#       fig5a_flag_vs_baseline_auc.tsv (AUC with gene-bootstrap 95% CI)
roc = rd('fig5a_roc_curves_u5.tsv'); ops = rd('fig5a_operating_points_u5.tsv'); auc = rd('fig5a_flag_vs_baseline_auc.tsv')
ax = ax_mm(fig, 17, 9, 42, 42)
ax.plot([0, 1], [0, 1], color='#BBBBBB', lw=0.6, ls=':', zorder=1)
for c in COHORT_ORDER:
    for model, col, lw, z in (('B1', GREY_B1, 0.7, 2), ('U', ORANGE, 0.9, 4)):
        r = roc[(roc.cohort == c) & (roc.model == model) & (roc.fpr >= 0)]
        ax.plot(r.fpr, r.tpr, color=col, lw=lw, ls='--' if c == 'keele_liver' else '-', zorder=z)
o6 = ops[ops.threshold == 0.6]
for _, r in o6.iterrows():   # one filled dot per cohort-tissue; Keele liver is identified by its dashed curve (open dots mean u_g < 0.7 in Fig. S8)
    ax.plot(r.fpr, r.tpr_recall, 'o', ms=4.2, mfc=ORANGE, mec='white', mew=0.8, zorder=6)
ax.set_xlim(0, 1); ax.set_ylim(0, 1.0); ax.set_aspect('equal')
ax.set_xlabel('False-positive rate'); ax.set_ylabel('True-positive rate (recall of tail genes)')
tk = [0, 0.25, 0.5, 0.75, 1]; tl_ = ['0', '0.25', '0.5', '0.75', '1']
ax.set_xticks(tk); ax.set_yticks(tk); ax.set_xticklabels(tl_); ax.set_yticklabels(tl_)
ax.text(0.98, 0.23, '● $u_g$ < 0.6', color=ORANGE, ha='right', va='bottom')
ax.text(0.98, 0.135, 'design flag $u_g$ (no fit)', color=ORANGE, ha='right', va='bottom')
ax.text(0.98, 0.04, 'gene-feature baseline', color=GREY_B1, ha='right', va='bottom')
letter(fig, 1, 5, 'a')

# AUC per cohort-tissue
ax2 = ax_mm(fig, 92, 9, 36, 42)
ys = {c: i for i, c in enumerate(COHORT_ORDER[::-1])}
for c in COHORT_ORDER:
    y = ys[c]
    for model, col, dy in (('B1', GREY_B1, -0.24), ('B2', GREY_B2, 0.0), ('U', ORANGE, 0.24)):
        r = auc[(auc.cohort == c) & (auc.model == model)].iloc[0]
        ax2.plot([r.auc_lo, r.auc_hi], [y + dy] * 2, color=col, lw=0.8, zorder=2)
        if model == 'U':
            ax2.plot(r.auc, y + dy, 'o', ms=4.4, mfc=ORANGE, mec='white', mew=0.6, zorder=4)
        else:
            ax2.plot(r.auc, y + dy, '|', ms=6.5, color=col, mew=1.5, zorder=3)
    if y > 0: ax2.axhline(y - 0.5, color='#EEEEEE', lw=0.5, zorder=0)
ax2.set_yticks([ys[c] for c in COHORT_ORDER]); ax2.set_yticklabels([COHORT_LABEL[c] for c in COHORT_ORDER])
ax2.set_ylim(-0.55, len(COHORT_ORDER) - 0.45); ax2.tick_params(axis='y', length=0)
ax2.set_xlim(0.5, 1.0); ax2.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0]); ax2.set_xticklabels(['0.5', '0.6', '0.7', '0.8', '0.9', '1'])
ax2.set_xlabel('AUC for the >2× tail (95% CI)')
mtext(fig, 92, 1.8, '● design flag $u_g$', color=ORANGE, va='center')
mtext(fig, 92, 4.5, '| baseline + multi-copy fraction', color=GREY_B2, va='center', fontweight='bold')
mtext(fig, 92, 7.2, '| gene-feature baseline', color=GREY_B1, va='center', fontweight='bold')

# figure-level key (top right)
lx, ly = 138, 7
mtext(fig, lx, ly, 'Pipelines', fontweight='bold', va='top')
for i, lab in enumerate(['PLUCK', 'Salmon-pan', 'Salmon-ref', 'STAR+fC', 'gene probes alone']):
    lax = ax_mm(fig, lx, ly + 3.0 + i * 3.5, 4, 3); lax.axis('off'); lax.set_xlim(0, 1); lax.set_ylim(0, 1)
    pmark(lax, 0.5, 0.5, lab, size=3.6)
    mtext(fig, lx + 5, ly + 4.5 + i * 3.5, lab, va='center')
gy = ly + 22.5
mtext(fig, lx, gy, 'Gene sets (b, c)', fontweight='bold', va='top')
for k, (fc, txt) in enumerate((('white', 'unflagged, $u_g$ ≥ 0.6 (left)'), (ORANGE, 'flagged, $u_g$ < 0.6 (right)'))):
    lax = ax_mm(fig, lx, gy + 3.0 + k * 3.5, 4, 3); lax.axis('off'); lax.set_xlim(0, 1); lax.set_ylim(0, 1)
    lax.add_patch(plt.Rectangle((0.15, 0.1), 0.7, 0.8, fc=fc, ec=ORANGE, lw=0.8, alpha=0.8 if fc != 'white' else 1))
    mtext(fig, lx + 5, gy + 4.5 + k * 3.5, txt, va='center')
cy = gy + 12
mtext(fig, lx, cy, 'Calibration (d)', fontweight='bold', va='top')
mtext(fig, lx, cy + 3.2, 'left, full size: inside the subset\n(as scored)\nright, smaller and lighter:\nall-gene constants', va='top', linespacing=1.15)
mtext(fig, lx, cy + 16.5, '† tuning cohort for the flag\nK, L, M in e: Keele kidney, Keele\nliver†, Munger founders', va='top', color='#333333', linespacing=1.15)

# =============================================================== b, c: per-gene ratio distributions, identical y axes
YL = (np.log2(0.25), np.log2(128)); YT = [0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128]
def half_violin(ax, x, vals, side, filled, width=0.40):
    v = np.clip(np.log2(np.maximum(vals, 1e-12)), YL[0], YL[1]); g = np.linspace(YL[0], YL[1], 400)
    k = gaussian_kde(v, bw_method=0.12)(g); k = k / k.max() * width
    sgn = -1 if side == 'L' else 1
    ax.fill_betweenx(g, x, x + sgn * k, color=ORANGE if filled else 'white', alpha=0.55 if filled else 1, lw=0, zorder=2)
    ax.plot(x + sgn * k, g, color=ORANGE, lw=0.6 if filled else 0.8, zorder=3)
    # every gene above the tail line as a dot, jittered across the half
    t = v[v > 1]
    ax.scatter(x + sgn * rng.uniform(0.03, width * 0.95, len(t)), t, s=0.6, color=ORANGE if filled else '#B04A00', lw=0, alpha=0.8, zorder=3.5, rasterized=False)
    m = np.median(v)
    ax.plot([x, x + sgn * width], [m, m], color='black', lw=0.9, zorder=4)
TISSUE6 = ['kidney\n50 bp', 'liver\n50 bp†', 'founders\n100 bp', 'liver\n101 bp', 'heart\n101 bp', 'kidney\n101 bp']
def group_brackets(ax, groups, y=-0.155):
    tr = ax.get_xaxis_transform()
    for a, b, lab in groups:
        ax.plot([a - 0.38, b + 0.38], [y, y], color='black', lw=0.6, transform=tr, clip_on=False)
        ax.text((a + b) / 2, y - 0.02, lab, ha='center', va='top', transform=tr)
def ratio_axis(ax, label):
    ax.set_ylim(*YL); ax.set_yticks(np.log2(YT)); ax.set_yticklabels([('%g' % t) for t in YT])
    ax.axhline(0, color='#888888', lw=0.6, zorder=1); ax.axhline(1, color='black', lw=0.7, ls='--', zorder=1)
    ax.set_ylabel(label, linespacing=1.1)

# b: data fig5b_pergene_resvar_ratio_u5.tsv (per-gene ratios), fig5a_operating_points_u5.tsv (tail rates and recall at u_g < 0.6)
pg = rd('fig5b_pergene_resvar_ratio_u5.tsv')
BX, BW = 33, 76
axb = ax_mm(fig, BX, 74, BW, 50)
for i, c in enumerate(COHORT_ORDER):
    d = pg[pg.cohort == c]
    half_violin(axb, i, d[d.flag == 0].ratio.values, 'L', False)
    half_violin(axb, i, d[d.flag == 1].ratio.values, 'R', True)
    if i: axb.axvline(i - 0.5, color='#EEEEEE', lw=0.5, zorder=0)
ratio_axis(axb, 'Residual variance per gene,\nPLUCK / Salmon-pan')
axb.set_xlim(-0.5, 5.5); axb.set_xticks(range(6)); axb.set_xticklabels(TISSUE6, linespacing=1.05)
axb.tick_params(axis='x', length=0)
group_brackets(axb, [(0, 1, 'Keele'), (2, 2, 'Munger'), (3, 5, 'GSE199702')])
axb.text(5.6, 1, 'tail:\n>2×', ha='left', va='center', linespacing=1.05)
axb.set_title('Real cohorts (ratio to Salmon-pan; no cohort truth used)', loc='left', pad=3)
o = o6.set_index('cohort')
rows = [('tail, all genes', 'tail_all'), ('tail, unflagged', 'tail_unflagged'), ('tail, flagged', 'tail_flagged'), ('genes flagged', 'flag_rate'), ('recall of tail', 'tpr_recall')]
for j, (lab, col) in enumerate(rows):
    yy = 141 + j * 3.1
    mtext(fig, BX - 1.5, yy, lab, ha='right', va='center', color='#333333')
    for i, c in enumerate(COHORT_ORDER):
        v = o.loc[c, col] * 100
        s = ('%.2f%%' % v) if v < 1 else ('%.1f%%' % v if (v < 10 or col == 'flag_rate') else '%.0f%%' % v)   # rounding as in the text
        mtext(fig, BX + BW * (i + 0.5) / 6, yy, s, ha='center', va='center', color=ORANGE if col in ('tail_flagged', 'tpr_recall') else 'black')
letter(fig, 1, 66, 'b')

# c: data fig5c_e2_pergene_errvar_u5.tsv (per-gene variances), fig5c_e2_pooled_ratio_ci_u5.tsv (pooled ratios, gene-bootstrap 95% CI)
e2 = rd('fig5c_e2_pergene_errvar_u5.tsv'); pool = rd('fig5c_e2_pooled_ratio_ci_u5.tsv')
LAY = ['50 bp SE', '101 bp SE', '101 bp PE']
CX, CW = 140, 38
axc = ax_mm(fig, CX, 74, CW, 50)
for i, L in enumerate(LAY):
    d = e2[e2.layout == L]
    for flag, side in ((0, 'L'), (1, 'R')):
        dd = d[d.flag == flag]
        half_violin(axc, i, (dd.var_hybrid / dd.var_salmon_pan).values, side, flag == 1)
        p = pool[(pool.layout == L) & (pool.subset == ('flagged' if flag else 'unflagged'))].iloc[0]
        sg = 1 if flag else -1
        for lab, col, lo, hi, val, off in (('PLUCK', ORANGE, p.ratio_hybrid_lo, p.ratio_hybrid_hi, p.ratio_hybrid, 0.24),
                                            ('Salmon-ref', SKY, p.ratio_salmon_ref_lo, p.ratio_salmon_ref_hi, p.ratio_salmon_ref, 0.36)):
            xp = i + sg * off
            axc.plot([xp, xp], np.log2([lo, hi]), color='white', lw=2.4, solid_capstyle='butt', zorder=5.9)   # halo so the interval shows over the violin fill
            axc.plot([xp, xp], np.log2([lo, hi]), color=col if lab != 'PLUCK' else '#8A3A00', lw=1.1, zorder=6)
            pmark(axc, xp, np.log2(val), lab, size=3.2, z=6)
    if i: axc.axvline(i - 0.5, color='#EEEEEE', lw=0.5, zorder=0)
ratio_axis(axc, 'Error variance per gene,\nPLUCK / Salmon-pan')
axc.set_xlim(-0.5, 2.5); axc.set_xticks(range(3)); axc.set_xticklabels(['Simulation\n50 bp SE', 'Simulation\n101 bp SE', 'Simulation\n101 bp PE'], linespacing=1.05)
axc.tick_params(axis='x', length=0)
axc.set_title('Simulation (expected-count truth)', loc='left', pad=3)
pl = pool.set_index(['layout', 'subset'])
mtext(fig, CX - 1.5, 141, 'pooled variance ratio:', ha='right', va='center', color='#333333')
rows = [('PLUCK, unflagged', 'unflagged', 'ratio_hybrid', 'black'), ('PLUCK, flagged', 'flagged', 'ratio_hybrid', ORANGE), ('Salmon-ref, flagged', 'flagged', 'ratio_salmon_ref', SKY)]
for j, (lab, sub, col, colr) in enumerate(rows):
    yy = 144.1 + j * 3.1
    mtext(fig, CX - 1.5, yy, lab, ha='right', va='center', color='#333333')
    for i, L in enumerate(LAY):
        mtext(fig, CX + CW * (i + 0.5) / 3, yy, half_up(pl.loc[(L, sub), col]), ha='center', va='center', color=colr)
letter(fig, 112, 66, 'c')

# =============================================================== d: flagged-gene shift and spurious rate, both calibrations
# data: fig5d_flagged_shift_both_calibrations.tsv (flagged rows, "subset (as scored)" and "all-gene"),
#       fig3_metrics_all_frames.tsv (scoring set A, flagged: delta_ci95 of the as-scored values)
fd = rd('fig5d_flagged_shift_both_calibrations.tsv'); fd = fd[fd.subset == 'flagged']
fr = rd('fig3_metrics_all_frames.tsv'); fr = fr[(fr.frame == 'A') & (fr.subset == 'flagged')]
METH = ['hybrid', 'salmon_pan', 'salmon_ref', 'star_fc']; FR_CODE = {'hybrid': 'vf'}
COH3 = ['keele_kidney', 'keele_liver', 'munger']
def dpanel(ax, col, pct):
    for i, c in enumerate(COH3):
        for j, m in enumerate(METH):
            lab = CODE2LABEL[m]; x0 = i + (j - 1.5) * 0.215
            a = fd[(fd.cohort == c) & (fd.method == m) & (fd.calibration == 'subset (as scored)')][col].iloc[0]
            b = fd[(fd.cohort == c) & (fd.method == m) & (fd.calibration == 'all-gene')][col].iloc[0]
            ya, yb = (a * 100, b * 100) if pct else (a, b)
            ax.plot([x0 - 0.05, x0 + 0.05], [ya, yb], color=PIPE[lab][0], lw=0.8, zorder=2)
            if not pct:
                lo, hi = [float(t) for t in fr[(fr.cohort == c) & (fr.method == FR_CODE.get(m, m))].delta_ci95.iloc[0].strip('[]').split(',')]
                ax.plot([x0 - 0.05] * 2, [lo, hi], color=PIPE[lab][0], lw=0.6, zorder=2)
            pmark(ax, x0 - 0.05, ya, lab, size=3.6)                          # as scored: full-size marker
            pmark(ax, x0 + 0.05, yb, lab, size=3.6, scale=0.75, alpha=0.55)  # all-gene: same marker, smaller and lighter (Fig. 3c/d convention)
        if i: ax.axvline(i - 0.5, color='#EEEEEE', lw=0.5, zorder=0)
    ax.set_xlim(-0.5, 2.5); ax.set_xticks(range(3)); ax.set_xticklabels([COHORT_SHORT[c] for c in COH3], linespacing=1.05)
    ax.tick_params(axis='x', length=0)
axd1 = ax_mm(fig, 22, 166, 48, 38)
dpanel(axd1, 'median_delta', False)
axd1.axhline(0, color='black', lw=0.6, zorder=1)
axd1.set_ylim(-0.08, 0.08); axd1.set_yticks([-0.08, -0.04, 0, 0.04, 0.08]); axd1.set_yticklabels(['−0.08', '−0.04', '0', '0.04', '0.08'])
axd1.set_ylabel('Median founder shift Δ (log$_2$),\nflagged genes', linespacing=1.1)
axd1.text(-0.42, 0.074, 'Munger, as scored → all-gene:\nPLUCK +0.058 → +0.016\nSTAR+fC +0.059 → −0.008\nSalmon-pan +0.001 → +0.006',
          va='top', ha='left', linespacing=1.15)
axd1.set_title('Flagged genes, both calibrations', loc='left', pad=3)
axd2 = ax_mm(fig, 89, 166, 37, 38)
dpanel(axd2, 'spurious_rate', True)
axd2.set_yscale('log'); axd2.set_ylim(1, 100); axd2.set_yticks([1, 3, 10, 30, 100]); axd2.set_yticklabels(['1%', '3%', '10%', '30%', '100%'])   # same log scale as Fig. 3c
axd2.minorticks_off(); axd2.axhline(5, color='black', lw=0.6, ls='--', zorder=1)
axd2.set_ylabel('Spurious strain-effect rate,\nflagged genes (log scale)', linespacing=1.1)
letter(fig, 1, 158, 'd')

# =============================================================== e: ablation, gene probes alone vs PLUCK on the same genes
# data: fig5e_ablation_arms.tsv (noise SD: arm1_vf26 vs arm3_hybrid_frozen__on_arm1_vf26_genes),
#       fig5e_ablation_tail_cis_scan.tsv (share of genes with >2x Salmon-pan residual variance, arm1 vs hybrid)
ab = rd('fig5e_ablation_arms.tsv'); tl = rd('fig5e_ablation_tail_cis_scan.tsv')
axe1 = ax_mm(fig, 141, 166, 22, 38)
xs = [0, 1, 2, 3.6, 4.6, 5.6]
for x, (s, c) in zip(xs, [(s, c) for s in ('all', 'flagged') for c in COH3]):
    a1 = ab[(ab.cohort == c) & (ab.arm == 'arm1_vf26') & (ab.subset == s)].noise_sd.iloc[0]
    a3 = ab[(ab.cohort == c) & (ab.arm == 'arm3_hybrid_frozen__on_arm1_vf26_genes') & (ab.subset == s)].noise_sd.iloc[0]
    axe1.plot([x, x], [a1, a3], color='#AAAAAA', lw=0.7, zorder=2)
    pmark(axe1, x, a1, 'gene probes alone', size=3.4); pmark(axe1, x, a3, 'PLUCK', size=3.2)
axe1.set_xticks(xs); axe1.set_xticklabels(['K', 'L†', 'M'] * 2); axe1.tick_params(axis='x', length=0)
axe1.set_xlim(-0.6, 6.2); axe1.set_ylim(0, 0.21); axe1.set_yticks([0, 0.05, 0.1, 0.15, 0.2]); axe1.set_yticklabels(['0', '0.05', '0.10', '0.15', '0.20'])
axe1.set_ylabel('Noise SD against cohort truth', labelpad=2)
axe1.annotate('', xy=(0, -0.11), xytext=(2, -0.11), xycoords=('data', 'axes fraction'), textcoords=('data', 'axes fraction'), arrowprops=dict(arrowstyle='-', lw=0.6))
axe1.text(1, -0.13, 'all', ha='center', va='top', transform=axe1.get_xaxis_transform())
axe1.annotate('', xy=(3.6, -0.11), xytext=(5.6, -0.11), xycoords=('data', 'axes fraction'), textcoords=('data', 'axes fraction'), arrowprops=dict(arrowstyle='-', lw=0.6))
axe1.text(4.6, -0.13, 'flagged', ha='center', va='top', transform=axe1.get_xaxis_transform())
axe1.set_title('EM branch ablation', loc='left', pad=3)
axe2 = ax_mm(fig, 170, 166, 8.5, 38)
for k, c in enumerate(['keele_kidney', 'keele_liver']):
    r = tl[(tl.tissue == c) & (tl.arm == 'arm1')].iloc[0]
    axe2.plot([k, k], [r.gt2pct_arm1, r.gt2pct_hybrid], color='#AAAAAA', lw=0.7, zorder=2)
    pmark(axe2, k, r.gt2pct_arm1, 'gene probes alone', size=3.4); pmark(axe2, k, r.gt2pct_hybrid, 'PLUCK', size=3.2)
axe2.set_xlim(-0.6, 1.6); axe2.set_xticks([0, 1]); axe2.set_xticklabels(['K', 'L†']); axe2.tick_params(axis='x', length=0)
axe2.set_ylim(0, 22); axe2.set_yticks([0, 5, 10, 15, 20]); axe2.set_yticklabels(['0', '5', '10', '15', '20'])
axe2.set_ylabel('% genes in >2× tail', labelpad=1)
letter(fig, 130, 158, 'e')

save(fig, 'Fig5')
print('wrote Fig5')
