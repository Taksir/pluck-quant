#!/usr/bin/env python3
"""Fig. S12: the flag-guided Salmon second pass pilot and the EM-ablation with a broader k-mer set.
Reads only figures/data/. Writes figures/FigS12.pdf and FigS12.svg."""
import numpy as np, pandas as pd
from u5_style import *

setup()
rd = lambda f: pd.read_csv(f'{DATA}/{f}', sep='\t', comment='#')
fig = plt.figure(figsize=(180 * MM, 120 * MM))
ORANGE = PIPE['PLUCK'][0]
A4NAME = 'flag-guided Salmon second pass'
A4COL = PIPE[A4NAME][0]
SAMPLES = ['Keele kidney CC001', 'Keele liver CC021', 'Munger A/J', 'Munger CAST']
SLAB = ['Keele kidney\nCC001, 50 bp', 'Keele liver\nCC021, 50 bp†', 'Munger A/J\n100 bp', 'Munger CAST\n100 bp']

def a4(ax, x, y):
    pmark(ax, x, y, A4NAME, size=6.0, z=4)

# a: data figS12a_arm4_accuracy_u5.tsv (flagged genes, per sample)
acc = rd('figS12a_arm4_accuracy_u5.tsv'); acc = acc[acc.subset == 'flagged']
for j, (col, lab) in enumerate((('median_abs_lr', 'Median |log$_2$ error|, flagged genes'), ('mad_sd_lr', 'MAD-SD of log$_2$ error, flagged genes'))):
    ax = ax_mm(fig, 16 + j * 48, 10, 36, 44)
    for i, smp in enumerate(SAMPLES):
        d = acc[acc['sample'] == smp].set_index('method')
        ys = [d.loc['hybrid_stored', col], d.loc['arm4', col], d.loc['salmon_pan_full', col]]
        ax.plot([i - 0.2, i, i + 0.2], ys, color='#BBBBBB', lw=0.6, zorder=1)
        pmark(ax, i - 0.2, ys[0], 'PLUCK', size=3.4); a4(ax, i, ys[1]); pmark(ax, i + 0.2, ys[2], 'Salmon-pan', size=3.4)
    ax.set_xticks(range(4)); ax.set_xticklabels(SLAB, rotation=0, linespacing=1.0, fontsize=7); ax.tick_params(axis='x', length=0)
    for t in ax.get_xticklabels(): t.set_fontsize(7)
    ax.set_xlim(-0.5, 3.5); ax.set_ylim(0, 1.3); ax.set_ylabel(lab)
    for t in ax.get_xticklabels(): t.set_rotation(35); t.set_ha('right')
letter(fig, 1, 6, 'a')

# b: data figS12b_arm4_timing_u5.tsv (median and range over 3 repetitions)
tt = rd('figS12b_arm4_timing_u5.tsv').set_index('sample')
ax = ax_mm(fig, 125, 10, 50, 44)
for i, smp in enumerate(SAMPLES):
    r = tt.loc[smp]
    ax.plot([i, i], [r.ratio_SR_over_A4_min, r.ratio_SR_over_A4_max], color=A4COL, lw=0.8)
    a4(ax, i + 0.0, r.ratio_SR_over_A4)
    pmark(ax, i + 0.0, r.ratio_SR_over_H, 'PLUCK', size=3.4)
    up = r.ratio_SR_over_H > 10   # PLUCK alone, same machine: label above the point when the ratio exceeds 10
    ax.text(i, r.ratio_SR_over_H * (1.13 if up else 0.89), '%.1f×' % r.ratio_SR_over_H, ha='center', va='bottom' if up else 'top', color=ORANGE)
ax.set_yscale('log'); ax.set_ylim(0.8, 30); ax.set_yticks([1, 2, 5, 10, 20]); ax.set_yticklabels(['1×', '2×', '5×', '10×', '20×']); ax.minorticks_off()
ax.set_xticks(range(4)); ax.set_xticklabels(SLAB, linespacing=1.0); ax.tick_params(axis='x', length=0)
for t in ax.get_xticklabels(): t.set_rotation(35); t.set_ha('right')
ax.set_xlim(-0.5, 3.5); ax.set_ylabel('Salmon-ref wall time / pipeline wall time\n(cloud virtual machine, 4 vCPU)', linespacing=1.1)
letter(fig, 112, 6, 'b')
# key for a and b
ky = 76
for k, (lab, kind) in enumerate((('PLUCK', 'p'), ('flag-guided Salmon second pass: PLUCK + Salmon-pan on reads with a flagged-gene probe', 'a'), ('Salmon-pan, all reads', 's'))):
    xx = [16, 34, 150][k]
    lax = ax_mm(fig, xx, ky - 1.5, 4, 3); lax.axis('off'); lax.set_xlim(0, 1); lax.set_ylim(0, 1)
    if kind == 'p': pmark(lax, 0.5, 0.5, 'PLUCK', size=3.4)
    elif kind == 's': pmark(lax, 0.5, 0.5, 'Salmon-pan', size=3.4)
    else: a4(lax, 0.5, 0.5)
    mtext(fig, xx + 5, ky, lab, va='center')

# c: data fig5e_ablation_arms.tsv (arm 2 variants and PLUCK on the same genes), fig5e_ablation_tail_cis_scan.tsv (tail share, detections)
ab = rd('fig5e_ablation_arms.tsv'); tl = rd('fig5e_ablation_tail_cis_scan.tsv')
ARMS = [('arm2_T1_S1_P', 'TSL1\nisoforms,\n1 gene'), ('arm2_TA_S1_P', 'all\ntranscripts,\n1 gene'), ('arm2_T1_S2_P', 'TSL1\nisoforms,\n≤2 genes')]
C2 = PIPE['all gene-unique k-mers'][0]
panels = [('noise_sd', 'Noise SD against cohort truth', None), ('gt2', '% genes in >2× tail', None), ('det', 'cis-eQTL detections\nrelative to PLUCK', None)]
for j, (kind, lab, _) in enumerate(panels):
    ax = ax_mm(fig, 16 + j * 56, 86, 40, 26)
    for i, (arm, al) in enumerate(ARMS):
        if kind == 'noise_sd':
            a = ab[(ab.cohort == 'keele_liver') & (ab.arm == arm) & (ab.subset == 'all')].noise_sd.iloc[0]
            h = ab[(ab.cohort == 'keele_liver') & (ab.arm == f'arm3_hybrid_frozen__on_{arm}_genes') & (ab.subset == 'all')].noise_sd.iloc[0]
        else:
            r = tl[(tl.tissue == 'keele_liver') & (tl.arm == arm)].iloc[0]
            if kind == 'gt2': a, h = r[f'gt2pct_{arm}'], r.gt2pct_hybrid
            else: a, h = r[f'det_{arm}'] / r.det_hybrid, 1.0
        ax.plot([i, i], [a, h], color='#BBBBBB', lw=0.7, zorder=1)
        pmark(ax, i, a, 'all gene-unique k-mers', size=5.0, z=3)
        pmark(ax, i, h, 'PLUCK', size=3.2)
    ax.set_xticks(range(3)); ax.set_xticklabels([a[1] for a in ARMS], linespacing=1.0); ax.tick_params(axis='x', length=0)
    ax.set_xlim(-0.5, 2.5); ax.set_ylabel(lab, linespacing=1.1)
    if kind == 'noise_sd': ax.set_ylim(0.07, 0.09)
    if kind == 'gt2': ax.set_ylim(8, 12)
    if kind == 'det': ax.set_ylim(0.95, 1.01); ax.axhline(1, color='#888888', lw=0.5, zorder=0)
    if j == 0: ax.set_title('Broader k-mer set (Keele liver†)', loc='left', pad=3)
letter(fig, 1, 82, 'c')
lax = ax_mm(fig, 112, 82.5, 4, 3); lax.axis('off'); lax.set_xlim(0, 1); lax.set_ylim(0, 1)
pmark(lax, 0.5, 0.5, 'all gene-unique k-mers', size=5.0)
mtext(fig, 117, 84, 'broader k-mer set', va='center')
lax = ax_mm(fig, 138, 82.5, 4, 3); lax.axis('off'); lax.set_xlim(0, 1); lax.set_ylim(0, 1)
pmark(lax, 0.5, 0.5, 'PLUCK', size=3.2)
mtext(fig, 143, 84, 'PLUCK, same genes', va='center')

save(fig, 'FigS12')
print('wrote FigS12')
