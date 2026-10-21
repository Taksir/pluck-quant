"""Fig. 6: PLUCK's place on the speed-accuracy frontier, priced in animals.

Panels
  a  cis-eQTL detections vs number of strains (strain-subsampling scan, Keele 50 bp), PLUCK and Salmon-pan on
     100/25/10% of reads; right: detection-matched extra animals in every cohort-tissue.
  b  strains needed to match PLUCK's 47-strain yield vs wall time per sample relative to PLUCK
     (4 threads, paired blocks on the main cluster, 50 bp samples).
  c  wall time per sample vs threads, three timed samples (main cluster, shared node), plus the
     second-cluster blocks for RNA-Skim drawn apart as a separate machine.
  d  speed-accuracy frontier: wall-time ratio to PLUCK at 4 threads vs within-gene Spearman and
     noise SD (scoring set A), Keele kidney (50 bp timing layout) and Munger founders (101 bp timing
     layout as the nearest timed layout; no 100 bp timing exists).
Run: python3 -I fig6.py  (reads figures/data only; writes figures/Fig6.pdf and .svg)
"""
import sys
# Folder constants below point to the original working folder; set them for your copy.
sys.path.insert(0, '/mnt/project-files/paper/figures/scripts')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, NullFormatter, FuncFormatter
from u6_common import ENC, CODE, MM, mkw, lkw, read, setup, panel_letter, save

setup()
fig = plt.figure(figsize=(180 * MM, 220 * MM))

# ---------------------------------------------------------------- layout (figure fractions)
L, R = 0.075, 0.985
ROW_A = (0.800, 0.955)
ROW_B = (0.588, 0.738)
ROW_C = (0.335, 0.478)
ROW_D1 = (0.170, 0.258)
ROW_D2 = (0.045, 0.133)


def ax_at(x0, x1, row):
    return fig.add_axes([x0, row[0], x1 - x0, row[1] - row[0]])


TISSUES = [('keele_kidney', 'Keele kidney, 50 bp'), ('keele_liver', 'Keele liver, 50 bp†')]
fmt_int = FuncFormatter(lambda v, p: f'{v:g}')

# ================================================================ panel a
# data: fig6a_r7_detection_curves.tsv (strain-subsampling scan), fig6a_r7_scan_summary.tsv (pan-equivalent strains),
#       fig6a_e6_extra_animals.tsv (extra animals, all cohort-tissues)
cur = read('fig6a_r7_detection_curves.tsv')
summ = read('fig6a_r7_scan_summary.tsv')
summ = summ[summ.variant == 'asis']
a_axes = [ax_at(L, 0.300, ROW_A), ax_at(0.440, 0.655, ROW_A)]
for ax, (tis, title) in zip(a_axes, TISSUES):
    d = cur[cur.tissue == tis]
    for code in ['pan_p10', 'pan_p25', 'salmon_pan', 'hybrid']:
        lab = CODE[code]
        s = d[d.method == code].sort_values('n_strains')
        ax.plot(s.n_strains, s.mean_detections, **lkw(lab), label=lab)
        last = s.iloc[-1]
        ax.plot([last.n_strains], [last.mean_detections], linestyle='none', **mkw(lab))
    hyb47 = d[(d.method == 'hybrid') & (d.n_strains == 47)].mean_detections.iloc[0]
    eq = summ[(summ.tissue == tis) & (summ.arm.str.startswith('hybrid'))].pan_equiv_n.iloc[0]
    ax.axhline(hyb47, color=ENC['PLUCK'][0], lw=0.5, ls=(0, (2, 2)), zorder=1)
    ax.plot([eq, eq], [0, hyb47], color=ENC['Salmon-pan'][0], lw=0.6, ls=(0, (2, 2)), zorder=1)
    top_a = 4900 if tis == 'keele_kidney' else 3150
    ax.annotate(f'Salmon-pan reaches PLUCK’s\n47-strain yield at {eq:.1f} strains',
                xy=(eq, hyb47), xytext=(8, top_a * 0.79), fontsize=7, va='bottom',
                arrowprops=dict(arrowstyle='-', lw=0.4, color='0.35'))
    # direct labels at line ends
    top = 4300 if tis == 'keele_kidney' else 2800
    for code, dy in [('hybrid', -0.035), ('salmon_pan', 0.035), ('pan_p25', 0.0), ('pan_p10', 0.0)]:
        v = d[(d.method == code) & (d.n_strains == 47)].mean_detections.iloc[0]
        ax.text(48.6, v + dy * top, CODE[code].replace(' 25%', '\n25%').replace(' 10%', '\n10%'), fontsize=7, va='center', ha='left', clip_on=False,
                color=ENC[CODE[code]][0] if code != 'pan_p10' else '#5E8FB0')
    ax.axvline(47, color='0.5', lw=0.5, ls=':', zorder=0)
    ax.set_xlim(6, 47.5)
    ax.set_ylim(0, 4900 if tis == 'keele_kidney' else 3150)
    ax.set_xlabel('Strains in the scan (one animal each)')
    ax.set_title(title, loc='left', fontsize=7)
    ax.xaxis.set_major_locator(FixedLocator([10, 20, 30, 40, 47]))
a_axes[0].set_ylabel('cis-eQTL detected (BH 5%)')
panel_letter(fig, a_axes[0], 'a', dx=-0.065, dy=0.012)

# a, right: extra animals across cohort-tissues
e6 = read('fig6a_e6_extra_animals.tsv')
axs = ax_at(0.795, R, ROW_A)
order = [('Keele kidney', 'Keele kidney, 50 bp', 'extra_n_detection_matched'),
         ('Keele liver', 'Keele liver, 50 bp†', 'extra_n_detection_matched'),
         ('Munger founders', 'Munger, 100 bp, strain test', 'noncentrality_ratio_strong_genes'),
         ('GSE199702 liver', 'GSE199702 liver, 101 bp', 'extra_n_detection_matched'),
         ('GSE199702 heart', 'GSE199702 heart, 101 bp', 'extra_n_detection_matched'),
         ('GSE199702 kidney', 'GSE199702 kidney, 101 bp', 'extra_n_detection_matched')]
for i, (coh, lab, col) in enumerate(order):
    v = e6[e6.cohort == coh][col].iloc[0]
    # Munger is a noncentrality ratio, not detection-matched extra animals: dashed stem
    axs.plot([1, v], [i, i], color=ENC['PLUCK'][0], lw=0.9, zorder=2,
             ls=(0, (2, 1.2)) if coh == 'Munger founders' else '-')
    axs.plot([v], [i], linestyle='none', **mkw('PLUCK'))
    axs.text(v + 0.003, i, f'{v:.3f}', fontsize=7, va='center')
    axs.text(1.002, i - 0.32, lab, fontsize=7, va='bottom', ha='left', color='0.15')
axs.axvline(1, color='black', lw=0.6)
axs.set_yticks([])
axs.spines['left'].set_visible(False)
axs.set_xlim(0.995, 1.12)
axs.set_xticks([1.0, 1.05, 1.10])
axs.set_xlabel('Animals vs Salmon-pan', fontsize=7)
axs.set_title('Animals for the same detections', loc='left', fontsize=7)
axs.set_ylim(len(order) - 0.5, -0.9)

# ================================================================ panel b
# data: fig6b_animals_vs_time.tsv (strain-scan values at 47 strains, wall ratios over the two 50 bp samples),
#       fig6a_r7_scan_summary.tsv (extrapolated strains, inference),
#       table2_timing_4thr.tsv (Salmon-pan 100% wall ratio for the 50 bp samples)
avt = read('fig6b_animals_vs_time.tsv')
t4 = read('table2_timing_4thr.tsv')
b_axes = [ax_at(L, 0.335, ROW_B), ax_at(0.425, 0.685, ROW_B)]
ARM_SUMMARY = {'Salmon-pan 25%': 'Salmon-pan, 25% of reads', 'Salmon-pan 10%': 'Salmon-pan, 10% of reads',
               'Salmon-ref 25%': 'Salmon-ref, 25% of reads', 'Salmon-ref 10%': 'Salmon-ref, 10% of reads'}
pan_full = t4[(t4.layout == 50) & (t4.tool == 'Salmon-pan')].ratio_median
for ax, (tis, title) in zip(b_axes, TISSUES):
    d = avt[avt.tissue == tis]
    s = summ[summ.tissue == tis]
    # reference lines: PLUCK's time and PLUCK's strains
    ax.axvline(1, color=ENC['PLUCK'][0], lw=0.5, ls=(0, (2, 2)), zorder=0)
    ax.axhline(47, color='0.5', lw=0.5, ls=':', zorder=0)
    # Salmon-pan on all reads: fewer strains, 74-81x the time
    eq = s[s.arm.str.startswith('hybrid')].pan_equiv_n.iloc[0]
    ax.plot([pan_full.min(), pan_full.max()], [eq, eq], **lkw('Salmon-pan'))
    for x in (pan_full.min(), pan_full.max()):
        ax.plot([x], [eq], linestyle='none', **mkw('Salmon-pan'))
    ax.text(pan_full.max() * 1.3, eq, 'Salmon-pan', fontsize=7, va='center', color=ENC['Salmon-pan'][0])
    for lab in ['Salmon-pan 25%', 'Salmon-pan 10%', 'Salmon-ref 25%', 'Salmon-ref 10%']:
        r = d[d.tool == lab].iloc[0]
        x0, x1 = r.wall_ratio_to_hybrid_4thr_median_min_over_samples, r.wall_ratio_to_hybrid_4thr_median_max_over_samples
        y = 47 * r.strains_ratio_to_hybrid_proportional
        ax.plot([x0, x1], [y, y], **lkw(lab))
        for x in (x0, x1):
            ax.plot([x], [y], linestyle='none', **mkw(lab))
        xm = np.sqrt(x0 * x1)
        row = s[s.arm == ARM_SUMMARY[lab]].iloc[0]
        if '25%' in lab:
            # extrapolated strains (inference): hollow marker, dashed line
            ytop = max(row.strains_needed_to_match_hybrid_extrap_lin, row.strains_needed_to_match_hybrid_extrap_loglog)
            ax.plot([xm, xm], [y, ytop], color=ENC[lab][0], lw=0.7, ls=(0, (2, 1.5)), zorder=2)
            ax.plot([xm], [ytop], linestyle='none', marker=ENC[lab][1], markersize=2.2,
                    markerfacecolor='white', markeredgecolor=ENC[lab][0], markeredgewidth=0.6, zorder=3)
            if 'pan' in lab:
                ax.text(xm * 1.12, ytop, 'inference', fontsize=7, ha='left', va='center', color='0.35')
            else:
                ax.text(xm / 1.12, ytop, 'inference', fontsize=7, ha='right', va='center', color='0.35')
        else:
            # proportional reading is a lower bound: arrow upward
            ax.annotate('', xy=(xm, y + 13), xytext=(xm, y + 1.5),
                        arrowprops=dict(arrowstyle='-|>', lw=0.6, color=ENC[lab][0], mutation_scale=5))
        if 'pan' in lab or (tis == 'keele_kidney' and lab == 'Salmon-ref 25%'):
            # kidney Salmon-ref 25% sits just above the 47-strain line, so its label goes to the right
            lx, ly, ha, va = x1 * 1.2, y, 'left', 'center'
        else:
            lx, ly, ha, va = xm, y - 3.4, 'center', 'top'
        ax.text(lx, ly, lab, fontsize=7, va=va, ha=ha,
                color=ENC[lab][0] if lab != 'Salmon-pan 10%' else '#5E8FB0')
    ax.plot([1], [47], linestyle='none', **mkw('PLUCK'))
    ax.text(1.12, 47 - 4.5, 'PLUCK', fontsize=7, color=ENC['PLUCK'][0], va='top')
    ax.set_xscale('log')
    ax.set_xlim(0.7, 300)
    ax.set_ylim(30, 112)
    ax.xaxis.set_major_locator(FixedLocator([1, 3, 10, 30, 100]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, p: f'{v:g}×'))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel('Wall time per sample relative to PLUCK (4 threads)')
    ax.set_title(title, loc='left', fontsize=7)
b_axes[0].set_ylabel('Strains needed for PLUCK’s\nyield at 47 strains')
# (panel b symbol note dropped for the 7 pt layout: the caption explains the symbols and each extrapolation has a direct 'inference' label)
panel_letter(fig, b_axes[0], 'b', dx=-0.065, dy=0.012)

# shared legend for every pipeline in the figure (right of row b)
from matplotlib.lines import Line2D
handles = []
for lab in ['PLUCK', 'Salmon-pan', 'Salmon-ref', 'STAR+fC', 'kallisto-pan', 'kallisto-ref', 'Kmerator+countTags',
            'RNA-Skim', 'Salmon-pan 25%', 'Salmon-pan 10%', 'Salmon-ref 25%', 'Salmon-ref 10%', 'GBRS']:
    kw = mkw(lab)
    kw.pop('zorder', None)
    ls = lkw(lab)
    if lab == 'GBRS':
        h = Line2D([], [], linestyle='none', **kw, label='GBRS (timing only)')
    else:
        h = Line2D([], [], linestyle=ls['linestyle'], linewidth=ls['linewidth'], **kw, label=lab)
    handles.append(h)
leg_ax = fig.add_axes([0.745, ROW_B[0] + 0.005, 0.24, ROW_B[1] - ROW_B[0] - 0.01])
leg_ax.axis('off')
leg_ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(0.745, 0.760), bbox_transform=fig.transFigure, frameon=False, handlelength=2.6, borderaxespad=0,
              labelspacing=0.28, fontsize=7)

# ================================================================ panel c
# data: fig6c_timing_ratios.tsv (main-cluster medians and min/max over 3 reps; kind=gbrs one rep),
#       fig6c_rnaskim_timing_second_cluster.tsv (second-cluster paired blocks, 4 threads)
tr = read('fig6c_timing_ratios.tsv')
reg = tr[tr.kind == 'regular']
gb = tr[(tr.kind == 'gbrs') & (tr.method == 'gbrs')]
SAMPLES = [('50', 'CC021', 'Keele liver CC021, 50 bp'), ('50', 'CC001', 'Keele liver CC001, 50 bp'),
           ('101', 'LiverCC001F', 'GSE199702 liver CC001 F, 101 bp')]
c_x = [(L, 0.285), (0.330, 0.540), (0.585, 0.795)]
c_axes = []
FULL = ['Kmerator+countTags', 'kallisto-ref', 'kallisto-pan', 'STAR+fC', 'Salmon-ref', 'Salmon-pan', 'PLUCK']
for (lay, smp, title), (x0, x1) in zip(SAMPLES, c_x):
    ax = ax_at(x0, x1, ROW_C)
    c_axes.append(ax)
    d = reg[(reg.layout.astype(str) == lay) & (reg['sample'] == smp)]
    for lab in FULL:
        code = [k for k, v in CODE.items() if v == lab]
        s = d[d.method.isin(code)].sort_values('T')
        if s.empty:
            continue
        ax.plot(s['T'], s.wall_median_s, **lkw(lab))
        for _, r in s.iterrows():
            kw = mkw(lab)
            ax.errorbar([r['T']], [r.wall_median_s],
                        yerr=[[r.wall_median_s - r.wall_min], [r.wall_max - r.wall_median_s]],
                        linestyle='none', elinewidth=0.6, capsize=1.2, capthick=0.6, ecolor=kw['color'], **kw)
    g = gb[(gb.layout.astype(str) == lay) & (gb['sample'] == smp)]
    for _, r in g.iterrows():
        ax.plot([r['T']], [r.wall_median_s], linestyle='none', **mkw('GBRS'))
    ax.set_xscale('log', base=2)
    ax.set_yscale('log')
    ax.set_xlim(0.8, 10)
    ax.set_ylim(3, 5000)
    ax.xaxis.set_major_locator(FixedLocator([1, 4, 8]))
    ax.xaxis.set_major_formatter(fmt_int)
    ax.xaxis.set_minor_locator(FixedLocator([]))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: f'{v:g}'))
    ax.set_xlabel('Threads')
    ax.set_title(title, loc='left', fontsize=7)
    if lay == '50':
        ax.text(2.4, 38, 'kallisto-pan not timed', fontsize=7, color=ENC['kallisto-pan'][0], ha='center')
c_axes[0].set_ylabel('Wall time per sample (s)')
for ax in c_axes[1:]:
    ax.tick_params(labelleft=False)
c_axes[0].text(1.02, 2585.7 * 1.25, 'GBRS: one rep.', fontsize=7, color='0.45', ha='left')
fig.text(L - 0.035, ROW_C[1] + 0.036, 'Compute node (AMD EPYC 9654), shared with other jobs; median of 3 repetitions, bars min to max (1 thread: one run)',
         fontsize=7, color='0.25', va='bottom')
panel_letter(fig, c_axes[0], 'c', dx=-0.065, dy=0.032)

# second-cluster blocks: separate machine, separated by a gap and a machine label
sy = read('fig6c_rnaskim_timing_second_cluster.tsv')
axy = ax_at(0.855, R, ROW_C)
xs = {'CC021': 0, 'CC001': 1, 'LiverCC001F': 2}
blocks = []
for smp in ['CC021', 'CC001', 'LiverCC001F']:
    s = sy[sy['sample'] == smp]
    rn = s.iloc[0]
    blocks.append((smp, 'RNA-Skim', rn.median_wall_rnaskim_s))
    for _, r in s.iterrows():
        blocks.append((smp, CODE[r.other], r.median_wall_other_s))
off = {'PLUCK': -0.24, 'Salmon-ref': -0.08, 'Kmerator+countTags': 0.08, 'RNA-Skim': 0.24}
for smp, lab, w in blocks:
    axy.plot([xs[smp] + off[lab]], [w], linestyle='none', **mkw(lab))
axy.set_yscale('log')
axy.set_ylim(3, 5000)
axy.set_xlim(-0.5, 2.5)
axy.set_xticks([0, 1, 2])
axy.set_xticklabels(['CC021, 50 bp', 'CC001, 50 bp', 'Liver CC001 F,\n101 bp'], fontsize=7, rotation=35, ha='right', rotation_mode='anchor')
axy.tick_params(labelleft=False)
axy.set_title('Second cluster,\n4 threads', loc='left', fontsize=7)
axy.text(2.3, 610.1 * 1.6, 'RNA-Skim', fontsize=7, ha='right', color=ENC['RNA-Skim'][0])

# ================================================================ panel d
# data: table2_accuracy_frameA.tsv (scoring set A, all core genes: within-gene Spearman, noise SD),
#       table2_timing_4thr.tsv (wall ratio to PLUCK, 4 threads, median of 3 reps per sample)
acc = read('table2_accuracy_frameA.tsv')
COLS = [('Keele kidney', '50', 'Keele kidney, 50 bp (timing: Keele liver 50 bp libraries)'),
        ('Munger founders', '101', 'Munger founders, 100 bp (timing: 101 bp sample)')]
D_LABS = ['PLUCK', 'Salmon-pan', 'Salmon-ref', 'STAR+fC', 'kallisto-pan', 'kallisto-ref', 'Kmerator+countTags']
d_axes = {}
d_x = [(L, 0.505), (0.555, R)]
for j, (coh, lay, title) in enumerate(COLS):
    for i, (metric, row, ylab) in enumerate([('within_gene_spearman', ROW_D1, 'Within-gene Spearman'),
                                             ('noise_sd', ROW_D2, 'Noise SD vs cohort truth')]):
        ax = ax_at(d_x[j][0], d_x[j][1], row)
        d_axes[(i, j)] = ax
        a = acc[acc.cohort == coh]
        t = t4[t4.layout.astype(str) == lay]
        for lab in reversed(D_LABS):
            aa = a[a.tool == ('hybrid' if lab == 'PLUCK' else lab)]
            tt = t[t.tool == ('hybrid' if lab == 'PLUCK' else lab)]
            if aa.empty or tt.empty:
                continue
            y = aa[metric].iloc[0]
            xr = tt.ratio_median.values
            if len(xr) > 1:
                ax.plot([xr.min(), xr.max()], [y, y], **lkw(lab, base_lw=0.8))
            for x in xr:
                ax.plot([x], [y], linestyle='none', **mkw(lab))
        # the gap to Salmon-pan, named
        yp = a[a.tool == 'hybrid'][metric].iloc[0]
        ys = a[a.tool == 'Salmon-pan'][metric].iloc[0]
        ax.plot([1.45, 1.45], [yp, ys], color='0.3', lw=0.6)
        for yy in (yp, ys):
            ax.plot([1.35, 1.45], [yy, yy], color='0.3', lw=0.6)
        if metric == 'within_gene_spearman':
            txt = f'{ys - yp:.4f} below\nSalmon-pan'
        else:
            txt = f'{yp / ys:.2f}× Salmon-pan'
        ax.text(1.6, (yp + ys) / 2, txt, fontsize=7, va='center', color='0.2')
        ax.set_xscale('log')
        ax.set_xlim(0.7, 130)
        ax.xaxis.set_major_locator(FixedLocator([1, 3, 10, 30, 100]))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, p: f'{v:g}×'))
        ax.xaxis.set_minor_formatter(NullFormatter())
        if metric == 'within_gene_spearman':
            ax.set_ylim(0.938, 0.999)
            ax.set_title(title, loc='left', fontsize=7)
            ax.tick_params(labelbottom=False)
        else:
            ax.set_ylim(0, 0.12)
            ax.set_xlabel('Wall time per sample relative to PLUCK (4 threads)')
        if j == 0:
            ax.set_ylabel(ylab)
        if coh == 'Keele kidney' and metric == 'within_gene_spearman':
            pass
        if coh == 'Munger founders' and metric == 'noise_sd':
            # data: figS12b_arm4_timing_u5.tsv (arm-4 pilot, cloud VM, two Munger samples)
            a4 = read('figS12b_arm4_timing_u5.tsv')
            mr = sorted(a4[a4.cohort == 'munger'].ratio_SR_over_H)
            ax.text(0.75, 0.004, 'Two Munger samples timed on a cloud virtual\n'
                    f'machine: Salmon-ref {mr[0]:.1f}× and {mr[-1]:.1f}× PLUCK (Fig. S12b)',
                    fontsize=7, color='0.25', ha='left', va='bottom')
        if coh == 'Keele kidney' and metric == 'within_gene_spearman':
            ax.text(67, 0.9405, 'kallisto-pan not run', fontsize=7, color=ENC['kallisto-pan'][0], ha='center')
panel_letter(fig, d_axes[(0, 0)], 'd', dx=-0.065, dy=0.012)

save(fig, 'Fig6')
fig.savefig('Fig6_preview.png', dpi=220)
print('done')
