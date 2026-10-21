"""Fig. S11: every timing block, repetition by repetition; CPU time; memory; cloud replication.

Panels
  a  wall-time ratio to PLUCK in the same paired block, every repetition (dots) and the median (bar),
     every pipeline including the reduced-read arms and GBRS; rows = threads (1, 4, 8),
     columns = the three timed samples (main-cluster node, shared with other jobs).
  b  CPU-time ratio to PLUCK at 4 threads, median of 3 repetitions, one point per sample.
  c  peak resident memory per sample at 4 threads (main-cluster medians; second-cluster block for RNA-Skim).
  d  cloud virtual machine replication (4 vCPU, paired runs): core tier and core + extended tier.
Run: python3 -I figS11.py
"""
import sys
# Folder constants below point to the original working folder; set them for your copy.
sys.path.insert(0, '/mnt/project-files/paper/figures/scripts')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
from u6_common import ENC, CODE, MM, mkw, lkw, read, setup, panel_letter, save

setup()
fig = plt.figure(figsize=(180 * MM, 225 * MM))
xfmt = FuncFormatter(lambda v, p: f'{v:g}×')

ORDER = ['Salmon-pan', 'Salmon-ref', 'STAR+fC', 'kallisto-pan', 'kallisto-ref', 'Kmerator+countTags',
         'Salmon-pan 25%', 'Salmon-pan 10%', 'Salmon-ref 25%', 'Salmon-ref 10%', 'GBRS']
SAMPLES = [('50', 'CC021', 'Keele liver CC021, 50 bp'), ('50', 'CC001', 'Keele liver CC001, 50 bp'),
           ('101', 'LiverCC001F', 'GSE199702 liver CC001 F, 101 bp')]

# ---------------------------------------------------------------- a: per-repetition ratios
# data: fig6c_timing_per_block.tsv (one row per repetition; kind=gbrs rows are one-repetition GBRS blocks),
#       fig6c_timing_ratios.tsv (medians)
pb = read('fig6c_timing_per_block.tsv')
pb['label'] = pb.method.map(CODE)
med = read('fig6c_timing_ratios.tsv')
med['label'] = med.method.map(CODE)
rows = [(1, 0.770, 0.920), (4, 0.580, 0.730), (8, 0.390, 0.540)]
cols = [(0.085, 0.355), (0.395, 0.665), (0.705, 0.975)]
axa = {}
for T, y0, y1 in rows:
    for (lay, smp, title), (x0, x1) in zip(SAMPLES, cols):
        ax = fig.add_axes([x0, y0, x1 - x0, y1 - y0])
        axa[(T, smp)] = ax
        d = pb[(pb.layout.astype(str) == lay) & (pb['sample'] == smp) & (pb['T'] == T)]
        m = med[(med.layout.astype(str) == lay) & (med['sample'] == smp) & (med['T'] == T)]
        # GBRS from its own one-repetition blocks; every other pipeline from the regular blocks only
        d = d[(d.kind == 'regular') | (d.label == 'GBRS')]
        m = m[(m.kind == 'regular') | (m.label == 'GBRS')]
        for i, lab in enumerate(ORDER):
            dd = d[d.label == lab]
            if dd.empty:
                txt = 'not timed' if lab != 'GBRS' else ''
                if txt:
                    ax.text(i, 1.6, txt, rotation=90, fontsize=7, ha='center', va='bottom', color='0.5')
                continue
            jit = np.linspace(-0.18, 0.18, len(dd)) if len(dd) > 1 else [0]
            kw = mkw(lab, ms_scale=0.65)
            ax.plot(i + np.asarray(jit), dd.ratio, linestyle='none', **kw)
            mm = m[m.label == lab]
            if not mm.empty:
                ax.plot([i - 0.3, i + 0.3], [mm.ratio_median.iloc[0]] * 2, color=ENC[lab][0], lw=1.1, zorder=5)
        ax.axhline(1, color=ENC['PLUCK'][0], lw=0.8)
        ax.set_yscale('log')
        ax.set_ylim(1, 1000)
        ax.yaxis.set_major_locator(FixedLocator([1, 3, 10, 30, 100, 300, 1000]))
        ax.yaxis.set_major_formatter(xfmt)
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.set_xlim(-0.6, len(ORDER) - 0.4)
        ax.set_xticks(range(len(ORDER)))
        if T == 8:
            ax.set_xticklabels(ORDER, rotation=60, ha='right', fontsize=7)
        else:
            ax.set_xticklabels([])
        if smp != 'CC021':
            ax.tick_params(labelleft=False)
        if T == 1:
            ax.set_title(title, loc='left', fontsize=7)
        ax.text(0.62, 0.98, f'{T} thread' + ('s' if T > 1 else '') + (' (one run)' if T == 1 else ''),
                transform=ax.transAxes, fontsize=7, ha='center', va='top', color='0.25')
for T, *_ in rows:
    axa[(T, 'CC021')].set_ylabel('Wall time / PLUCK\n(same block)')
# the 8-thread Salmon-ref ratios
for lay, smp, _ in SAMPLES:
    m = med[(med.kind == 'regular') & (med.layout.astype(str) == lay) & (med['sample'] == smp) & (med['T'] == 8) & (med.label == 'Salmon-ref')]
    v = m.ratio_median.iloc[0]
    axa[(8, smp)].annotate(f'{v:.1f}×', xy=(1, v), xytext=(1.6, 2.2), fontsize=7, color='0.2',
                           arrowprops=dict(arrowstyle='-', lw=0.4, color='0.4'))
fig.text(0.085, 0.945, 'Compute node (AMD EPYC 9654), shared with other jobs. Dots: repetitions; bar: median; '
         'GBRS: one repetition, timing only, in its own blocks', fontsize=7, color='0.25')
panel_letter(fig, axa[(1, 'CC021')], 'a', dx=-0.075, dy=0.03)

# ---------------------------------------------------------------- b: CPU-time ratio
# data: fig6c_timing_ratios.tsv (cpu_ratio_median, T = 4, kind = regular)
axb = fig.add_axes([0.16, 0.05, 0.205, 0.19])
reg4 = med[(med.kind == 'regular') & (med['T'] == 4)]
LAB_B = ['Salmon-pan', 'Salmon-ref', 'STAR+fC', 'kallisto-pan', 'kallisto-ref', 'Kmerator+countTags',
         'Salmon-pan 25%', 'Salmon-pan 10%', 'Salmon-ref 25%', 'Salmon-ref 10%']
for i, lab in enumerate(LAB_B):
    s = reg4[reg4.label == lab]
    axb.plot(s.cpu_ratio_median, [i] * len(s), linestyle='none', **mkw(lab, ms_scale=0.8))
axb.axvline(1, color=ENC['PLUCK'][0], lw=0.8)
axb.set_xscale('log')
axb.set_xlim(0.8, 150)
axb.xaxis.set_major_locator(FixedLocator([1, 3, 10, 30, 100]))
axb.xaxis.set_major_formatter(xfmt)
axb.xaxis.set_minor_formatter(NullFormatter())
axb.set_yticks(range(len(LAB_B)))
axb.set_yticklabels(LAB_B, fontsize=7)
axb.invert_yaxis()
axb.set_xlabel('CPU time / PLUCK, 4 threads')
axb.set_title('CPU time (3 samples)', loc='left', fontsize=7)
panel_letter(fig, axb, 'b', dx=-0.125, dy=0.012)

# ---------------------------------------------------------------- c: peak memory
# data: table2_timing_4thr.tsv (rss_mb_median, main cluster), fig6c_rnaskim_peak_rss.tsv (second cluster)
t4 = read('table2_timing_4thr.tsv')
rss = read('fig6c_rnaskim_peak_rss.tsv')
axc = fig.add_axes([0.475, 0.05, 0.20, 0.19])
LAB_C = ['PLUCK', 'Salmon-pan', 'Salmon-ref', 'STAR+fC', 'kallisto-pan', 'kallisto-ref', 'Kmerator+countTags']
for i, lab in enumerate(LAB_C):
    s = t4[t4.tool == ('hybrid' if lab == 'PLUCK' else lab)]
    axc.plot(s.rss_mb_median / 1000, [i] * len(s), linestyle='none', **mkw(lab, ms_scale=0.8))
i = len(LAB_C) + 0.6
s = rss[rss.method == 'rnaskim']
axc.plot(s.rss_mb / 1000, [i] * len(s), linestyle='none', **mkw('RNA-Skim', ms_scale=0.8))
axc.axhline(len(LAB_C) - 0.2, color='0.6', lw=0.5, ls=':')
axc.set_yticks(list(range(len(LAB_C))) + [i])
axc.set_yticklabels(LAB_C + ['RNA-Skim (second cluster)'], fontsize=7)
axc.invert_yaxis()
axc.set_xscale('log')
axc.set_xlim(0.7, 40)
axc.xaxis.set_major_locator(FixedLocator([1, 3, 10, 30]))
axc.xaxis.set_major_formatter(FuncFormatter(lambda v, p: f'{v:g}'))
axc.xaxis.set_minor_formatter(NullFormatter())
axc.set_xlabel('Peak memory (GB), 4 threads')
axc.set_title('Memory', loc='left', fontsize=7)
panel_letter(fig, axc, 'c', dx=-0.155, dy=0.012)

# ---------------------------------------------------------------- d: cloud replication
# data: fig6c_e3_cloud_paired_speedups.tsv (paired repetitions; v1 = core tier, v2 = core + extended tier)
e3 = read('fig6c_e3_cloud_paired_speedups.tsv')
axd = fig.add_axes([0.775, 0.05, 0.21, 0.19])
groups = [('50', 'CC021', 'CC021\n50 bp'), ('50', 'CC001', 'CC001\n50 bp'), ('101', 'liverCC001F', 'Liver\nCC001 F\n101 bp')]
for gi, (lay, smp, _) in enumerate(groups):
    d = e3[(e3.layout.astype(str) == lay) & (e3['sample'] == smp)]
    for col, dx, fill, comp in [('ref_v1', -0.27, True, 'Salmon-ref'), ('ref_v2', -0.09, False, 'Salmon-ref'),
                                ('pan_v1', 0.09, True, 'Salmon-pan'), ('pan_v2', 0.27, False, 'Salmon-pan')]:
        v = d[col].dropna()
        kw = mkw('PLUCK', ms_scale=0.55)
        kw.update(markeredgecolor=ENC[comp][0], markerfacecolor=ENC[comp][0] if fill else 'white')
        kw['color'] = ENC[comp][0]
        axd.plot([gi + dx] * len(v), v, linestyle='none', **kw)
axd.set_yscale('log')
axd.set_ylim(2.6, 140)
axd.yaxis.set_major_locator(FixedLocator([5, 10, 25, 50]))
axd.yaxis.set_major_formatter(xfmt)
axd.yaxis.set_minor_formatter(NullFormatter())
axd.set_xticks(range(3))
axd.set_xticklabels([g[2] for g in groups], fontsize=7)
axd.set_xlim(-0.5, 2.5)
axd.set_ylabel('Salmon time / PLUCK', fontsize=7)
axd.set_title('Cloud virtual machine, 4 vCPU', loc='left', fontsize=7)
from matplotlib.lines import Line2D
h = [Line2D([], [], linestyle='none', marker='o', markersize=3.2, markerfacecolor='0.3', markeredgecolor='0.3', label='core tier'),
     Line2D([], [], linestyle='none', marker='o', markersize=3.2, markerfacecolor='white', markeredgecolor='0.3', label='core + extended tier')]
axd.legend(handles=h, loc='upper left', frameon=False, fontsize=7, handletextpad=0.2, borderaxespad=0.1)
panel_letter(fig, axd, 'd', dx=-0.07, dy=0.012)

save(fig, 'FigS11')
fig.savefig('FigS11_preview.png', dpi=200)
print('done')
