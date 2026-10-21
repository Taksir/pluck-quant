"""Fig. S10: the power model in every cohort (redraw of e6_power/out/e6_power.png).

Panels
  a  cis-eQTL detections vs animals in the truth-free scan, five CC cohort-tissues (power-model grid).
  b  excess residual variance of PLUCK over Salmon-pan vs gene reads per sample, with the fit
     D(T) = delta/T + gamma, six cohort-tissues.
  c  read factor for equal per-gene noise (multiplier on reads, not an SD ratio), by source.
  d  heavy tail: share of genes with more than 2x (and 4x) Salmon-pan's residual variance;
     Salmon-ref beside it as the same-algorithm, different-index baseline.
  e  Munger founders: observed vs predicted within-strain variance ratio by expression bin.
  f  animals needed relative to Salmon-pan for one gene vs its reads per sample (calculator).
Run: python3 -I figS10.py
"""
import sys
# Folder constants below point to the original working folder; set them for your copy.
sys.path.insert(0, '/mnt/project-files/paper/figures/scripts')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
from u6_common import ENC, CODE, MM, mkw, lkw, read, setup, panel_letter, save

setup()
fig = plt.figure(figsize=(180 * MM, 220 * MM))
L, R = 0.075, 0.985
ROW_A = (0.785, 0.945)
ROW_B = (0.545, 0.690)
ROW_C = (0.285, 0.445)
ROW_D = (0.050, 0.200)
gfmt = FuncFormatter(lambda v, p: f'{v:g}')

COH = [('Keele kidney', 50, 'Keele kidney, 50 bp'), ('Keele liver', 50, 'Keele liver, 50 bp†'),
       ('GSE199702 liver', 101, 'GSE199702 liver, 101 bp'), ('GSE199702 heart', 101, 'GSE199702 heart, 101 bp'),
       ('GSE199702 kidney', 101, 'GSE199702 kidney, 101 bp')]


def row_axes(row, n, gap=0.035, x0=L, x1=R):
    w = (x1 - x0 - gap * (n - 1)) / n
    return [fig.add_axes([x0 + i * (w + gap), row[0], w, row[1] - row[0]]) for i in range(n)]


# ---------------------------------------------------------------- a: detection curves
# data: figS10_e6_detection_curves.tsv
dc = read('figS10_e6_detection_curves.tsv')
axa = row_axes(ROW_A, 5, gap=0.045)
for ax, (coh, rl, title) in zip(axa, COH):
    d = dc[dc.cohort == coh]
    for code in ['ref_unique', 'salmon_ref', 'salmon_pan', 'hybrid']:
        lab = CODE[code]
        s = d[d.method == code].sort_values('n_animals')
        ax.plot(s.n_animals, s.mean_detections, **lkw(lab))
        ax.plot(s.n_animals, s.mean_detections, linestyle='none', **mkw(lab, ms_scale=0.6))
    ax.set_title(title.replace(', ', ',\n'), loc='left', fontsize=7)
    ax.set_xlabel('Animals in the scan')
    ax.set_ylim(0, None)
    ax.yaxis.set_major_formatter(gfmt)
axa[0].set_ylabel('cis-eQTL detected (BH 5%)')
from matplotlib.lines import Line2D
hs = []
for lab in ['PLUCK', 'Salmon-pan', 'Salmon-ref', 'ref-unique probes']:
    kw = mkw(lab, ms_scale=0.8); kw.pop('zorder')
    l = lkw(lab)
    hs.append(Line2D([], [], linestyle=l['linestyle'], linewidth=l['linewidth'], label=lab, **kw))
fig.legend(handles=hs, loc='upper right', bbox_to_anchor=(R, 0.998), ncol=4, frameon=False, handlelength=2.4)
panel_letter(fig, axa[0], 'a', dx=-0.065, dy=0.015)

# ---------------------------------------------------------------- b: excess-variance fits
# data: figS10_e6_excess_bins.tsv
ex = read('figS10_e6_excess_bins.tsv')
COHB = [('Keele kidney', 'Keele kidney, 50 bp'), ('Keele liver', 'Keele liver, 50 bp†'),
        ('Munger founders', 'Munger founders, 100 bp'), ('GSE199702 liver', 'GSE199702 liver, 101 bp'),
        ('GSE199702 heart', 'GSE199702 heart, 101 bp'), ('GSE199702 kidney', 'GSE199702 kidney, 101 bp')]
axb = row_axes(ROW_B, 6, gap=0.022, x0=L + 0.02)
for ax, (coh, title) in zip(axb, COHB):
    d = ex[ex.cohort == coh]
    de, ga = d.delta_median.iloc[0], d.gamma_median.iloc[0]
    T = np.logspace(np.log10(40), np.log10(8000), 100)
    ax.plot(T, de / T + ga, color=ENC['PLUCK'][0], lw=1.0)
    ax.plot(d.bin_T_median_reads, d.D_median, linestyle='none', **mkw('PLUCK', ms_scale=0.6))
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlim(40, 8000); ax.set_ylim(2e-4, 6e-2)
    ax.xaxis.set_major_locator(FixedLocator([100, 1000]))
    ax.xaxis.set_major_formatter(gfmt)
    # plain decimals on the log y axis: no mathtext superscripts (Elsevier >= 6 pt rule)
    ax.yaxis.set_major_locator(FixedLocator([1e-3, 1e-2]))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: f'{v:g}'))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_title(title.replace(', ', ',\n'), loc='left', fontsize=7)
    ax.text(0.04, 0.05, f'δ = {de:.2f}\nγ = {ga:.5f}', transform=ax.transAxes, fontsize=7)
    if ax is not axb[0]:
        ax.tick_params(labelleft=False)
axb[0].set_ylabel('Excess residual variance\nof PLUCK over Salmon-pan')
# one shared x label for the six facets (per-facet labels run into each other at 7 pt)
_b0, _b1 = axb[0].get_position(), axb[-1].get_position()
fig.text((_b0.x0 + _b1.x1) / 2, _b0.y0 - 0.024, 'Gene reads per sample (T)', ha='center', va='top', fontsize=7)
panel_letter(fig, axb[0], 'b', dx=-0.065, dy=0.018)

# ---------------------------------------------------------------- c: read factor
# data: figS10_e6_read_factor.tsv
rf = read('figS10_e6_read_factor.tsv')
axc = fig.add_axes([0.205, ROW_C[0], 0.27, ROW_C[1] - ROW_C[0]])
axd = fig.add_axes([0.695, ROW_C[0], 0.29, ROW_C[1] - ROW_C[0]])
axe = fig.add_axes([0.105, ROW_D[0], 0.37, ROW_D[1] - ROW_D[0]])
axf = fig.add_axes([0.615, ROW_D[0], 0.37, ROW_D[1] - ROW_D[0]])
for i, r in rf.iterrows():
    sim = r.source.startswith('E2')
    col = ENC['PLUCK'][0]
    if r.read_factor_high > r.read_factor_low:
        axc.plot([r.read_factor_low, r.read_factor_high], [i, i], color=col, lw=1.4, solid_capstyle='butt')
    kw = mkw('PLUCK', ms_scale=0.8)
    if sim:
        # simulation: a neutral bar marker that no pipeline uses
        kw = dict(marker='|', markersize=7, markeredgewidth=1.3, color=ENC['PLUCK'][0])
    axc.plot([r.read_factor_low], [i], linestyle='none', **kw)
    axc.text(r.read_factor_high + 0.04, i, r.read_factor_text.replace('-', '–').replace('<=', '≤'),
             fontsize=7, va='center')
axc.set_yticks(range(len(rf)))
axc.set_yticklabels([f"{s.replace('E2 simulation', 'Simulation')}, {l}" + ('†' if s == 'Keele liver' else '') for s, l in zip(rf.source, rf.layout)], fontsize=7)
axc.invert_yaxis()
axc.axvline(1, color='black', lw=0.6)
axc.set_xlim(0.95, 2.6)
axc.set_xlabel('Read factor for equal noise\n(× Salmon-pan’s reads)')
axc.set_title('Reads for equal noise', loc='left', fontsize=7)
panel_letter(fig, axc, 'c', dx=-0.195, dy=0.012)

# ---------------------------------------------------------------- d: tail shares
# data: figS10_e6_tail.tsv
tl = read('figS10_e6_tail.tsv')
ordd = ['Keele kidney', 'Keele liver', 'Munger founders', 'GSE199702 liver', 'GSE199702 heart', 'GSE199702 kidney']
for i, coh in enumerate(ordd):
    r = tl[tl.cohort == coh].iloc[0]
    axd.barh(i - 0.18, 100 * r.pluck_frac_gt2, height=0.32, color=ENC['PLUCK'][0])
    axd.barh(i - 0.18, 100 * r.pluck_frac_gt4, height=0.32, color='#8C2D04')
    axd.barh(i + 0.18, 100 * r.salmon_ref_frac_gt2, height=0.32, color=ENC['Salmon-ref'][0])
    axd.text(100 * r.pluck_frac_gt2 + 0.2, i - 0.18, f'{100 * r.pluck_frac_gt2:.1f}%', fontsize=7, va='center')
    axd.text(100 * r.salmon_ref_frac_gt2 + 0.2, i + 0.18, f'{100 * r.salmon_ref_frac_gt2:.1f}%', fontsize=7,
             va='center', color='0.3')
axd.set_yticks(range(len(ordd)))
axd.set_yticklabels(['Keele kidney, 50 bp', 'Keele liver, 50 bp†', 'Munger founders, 100 bp', 'GSE199702 liver, 101 bp',
                     'GSE199702 heart, 101 bp', 'GSE199702 kidney, 101 bp'], fontsize=7)
axd.invert_yaxis()
axd.set_xlim(0, 12)
axd.set_xlabel('Genes above the ratio (%)')
axd.set_title('Heavy tail', loc='left', fontsize=7)
from matplotlib.patches import Patch
axd.legend(handles=[Patch(color=ENC['PLUCK'][0], label='PLUCK >2×'), Patch(color='#8C2D04', label='PLUCK >4×'),
                    Patch(color=ENC['Salmon-ref'][0], label='Salmon-ref >2×')],
           loc='lower right', frameon=False, fontsize=7, handlelength=1.2)
panel_letter(fig, axd, 'd', dx=-0.185, dy=0.012)

# ---------------------------------------------------------------- e: Munger model check
# data: figS10_e6_munger_model_check.tsv
mc = read('figS10_e6_munger_model_check.tsv')
axe.plot(mc.T_median, mc.hybrid_ratio_pred_med, color=ENC['PLUCK'][0], lw=0.9, ls='--', label='PLUCK, predicted (simulation α)')
axe.plot(mc.T_median, mc.hybrid_ratio_obs_med, linestyle='-', lw=0.9, label='PLUCK, observed', **{k: v for k, v in mkw('PLUCK', 0.7).items()})
axe.plot(mc.T_median, mc.salmon_ref_ratio_obs_med, linestyle='-', lw=0.8, label='Salmon-ref, observed', **mkw('Salmon-ref', 0.8))
axe.axhline(1, color='black', lw=0.6)
axe.set_xscale('log')
axe.xaxis.set_major_formatter(gfmt)
axe.set_ylim(0.98, 1.12)
axe.set_xlabel('Gene reads per sample (bin median)')
axe.set_ylabel('Within-strain variance\nratio to Salmon-pan')
axe.set_title('Munger founders, 100 bp', loc='left', fontsize=7)
axe.legend(frameon=False, fontsize=7, loc='upper right', handlelength=2.2)
panel_letter(fig, axe, 'e', dx=-0.095, dy=0.012)

# ---------------------------------------------------------------- f: calculator
# data: figS10_e6_extra_n_grid.tsv (resid_sd 0.25 log2 shown; 0.15 and 0.35 in Table S11)
g = read('figS10_e6_extra_n_grid.tsv')
g = g[g.resid_sd == 0.25]
# PLUCK only, layouts told apart by line style, width and direct legend labels (no pipeline marker shapes)
styles = [('50 SE', 'extra_n_real_median_gene', '50 bp SE, Keele liver† real fit', '-', 1.3),
          ('50 SE', 'extra_n_sim', '50 bp SE, simulation', '--', 0.9),
          ('101 SE', 'extra_n_real_median_gene', '100 bp SE, Munger real fit', '-', 0.7),
          ('101 PE', 'extra_n_sim', '101 bp PE, simulation', ':', 0.9)]
for lay, col, lab, ls, lw in styles:
    s = g[g.layout == lay].sort_values('T')
    axf.plot(s['T'], s[col], color=ENC['PLUCK'][0], ls=ls, lw=lw, label=lab)
axf.axhline(1, color='black', lw=0.6)
axf.set_xscale('log')
axf.xaxis.set_major_formatter(gfmt)
axf.set_xlabel('Gene reads per sample (T)')
axf.set_ylabel('Animals relative to Salmon-pan')
axf.set_title('Calculator: one gene, residual SD 0.25 log2', loc='left', fontsize=7)
axf.legend(frameon=False, fontsize=7, loc='upper right', handlelength=2.2)
panel_letter(fig, axf, 'f', dx=-0.095, dy=0.012)

save(fig, 'FigS10')
fig.savefig('FigS10_preview.png', dpi=200)
print('done')
