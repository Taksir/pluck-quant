#!/usr/bin/env python3
"""Fig. 2 (simulations with known truth). Vector output figures/Fig2.pdf and Fig2.svg.

Reads only figures/data/. The data file of each panel is named in the comment above the panel.
Run: python3 figures/scripts/fig2.py   (after fig2_recompute_e2.py and fig2_recompute_e2_delta.py)
Text is 7 pt (annotations 6.5 pt); panels c and d each take a full-width row.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.transforms import offset_copy
from fig2_style import P, LAYOUT, FOUNDER, FOUNDERS, GREY, MM, DATA, FS_ANN, setup, mstyle, panel_letter, save

setup()
rd = lambda f: pd.read_csv(f'{DATA}/{f}', sep='\t', comment='#')
LAYS = ['se50', 'se101', 'pe']
MAIN = ['hybrid', 'salmon_pan', 'salmon_ref', 'ref']          # plotting order

fig = plt.figure(figsize=(180 * MM, 222 * MM))
gs0 = fig.add_gridspec(6, 1, height_ratios=[1.1, 0.62, 0.6, 0.92, 1.1, 1.45], hspace=0.55,
                       left=0.1, right=0.98, top=0.965, bottom=0.045)

# ---------------------------------------------------------------- legend shared by a-c
handles = [Line2D([], [], linestyle='none', label=P[c][0], marker=P[c][2], markersize=5 * (1.25 if c == 'hybrid' else 1),
                  markerfacecolor=P[c][1] if P[c][3] else 'white', markeredgecolor='white' if c == 'hybrid' else P[c][1],
                  markeredgewidth=0.6 if c == 'hybrid' else 0.9) for c in MAIN]
fig.legend(handles=handles, loc='upper center', ncol=4, frameon=False, bbox_to_anchor=(0.53, 0.997), handletextpad=0.3,
           columnspacing=1.8)

# ---------------------------------------------------------------- panel a: per-sample error, clean simulated reads, expected truth
# data: fig2a_e2_per_sample_error.tsv (recomputed by fig2_recompute_e2.py); ratios and summary values: fig2a_e2_metrics_all.tsv
ps = rd('fig2a_e2_per_sample_error.tsv'); ps = ps[(ps.cond == 'clean') & (ps.tier == 'core')]
mt = rd('fig2a_e2_metrics_all.tsv'); mt = mt[(mt.cond == 'clean') & (mt.truth == 'expected') & (mt.tier == 'core')]
ga = gs0[0].subgridspec(1, 3, wspace=0.1)
gstrip = gs0[1].subgridspec(1, 3, wspace=0.1)
rng = np.random.default_rng(1)
axa = []
for j, lay in enumerate(LAYS):
    ax = fig.add_subplot(ga[j], sharey=axa[0] if axa else None); axa.append(ax)
    for i, c in enumerate(MAIN):
        v = ps[(ps.layout == lay) & (ps.method == c)].sort_values('sample').err_B6.values
        x = i + (rng.random(len(v)) - 0.5) * 0.36
        ax.scatter(x, v, **mstyle(c, 11))
        ax.plot([i - 0.28, i + 0.28], [np.median(v)] * 2, color=P[c][1], lw=1.4 if c == 'hybrid' else 1.0,
                zorder=6 if c == 'hybrid' else 4)
    r = mt[(mt.layout == lay) & (mt.method == 'hybrid')].iloc[0]
    hv = ps[(ps.layout == lay) & (ps.method == 'hybrid')].err_B6.max()
    ax.text(-0.35, hv + 0.0025, f'{r.M2_ratio_to_pan:.2f}× Salmon-pan', ha='left', va='bottom', fontsize=FS_ANN,
            color=P['hybrid'][1])
    ax.set_xticks(range(len(MAIN))); ax.set_xticklabels(['PLUCK', 'Salmon-\npan', 'Salmon-\nref', 'ref-unique\nprobes'])
    ax.set_xlim(-0.6, len(MAIN) - 0.4)
    ax.set_title(LAYOUT[lay])
    if j == 0:
        ax.set_ylabel('Median |log$_2$ error|\nper sample (B6 genes)')
    else:
        plt.setp(ax.get_yticklabels(), visible=False)
axa[0].set_ylim(0.04, 0.1)
panel_letter(axa[0], 'a', x=-0.3)

# value strip under each facet (values from fig2a_e2_metrics_all.tsv)
for j, lay in enumerate(LAYS):
    ax = fig.add_subplot(gstrip[j]); ax.axis('off')
    h = mt[(mt.layout == lay) & (mt.method == 'hybrid')].iloc[0]
    p = mt[(mt.layout == lay) & (mt.method == 'salmon_pan')].iloc[0]
    nh, npn, n = round(h.M4_spurious * h.n_M4), round(p.M4_spurious * p.n_M4), int(h.n_M4)
    lines = [('Within-gene r²', f'{h.M1_within_r2:.4f}'),
             ('Error vs Salmon-pan', f'{h.M2_ratio_to_pan:.2f}×'),
             ('Founder |Δ|, maximum', f'{h.M3_max_abs_delta:.3f}'),
             ('Spurious genes, vs Salmon-pan', f'{nh} vs {npn}')]   # gene totals (340, 224) in the caption and Table S9
    row = lambda k: offset_copy(ax.transAxes, fig=fig, y=-8.6 * k, units='points')   # fixed 8.6 pt line pitch at 7 pt text
    for k, (name, val) in enumerate(lines):
        ax.text(0.0, 1.0, name, fontsize=FS_ANN, va='top', transform=row(k))
        ax.text(1.0, 1.0, val, fontsize=FS_ANN, va='top', ha='right', transform=row(k))
    ax.add_line(Line2D([0, 1], [1.06, 1.06], transform=ax.transAxes, color=GREY, lw=0.4))

# ---------------------------------------------------------------- panel b: per-founder median shift with gene-bootstrap 95% CI
# data: fig2b_e2_delta_founder_ci.tsv (recomputed by fig2_recompute_e2_delta.py; medians equal fig2b_e2_delta_founder.tsv)
# top row: PLUCK and Salmon-pan only, zoomed to +-0.02; bottom row: all four pipelines on the full scale
dl = rd('fig2b_e2_delta_founder_ci.tsv')
gbz = gs0[2].subgridspec(1, 3, wspace=0.1)
gb = gs0[3].subgridspec(1, 3, wspace=0.1)
axz, axb = [], []
for j, lay in enumerate(LAYS):
    for row, codes, offs in (('zoom', ['hybrid', 'salmon_pan'], [-0.14, 0.14]), ('full', MAIN, [-0.27, -0.09, 0.09, 0.27])):
        lst = axz if row == 'zoom' else axb
        spec = gbz[j] if row == 'zoom' else gb[j]
        ax = fig.add_subplot(spec, sharey=lst[0] if lst else None, sharex=axz[j] if row == 'full' else None); lst.append(ax)
        ax.axhline(0, color='black', lw=0.5, zorder=1)
        off = dict(zip(codes, offs))
        for c in codes[::-1]:
            d = dl[(dl.layout == lay) & (dl.method == c)].set_index('founder').loc[FOUNDERS]
            x = np.arange(7) + off[c]
            ax.errorbar(x, d.median_delta, yerr=[d.median_delta - d.ci_lo, d.ci_hi - d.median_delta], fmt='none',
                        ecolor=P[c][1], elinewidth=0.7, capsize=0, zorder=4 if c == 'hybrid' else 2)
            ax.scatter(x, d.median_delta, **mstyle(c, 12))
        ax.set_xlim(-0.6, 6.6)
        if row == 'zoom':
            ax.set_title(LAYOUT[lay])
            plt.setp(ax.get_xticklabels(), visible=False)
            ax.tick_params(axis='x', length=0)
        else:
            ax.set_xticks(range(7)); ax.set_xticklabels([FOUNDER[f] for f in FOUNDERS], rotation=45, ha='right',
                                                         rotation_mode='anchor')
        if j > 0:
            plt.setp(ax.get_yticklabels(), visible=False)
axz[0].set_ylim(-0.022, 0.022); axz[0].set_yticks([-0.02, 0, 0.02])
axz[0].text(-0.45, 0.0195, 'zoom: PLUCK and Salmon-pan', fontsize=FS_ANN, color=GREY, ha='left', va='top')
axb[0].set_ylim(-0.225, 0.035)
bb0, bb1 = axz[0].get_position(), axb[0].get_position()
fig.text(bb0.x0 - 0.062, (bb0.y1 + bb1.y0) / 2, 'Median founder-haplotype shift Δ\n(log$_2$, hit genes, 95% CI)', rotation=90,
         ha='center', va='center', fontsize=7)
panel_letter(axz[0], 'b', x=-0.3, y=1.12)

# ---------------------------------------------------------------- panel c: planted cis-eQTL simulation
# data: fig2c_idea1_gates.tsv (G1a ratios and 95% CIs); fig2c_idea1_rates.tsv (false-cis rate of each pipeline and of the
#       expected-count truth, CC; DO has no truth row); fig2c_idea1_power.tsv (detections ratios; DO reported only).
g1 = rd('fig2c_idea1_gates.tsv')
dor = rd('fig2c_idea1_rates.tsv')
gc = gs0[4].subgridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.42)
ax = fig.add_subplot(gc[0])
offc = dict(zip(MAIN, [-0.27, -0.09, 0.09, 0.27]))
groups = [('cc', 'se50', 'CC, 50 bp'), ('cc', 'se101', 'CC, 101 bp'), ('do', 'se50', 'DO, 50 bp'), ('do', 'se101', 'DO, 101 bp')]
for gi, (pop, lay, _) in enumerate(groups):
    for c in MAIN[::-1]:
        v = dor[(dor.population == pop) & (dor.layout == lay) & (dor.method == c)].fpr.iloc[0]
        ax.scatter(gi + offc[c], 100 * v, **mstyle(c, 13))
    t = dor[(dor.population == pop) & (dor.layout == lay) & (dor.method == 'truth')]
    if len(t):
        ax.plot([gi - 0.4, gi + 0.4], [100 * t.fpr.iloc[0]] * 2, color='black', lw=0.6, ls=(0, (1, 1)), zorder=1)
ax.text(0.6, 100 * dor[(dor.population == 'cc') & (dor.method == 'truth')].fpr.iloc[0] * 0.62,
        'dotted: test on the\nexpected counts', fontsize=FS_ANN, ha='center', va='top')
ax.set_yscale('log'); ax.set_ylim(0.1, 40)
ax.set_yticks([0.1, 1, 10]); ax.set_yticklabels(['0.1%', '1%', '10%'])
ax.set_xticks(range(4)); ax.set_xticklabels([g[2].replace(', ', '\n') for g in groups])
ax.set_xlim(-0.55, 3.55)
ax.set_ylabel('False cis-eQTL rate\n(null genes called at BH 5%)')
ax.set_title('Planted cis-eQTL: null genes called')
panel_letter(ax, 'c', x=-0.27)

# forest: PLUCK / Salmon-pan ratios
ax = fig.add_subplot(gc[1])
rows = []
# DO uses the plain ratio (fig2c_do_fpr_plain.py): the scorer's 0.5% floor on Salmon-pan's rate binds in DO only
g1do = rd('fig2c_do_fpr_plain.tsv')
for pop, gate, lab in (('cc', 'G1a_fpr_cc', 'CC'), ('do', 'G1a_fpr_do', 'DO')):
    for lay in ('se50', 'se101'):
        src = g1do if pop == 'do' else g1
        r = src[(src.layout == lay) & (src.gate == gate)].iloc[0]
        rows.append((f'{lab}, {lay[2:]} bp', r.ratio, r.ci_lo, r.ci_hi))
pw = rd('fig2c_idea1_power.tsv')
for pop in ('cc', 'do'):
    for lay in ('se50', 'se101'):
        r = pw[(pw.population == pop) & (pw.layout == lay)].iloc[0]
        rows.append((f'{pop.upper()}, {lay[2:]} bp', r.ratio, np.nan, np.nan))
yy = np.array([8.4, 7.4, 6.4, 5.4, 3.0, 2.0, 1.0, 0.0])
ax.axvline(1, color='black', lw=0.5)
for y, (lab, v, lo, hi) in zip(yy, rows):
    if np.isfinite(lo):
        ax.plot([lo, hi], [y, y], color=P['hybrid'][1], lw=0.9, zorder=3)
    ax.scatter([v], [y], **mstyle('hybrid', 13))
r50 = rows[0]
ax.text(r50[2] - 0.03, yy[0], f'{r50[1]:.2f} [{r50[2]:.2f}, {r50[3]:.2f}]', fontsize=FS_ANN, va='center', ha='right')
ax.set_yticks(yy); ax.set_yticklabels([r[0] for r in rows])
ax.set_xlim(0, 1.75); ax.set_ylim(-0.7, 9.9)
ax.axhline(4.3, color=GREY, lw=0.4)
ax.text(0.03, 9.4, 'False-cis rate', fontsize=FS_ANN, va='center')
ax.text(0.03, 3.75, 'Detections', fontsize=FS_ANN, va='center')
ax.set_xlabel('PLUCK / Salmon-pan')
ax.set_title('Ratios to Salmon-pan')

# ---------------------------------------------------------------- panel d: isoform switches that follow genotype (CC, 50 bp)
# data: fig2d_g5_fixed_unit.tsv (PLUCK, branch, Salmon and truth fixed-unit rows);
#       fig2c_idea1_gates.tsv G5 rows (Salmon-pan in the matched unit)
g5 = rd('fig2d_g5_fixed_unit.tsv'); g5 = g5[g5.layout == 'se50']
cr = g5[g5.analysis.str.startswith('call rate')]
get = lambda unit, code: cr[(cr.unit == unit) & (cr.tool_code == code)].iloc[0]
gd = gs0[5].subgridspec(1, 2, width_ratios=[1.15, 0.75], wspace=0.62)
ax = fig.add_subplot(gd[0])
vf, emg = get('reads', 'vf'), get('reads', 'emg')
ax.plot([100 * vf.swM_call_rate, 100 * emg.swM_call_rate], [100 * vf.swR_call_rate, 100 * emg.swR_call_rate],
        color=P['hybrid'][1], lw=0.6, ls=':', zorder=2)
# label positions chosen by hand; leader lines for the crowded corner near (33, 0)
pts = [('tpm', 'salmon_ref', 'salmon_ref', 'Salmon-ref, TPM', (6, 36.5), 'left'),
       ('tpm', 'salmon_pan', 'salmon_pan', 'Salmon-pan, TPM', (6, 32.5), 'left'),
       ('reads', 'vf', 'vf', 'gene probes alone', (6, 28.5), 'left'),
       ('reads', 'hybrid', 'hybrid', 'PLUCK', (24, 3.5), 'right'),
       ('reads', 'emg', 'emg', 'EM alone', (44, 11), 'right'),
       ('reads', 'salmon_ref', 'salmon_ref', 'Salmon-ref, reads', (44, 23), 'right'),
       ('reads', 'salmon_pan', 'salmon_pan', 'Salmon-pan, reads', (44, 19), 'right')]
for unit, code, sty, lab, (tx, ty), ha in pts:
    r = get(unit, code); x, y = 100 * r.swM_call_rate, 100 * r.swR_call_rate
    ax.scatter([x], [y], **mstyle(sty, 16))
    ax.annotate(lab, xy=(x, y), xytext=(tx, ty), fontsize=FS_ANN, ha=ha, va='center',
                color=P[sty][1] if sty in ('hybrid', 'vf', 'emg') else 'black', fontweight='bold' if sty == 'hybrid' else 'normal',
                arrowprops=dict(arrowstyle='-', color='#AAAAAA', lw=0.35, shrinkA=1, shrinkB=2))
t = get('reads', 'truth')
ax.scatter([100 * t.swM_call_rate], [100 * t.swR_call_rate], marker='*', s=26, facecolors='black', edgecolors='black', lw=0.3, zorder=6)
ax.annotate('truth, reads', xy=(100 * t.swM_call_rate, 100 * t.swR_call_rate), xytext=(44, 15), fontsize=FS_ANN, ha='right',
            va='center', arrowprops=dict(arrowstyle='-', color='#AAAAAA', lw=0.35, shrinkA=1, shrinkB=2))
gm = g1[(g1.layout == 'se50') & (g1.gate == 'G5_swM')].iloc[0]
gr = g1[(g1.layout == 'se50') & (g1.gate == 'G5_swR')].iloc[0]
ax.scatter([100 * gm.salmon_pan], [100 * gr.salmon_pan], marker='s', s=10, facecolors=P['salmon_pan'][1], edgecolors='white',
           lw=0.3, zorder=5)
ax.annotate('Salmon-pan, matched\nunit (oracle)', xy=(100 * gm.salmon_pan, 100 * gr.salmon_pan), xytext=(6, 8.5),
            fontsize=FS_ANN, ha='left', va='center', arrowprops=dict(arrowstyle='-', color='#888888', lw=0.45, shrinkA=1, shrinkB=1,
                                                                    relpos=(0, 0.5)), zorder=7)
ax.text(44, 49.5, f'PLUCK: {100 * gm.hybrid:.1f}% (copy-pres.),\n{100 * gr.hybrid:.1f}% (read-pres.)', fontsize=FS_ANN, ha='right', va='top', fontweight='bold', color=P['hybrid'][1])
ax.text(0, 49.5, 'Salmon and truth in\nfixed units', fontsize=FS_ANN, ha='left', va='top', color=GREY)
ax.set_xlim(-1.5, 45); ax.set_ylim(-1.5, 50); ax.set_yticks([0, 10, 20, 30, 40])
ax.set_xlabel('Copy-preserving switch genes called (%)'); ax.set_ylabel('Read-preserving switch genes called (%)')
ax.set_title('Isoform switches, CC 50 bp')
panel_letter(ax, 'd', x=-0.2)

# FDP with switch genes in the BH family, and after the ratio-test exclusion
ax = fig.add_subplot(gd[1])
bh = g5[g5.analysis.str.startswith('BH 5%')]; ex = g5[g5.analysis.str.startswith('after excluding')]
fd = [('PLUCK', 'reads', 'hybrid', 'hybrid'), ('Salmon-pan, reads', 'reads', 'salmon_pan', 'salmon_pan'),
      ('Salmon-pan, TPM', 'tpm', 'salmon_pan', 'salmon_pan'), ('Salmon-ref, reads', 'reads', 'salmon_ref', 'salmon_ref'),
      ('truth, reads', 'reads', 'truth', None), ('Salmon-pan, matched\nunit (oracle)', 'oracle', 'salmon_pan', 'salmon_pan')]
yy = np.arange(len(fd))[::-1]
ax.axvline(5, color=GREY, lw=0.6, ls=(0, (3, 2)))
for y, (lab, unit, code, sty) in zip(yy, fd):
    v = 100 * bh[(bh.unit == unit) & (bh.tool_code == code)].FDP_null_plus_switch.iloc[0]
    e = ex[(ex.unit == unit) & (ex.tool_code == code)]
    if len(e):
        w = 100 * e.FDP_null_plus_switch.iloc[0]
        ax.annotate('', xy=(w, y), xytext=(v, y), arrowprops=dict(arrowstyle='-|>', color=P[sty][1], lw=0.6, mutation_scale=5,
                                                                    shrinkA=3, shrinkB=0))
        ax.scatter([w], [y], marker='|', s=20, c=P[sty][1], lw=0.9, zorder=4)
    if sty:
        ax.scatter([v], [y], **mstyle(sty, 15))
    else:
        ax.scatter([v], [y], marker='*', s=26, c='black', zorder=4)
ax.set_yticks(yy); ax.set_yticklabels([f[0] for f in fd])
ax.set_xlim(0, 25); ax.set_ylim(-0.6, len(fd) + 0.4)
ax.scatter([8.0], [len(fd) - 0.15], marker='|', s=20, c=GREY, lw=0.9)
ax.text(9.0, len(fd) - 0.15, 'after ratio screen', fontsize=FS_ANN, va='center', color=GREY)
ax.set_xticks([0, 5, 10, 15, 20, 25]); ax.set_xticklabels(['0', '5', '10', '15', '20', '25'])
ax.set_xlabel('FDP at BH 5%, switch genes\nin the tested family (%)')
ax.set_title('FDP with switch genes, CC 50 bp')

fig.canvas.draw()
r = fig.canvas.get_renderer()
for a in fig.axes:
    bb = a.get_tightbbox(r)
    if bb.x0 < 0 or bb.x1 > fig.bbox.x1:
        print('outside', a.get_title(), round(bb.x0), round(bb.x1), round(fig.bbox.x1))
save(fig, 'Fig2')
print('saved Fig2')
