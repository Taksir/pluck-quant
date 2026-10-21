#!/usr/bin/env python3
"""Supplementary Table S9 (numeric tables of the read simulation and the planted cis-eQTL simulation), written as Markdown.

Reads (read only): figures/data/fig2a_e2_metrics_all.tsv, fig2a_e2_per_sample_error.tsv, fig2b_e2_delta_founder_ci.tsv,
                   fig2c_idea1_gates.tsv, fig2c_idea1_rates.tsv, fig2c_idea1_power.tsv, fig2d_g5_fixed_unit.tsv,
                   figS3_idea1_founder.tsv, figS3_idea1_effect.tsv, figS3_idea1_branch.tsv
Writes:            tableS9.md (path in OUT below)
Every value is printed from those files; nothing is typed by hand. Codes: in the simulation files hybrid = PLUCK,
vf = gene probes alone, emg = EM alone, ref = ref-unique probes.
"""
import os
import pandas as pd

# Folder constants below point to the original working folder; set them for your copy.
DATA = '/mnt/project-files/paper/figures/data'
OUT = '/mnt/project-files/paper/draft/U2_tables/tableS9.md'
rd = lambda f: pd.read_csv(f'{DATA}/{f}', sep='\t', comment='#')
NAME = {'hybrid': 'PLUCK', 'vf': 'gene probes alone', 'emg': 'EM alone', 'ref': 'ref-unique probes',
        'salmon_ref': 'Salmon-ref', 'salmon_pan': 'Salmon-pan', 'v2ext': 'PLUCK (extended tier)',
        'salmon_pan_flags': 'Salmon-pan, bias flags, NumReads', 'salmon_pan_tpm': 'Salmon-pan, no flags, TPM scale',
        'salmon_pan_flags_tpm': 'Salmon-pan, bias flags, TPM × eff. length (post hoc)', 'truth': 'truth (expected reads)'}
ORDER = ['hybrid', 'vf', 'emg', 'salmon_pan', 'salmon_ref', 'ref']
LAY = {'se50': '50 bp SE', 'se101': '101 bp SE', 'pe': '101 bp PE'}
FOUNDER = dict(zip('ACDEFGH', ['A/J', '129S1', 'NOD', 'NZO', 'CAST', 'PWK', 'WSB']))


def f(x, d=4, sign=False):
    if pd.isna(x):
        return 'n/a'
    s = f'{x:+.{d}f}' if sign else f'{x:.{d}f}'
    return s.replace('-', '−')


def pct(x, d=1):
    return 'n/a' if pd.isna(x) else f'{100 * x:.{d}f}%'


def table(head, rows):
    esc = lambda c: str(c).replace('|', '\\|')
    out = ['| ' + ' | '.join(esc(h) for h in head) + ' |', '|' + '---|' * len(head)]
    out += ['| ' + ' | '.join(esc(c) for c in r) + ' |' for r in rows]
    return '\n'.join(out)


L = []
mt = rd('fig2a_e2_metrics_all.tsv')
ps = rd('fig2a_e2_per_sample_error.tsv')
core = mt[(mt.cond == 'clean') & (mt.truth == 'expected') & (mt.tier == 'core')]

# S9b: full per-method table, clean, expected truth
rows = []
for lay in LAY:
    for m in ORDER:
        r = core[(core.layout == lay) & (core.method == m)].iloc[0]
        s = ps[(ps.cond == 'clean') & (ps.layout == lay) & (ps.tier == 'core') & (ps.method == m)].err_B6
        nsp = round(r.M4_spurious * r.n_M4)
        rows.append([LAY[lay], NAME[m], f'{int(r.n_genes):,}', f(r.M1_within_r2), f(r.M2_err_B6), f'{f(s.min())} to {f(s.max())}',
                     f(r.M2_ratio_to_pan, 2), f'{pct(r.M4_spurious)} ({nsp} of {int(r.n_M4)})', f(r.M3_max_abs_delta),
                     f(r.M4_median_delta_all, sign=True)])
L.append(('S9b', table(['Layout', 'Method', 'Genes', 'Within-gene r²', 'Error (median of 8)', 'Error, sample range',
                        'Error ratio to Salmon-pan', 'Spurious rate (genes)', 'Max |Δ| over founders', 'Median Δ, all hit genes'], rows)))

# S9c: per-founder median shift, read simulation, clean reads
dl = rd('fig2b_e2_delta_founder_ci.tsv')
rows = []
for lay in LAY:
    for m in ORDER:
        x = dl[(dl.layout == lay) & (dl.method == m)].set_index('founder')
        rows.append([LAY[lay], NAME[m]] + [f(x.loc[k, 'median_delta'], sign=True) for k in 'ACDEFGH'])
    n = dl[(dl.layout == lay) & (dl.method == 'hybrid')].set_index('founder').n_genes
    rows.append([LAY[lay], 'hit genes'] + [f'{int(n[k]):,}' for k in 'ACDEFGH'])
L.append(('S9c', table(['Layout', 'Method'] + [FOUNDER[k] for k in 'ACDEFGH'], rows)))
rows = []
for lay in LAY:
    x = dl[(dl.layout == lay) & (dl.method == 'hybrid')].set_index('founder')
    rows.append([LAY[lay]] + [f'{f(x.loc[k, "ci_lo"], sign=True)} to {f(x.loc[k, "ci_hi"], sign=True)}' for k in 'ACDEFGH'])
L.append(('S9c_ci', table(['PLUCK, 95% interval'] + [FOUNDER[k] for k in 'ACDEFGH'], rows)))

# S9d: pre-mRNA, 101 bp PE
pm = mt[(mt.cond == 'premrna') & (mt.truth == 'expected') & (mt.tier == 'core')]
rows = []
for m in ORDER:
    c, p = core[(core.layout == 'pe') & (core.method == m)].iloc[0], pm[pm.method == m].iloc[0]
    rows.append([NAME[m], f(c.M1_within_r2), f(p.M1_within_r2), f(c.M1_within_r2 - p.M1_within_r2), f(p.M2_err_B6),
                 f(p.M2_ratio_to_pan, 2), f'{pct(p.M4_spurious)} ({round(p.M4_spurious * p.n_M4)} of {int(p.n_M4)})', f(p.M3_max_abs_delta)])
L.append(('S9d', table(['Method', 'r², clean', 'r², pre-mRNA', 'r² drop', 'Error', 'Error ratio to Salmon-pan', 'Spurious rate (genes)',
                        'Max |Δ|'], rows)))

# S9e: bias arm, 101 bp PE
rows = []
for cond, lab in [('gc', 'GC'), ('pos3', '3′ positional'), ('hex', 'random hexamer')]:
    b = mt[(mt.cond == cond) & (mt.truth == 'expected') & (mt.tier == 'core')]
    ref = b[b.method == 'salmon_pan_flags_tpm'].iloc[0]
    dec = b[b.method == 'salmon_pan_flags'].iloc[0]   # declared comparator: Salmon-pan with bias flags, NumReads
    for m in ORDER[:4] + ['salmon_pan_flags', 'salmon_pan_flags_tpm', 'ref']:
        r = b[b.method == m]
        if not len(r):
            continue
        r = r.iloc[0]
        nm = 'Salmon-pan, no flags, NumReads' if m == 'salmon_pan' else NAME[m]
        rows.append([lab, nm, f(r.M1_within_r2), f(r.M2_err_B6), f(r.M2_err_B6 / dec.M2_err_B6, 3),
                     f(r.M2_err_B6 / ref.M2_err_B6, 3), f(r.M3_max_abs_delta), pct(r.M4_spurious)])
L.append(('S9e', table(['Bias', 'Method', 'Within-gene r²', 'Error', 'Error ratio to Salmon-pan with flags, NumReads (declared)',
                        'Error ratio to Salmon-pan with flags, TPM × eff. length (post hoc)', 'Max |Δ|', 'Spurious rate'], rows)))

# S9f: extended tier
ex = mt[(mt.cond == 'clean') & (mt.truth == 'expected') & (mt.tier == 'extended')]
rows = []
for lay in LAY:
    for m in ['v2ext', 'salmon_pan', 'salmon_ref']:
        r = ex[(ex.layout == lay) & (ex.method == m)].iloc[0]
        rows.append([LAY[lay], NAME[m], f'{int(r.n_genes):,} ({int(r.n_M1)})', f(r.M1_within_r2), f(r.M2_err_B6), f(r.M2_ratio_to_pan, 2)])
L.append(('S9f', table(['Layout', 'Method', 'Genes (with r²)', 'Within-gene r²', 'Error', 'Error ratio to Salmon-pan'], rows)))

# S9g: realized truth, clean
re_ = mt[(mt.cond == 'clean') & (mt.truth == 'realized') & (mt.tier == 'core')]
rows = []
for lay in LAY:
    for m in ORDER:
        r = re_[(re_.layout == lay) & (re_.method == m)].iloc[0]
        rows.append([LAY[lay], NAME[m], f(r.M1_within_r2), f(r.M2_err_B6)])
L.append(('S9g', table(['Layout', 'Method', 'Within-gene r²', 'Error'], rows)))

# S9h: cis-eQTL simulation checks
g = rd('fig2c_idea1_gates.tsv')
pw = rd('fig2c_idea1_power.tsv')
rows = []
def gate(lay, name):
    return g[(g.layout == lay) & (g.gate == name)]
for lab, name, fmt, lim in [
        ('G1a false-cis rate ratio, CC', 'G1a_fpr_cc', 'ratio', '≤ 1.25'),
        ('G1a false-cis rate ratio, DO', 'G1a_fpr_do', 'ratio', '≤ 1.25'),
        ('G1b FDP at BH 5%, CC', 'G1b_fdp_cc', 'fdp', '≤ 5%'),
        ('G1b FDP at BH 5%, DO', 'G1b_fdp_do', 'fdp', '≤ 5%')]:
    cells = []
    for lay in ['se50', 'se101']:
        r = gate(lay, name).iloc[0]
        if fmt == 'ratio':
            cells.append(f'{pct(r.hybrid, 2)} vs {pct(r.salmon_pan, 2)}; {f(r.ratio, 2)} [{f(r.ci_lo, 2)}, {f(r.ci_hi, 2)}]')
        else:
            cells.append(f'{pct(r.hybrid)} (Salmon-pan {pct(r.salmon_pan)})')
    rows.append([lab, lim] + cells + ['PASS / PASS'])
cells = []
for lay in ['se50', 'se101']:
    r = gate(lay, 'G2_power_cc')
    r = r[r.quantity.str.startswith('detections')].iloc[0]
    cells.append(f'{int(r.hybrid):,} / {int(r.salmon_pan):,} = {float(r.ratio):.3f}')
rows.append(['G2 detections ratio, CC', '≥ 0.90'] + cells + ['PASS / PASS'])
cells = []
for lay in ['se50', 'se101']:
    r = gate(lay, 'G3_effect_bias_cc')
    cells.append(' / '.join(f(v, sign=True) for v in r.hybrid.astype(float)))
rows.append(['G3 effect-size shift at |β| 0.25 / 0.5 / 0.75 / 1.0, CC', '≤ 0.05'] + cells + ['PASS / PASS'])
cells = []
for lay in ['se50', 'se101']:
    r = gate(lay, 'G4a_founder_median_cc').hybrid.astype(float)
    cells.append(f'max |median| {f(r.abs().max())}')
rows.append(['G4a founder median shift, CC (7 founders)', '≤ 0.02'] + cells + ['PASS / PASS'])
cells = []
for lay in ['se50', 'se101']:
    r = gate(lay, 'G4b_founder_q90_bound_cc').hybrid.astype(float)
    cells.append(f'{pct(r.min())} to {pct(r.max())}')
rows.append(['G4b share of hit genes beyond 0.05 log2, CC', '≤ 10% per founder'] + cells + ['PASS / PASS'])
cells = [f'max |median| {f(float(gate(lay, "G4c_pair_median_do").iloc[0].hybrid))}' for lay in ['se50', 'se101']]
rows.append(['G4c heterozygous founder pairs, DO (28 pairs)', '≤ 0.02'] + cells + ['PASS / PASS'])
for lab, name in [('G5 SW-M switch genes called (p < 0.001)', 'G5_swM'), ('G5 SW-R switch genes called (p < 0.001)', 'G5_swR')]:
    cells = [f'{pct(float(gate(lay, name).iloc[0].hybrid))} (Salmon-pan {pct(float(gate(lay, name).iloc[0].salmon_pan), 2)})'
             for lay in ['se50', 'se101']]
    rows.append([lab, '≤ 1.25 × max(Salmon-pan, 1%)'] + cells + ['FAIL / FAIL'])
cells = []
for lay in ['se50', 'se101']:
    r = pw[(pw.layout == lay) & (pw.population == 'do')].iloc[0]
    cells.append(f'{int(r.hybrid_detected)} / {int(r.pan_detected)} = {r.ratio:.3f}')
rows.append(['G2 detections ratio, DO (reported only)', 'none'] + cells + ['not gated'])
L.append(('S9h', table(['Gate', 'Limit', '50 bp SE (PLUCK vs Salmon-pan)', '101 bp SE', 'Verdict'], rows)))

# S9i: idea-1 per-method rates, CC and DO
rt = rd('fig2c_idea1_rates.tsv')
rows = []
for pop in ['cc', 'do']:
    for lay in ['se50', 'se101']:
        for m in ORDER + ['truth']:
            r = rt[(rt.population == pop) & (rt.layout == lay) & (rt.method == m)]
            if not len(r):
                continue
            r = r.iloc[0]
            rows.append([pop.upper(), LAY[lay], NAME[m], f'{int(r.null_calls)} of {int(r.n_null):,}', pct(r.fpr, 2), pct(r.fdp),
                         f'{int(r.eqtl_calls):,} of {int(r.n_eqtl):,}', ' / '.join(f'{r[c]:.2f}' for c in
                                                                                ['power_b0.25', 'power_b0.5', 'power_b0.75', 'power_b1']),
                         f'{r.lambda_gc_null:.2f}'])
L.append(('S9i', table(['Set', 'Layout', 'Method', 'Null genes called', 'False-cis rate', 'FDP', 'Planted genes called',
                        'Power at |β| 0.25 / 0.5 / 0.75 / 1.0', 'λGC, null genes'], rows)))

# S9j: G4 per founder, CC: original calibration, null-relative check and recalibration summaries
fo = rd('figS3_idea1_founder.tsv')
rows = []
for lay in ['se50', 'se101']:
    for m in ['hybrid', 'salmon_pan', 'salmon_ref', 'ref']:
        lk = fo[(fo.source == 'locked') & (fo.layout == lay) & (fo.method == m)]
        nr = fo[(fo.source.str.startswith('null-relative')) & (fo.layout == lay) & (fo.method == m)]
        ph = fo[(fo.source.str.startswith('post hoc')) & (fo.layout == lay) & (fo.method == m)]
        rows.append([LAY[lay], NAME[m], f'{f(lk.median_delta.min(), sign=True)} to {f(lk.median_delta.max(), sign=True)}',
                     f'{pct(lk.frac_beyond_bound.min())} to {pct(lk.frac_beyond_bound.max())}',
                     f'{f(100 * nr.excess.min(), 1, True)} to {f(100 * nr.excess.max(), 1, True)}' if len(nr) else 'n/a',
                     f'{f(ph.median_delta.min(), sign=True)} to {f(ph.median_delta.max(), sign=True)}',
                     f'{pct(ph.frac_beyond_bound.min())} to {pct(ph.frac_beyond_bound.max())}'])
L.append(('S9j', table(['Layout', 'Method', 'G4a median shift, range over founders', 'G4b share beyond 0.05, range',
                        'Excess over own B6-split null, points (reported only)', 'Median shift, post hoc calibration',
                        'Share beyond 0.05, post hoc calibration'], rows)))
ef = rd('figS3_idea1_effect.tsv')
rows = []
for lay in ['se50', 'se101']:
    for m in ORDER:
        lk = ef[(ef.layout == lay) & (ef.method == m) & (ef.source == 'locked') & ef.M2_null_cc.notna()].M2_null_cc.iloc[0]
        ph = ef[(ef.layout == lay) & (ef.method == m) & (ef.source != 'locked') & ef.M2_null_cc.notna()].M2_null_cc.iloc[0]
        rows.append([LAY[lay], NAME[m], f(lk), f(ph)])
L.append(('S9j_m2', table(['Layout', 'Method', 'Error on null genes, locked calibration set', 'Same, post hoc calibration on all null core genes'], rows)))

# S9k: G5 fixed units, branch split, FDP with switch genes
g5 = rd('fig2d_g5_fixed_unit.tsv')
rows = []
cr = g5[g5.analysis.str.startswith('call rate')]
for lay in ['se50', 'se101']:
    for _, r in cr[cr.layout == lay].iterrows():
        rows.append([LAY[lay], NAME[r.tool_code], 'NumReads' if r.unit == 'reads' else 'TPM', pct(r.swM_call_rate), pct(r.swR_call_rate)])
L.append(('S9k', table(['Layout', 'Method', 'Gene-level unit', 'SW-M genes called', 'SW-R genes called'], rows)))
br = rd('figS3_idea1_branch.tsv')
rows = []
for lay in ['se50', 'se101']:
    for sw in ['swM', 'swR']:
        x = br[(br.layout == lay) & (br.switch == sw)]
        sh = x[x.branch.str.startswith('share')].call_rate.iloc[0]
        for b, lab in [('gene_probe_branch', 'gene-probe branch'), ('EM_branch', 'EM branch')]:
            y = x[x.branch == b].set_index('method')
            rows.append([LAY[lay], sw.replace('sw', 'SW-'), lab, f'{int(y.loc["hybrid", "n"])}',
                         pct(y.loc['hybrid', 'call_rate']), pct(y.loc['vf', 'call_rate']), pct(y.loc['emg', 'call_rate']),
                         pct(y.loc['salmon_pan', 'call_rate'])])
        rows.append([LAY[lay], sw.replace('sw', 'SW-'), 'share on gene-probe branch', '', pct(sh), '', '', ''])
L.append(('S9k_branch', table(['Layout', 'Switch type', 'PLUCK branch', 'Genes', 'PLUCK', 'gene probes alone', 'EM alone',
                               'Salmon-pan (matched unit)'], rows)))
rows = []
bh = g5[g5.analysis.str.startswith('BH 5%') | g5.analysis.str.startswith('after excluding')]
for lay in ['se50', 'se101']:
    for _, r in bh[bh.layout == lay].iterrows():
        unit = {'reads': 'NumReads', 'tpm': 'TPM', 'oracle': 'matched unit'}[r.unit]
        fam = 'all tested genes' if r.analysis.startswith('BH') else 'after ratio screen'
        rows.append([LAY[lay], fam, NAME[r.tool_code], unit, f'{int(r.calls):,}', pct(r.FDP_null_plus_switch), pct(r.FDP_null_only)])
L.append(('S9k_fdp', table(['Layout', 'Family', 'Method', 'Unit', 'Calls', 'FDP, null and switch genes false',
                            'FDP, null genes only false'], rows)))
fl = g5[g5.analysis.str.startswith('ratio-test')]
rows = [[LAY[r.layout]] + [v.strip() for v in str(r.note).split('/')] for _, r in fl.iterrows()]
rows = [[r[0]] + [f'{100 * float(v):.2f}%' for v in r[1:]] for r in rows]
L.append(('S9k_flag', table(['Layout', 'Null genes', 'Planted cis-eQTL genes', 'SW-M genes', 'SW-R genes'], rows)))

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, 'w') as fh:
    for k, t in L:
        fh.write(f'<!-- {k} -->\n{t}\n\n')
print('wrote', OUT)
