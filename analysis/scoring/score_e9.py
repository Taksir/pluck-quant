#!/usr/bin/env python3
"""Score per-gene probe estimates against truth counts, per tissue, with fixed thresholds.

    python score_e9.py RESULTS_DIR [--log LOG.jsonl] [--extra extra.tsv] [--timing timing.tsv] [--out OUTDIR]

RESULTS_DIR holds per_gene.tsv (required) and summary.tsv (optional) written by e9.py score.
--log     log.jsonl from the probe run (default: RESULTS_DIR/../log.jsonl if present); gives probe query times.
--extra   other methods' gene estimates, long format: srr, gene, method, estimate (e.g. salmon_ref, salmon_pan).
--timing  per-sample wall time of other methods: srr, method, seconds, threads
          (method salmon_pan; optionally probe_stream = streaming FASTQ probe counter).
Writes OUTDIR/report.md, metrics.json and delta_per_gene_<tissue>.tsv (default OUTDIR = RESULTS_DIR/score_out).
The thresholds are constants at the top of this file.
"""
import argparse
import json
import os
from collections import defaultdict

import numpy as np
from scipy.stats import mannwhitneyu, spearmanr

# ---- thresholds ----
MIN_TRUTH = 50          # reads; per sample-gene
MIN_HIT, MIN_B6 = 5, 3  # samples per gene in the non-B6-hit and B6 groups
MIN_SPEARMAN_N = 10     # samples per gene for within-gene Spearman
VF_BAND, REF_MAX = 0.03, -0.08   # C1: |median delta| of vf <= 0.03 and median delta of ref <= -0.08
SPUR_EFFECT, SPUR_P, SPUR_RATIO = 0.10, 0.05, 2.0  # C2: spurious gene = |delta| > 0.10 and p < 0.05; ref rate >= 2x vf rate
MATCH_RATE = 1.25       # C3: competitor "matches" vf if its spurious rate <= 1.25x vf's
SPEED_RATIO = 5.0       # C3: minimum ratio of salmon_pan time to probe FASTQ-to-counts time
CTRL_BAND = 0.03        # negative control on B6 sample-genes
MIN_SAMPLES, MIN_GENES = 30, 200   # minimum samples and genes to evaluate C1-C3
COHORT_EXTRA = 30       # an extra method needs this many samples for cohort-scale C3
WILD, DOMESTIC = set('FGH'), set('ACDE')  # CAST, PWK, WSB vs A/J, 129, NOD, NZO
NBOOT, SEED = 2000, 20261004
PROBE_METHODS = ('ref', 'vf')


def read_tsv(path):
    with open(path) as f:
        hdr = f.readline().rstrip('\n').split('\t')
        return [dict(zip(hdr, l.rstrip('\n').split('\t'))) for l in f if l.strip()]


def load(a):
    """Return {tissue: {(srr, gene): rec}} with rec = dict(strain, stratum, founders, hit, truth, est{method})."""
    data = defaultdict(dict)
    for r in read_tsv(os.path.join(a.results, 'per_gene.tsv')):
        data[r['tissue']][(r['srr'], r['gene'])] = dict(
            strain=r['strain'], stratum=r['stratum'], founders=r['founders'], hit=r['hit'] == 'True',
            truth=float(r['truth']), est={'ref': float(r['est_ref']), 'vf': float(r['est_vf'])})
    extra = set()
    if a.extra:
        for r in read_tsv(a.extra):
            for t in data.values():
                rec = t.get((r['srr'], r['gene']))
                if rec is not None:
                    rec['est'][r['method']] = float(r['estimate'])
                    extra.add(r['method'])
    return data, sorted(extra)


def calibrate(recs, methods):
    """Per sample and method: C = median(truth/est) over B6 sample-genes with truth >= MIN_TRUTH (as e9.py)."""
    by = defaultdict(lambda: defaultdict(list))
    for (srr, _), r in recs.items():
        if r['stratum'] == 'B6' and r['truth'] >= MIN_TRUTH:
            for m in methods:
                e = r['est'].get(m)
                if e is not None and e > 0:
                    by[(srr, m)]['x'].append(r['truth'] / e)
    return {k: float(np.median(v['x'])) for k, v in by.items()}


def boot_median(x, rng):
    x = np.asarray(x)
    if len(x) < 2:
        return [float('nan')] * 2
    b = np.median(x[rng.integers(0, len(x), (NBOOT, len(x)))], axis=1)
    return [round(float(v), 4) for v in np.percentile(b, [2.5, 97.5])]


def score_tissue(recs, methods, rng, samples=None):
    """Metrics for each method, restricted to `samples` (None = all) and to sample-genes every method covers."""
    C = calibrate(recs, methods)
    genes = defaultdict(list)
    for (srr, g), r in recs.items():
        if samples is not None and srr not in samples:
            continue
        if r['truth'] < MIN_TRUTH or r['stratum'] not in ('B6', 'nonB6'):
            continue
        if any(r['est'].get(m) is None or (srr, m) not in C for m in methods):
            continue
        lr = {m: float(np.log2((C[(srr, m)] * r['est'][m] + 1) / (r['truth'] + 1))) for m in methods}
        lt = float(np.log2(r['truth'] + 1))
        le = {m: float(np.log2(C[(srr, m)] * r['est'][m] + 1)) for m in methods}
        grp = 'B6' if r['stratum'] == 'B6' else ('hit' if r['hit'] else 'nohit')
        f = set(r['founders'])
        fam = 'wild' if f and f <= WILD else ('domestic' if f and f <= DOMESTIC else None)
        genes[g].append((grp, fam, lr, lt, le))
    nsamp = len({s for (s, _) in recs if samples is None or s in samples})
    out = {m: dict(delta=[], spur=[], fam={'wild': [], 'domestic': []}, rho=[], pooled={'B6': [], 'hit': []})
           for m in methods}
    rows = []
    for g, ps in genes.items():
        for m in methods:
            for grp, _, lr, _, _ in ps:
                if grp in ('B6', 'hit'):
                    out[m]['pooled'][grp].append(lr[m])
        if len(ps) >= MIN_SPEARMAN_N:
            lt = np.array([p[3] for p in ps])
            for m in methods:
                le = np.array([p[4][m] for p in ps])
                if lt.std() > 0 and le.std() > 0:
                    out[m]['rho'].append(float(spearmanr(le, lt).statistic))
        b6 = [p for p in ps if p[0] == 'B6']
        hit = [p for p in ps if p[0] == 'hit']
        if len(b6) < MIN_B6 or len(hit) < MIN_HIT:
            continue
        row = [g, len(hit), len(b6)]
        for m in methods:
            xb, xh = [p[2][m] for p in b6], [p[2][m] for p in hit]
            d = float(np.median(xh) - np.median(xb))
            p = float(mannwhitneyu(xh, xb, alternative='two-sided').pvalue) if len(set(xh + xb)) > 1 else 1.0
            out[m]['delta'].append(d)
            out[m]['spur'].append(abs(d) > SPUR_EFFECT and p < SPUR_P)
            row += [round(d, 4), round(p, 4)]
            for fam in ('wild', 'domestic'):
                xf = [p_[2][m] for p_ in hit if p_[1] == fam]
                if len(xf) >= 3:
                    out[m]['fam'][fam].append(float(np.median(xf) - np.median(xb)))
        rows.append(row)
    met = {}
    for m in methods:
        o = out[m]
        met[m] = dict(
            n_samples=nsamp, n_genes=len(o['delta']),
            median_delta=round(float(np.median(o['delta'])), 4) if o['delta'] else None,
            delta_ci95=boot_median(o['delta'], rng),
            spurious_rate=round(float(np.mean(o['spur'])), 4) if o['spur'] else None,
            pooled_hit_med_log2=round(float(np.median(o['pooled']['hit'])), 4) if o['pooled']['hit'] else None,
            pooled_B6_med_log2=round(float(np.median(o['pooled']['B6'])), 4) if o['pooled']['B6'] else None,
            within_gene_spearman=round(float(np.median(o['rho'])), 4) if o['rho'] else None,
            n_genes_spearman=len(o['rho']),
            delta_wild=round(float(np.median(o['fam']['wild'])), 4) if o['fam']['wild'] else None,
            delta_domestic=round(float(np.median(o['fam']['domestic'])), 4) if o['fam']['domestic'] else None,
            C_median=round(float(np.median([v for (s, mm), v in C.items() if mm == m and (samples is None or s in samples)])), 4))
    return met, rows


def probe_seconds(log_path):
    if not log_path or not os.path.exists(log_path):
        return None
    q, b = [], []
    for line in open(log_path):
        r = json.loads(line)
        if r.get('step') == 'sample':
            q.append(r['query_s']); b.append(r['build_s'] + r['query_s'])
    return dict(n=len(q), query_s=float(np.median(q)), build_plus_query_s=float(np.median(b))) if q else None


def other_seconds(path, method):
    if not path:
        return None
    x = [(float(r['seconds']), r.get('threads', '')) for r in read_tsv(path) if r['method'] == method]
    return dict(n=len(x), seconds=float(np.median([v for v, _ in x])), threads=sorted({t for _, t in x})) if x else None


def gate(met, cmp_met, speed):
    r, v = met['ref'], met['vf']
    res = {}
    ok = r['n_samples'] >= MIN_SAMPLES and v['n_genes'] >= MIN_GENES
    ctrl = all(met[m]['pooled_B6_med_log2'] is not None and abs(met[m]['pooled_B6_med_log2']) <= CTRL_BAND
               for m in PROBE_METHODS)
    res['valid'] = dict(pass_=ok, n_samples=r['n_samples'], n_genes=v['n_genes'])
    res['control'] = dict(pass_=ctrl, ref=r['pooled_B6_med_log2'], vf=v['pooled_B6_med_log2'])
    if not ok:
        return res, 'too few samples to score'
    c1 = abs(v['median_delta']) <= VF_BAND and r['median_delta'] <= REF_MAX
    res['C1'] = dict(pass_=c1, vf=v['median_delta'], ref=r['median_delta'])
    c2 = r['spurious_rate'] > 0 and r['spurious_rate'] >= SPUR_RATIO * v['spurious_rate']
    res['C2'] = dict(pass_=c2, ref=r['spurious_rate'], vf=v['spurious_rate'],
                     ratio=round(r['spurious_rate'] / v['spurious_rate'], 2) if v['spurious_rate'] else None)
    c3 = None
    if cmp_met is not None:
        p, vv, scale = cmp_met['salmon_pan'], cmp_met['vf'], cmp_met['_scale']
        if scale == 'cohort':
            match = (p['median_delta'] is not None and abs(p['median_delta']) <= VF_BAND
                     and p['spurious_rate'] <= MATCH_RATE * vv['spurious_rate'])
        else:  # pilot scale: pooled hit-gene bias only
            match = abs(p['pooled_hit_med_log2']) <= max(VF_BAND, abs(vv['pooled_hit_med_log2']) + 0.01)
        ratio = speed['ratio_gate_fastq'] if speed else None
        c3 = (not match) or (ratio is not None and ratio >= SPEED_RATIO)
        res['C3'] = dict(pass_=c3, scale=scale, pan_matches_vf=match, speed=speed,
                         provisional=scale != 'cohort')
    else:
        res['C3'] = dict(pass_=None, note='no salmon_pan estimates supplied')
    verdict = f'conditions met: C1={c1}, C2={c2}, C3={c3}' + (' (C3 at pilot scale)' if c3 and res['C3']['provisional'] else '')
    if not ctrl:
        verdict += ' (the B6 control check did not hold; read the metrics with care)'
    return res, verdict


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('results')
    ap.add_argument('--log')
    ap.add_argument('--extra')
    ap.add_argument('--timing')
    ap.add_argument('--out')
    a = ap.parse_args()
    log_path = a.log or os.path.join(a.results, '..', 'log.jsonl')
    out_dir = a.out or os.path.join(a.results, 'score_out')
    os.makedirs(out_dir, exist_ok=True)
    data, extra = load(a)
    psec = probe_seconds(log_path)
    report, allmet = ['# Scoring report', ''], {}
    for tissue, recs in sorted(data.items()):
        rng = np.random.default_rng(SEED)
        met, rows = score_tissue(recs, list(PROBE_METHODS), rng)
        cmp_met = None
        if 'salmon_pan' in extra:  # C3 on the samples and sample-genes every method covers
            cov = {s for (s, _), r in recs.items() if 'salmon_pan' in r['est']}
            ms = list(PROBE_METHODS) + [m for m in extra if m in ('salmon_pan', 'salmon_ref')]
            cmp_met, _ = score_tissue(recs, ms, np.random.default_rng(SEED), samples=cov)
            cmp_met['_scale'] = 'cohort' if len(cov) >= COHORT_EXTRA else 'pilot'
        speed = None
        osec, stream = other_seconds(a.timing, 'salmon_pan'), other_seconds(a.timing, 'probe_stream')
        if psec and osec:  # speed ratio: FASTQ-to-counts time for both methods; the amortized ratio is only reported
            fq = min(psec['build_plus_query_s'], stream['seconds'] if stream else float('inf'))
            speed = dict(probe_query_s=psec['query_s'], probe_build_plus_query_s=psec['build_plus_query_s'],
                         probe_stream_s=stream and stream['seconds'], salmon_pan_s=osec['seconds'],
                         salmon_threads=osec['threads'], ratio_gate_fastq=round(osec['seconds'] / fq, 2),
                         ratio_amortized=round(osec['seconds'] / psec['query_s'], 2))
        res, verdict = gate(met, cmp_met, speed)
        allmet[tissue] = dict(metrics=met, compare=cmp_met, gate=res, verdict=verdict, probe_timing=psec)
        with open(os.path.join(out_dir, f'delta_per_gene_{tissue}.tsv'), 'w') as f:
            f.write('gene\tn_hit\tn_B6\t' + '\t'.join(f'delta_{m}\tp_{m}' for m in PROBE_METHODS) + '\n')
            f.writelines('\t'.join(map(str, r)) + '\n' for r in rows)
        report += [f'## {tissue}: **{verdict}**', '',
                   '| method | samples | genes | median Δ (95% CI) | spurious rate | pooled hit log2 | pooled B6 log2 '
                   '| within-gene Spearman | Δ wild | Δ domestic | C |', '|' + '---|' * 11]
        for src in (met, cmp_met or {}):
            for m, x in src.items():
                if m.startswith('_') or (src is cmp_met and m in PROBE_METHODS and cmp_met['_scale'] == 'cohort'
                                         and x['n_samples'] == met[m]['n_samples']):
                    continue
                tag = m if src is met else f"{m} (C3 subset, {cmp_met['_scale']})"
                report.append(f"| {tag} | {x['n_samples']} | {x['n_genes']} | {x['median_delta']} {x['delta_ci95']} "
                              f"| {x['spurious_rate']} | {x['pooled_hit_med_log2']} | {x['pooled_B6_med_log2']} "
                              f"| {x['within_gene_spearman']} (n={x['n_genes_spearman']}) | {x['delta_wild']} "
                              f"| {x['delta_domestic']} | {x['C_median']} |")
        report += ['', '| check | pass | detail |', '|---|---|---|']
        for k, v in res.items():
            report.append(f"| {k} | {v.get('pass_')} | {json.dumps({kk: vv for kk, vv in v.items() if kk != 'pass_'})} |")
        report.append('')
    if os.path.exists(os.path.join(a.results, 'summary.tsv')):  # secondary: across-gene Spearman, e9.py's own C
        s = read_tsv(os.path.join(a.results, 'summary.tsv'))
        report += ['## Secondary (from summary.tsv, per-sample medians)', '', '| set | stratum | samples | median Spearman | median log2 | median C |', '|---|---|---|---|---|---|']
        for st in ('all', 'B6', 'nonB6_hit'):
            for m in PROBE_METHODS:
                x = [r for r in s if r['set'] == m and r['stratum'] == st and r['spearman'] != 'nan']
                if x:
                    md = [float(r['med_log2']) for r in x if r['med_log2'] not in ('None', '')]
                    report.append(f"| {m} | {st} | {len(x)} | {np.median([float(r['spearman']) for r in x]):.4f} "
                                  f"| {np.median(md) if md else 'NA'} | {np.median([float(r['C']) for r in x]):.3f} |")
    with open(os.path.join(out_dir, 'report.md'), 'w') as f:
        f.write('\n'.join(report) + '\n')
    with open(os.path.join(out_dir, 'metrics.json'), 'w') as f:
        json.dump(allmet, f, indent=1, default=str)
    print('\n'.join(report))


if __name__ == '__main__':
    main()
