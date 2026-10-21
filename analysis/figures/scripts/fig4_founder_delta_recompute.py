#!/usr/bin/env python3
"""Fig. 4a source data: median founder-haplotype shift per founder in the Munger founder livers (100 bp),
for PLUCK (scorer code vf), ref-unique probes, Salmon-ref, Salmon-pan and Kmerator+countTags, on one scoring set.

Why this exists: e5_munger/metrics.json has per-founder values for four methods only, and the scored k-mer-tool results
(kmer_tools/results/score/results.json) carry only the pooled wild/domestic split. This script recomputes the per-founder
split for Kmerator+countTags with the scoring code, unchanged:
  - scoring/score_e9.py, score_tissue, called as the Munger scoring does for per_founder_core
    (records of B6 samples plus the samples of one founder; the Munger non-hit calibration, set by fc_gather_score.set_calibration);
  - scoring/fc_gather_score.py: load_cohort, read_tool, covered_samples.
Scoring set = P_kmerator (ref, vf, salmon_ref, salmon_pan, kmerator), the set of acc_scores/results.json in which the
Kmerator numbers were scored. RNA-Skim per-sample tables are not in the shared folder, so RNA-Skim
appears in Fig. 4a only through its pooled domestic/wild values from kmer_tools/results/score/results.json (set K).

Checks written to the output header: the all-founder run must reproduce acc_scores P_kmerator Munger 'all' median deltas
exactly (4 decimals); the per-founder values of the four methods of the Munger scoring are compared with e5_munger/metrics.json.

Reads data only from revision/ in the original working folder; writes figures/data/fig4a_munger_founder_delta_frameP.tsv.
"""
import json, os, sys, time
import numpy as np

# Folder constants below point to the original working folder; set them for your copy.
REV = '/mnt/project-files/revision'
OUT = '/mnt/project-files/paper/figures/data/fig4a_munger_founder_delta_frameP.tsv'
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'scoring'))
import fc_gather_score as G  # noqa: E402

METHODS = ['ref', 'vf', 'salmon_ref', 'salmon_pan', 'kmerator']
LABEL = {'ref': 'ref-unique probes', 'vf': 'hybrid', 'salmon_ref': 'Salmon-ref', 'salmon_pan': 'Salmon-pan', 'kmerator': 'Kmerator+countTags'}
FOUNDERS = ['A', 'C', 'D', 'E', 'H', 'F', 'G']   # A/J 129S1 NOD NZO WSB CAST PWK (domestic, then wild-derived)
NAME = dict(A='A/J', C='129S1', D='NOD', E='NZO', F='CAST', G='PWK', H='WSB')


def main():
    t0 = time.time()
    G.load_scorer()
    S9 = G.S9
    c = G.load_cohort('munger', f'{REV}/final_comparison/frames'); G._CTX['munger'] = c
    c.V, c.missing = {}, {}
    c.V['kmerator'], c.missing['kmerator'] = G.read_tool(c, f'{REV}/final_comparison/results_ll/acc', 'kmerator')
    for (srr, g), r in c.recs.items():
        v = c.V['kmerator'][c.si[srr], c.gi[g]] if g in c.gi else np.nan
        if v == v: r['est']['kmerator'] = float(v)
    G.set_calibration('munger')
    cov = G.covered_samples(c, METHODS)
    print(f'records {len(c.recs)}, samples {len(c.srrs)}, covered {len(cov)}, kmerator missing {len(c.missing["kmerator"])}', flush=True)

    rows, notes = [], []
    # check 1: all founders together reproduces the scored P_kmerator set
    met, _ = S9.score_tissue(c.recs, METHODS, np.random.default_rng(S9.SEED), samples=cov)
    pub = json.load(open(f'{REV}/final_comparison/scores/acc_scores/results.json'))['cohorts']['munger']['frames']['P_kmerator']['all']['metrics']
    ok = all(met[m]['median_delta'] == pub[m]['median_delta'] and met[m]['spurious_rate'] == pub[m]['spurious_rate'] for m in METHODS)
    notes.append('check 1 (all founders vs acc_scores P_kmerator munger all, median_delta and spurious_rate): ' + ('match' if ok else 'MISMATCH') + ' '
                 + '; '.join(f'{m} {met[m]["median_delta"]}/{pub[m]["median_delta"]}' for m in METHODS))
    for m in METHODS:
        rows.append(('all', 'all non-B6 founders', m, LABEL[m], met[m]['median_delta'], json.dumps(met[m]['delta_ci95']), met[m]['n_genes'], met[m]['spurious_rate']))
    print(notes[-1], flush=True)

    e5 = json.load(open(f'{REV}/e5_munger/metrics.json'))['per_founder_core']
    diffs = []
    for L in FOUNDERS:
        sub = {k: x for k, x in c.recs.items() if x['stratum'] == 'B6' or x['founders'] == L}
        m_, _ = S9.score_tissue(sub, METHODS, np.random.default_rng(S9.SEED), samples=cov)
        for m in METHODS:
            rows.append((L, NAME[L], m, LABEL[m], m_[m]['median_delta'], json.dumps(m_[m]['delta_ci95']), m_[m]['n_genes'], m_[m]['spurious_rate']))
            if m in e5[L]:
                diffs.append(abs(m_[m]['median_delta'] - e5[L][m]['delta']))
        print(L, {m: m_[m]['median_delta'] for m in METHODS}, f'{time.time() - t0:.0f}s', flush=True)
    notes.append(f'check 2 (per-founder values of ref, vf, salmon_ref, salmon_pan vs e5_munger/metrics.json per_founder_core): max |difference| {max(diffs):.4f} over {len(diffs)} values')
    print(notes[-1], flush=True)

    with open(OUT, 'w') as f:
        f.write('# file: fig4a_munger_founder_delta_frameP.tsv\n')
        f.write('# source: recomputed by fig4_founder_delta_recompute.py from final_comparison/frames/munger and results_ll/acc/munger/*/kmerator.genes.tsv.gz with score_e9.py and fc_gather_score.py, unchanged\n')
        f.write(f'# written: {time.strftime("%Y-%m-%d", time.gmtime())} by fig4_founder_delta_recompute.py; derived (per-founder split of set P_kmerator, Munger non-hit calibration, as in per_founder_core)\n')
        for n in notes: f.write(f'# {n}\n')
        f.write('# note: founder letters A C D E F G H = A/J 129S1 NOD NZO CAST PWK WSB; tool "hybrid" is plotted as PLUCK; spurious_rate at the nominal 5% test level\n')
        f.write('founder\tfounder_name\tmethod\ttool\tmedian_delta\tdelta_ci95\tn_genes\tspurious_rate\n')
        for r in rows: f.write('\t'.join(str(x) for x in r) + '\n')
    print('wrote', OUT, f'{time.time() - t0:.0f}s')


if __name__ == '__main__':
    main()
