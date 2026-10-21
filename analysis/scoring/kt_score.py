#!/usr/bin/env python3
"""kt_score.py : scoring wrapper for the RNA-Skim comparison.

The scorer (score_e9.py) and the gather script (fc_gather_score.py) are imported and called unchanged (the gather script's sha256 is
recorded in the output).  Nothing here defines a metric or a threshold; it only registers the method `rnaskim`, builds the scoring sets
below, applies the claim rules (RK1-RK3) and the integrity checks (I1-I4), and summarises the timing blocks (RK4).

Sub-commands
  score   --frames FRAMES --rs-root DIR --old-root DIR --out DIR [--cohorts ...] [--procs 2] [--skip-repro-gate] [--resume]
          FRAMES   = final_comparison/frames (plain input files, frames.sha256 checked)
          RS-ROOT  = RNA-Skim outputs:  <RS-ROOT>/<cohort>/<SRR>/rnaskim.genes.tsv.gz (+ integrity.json, DONE or FAILED)   [= kmer_tools/out/acc]
          OLD-ROOT = earlier comparison outputs (for the Kmerator arm of set 'K'): <OLD-ROOT>/<cohort>/<SRR>/kmerator.genes.tsv.gz  [= results_ll/acc]
  timing  --raw out/timing_raw.tsv --out DIR      RK4: paired wall-time ratios RNA-Skim / other tool per sample at 4 threads (blocks are paired by sample, threads, rep)

Scoring sets: P_rnaskim = [ref, vf, salmon_ref, salmon_pan, rnaskim]; K = [ref, vf, salmon_ref, salmon_pan, kmerator, rnaskim] (descriptive; dropped when the
scorer's validity rule fails on the all-core-gene rows); P_rnaskim_z = P_rnaskim with 0 instead of missing for genes without a quantifiable transcript
(sensitivity, descriptive); each for subsets all / unflagged / flagged.  The four existing methods are scored first and compared with the published
values; on a mismatch nothing else is scored unless --skip-repro-gate is given (fake-truth plumbing test only).
"""
import argparse, collections, json, os, sys, time
import numpy as np
import pandas as pd

# Folder constants below point to the original working folder; set them for your copy.
GATHER_DIR = os.environ.get('KT_GATHER_DIR', os.path.dirname(os.path.abspath(__file__)))   # KT_* environment variables override paths only
sys.path.insert(0, GATHER_DIR)
import fc_gather_score as G                                    # noqa: E402
for _env, _attr in (('KT_GENES_JSON', 'GENES_JSON'), ('KT_FLAG_TABLE', 'FLAG_TABLE')):
    if os.environ.get(_env): setattr(G, _attr, os.environ[_env])

PRIOR_RESULTS = os.environ.get('KT_PRIOR_RESULTS', '/mnt/project-files/revision/final_comparison/scores/acc_scores/results.json')    # P_kmerator scores from the earlier comparison run (used by RK3)
RS_FILE = 'rnaskim.genes.tsv.gz'
EXISTING = G.EXISTING
FRAMES = collections.OrderedDict([('P_rnaskim', EXISTING + ['rnaskim']), ('K', EXISTING + ['kmerator', 'rnaskim']), ('P_rnaskim_z', EXISTING + ['rnaskim_z'])])
STOP_FAILED = G.STOP_FAILED
RATIO_LO, RATIO_HI = G.RATIO_LO, G.RATIO_HI


# ------------------------------------------------------------------ reading RNA-Skim outputs
def read_rnaskim(c, rs_root):
    """Matrices (samples x genes): rnaskim (NaN = no estimate: gene without a quantifiable transcript) and rnaskim_z (same, 0 instead of NaN for such genes).
    Genes outside the RNA-Skim index (not protein-coding) are absent from the file and stay NaN in both.  Returns (V, Vz, missing samples)."""
    V = np.full((len(c.srrs), len(c.genes)), np.nan); Vz = V.copy(); missing = []
    for s in c.srrs:
        p = f'{rs_root}/{c.name}/{s}/{RS_FILE}'
        if not os.path.exists(p) or os.path.exists(f'{rs_root}/{c.name}/{s}/FAILED') or not os.path.exists(f'{rs_root}/{c.name}/{s}/DONE'):
            missing.append(s); continue
        try:
            d = pd.read_csv(p, sep='\t', usecols=['gene', 'est_reads', 'n_tx_quant'], dtype={'gene': str})
        except Exception as e:
            missing.append(s); print(f'  unreadable {p}: {e}', file=sys.stderr); continue
        g = d['gene'].map(G.strip_ver)
        if c.gene_of_id is not None: g = g.map(c.gene_of_id)
        idx = g.map(c.gi); ok = idx.notna().to_numpy()
        i = c.si[s]; cols = idx[ok].astype(int).to_numpy()
        est = pd.to_numeric(d['est_reads'], errors='coerce').to_numpy()[ok]
        nq = pd.to_numeric(d['n_tx_quant'], errors='coerce').to_numpy()[ok]
        Vz[i, cols] = est
        V[i, cols] = np.where(nq > 0, est, np.nan)
    return V, Vz, missing


def attach(c, rs_root, old_root):
    c.V, c.missing = {}, {}
    c.V['rnaskim'], c.V['rnaskim_z'], c.missing['rnaskim'] = read_rnaskim(c, rs_root)
    c.missing['rnaskim_z'] = c.missing['rnaskim']
    c.V['kmerator'], c.missing['kmerator'] = G.read_tool(c, old_root, 'kmerator')
    for m in ('rnaskim', 'rnaskim_z', 'kmerator'):
        for (srr, g), r in c.recs.items():
            v = c.V[m][c.si[srr], c.gi[g]] if g in c.gi else np.nan
            if v == v: r['est'][m] = float(v)


# ------------------------------------------------------------------ integrity checks I1-I4
def integrity(c, rs_root):
    rows = []
    ena = dict(zip(c.samples.srr, c.samples.ena_read_count))
    for s in c.srrs:
        d = f'{rs_root}/{c.name}/{s}'; r = dict(srr=s, ena_read_count=ena[s], dir_present=os.path.isdir(d), DONE=os.path.exists(f'{d}/DONE'), FAILED=os.path.exists(f'{d}/FAILED'))
        try: ij = json.load(open(f'{d}/integrity.json'))
        except Exception: ij = {}
        r['fastq_reads_counted'] = ij.get('fastq_reads_counted'); r['estimated_total_reads'] = ij.get('estimated_total_reads'); r['read_file_source'] = ij.get('read_file_source')
        r['rs_hits'] = ij.get('rs_hits'); r['rs_misses'] = ij.get('rs_misses')
        r['est_reads_per_read'] = (ij['estimated_total_reads'] / ena[s]) if ij.get('estimated_total_reads') is not None and ena[s] else None
        try:
            t = pd.read_csv(f'{d}/{RS_FILE}', sep='\t', usecols=['est_reads', 'n_tx_quant'])
            r['n_genes_with_estimate'] = int(((t.est_reads > 0) & (t.n_tx_quant > 0)).sum()); r['n_genes_unquantifiable'] = int((t.n_tx_quant == 0).sum())
        except Exception:
            r['n_genes_with_estimate'] = None; r['n_genes_unquantifiable'] = None
        rows.append(r)
    df = pd.DataFrame(rows); flags = []
    med = {k: float(df[k].astype(float).median()) for k in ('est_reads_per_read', 'n_genes_with_estimate') if df[k].notna().any()}
    for _, r in df.iterrows():
        s = r['srr']
        if not r['dir_present']: flags.append(dict(srr=s, check='output_dir_missing', value=None, expected=None, ratio=None, flag=True)); continue
        if r['FAILED']: flags.append(dict(srr=s, check='FAILED_marker', value=None, expected=None, ratio=None, flag=True))
        if not r['DONE']: flags.append(dict(srr=s, check='DONE_marker_missing', value=None, expected=None, ratio=None, flag=True))
        v, e = r['fastq_reads_counted'], r['ena_read_count']      # I1 (I4, the loader message, makes kt_sample.sh fail the sample, so it shows up as FAILED)
        flags.append(dict(srr=s, check='I1_reads_counted_vs_ena', value=v, expected=int(e), ratio=None, flag=bool(v is None or v != v or int(v) != int(e))))
        for k, name in (('est_reads_per_read', 'I2_est_total_reads_per_read_ratio_to_median'), ('n_genes_with_estimate', 'I3_genes_with_estimate_ratio_to_median')):
            val = r[k]
            if val is None or val != val: flags.append(dict(srr=s, check=name, value=None, expected=med.get(k), ratio=None, flag=True, note='value missing')); continue
            ratio = float(val) / med[k] if med.get(k) else None
            flags.append(dict(srr=s, check=name, value=float(val), expected=med.get(k), ratio=ratio, flag=bool(ratio is None or ratio < RATIO_LO or ratio > RATIO_HI)))
    fl = pd.DataFrame(flags)
    failed = [s for s in c.srrs if s in c.missing['rnaskim']]
    summary = dict(medians=med, n_samples=len(c.srrs), n_flagged_checks=int(fl.flag.sum()) if len(fl) else 0,
                   flagged_samples=sorted(fl[fl.flag].srr.unique().tolist()) if len(fl) else [], failed_or_missing_samples=failed,
                   stop_rule_triggered=bool(len(failed) > STOP_FAILED))
    return df, fl, summary


# ------------------------------------------------------------------ claim rules RK1-RK3 (numbers and booleans only)
def get(R, coh, frame, subset='all'):
    return (R.get(coh, {}).get(frame, {}) or {}).get(subset)


def prior_kmerator():
    """P_kmerator all-core-gene rows exactly as already scored (scores/acc_scores/results.json)."""
    try:
        r = json.load(open(PRIOR_RESULTS)); out = {}
        for coh, v in r['cohorts'].items():
            f = (v['frames'].get('P_kmerator') or {}).get('all')
            if f: out[coh] = dict(vf=f['metrics']['vf'], kmerator=f['metrics']['kmerator'])
        return out, G.sha256(PRIOR_RESULTS)
    except Exception as e:
        return {}, f'unreadable: {e}'


def claim_rules(R, held):
    out = dict(RK1={}, RK2={}, RK3={})
    r1 = {}; r2 = {}
    for coh in G.COHORTS:
        f = get(R, coh, 'P_rnaskim')
        if not f: continue
        h, m = f['metrics']['vf'], f['metrics']['rnaskim']
        hs, ms = h['within_gene_spearman'], m['within_gene_spearman']
        gap = round(ms - hs, 4) if hs is not None and ms is not None else None
        r1[coh] = dict(gap_rnaskim_minus_hybrid=gap, threshold=0.02, passes=None if gap is None else bool(gap <= 0.02), counts_for_rule=G.COHORTS[coh]['heldout'], hybrid_spearman=hs, rnaskim_spearman=ms)
        sh, sm = h['spurious_rate'], m['spurious_rate']
        r2[coh] = dict(hybrid_spurious=sh, rnaskim_spurious=sm, ratio_hybrid_over_rnaskim=(round(sh / sm, 3) if sm else None), threshold=G.S9.MATCH_RATE,
                       passes=None if sh is None or sm is None else bool(sh <= G.S9.MATCH_RATE * sm), counts_for_rule=G.COHORTS[coh]['heldout'])
    for name, d in (('R1', r1), ('R2', r2)):
        vals = [d[k]['passes'] for k in held if k in d]
        out['RK1'][name] = dict(per_cohort=d, holds_in_both_heldout=(bool(all(vals)) if len(vals) == len(held) and None not in vals else None))
    for coh in G.COHORTS:
        for fr in ('P_rnaskim', 'P_rnaskim_z'):
            f = get(R, coh, fr)
            if f:
                mm = f['metrics'][fr[2:]]
                out['RK2'].setdefault(fr, {})[coh] = dict(median_delta=mm['median_delta'], ci95=mm['delta_ci95'], cls=G.classify_r3(mm['median_delta'], mm['delta_ci95']))
    # RK3: hybrid ahead of X on rank accuracy (Spearman gap > 0.02 in P(X)) and on noise (lower noise SD), both held-out cohorts; X = kmerator (as already scored), rnaskim
    pk, pk_sha = prior_kmerator()
    out['RK3']['prior_kmerator_results_sha256'] = pk_sha
    for X in ('kmerator', 'rnaskim'):
        per = {}
        for coh in held:
            if X == 'kmerator': src = pk.get(coh)
            else:
                f = get(R, coh, 'P_rnaskim'); src = dict(vf=f['metrics']['vf'], rnaskim=f['metrics']['rnaskim']) if f else None
            if not src: continue
            hv, xv = src['vf'], src[X]
            hs, xs = hv['within_gene_spearman'], xv['within_gene_spearman']; hn, xn = hv.get('noise_sd_raw'), xv.get('noise_sd_raw')
            per[coh] = dict(spearman_hybrid=hs, spearman_X=xs, gap_hybrid_minus_X=(round(hs - xs, 4) if hs is not None and xs is not None else None),
                            ahead_rank=(bool(hs - xs > 0.02) if hs is not None and xs is not None else None), noise_hybrid=hn, noise_X=xn,
                            ahead_noise=(bool(hn < xn) if hn is not None and xn is not None else None))
        ar = [per[c]['ahead_rank'] for c in held if c in per]; an = [per[c]['ahead_noise'] for c in held if c in per]
        out['RK3'][X] = dict(per_cohort=per, ahead_rank_both_heldout=(bool(all(ar)) if len(ar) == len(held) and None not in ar else None),
                             ahead_noise_both_heldout=(bool(all(an)) if len(an) == len(held) and None not in an else None))
    ok = [out['RK3'][X].get('ahead_rank_both_heldout') is True and out['RK3'][X].get('ahead_noise_both_heldout') is True for X in ('kmerator', 'rnaskim')]
    out['RK3']['sentence_more_accurate_than_the_kmer_set_quantifiers_we_benchmarked_allowed'] = bool(all(ok))
    return out


# ------------------------------------------------------------------ score
def cmd_score(a):
    sha = G.load_scorer(a.scorer, a.allow_scorer_mismatch)
    if G.sha256(G.FLAG_TABLE) != G.FLAG_SHA: sys.exit('flag_table.tsv hash differs from the expected value')
    for _m in ('rnaskim', 'rnaskim_z'):                                                      # coverage_table walks NEW_ALL
        if _m not in G.NEW_ALL: G.NEW_ALL.append(_m)
    cohorts = a.cohorts or list(G.COHORTS); os.makedirs(a.out, exist_ok=True)
    held = [k for k in cohorts if G.COHORTS[k]['heldout']]
    meta = dict(kt_score_sha256=G.sha256(os.path.abspath(__file__)), gather_script_sha256=G.sha256(G.__file__), scorer_sha256=sha, flag_table_sha256=G.FLAG_SHA, seed=G.S9.SEED,
                started=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), rs_root=a.rs_root, old_root=a.old_root, frames=a.frames, cohorts=cohorts, label=a.label,
                frames_defined={k: v for k, v in FRAMES.items()}, python=sys.version.split()[0], numpy=np.__version__, pandas=pd.__version__)
    out = dict(meta=meta, cohorts={}, reproduction={}, claims=None)
    repro_ok = True
    for name in cohorts:
        c = G.load_cohort(name, a.frames); G._CTX[name] = c; attach(c, a.rs_root, a.old_root)
        print(f'[{name}] records {len(c.recs)}, samples {len(c.srrs)}; missing: rnaskim {len(c.missing["rnaskim"])}, kmerator {len(c.missing["kmerator"])}', flush=True)
        ok, rows, nrows, _ = G.repro_one(c)
        out['reproduction'][name] = dict(pass_=ok, compare_block=rows, noise_block_informational=nrows); G.write_tsv(f'{a.out}/reproduction_{name}.tsv', pd.DataFrame(rows + nrows))
        print(f'[{name}] reproduction of published values: {"match" if ok else "MISMATCH"}', flush=True); repro_ok &= ok
    if not repro_ok and not a.skip_repro_gate:
        json.dump(G.jclean(out), open(f'{a.out}/results.json', 'w'), indent=1)
        sys.exit('Reproduction of the published values failed; RNA-Skim is not scored')
    meta['repro_gate'] = 'passed' if repro_ok else 'SKIPPED (plumbing test)'
    # stop rule before scoring anything: more than 2 failed samples in one cohort stops the run
    ints = {}
    for name in cohorts:
        ints[name] = integrity(G._CTX[name], a.rs_root)
        if ints[name][2]['stop_rule_triggered'] and not a.skip_repro_gate:
            json.dump(G.jclean(out), open(f'{a.out}/results.json', 'w'), indent=1)
            sys.exit(f'STOP RULE: {name} has {len(ints[name][2]["failed_or_missing_samples"])} failed/missing RNA-Skim samples (> {STOP_FAILED})')
    cache = f'{a.out}/task_cache' if a.resume else None
    tasks, subsets_by, flaginfo = [], {}, {}
    for name in cohorts:
        c = G._CTX[name]; subsets_by[name], flaginfo[name] = G.flag_sets(c)
        for sub, genes in subsets_by[name].items():
            for fname, methods in FRAMES.items(): tasks.append((name, sub, fname, methods, genes))
    results = G.run_tasks(tasks, a.procs, cache)
    R = {}
    for r in results: R.setdefault(r['cohort'], {}).setdefault(r['frame'], {})[r['subset']] = r
    for name in cohorts:                         # set 'K': marked dropped when the scorer's validity rule fails on the all-core-gene rows
        k = get(R, name, 'K')
        if k and not k['valid']: k['dropped_for_validity'] = True
    for name in cohorts:
        c = G._CTX[name]; d = f'{a.out}/{name}'; os.makedirs(d, exist_ok=True)
        mr = G.metric_rows([r for r in results if r['cohort'] == name]); cov = G.coverage_table(c, subsets_by[name]); cov_df = pd.DataFrame(cov)
        if len(mr):
            mr = mr.merge(cov_df[['subset', 'method', 'coverage', 'n_covered', 'n_core_genes_truth50_in10']], on=['subset', 'method'], how='left')
            mr['delta_ci95'] = mr.delta_ci95.map(lambda x: json.dumps(x) if x is not None else '')
        G.write_tsv(f'{d}/metrics.tsv', mr); G.write_tsv(f'{d}/coverage.tsv', cov_df)
        idf, ifl, isum = ints[name]
        G.write_tsv(f'{d}/integrity_samples.tsv', idf); G.write_tsv(f'{d}/integrity_flags.tsv', ifl); G.write_tsv(f'{d}/integrity_flagged_only.tsv', ifl[ifl.flag] if len(ifl) else ifl)
        out['cohorts'][name] = dict(n_samples=len(c.srrs), subsets=flaginfo[name], rnaskim_missing_samples=c.missing['rnaskim'], frames={f: dict(v) for f, v in R.get(name, {}).items()},
                                    coverage=cov, integrity=isum)
        print(f'[{name}] wrote {d}', flush=True)
    out['claims'] = claim_rules(R, held); cl = []
    for rule, body in out['claims'].items(): G.flat(body, [rule], cl)
    G.write_tsv(f'{a.out}/claims.tsv', pd.DataFrame(cl, columns=['path', 'value']))
    meta['decisions'] = G.DECISIONS + ['K1 rnaskim = sum of rs_estimate transcript estimates per gene (kt_gsum.py); missing where no transcript carries a selected sig-mer (n_tx_quant == 0); genes outside the protein-coding index are absent (missing).',
                                         'K2 rnaskim_z = same file with 0 instead of missing for genes in the index without a quantifiable transcript (sensitivity, descriptive).',
                                         'K3 Frame K dropped (reported, not used for claims) when the scorer validity rule fails on the all-core-gene rows.',
                                         'K4 Integrity I2 = rs_estimate "estimated total reads" / ENA read count, I3 = genes with est_reads > 0 and n_tx_quant > 0; both within 0.5-2x of the cohort-tissue median (flagged samples included in the median, as D8).']
    json.dump(G.jclean(out), open(f'{a.out}/results.json', 'w'), indent=1)
    print('results.json written')


# ------------------------------------------------------------------ timing (RK4)
def cmd_timing(a):
    os.makedirs(a.out, exist_ok=True)
    cols = ['layout', 'sample', 'method', 'threads', 'rep', 'wall_s', 'cpu_s', 'rss_mb', 'rc']
    t = pd.read_csv(a.raw, sep='\t', header=None, names=cols)
    t = t[t.rc == 0]
    G.write_tsv(f'{a.out}/timing_per_run.tsv', t)
    rows = []
    t = t[(t.threads == 4) | (t.method == 'kmerator_counttags')]       # 4-thread blocks only (RK4); Kmerator+countTags is always single-threaded inside them
    p = t.pivot_table(index=['layout', 'sample', 'rep'], columns='method', values='wall_s')
    for other in ('hybrid', 'salmon_ref', 'kmerator_counttags'):
        if other not in p or 'rnaskim' not in p: continue
        r = (p['rnaskim'] / p[other]).dropna()
        for (lay, smp), v in r.groupby(level=[0, 1]):
            rows.append(dict(layout=lay, sample=smp, other=other, n_reps=len(v), median_ratio_rnaskim_over_other=round(float(v.median()), 2), min_ratio=round(float(v.min()), 2), max_ratio=round(float(v.max()), 2),
                             median_wall_rnaskim_s=round(float(p.loc[(lay, smp), 'rnaskim'].median()), 1), median_wall_other_s=round(float(p.loc[(lay, smp), other].median()), 2)))
    G.write_tsv(f'{a.out}/timing_ratios_4thr.tsv', pd.DataFrame(rows))
    rss = t.groupby(['sample', 'method']).rss_mb.max().reset_index()
    G.write_tsv(f'{a.out}/timing_peak_rss.tsv', rss)
    open(f'{a.out}/TIMING_NOTE.txt', 'w').write('RNA-Skim timing blocks were run on AMD EPYC nodes. Ratios are valid within their own paired blocks (hybrid, Salmon-ref, Kmerator+countTags and RNA-Skim in one block on one node). '
                                                'They should not be compared with the timing of the main comparison (different machine). RK4: ratios at 4 threads only, stated with sample, read length and thread count.\n')
    print(pd.DataFrame(rows).to_string(index=False))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter); sub = ap.add_subparsers(dest='cmd', required=True)
    s = sub.add_parser('score'); s.add_argument('--frames', required=True); s.add_argument('--rs-root', required=True); s.add_argument('--old-root', required=True); s.add_argument('--out', required=True)
    s.add_argument('--scorer', default=G.SCORER_PATH); s.add_argument('--cohorts', nargs='*', choices=list(G.COHORTS)); s.add_argument('--procs', type=int, default=2)
    s.add_argument('--skip-repro-gate', action='store_true'); s.add_argument('--allow-scorer-mismatch', action='store_true'); s.add_argument('--label', default=''); s.add_argument('--resume', action='store_true')
    t = sub.add_parser('timing'); t.add_argument('--raw', required=True); t.add_argument('--out', required=True)
    a = ap.parse_args()
    cmd_score(a) if a.cmd == 'score' else cmd_timing(a)


if __name__ == '__main__':
    main()
