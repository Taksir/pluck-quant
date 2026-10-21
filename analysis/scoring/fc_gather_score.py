#!/usr/bin/env python3
"""Gather the outputs of all quantifiers and score them with score_e9.py.  The scorer (score_e9.py, in this folder) is
imported and its functions are called as they are; this file only builds records, calls it, and writes tables.  It sets no
new rule or threshold.

Sub-commands
  frames   Turn the rebuilt Keele and Munger scorer inputs into plain input files (records, universe, samples) under --frames.
           Keele: needs the make_tables.py output dir of that tissue (per_gene.tsv, extra.tsv, core_all.tsv.gz).
           Munger: imports the Munger scoring script score_e5_post_lock_fix.py (needs the Munger work tree in the original working folder and
           out/est.npz, out/salmon.pkl) and calls its own truth()/build_recs() unchanged.
  repro    Score the four existing methods (no new tool) per cohort-tissue and compare with the published values (4 decimals).
  score    Reproduction check (unless --skip-repro-gate), then every scoring set with new tools, integrity, coverage and claim rules.
"""
import argparse, collections, csv, gzip, hashlib, importlib.util, json, math, os, pickle, sys, time
import numpy as np
import pandas as pd

# Folder constants below point to the original working folder; set them for your copy.
REV = '/mnt/project-files/revision'
SCORER_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'score_e9.py')
SCORER_SHA = '2662e807f2cd5b105c02ee0aab17720e3ba5e6bc70b23342c740b600bfcc0c9b'
GENES_JSON = f'{REV}/fast_competitors/e9_prep_cache/ref/genes.json'
FLAG_TABLE = f'{REV}/reliability_flag/out/flag_table.tsv'
FLAG_SHA = 'fff978bf6a1199c85e34bee564d4f5b40d87cb90639d7036cf7c333de757d506'
OUT_ROOT = f'{REV}/final_comparison/out'
RUNINFO = f'{REV}/fast_competitors/e9_prep_cache/ref/runinfo.tsv'

COHORTS = collections.OrderedDict([
    ('keele_kidney', dict(kind='keele', tissue='kidney', flag_col='flag_50bp', heldout=True)),
    ('keele_liver', dict(kind='keele', tissue='liver', flag_col='flag_50bp', heldout=False)),
    ('munger', dict(kind='munger', tissue='liver', flag_col='flag_100bp', heldout=True)),
])
EXISTING = ['ref', 'vf', 'salmon_ref', 'salmon_pan']             # the four existing methods: scored alone (set 'B') and as the base of every 'P_<tool>' set
NEW_ACC = ['kallisto_ref', 'kallisto_pan', 'kmerator', 'star_fc']     # set 'A' = EXISTING + NEW_ACC
# kallisto-pan was run on Munger only, not on the 50 bp Keele cohorts; there set 'A' uses NEW_ACC without it.
POLICY_SKIPPED = {'kallisto_pan': ('keele_kidney', 'keele_liver')}
def acc_for(cohort): return [m for m in NEW_ACC if cohort not in POLICY_SKIPPED.get(m, ())]
NEW_RED = ['salmon_ref_p10', 'salmon_ref_p25', 'salmon_pan_p10', 'salmon_pan_p25']
FRAME_S = ['vf', 'salmon_ref', 'salmon_pan'] + NEW_RED               # set 'S': vf, full-read Salmon and the reduced-read Salmon arms
# method -> (file, value column)
TOOL_FILES = {'kallisto_ref': ('kallisto_ref.genes.tsv.gz', 'est_counts'), 'kallisto_pan': ('kallisto_pan.genes.tsv.gz', 'est_counts'),
              'kmerator': ('kmerator.genes.tsv.gz', 'mean'), 'star_fc': ('featurecounts.genes.tsv.gz', 'count')}
for _s in ('ref', 'pan'):
    for _p in (10, 25):
        TOOL_FILES[f'salmon_{_s}_p{_p}'] = (f'salmon_{_s}_p{_p}.genes.tsv.gz', 'est_counts')
NEW_ALL = NEW_ACC + NEW_RED
RATIO_LO, RATIO_HI = 0.5, 2.0            # allowed range of a sample's value relative to the cohort median
SUBSAMPLE_TOL = 0.01                     # subsample fraction tolerance, relative to nominal (D6 in DECISIONS)
STOP_FAILED = 2                          # more than 2 failed samples of one cohort for one tool -> listed in stop_rule_tools
SUBSETS = ['all', 'unflagged', 'flagged']
S9 = None


# ---------------------------------------------------------------- utilities
def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def load_scorer(path=SCORER_PATH, allow_mismatch=False):
    global S9
    got = sha256(path)
    if got != SCORER_SHA and not allow_mismatch:
        sys.exit(f'scorer hash mismatch: {path} {got} != {SCORER_SHA}')
    spec = importlib.util.spec_from_file_location('score_e9', path)
    S9 = importlib.util.module_from_spec(spec); spec.loader.exec_module(S9)
    S9._orig_calibrate = S9.calibrate
    return got


def calibrate_nonhit(recs, methods):
    """Munger calibration, same logic as score_e5_post_lock_fix.py: C per sample = median(truth/est) over the
    sample's NON-HIT sample-genes with truth >= 50."""
    by = collections.defaultdict(list)
    for (srr, _), r in recs.items():
        if not r['hit'] and r['truth'] >= S9.MIN_TRUTH:
            for m in methods:
                e = r['est'].get(m)
                if e is not None and e > 0:
                    by[(srr, m)].append(r['truth'] / e)
    return {k: float(np.median(v)) for k, v in by.items()}


def set_calibration(cohort):
    """Keele = the B6 rule (scorer default); Munger = the non-hit rule (score_e5_post_lock_fix.py sets it the same way, by monkey-patching)."""
    S9.calibrate = calibrate_nonhit if COHORTS[cohort]['kind'] == 'munger' else S9._orig_calibrate


def num(x):
    if x is None:
        return None
    try:
        x = float(x)
    except (TypeError, ValueError):
        return None
    return None if (math.isnan(x) or math.isinf(x)) else x


def jclean(o):
    if isinstance(o, dict):
        return {str(k): jclean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jclean(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return None if (math.isnan(o) or math.isinf(o)) else float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return o


def strip_ver(g):
    return g.split('.')[0] if g.startswith('ENSMUSG') else g


def sym_core_map():
    """gene_id -> symbol for core genes with a unique symbol in genes.json (make_tables.py: sym_core)."""
    G = json.load(open(GENES_JSON))
    names = np.array([g['gene_name'] for g in G])
    dup = {n for n, c in zip(*np.unique(names, return_counts=True)) if c > 1}
    return {g['gene_id']: g['gene_name'] for g in G if g['gene_name'] not in dup}


# ---------------------------------------------------------------- frames (rebuilt inputs -> plain files)
REC_COLS = ['srr', 'gene', 'strain', 'stratum', 'founders', 'hit', 'truth', 'est_ref', 'est_vf', 'salmon_ref', 'salmon_pan']


def write_frame(fdir, rec, universe, samples):
    os.makedirs(fdir, exist_ok=True)
    rec[REC_COLS].to_csv(f'{fdir}/records.tsv.gz', sep='\t', index=False)
    universe.to_csv(f'{fdir}/universe.tsv.gz', sep='\t', index=False)
    samples.to_csv(f'{fdir}/samples.tsv', sep='\t', index=False)
    print(f'frame {fdir}: records {len(rec)}, samples {len(samples)}, universe rows {len(universe)}')


def cmd_frames(a):
    ri = pd.read_csv(RUNINFO, sep='\t', usecols=['run_accession', 'sample_title', 'read_count']).set_index('run_accession')
    if a.cohort in ('keele_kidney', 'keele_liver'):
        t = a.tabdir
        pg = pd.read_csv(f'{t}/per_gene.tsv', sep='\t', dtype={'strain': str, 'founders': str}, float_precision='round_trip')
        ex = pd.read_csv(f'{t}/extra.tsv', sep='\t', float_precision='round_trip').drop_duplicates(['srr', 'gene', 'method'], keep='last')   # the scorer's loader keeps the last row
        sal = ex.pivot(index=['srr', 'gene'], columns='method', values='estimate').reset_index()
        for m in ('salmon_ref', 'salmon_pan'):
            if m not in sal: sal[m] = np.nan
        rec = pg.merge(sal[['srr', 'gene', 'salmon_ref', 'salmon_pan']], on=['srr', 'gene'], how='left')
        rec['founders'] = rec['founders'].fillna('')
        rec['hit'] = rec['hit'].astype(str)
        rec = rec[REC_COLS]
        ca = pd.read_csv(f'{t}/core_all.tsv.gz', sep='\t')
        key = pd.MultiIndex.from_frame(ca[['srr', 'gene']])
        uni = pd.DataFrame({'srr': ca.srr, 'gene': ca.gene, 'truth': ca.truth, 'has_vf': ca.hyb.notna().astype(int), 'has_ref': ca.ref26.notna().astype(int)})
        for m in ('salmon_ref', 'salmon_pan'):
            k = ex[ex.method == m]
            uni['has_' + m] = key.isin(pd.MultiIndex.from_frame(k[['srr', 'gene']])).astype(int)
        srrs = sorted(pg.srr.unique())
        samp = pd.DataFrame({'srr': srrs, 'title': [ri.sample_title[s] for s in srrs], 'ena_read_count': [int(ri.read_count[s]) for s in srrs]})
        if a.cohort == 'keele_kidney' and a.sample_table:     # Keele kidney sample table: must agree with runinfo
            st = pd.read_csv(a.sample_table, sep='\t').set_index('srr')
            bad = [s for s in srrs if int(st.read_count[s]) != int(ri.read_count[s])]
            if bad: sys.exit(f'read_count disagreement between samples.tsv and runinfo.tsv: {bad}')
        write_frame(f'{a.frames}/{a.cohort}', rec, uni, samp)
    else:   # munger: the Munger scoring script's own functions, unchanged
        spec = importlib.util.spec_from_file_location('score_e5_post_lock_fix', f'{REV}/e5_munger/scripts/score_e5_post_lock_fix.py')
        s5 = importlib.util.module_from_spec(spec); spec.loader.exec_module(s5)
        E = dict(np.load(f'{s5.O}/est.npz', allow_pickle=True)); sal = pickle.load(open(f'{s5.O}/salmon.pkl', 'rb')); T = s5.truth()
        core_idx = np.nonzero(E['core'])[0]
        recs = s5.build_recs(E, sal, T, core_idx, dict(vf=E['hybrid'], ref=E['ref'], probes_only=E['vf'], emg=E['emg']))
        rows = [(s, g, r['strain'], r['stratum'], r['founders'], str(r['hit']), r['truth'], r['est'].get('ref', np.nan), r['est'].get('vf', np.nan),
                 r['est'].get('salmon_ref', np.nan), r['est'].get('salmon_pan', np.nan)) for (s, g), r in recs.items()]
        rec = pd.DataFrame(rows, columns=REC_COLS)
        uni = pd.DataFrame({'srr': rec.srr, 'gene': rec.gene, 'truth': rec.truth, 'has_vf': rec.est_vf.notna().astype(int), 'has_ref': rec.est_ref.notna().astype(int),
                            'has_salmon_ref': rec.salmon_ref.notna().astype(int), 'has_salmon_pan': rec.salmon_pan.notna().astype(int)})
        title = {r['srr']: r['title'] for r in s5.SAMP}; reads = {r['srr']: int(r['reads']) for r in s5.SAMP}
        srrs = [str(s) for s in E['srrs']]
        samp = pd.DataFrame({'srr': srrs, 'title': [title[s] for s in srrs], 'ena_read_count': [reads[s] for s in srrs]})
        write_frame(f'{a.frames}/munger', rec, uni, samp)


# ---------------------------------------------------------------- loading a cohort-tissue
class Cohort:
    pass


def load_cohort(name, frames):
    c = Cohort(); c.name = name; fd = f'{frames}/{name}'
    c.samples = pd.read_csv(f'{fd}/samples.tsv', sep='\t')
    c.srrs = list(c.samples.srr); c.si = {s: i for i, s in enumerate(c.srrs)}
    rec = pd.read_csv(fd + '/records.tsv.gz', sep='\t', dtype={'strain': str, 'founders': str, 'hit': str}, keep_default_na=False, na_values=[''], float_precision='round_trip')
    c.recs = {}
    for r in rec.itertuples(index=False):
        est = {}
        for m, v in (('ref', r.est_ref), ('vf', r.est_vf), ('salmon_ref', r.salmon_ref), ('salmon_pan', r.salmon_pan)):
            if v == v: est[m] = float(v)
        c.recs[(r.srr, r.gene)] = dict(strain=r.strain, stratum=r.stratum, founders='' if r.founders != r.founders else r.founders,
                                       hit=r.hit == 'True', truth=float(r.truth), est=est)
    c.universe = pd.read_csv(fd + '/universe.tsv.gz', sep='\t')
    if COHORTS[name]['kind'] == 'keele':
        s2s = sym_core_map(); c.gene_of_id = {gid: sym for gid, sym in s2s.items()}
        c.genes = sorted(set(s2s.values()))
    else:
        c.gene_of_id = None
        c.genes = sorted(set(rec.gene) | set(c.universe.gene))
    c.gi = {g: i for i, g in enumerate(c.genes)}
    return c


def flag_sets(c):
    ft = pd.read_csv(FLAG_TABLE, sep='\t', dtype={'chrom': str})
    col = COHORTS[c.name]['flag_col']
    fl = dict(zip(ft.gene_id, ft[col].astype(str) == 'True'))
    if c.gene_of_id is not None:
        keyflag = {c.gene_of_id[g]: v for g, v in fl.items() if g in c.gene_of_id}
    else:
        keyflag = dict(fl)
    allg = {g for (_, g) in c.recs} | set(c.universe.gene)
    flagged = {g for g in allg if keyflag.get(g) is True}
    unflagged = {g for g in allg if keyflag.get(g) is False}
    return dict(all=None, flagged=flagged, unflagged=unflagged), dict(n_genes_in_frame=len(allg), n_flagged=len(flagged), n_unflagged=len(unflagged),
                                                                     n_not_in_flag_table=len(allg) - len(flagged) - len(unflagged))


# ---------------------------------------------------------------- new-tool outputs
def read_tool(c, out_root, method):
    """Returns matrix (n_samples x n_genes) of raw estimates (NaN = no estimate) and the list of samples with a missing file."""
    fn, col = TOOL_FILES[method]
    V = np.full((len(c.srrs), len(c.genes)), np.nan); missing = []
    for s in c.srrs:
        p = f'{out_root}/{c.name}/{s}/{fn}'
        if not os.path.exists(p):
            missing.append(s); continue
        try:
            d = pd.read_csv(p, sep='\t', usecols=['gene', col], dtype={'gene': str})
        except Exception as e:           # unreadable file counts as a technical failure
            missing.append(s); print(f'  unreadable {p}: {e}', file=sys.stderr); continue
        g = d['gene'].map(strip_ver)
        if c.gene_of_id is not None: g = g.map(c.gene_of_id)       # Ensembl id -> symbol (unique-symbol core genes only)
        idx = g.map(c.gi)
        ok = idx.notna().to_numpy()
        V[c.si[s], idx[ok].astype(int).to_numpy()] = pd.to_numeric(d[col], errors='coerce').to_numpy()[ok]
    return V, missing


def attach_tools(c, out_root, methods):
    c.V, c.missing = {}, {}
    for m in methods:
        c.V[m], c.missing[m] = read_tool(c, out_root, m)
        for (srr, g), r in c.recs.items():
            v = c.V[m][c.si[srr], c.gi[g]] if g in c.gi else np.nan
            if v == v: r['est'][m] = float(v)


def covered_samples(c, methods):
    cov = set(c.srrs)
    for m in methods:
        if m in c.V:
            cov &= {s for s in c.srrs if np.isfinite(c.V[m][c.si[s]]).any()}
    return cov


# ---------------------------------------------------------------- scoring one frame
def noise_stats(recs, methods, samples):
    """score_noise.py definition on the frame's sample-genes (the same ones score_tissue uses): lr = log2((C est + 1)/(truth + 1)),
    within-gene SD = 1.4826 MAD over samples (genes with >= 10 sample-genes), reported = median over genes; median |lr| over sample-genes."""
    C = S9.calibrate(recs, methods)
    per = {m: collections.defaultdict(list) for m in methods}; ab = {m: [] for m in methods}
    for (srr, g), r in recs.items():
        if samples is not None and srr not in samples: continue
        if r['truth'] < S9.MIN_TRUTH or r['stratum'] not in ('B6', 'nonB6'): continue
        if any(r['est'].get(m) is None or (srr, m) not in C for m in methods): continue
        for m in methods:
            v = float(np.log2((C[(srr, m)] * r['est'][m] + 1) / (r['truth'] + 1)))
            per[m][g].append(v); ab[m].append(abs(v))
    out = {}
    for m in methods:
        sd = [1.4826 * float(np.median(np.abs(np.array(v) - np.median(v)))) for v in per[m].values() if len(v) >= 10]
        out[m] = dict(noise_sd_raw=float(np.median(sd)) if sd else None, noise_sd=round(float(np.median(sd)), 4) if sd else None,
                      median_abs_err=round(float(np.median(ab[m])), 4) if ab[m] else None, n_genes_noise=len(sd), n_sample_genes=len(ab[m]))
    return out


_CTX = {}      # cohort name -> Cohort (inherited by forked workers)


def score_frame(task):
    cname, subset, fname, methods, genes = task
    c = _CTX[cname]; set_calibration(cname)
    recs = c.recs if genes is None else {k: r for k, r in c.recs.items() if k[1] in genes}
    cov = covered_samples(c, methods)
    t0 = time.time()
    met, _ = S9.score_tissue(recs, methods, np.random.default_rng(S9.SEED), samples=cov)
    nz = noise_stats(recs, methods, cov)
    ng = met[methods[0]]['n_genes']; ns = met[methods[0]]['n_samples']
    res = dict(cohort=cname, subset=subset, frame=fname, methods=methods, n_samples=ns, n_delta_genes=ng,
               valid=bool(ns >= S9.MIN_SAMPLES and ng >= S9.MIN_GENES),
               below_validity_floor=bool(ns < S9.MIN_SAMPLES or ng < S9.MIN_GENES), metrics={}, seconds=round(time.time() - t0, 1))
    for m in methods:
        x = dict(met[m]); x.update(nz[m])
        if m == 'kmerator': x['median_abs_err'] = None            # not reported for Kmerator
        res['metrics'][m] = x
    return res


def coverage_table(c, subsets):
    """Among core genes (universe) with truth >= 50 in >= 10 samples, share with an estimate in >= 10 of those samples."""
    U = c.universe; U = U[U.truth >= S9.MIN_TRUTH].copy()
    gi = U.gene.map(c.gi); si = U.srr.map(c.si)
    for m in NEW_ALL:
        if m in c.V:
            ok = gi.notna() & si.notna()
            v = np.zeros(len(U), bool)
            v[ok.to_numpy()] = np.isfinite(c.V[m][si[ok].astype(int).to_numpy(), gi[ok].astype(int).to_numpy()])
            U['has_' + m] = v.astype(int)
    out = []
    for sub, genes in subsets.items():
        W = U if genes is None else U[U.gene.isin(genes)]
        n = W.groupby('gene').size(); q = n[n >= 10].index; Wq = W[W.gene.isin(q)]
        for col in [x for x in Wq.columns if x.startswith('has_')]:
            k = Wq.groupby('gene')[col].sum()
            out.append(dict(cohort=c.name, subset=sub, method={'has_vf': 'vf', 'has_ref': 'ref'}.get(col, col[4:]), n_core_genes_truth50_in10=int(len(q)),
                            n_covered=int((k >= 10).sum()), coverage=round(float((k >= 10).mean()), 4) if len(q) else None))
    return out


# ---------------------------------------------------------------- integrity checks (new tools)
def parse_kv(path, key, sep='|'):
    for l in open(path):
        if key in l:
            return int(l.split(sep)[-1].strip().replace(',', ''))
    return None


def read_json(p):
    try:
        return json.load(open(p))
    except Exception:
        return None


def sample_integrity(d, ena):
    """All per-sample raw numbers needed for the integrity checks. Prefers the raw tool files, falls back to integrity.json."""
    v = {}; ij = read_json(f'{d}/integrity.json') or {}
    v['integrity_json_ena_read_count'] = ij.get('ena_read_count')
    for s in ('ref', 'pan'):
        ri = read_json(f'{d}/kallisto_{s}_run_info.json') or {}
        v[f'kallisto_{s}_processed'] = ri.get('n_processed', ij.get(f'kallisto_{s}_n_processed'))
        pa = ri.get('n_pseudoaligned', ij.get(f'kallisto_{s}_n_pseudoaligned'))
        v[f'kallisto_{s}_pseudo_frac'] = pa / v[f'kallisto_{s}_processed'] if pa is not None and v[f'kallisto_{s}_processed'] else None
    sl = f'{d}/star_Log.final.out'
    if os.path.exists(sl):
        v['star_processed'] = parse_kv(sl, 'Number of input reads'); uq = parse_kv(sl, 'Uniquely mapped reads number')
    else:
        v['star_processed'] = ij.get('star_input_reads'); uq = ij.get('star_unique')
    v['star_unique_frac'] = uq / v['star_processed'] if uq is not None and v['star_processed'] else None
    ct = f'{d}/counttags_summary.txt'; tr = None
    if os.path.exists(ct):
        for l in open(ct):
            if l.startswith('total_reads'):
                tr = int(l.split('\t')[-1].strip())
    if tr is None and ij.get('counttags_total_reads_line'):
        tr = int(ij['counttags_total_reads_line'].split()[-1])
    v['counttags_processed'] = tr
    try:
        km = pd.read_csv(f'{d}/kmerator.genes.tsv.gz', sep='\t', usecols=['sum'])
        v['counttags_hits_per_read'] = float(km['sum'].sum()) / tr if tr else None
    except Exception:
        v['counttags_hits_per_read'] = None
    fs = f'{d}/featurecounts.summary.txt'
    if os.path.exists(fs):
        st = {}
        for l in open(fs):
            x = l.rstrip('\n').split('\t')
            if len(x) >= 2 and x[0] != 'Status':
                st[x[0]] = int(x[1])
        tot = sum(st.values()); v['featurecounts_assigned_frac'] = st.get('Assigned', 0) / tot if tot else None
    else:
        v['featurecounts_assigned_frac'] = None
    sub = {}
    sp = f'{d}/subsample_reads.txt'
    if os.path.exists(sp):
        for l in open(sp):
            x = l.split()
            if len(x) >= 3: sub[int(x[0])] = (int(x[1]), int(x[2]))
    for pc in (10, 25):
        n_sub, n_full = sub.get(pc, (None, None))
        v[f'sub_p{pc}_reads'] = n_sub; v[f'sub_p{pc}_full_reads_file'] = n_full
        v[f'sub_p{pc}_frac'] = n_sub / ena if n_sub is not None and ena else None
        for s in ('ref', 'pan'):
            mi = read_json(f'{d}/salmon_{s}_p{pc}_meta_info.json') or {}
            pr = mi.get('num_processed', ij.get(f'salmon_{s}_p{pc}_processed')); mp = mi.get('num_mapped', ij.get(f'salmon_{s}_p{pc}_mapped'))
            v[f'salmon_{s}_p{pc}_processed'] = pr
            v[f'salmon_{s}_p{pc}_mapped_frac'] = mp / pr if mp is not None and pr else None
    return v


EQUAL_CHECKS = [('kallisto_ref_processed', 'ena'), ('kallisto_pan_processed', 'ena'), ('star_processed', 'ena'), ('counttags_processed', 'ena')] + \
               [(f'salmon_{s}_p{p}_processed', f'sub_p{p}_reads') for s in ('ref', 'pan') for p in (10, 25)]
RATIO_CHECKS = ['kallisto_ref_pseudo_frac', 'kallisto_pan_pseudo_frac', 'star_unique_frac', 'featurecounts_assigned_frac', 'counttags_hits_per_read'] + \
               [f'salmon_{s}_p{p}_mapped_frac' for s in ('ref', 'pan') for p in (10, 25)]
EXPECTED_FILES = {m: TOOL_FILES[m][0] for m in NEW_ALL}


def integrity(c, out_root):
    rows, flags, miss_by_tool = [], [], {}
    ena = dict(zip(c.samples.srr, c.samples.ena_read_count))
    for s in c.srrs:
        d = f'{out_root}/{c.name}/{s}'
        v = sample_integrity(d, ena[s]) if os.path.isdir(d) else {}
        v['srr'] = s; v['ena_read_count'] = ena[s]; v['dir_present'] = os.path.isdir(d)
        v['FAILED_marker'] = os.path.exists(f'{d}/FAILED'); v['DONE_marker'] = os.path.exists(f'{d}/DONE')
        for m, fn in EXPECTED_FILES.items():
            if c.name in POLICY_SKIPPED.get(m, ()): continue
            if not os.path.exists(f'{d}/{fn}'): miss_by_tool.setdefault(m, []).append(s)
        rows.append(v)
    df = pd.DataFrame(rows)
    med = {k: float(df[k].astype(float).median()) for k in RATIO_CHECKS if k in df and df[k].notna().any()}
    for _, r in df.iterrows():
        s = r['srr']
        if not r['dir_present']:
            flags.append(dict(srr=s, check='output_dir_missing', value=None, expected=None, ratio=None, flag=True)); continue
        if r['FAILED_marker']:
            flags.append(dict(srr=s, check='FAILED_marker', value=None, expected=None, ratio=None, flag=True))
        if r.get('integrity_json_ena_read_count') is not None and r['integrity_json_ena_read_count'] == r['integrity_json_ena_read_count'] and int(r['integrity_json_ena_read_count']) != int(r['ena_read_count']):
            flags.append(dict(srr=s, check='integrity_json_ena_vs_sample_table', value=int(r['integrity_json_ena_read_count']), expected=int(r['ena_read_count']), ratio=None, flag=True))
        for chk, ref in EQUAL_CHECKS:
            if chk.startswith('kallisto_pan') and c.name in POLICY_SKIPPED['kallisto_pan']: continue
            val = r.get(chk); exp = r['ena_read_count'] if ref == 'ena' else r.get(ref)
            if val is None or val != val or exp is None or exp != exp:
                flags.append(dict(srr=s, check=chk + '_vs_' + ref, value=None, expected=None, ratio=None, flag=True, note='value missing')); continue
            flags.append(dict(srr=s, check=chk + '_vs_' + ref, value=int(val), expected=int(exp), ratio=None, flag=bool(int(val) != int(exp))))
        for k in RATIO_CHECKS:
            if k.startswith('kallisto_pan') and c.name in POLICY_SKIPPED['kallisto_pan']: continue
            val = r.get(k)
            if val is None or val != val:
                flags.append(dict(srr=s, check=k + '_ratio_to_median', value=None, expected=med.get(k), ratio=None, flag=True, note='value missing')); continue
            ratio = float(val) / med[k] if med.get(k) else None
            flags.append(dict(srr=s, check=k + '_ratio_to_median', value=float(val), expected=med.get(k), ratio=ratio,
                              flag=bool(ratio is None or ratio < RATIO_LO or ratio > RATIO_HI)))
        for p, nom in ((25, 0.25), (10, 0.10)):
            fr = r.get(f'sub_p{p}_frac')
            if fr is None or fr != fr:
                flags.append(dict(srr=s, check=f'subsample_p{p}_realised_fraction', value=None, expected=nom, ratio=None, flag=True, note='value missing')); continue
            flags.append(dict(srr=s, check=f'subsample_p{p}_realised_fraction', value=float(fr), expected=nom, ratio=float(fr) / nom,
                              flag=bool(abs(float(fr) / nom - 1) > SUBSAMPLE_TOL)))
    fl = pd.DataFrame(flags)
    summary = dict(medians=med, n_samples=len(c.srrs), n_flagged_checks=int(fl.flag.sum()) if len(fl) else 0,
                   flagged_samples=sorted(fl[fl.flag].srr.unique().tolist()) if len(fl) else [],
                   missing_outputs_by_tool={m: v for m, v in miss_by_tool.items()},
                   stop_rule_tools=[m for m, v in miss_by_tool.items() if len(v) > STOP_FAILED])
    return df, fl, summary


# ---------------------------------------------------------------- claim rules R1-R6, numbers and booleans only
def classify_r3(med, ci):
    if med is None: return None
    if abs(med) <= S9.VF_BAND: return 'unbiased'
    if med < -S9.VF_BAND and ci and ci[1] is not None and ci[1] < 0: return 'biased'
    return 'inconclusive'


def claim_rules(R):
    """R[cohort][frame]['all'] = score_frame result."""
    def g(coh, frame):
        return (R.get(coh, {}).get(frame, {}) or {}).get('all')
    held = [k for k, v in COHORTS.items() if v['heldout']]
    out = dict(R1={}, R2={}, R3={}, R4={}, R5={}, R6={}, hybrid_spurious_rate={})
    for coh in COHORTS:
        f = g(coh, 'P_' + NEW_ACC[0]) or g(coh, 'B')
        if f: out['hybrid_spurious_rate'][coh] = f['metrics']['vf']['spurious_rate']
    for M in NEW_ACC:
        r1 = {}; r2 = {}
        for coh in COHORTS:
            f = g(coh, 'P_' + M)
            if not f: continue
            h, m = f['metrics']['vf'], f['metrics'][M]
            gap = round(m['within_gene_spearman'] - h['within_gene_spearman'], 4) if m['within_gene_spearman'] is not None and h['within_gene_spearman'] is not None else None
            r1[coh] = dict(gap_M_minus_hybrid=gap, threshold=0.02, passes=None if gap is None else bool(gap <= 0.02), counts_for_rule=COHORTS[coh]['heldout'],
                           hybrid_spearman=h['within_gene_spearman'], M_spearman=m['within_gene_spearman'])
            hs, ms = h['spurious_rate'], m['spurious_rate']
            r2[coh] = dict(hybrid_spurious=hs, M_spurious=ms, ratio_hybrid_over_M=(round(hs / ms, 3) if ms else None), threshold=S9.MATCH_RATE,
                           passes=None if hs is None or ms is None else bool(hs <= S9.MATCH_RATE * ms), counts_for_rule=COHORTS[coh]['heldout'])
        for name, d in (('R1', r1), ('R2', r2)):
            vals = [d[k]['passes'] for k in held if k in d]
            d_all = dict(per_cohort=d, holds_in_both_heldout=(bool(all(vals)) if len(vals) == len(held) and None not in vals else None))
            out[name][M] = d_all
    for M in NEW_ACC + NEW_RED:
        r3 = {}
        for coh in COHORTS:
            f = g(coh, 'P_' + M)
            if not f: continue
            m = f['metrics'][M]
            r3[coh] = dict(median_delta=m['median_delta'], ci95=m['delta_ci95'], cls=classify_r3(m['median_delta'], m['delta_ci95']))
        out['R3'][M] = r3
    k3 = out['R3'].get('kmerator', {})
    if len(k3) == len(COHORTS):
        le = {c: bool(k3[c]['median_delta'] is not None and k3[c]['median_delta'] <= -0.08) for c in k3}
        bi = {c: k3[c]['cls'] == 'biased' for c in k3}
        out['R4'] = dict(median_delta=k3 and {c: k3[c]['median_delta'] for c in k3}, le_minus_0_08=le, r3_biased=bi,
                         outcome='confirmed' if all(le.values()) else ('partial' if all(bi.values()) else 'not confirmed'))
    # R5: noise SD of the reduced-read arm <= hybrid noise SD (same 'P_<arm>' set) in both held-out cohorts
    for base in ('salmon_ref', 'salmon_pan'):
        per_arm = {}
        for pc in (10, 25):
            arm = f'{base}_p{pc}'; d = {}
            for coh in held:
                f = g(coh, 'P_' + arm)
                if not f: continue
                a_, h_ = f['metrics'][arm]['noise_sd_raw'], f['metrics']['vf']['noise_sd_raw']
                d[coh] = dict(noise_sd_arm=f['metrics'][arm]['noise_sd'], noise_sd_hybrid=f['metrics']['vf']['noise_sd'],
                              precision_matched=None if a_ is None or h_ is None else bool(a_ <= h_),
                              spearman_arm=f['metrics'][arm]['within_gene_spearman'], spearman_hybrid=f['metrics']['vf']['within_gene_spearman'])
            ok = [d[c]['precision_matched'] for c in held if c in d]
            per_arm[arm] = dict(per_cohort=d, precision_matched_both_heldout=(bool(all(ok)) if len(ok) == len(held) and None not in ok else None))
            lv = g('keele_liver', 'P_' + arm)
            if lv: per_arm[arm]['keele_liver_reported'] = dict(noise_sd_arm=lv['metrics'][arm]['noise_sd'], noise_sd_hybrid=lv['metrics']['vf']['noise_sd'])
        fstar = None
        for pc in (10, 25):
            v = per_arm[f'{base}_p{pc}']['precision_matched_both_heldout']
            if v is True: fstar = f'{pc}%'; break
        complete = all(per_arm[f'{base}_p{pc}']['precision_matched_both_heldout'] is not None for pc in (10, 25))
        out['R5'][base] = dict(arms=per_arm, f_star=(fstar if fstar else ('none below 100%' if complete else None)))
    # R6: R3 class per reduced-read arm and cohort-tissue; bias_sentence_allowed is True only where the class is 'biased'
    for arm in NEW_RED:
        out['R6'][arm] = {coh: dict(r3=v['cls'], bias_sentence_allowed=bool(v['cls'] == 'biased')) for coh, v in out['R3'].get(arm, {}).items()}
    return out


# ---------------------------------------------------------------- repro (the four existing methods vs published compare blocks)
PUBLISHED = dict(keele_kidney=(f'{REV}/e4_kidney/kidney_frame/e4_gates.json', ('compare',)),
                 keele_liver=(f'{REV}/e4_kidney/liver_frame/e4_gates.json', ('compare',)),
                 munger=(f'{REV}/e5_munger/metrics.json', ('core', 'compare')))
PUB_NOISE = dict(keele_kidney=(f'{REV}/e4_kidney/kidney_frame/e4_gates.json', ('noise',)), keele_liver=(f'{REV}/e4_kidney/liver_frame/e4_gates.json', ('noise',)))
REPRO_KEYS = ['n_samples', 'n_genes', 'within_gene_spearman', 'spurious_rate', 'median_delta']


def dig(d, path):
    for p in path: d = d[p]
    return d


def repro_one(c, published_override=None):
    res = score_frame((c.name, 'all', 'B', EXISTING, None))
    pub = dig(*[json.load(open(PUBLISHED[c.name][0])), PUBLISHED[c.name][1]])
    rows, ok = [], True
    for m in EXISTING:
        for k in REPRO_KEYS:
            a, b = res['metrics'][m][k], pub[m][k]
            same = (a == b) if k in ('n_samples', 'n_genes') else (a is not None and b is not None and f'{a:.4f}' == f'{b:.4f}')
            ok &= bool(same); rows.append(dict(cohort=c.name, method=m, key=k, rebuilt=a, published=b, equal=bool(same)))
    noise_rows = []
    if c.name in PUB_NOISE:       # informational: noise block of e4_gates.json (hybrid=vf, ref, salmon_ref, salmon_pan)
        pn = dig(json.load(open(PUB_NOISE[c.name][0])), PUB_NOISE[c.name][1])
        for m, pm in (('vf', 'hybrid'), ('ref', 'ref'), ('salmon_ref', 'salmon_ref'), ('salmon_pan', 'salmon_pan')):
            for k, kk in (('noise_sd', 'within_gene_sd'), ('median_abs_err', 'median_abs_err'), ('n_genes_noise', 'n_genes')):
                a, b = res['metrics'][m][k], pn[pm][kk]
                noise_rows.append(dict(cohort=c.name, method=m, key=k, rebuilt=a, published=b, equal=bool(a is not None and f'{a:.4f}' == f'{b:.4f}')))
    elif c.name == 'munger':
        pn = json.load(open(PUBLISHED['munger'][0]))['noise_core']
        for m in EXISTING:
            a, b = res['metrics'][m]['noise_sd'], pn[m]['median_within_gene_sd']
            noise_rows.append(dict(cohort=c.name, method=m, key='noise_sd', rebuilt=a, published=b, equal=bool(a is not None and f'{a:.4f}' == f'{b:.4f}')))
    return ok, rows, noise_rows, res


# ---------------------------------------------------------------- score
def write_tsv(path, df):
    df.to_csv(path, sep='\t', index=False)


def validity_note(r):
    if not r['below_validity_floor']: return ''
    why = []
    if r['n_delta_genes'] < S9.MIN_GENES: why.append(f"n_delta_genes={r['n_delta_genes']} < {S9.MIN_GENES}")
    if r['n_samples'] < S9.MIN_SAMPLES: why.append(f"n_samples={r['n_samples']} < {S9.MIN_SAMPLES}")
    return "below the scorer's validity floor (" + '; '.join(why) + ')'


def metric_rows(results):
    rows = []
    for r in results:
        for m, x in r['metrics'].items():
            rows.append(dict(cohort=r['cohort'], frame=r['frame'], subset=r['subset'], method=m, n_samples=r['n_samples'], n_delta_genes=r['n_delta_genes'],
                             **{k: x.get(k) for k in ('median_delta', 'delta_ci95', 'spurious_rate', 'within_gene_spearman', 'n_genes_spearman', 'delta_wild',
                                                      'delta_domestic', 'pooled_hit_med_log2', 'pooled_B6_med_log2', 'C_median', 'noise_sd', 'median_abs_err', 'n_genes_noise')},
                             valid=r['valid'], note=validity_note(r)))
    return pd.DataFrame(rows)


def cmd_score(a):
    sha = load_scorer(a.scorer, a.allow_scorer_mismatch)
    if sha256(FLAG_TABLE) != FLAG_SHA: sys.exit('flag_table.tsv hash differs from the expected value')
    cohorts = a.cohorts or list(COHORTS)
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time(); meta = dict(script_sha256=sha256(os.path.abspath(__file__)), scorer_sha256=sha, flag_table_sha256=FLAG_SHA, seed=S9.SEED,
                                  started=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), out_root=a.out_root, frames=a.frames, cohorts=cohorts,
                                  label=a.label, python=sys.version.split()[0], numpy=np.__version__, pandas=pd.__version__,
                                  methods=dict(existing=EXISTING, new_accuracy=NEW_ACC, new_reduced=NEW_RED))
    out = dict(meta=meta, cohorts={}, reproduction={}, claims=None)
    # --- load + reproduction check
    repro_ok = True
    for name in cohorts:
        c = load_cohort(name, a.frames); _CTX[name] = c; attach_tools(c, a.out_root, NEW_ALL)
        print(f'[{name}] records {len(c.recs)}, samples {len(c.srrs)}; tools missing: {{{", ".join(m + ":" + str(len(v)) for m, v in c.missing.items())}}}', flush=True)
        ok, rows, nrows, _ = repro_one(c)
        out['reproduction'][name] = dict(pass_=ok, compare_block=rows, noise_block_informational=nrows)
        write_tsv(f'{a.out}/reproduction_{name}.tsv', pd.DataFrame(rows + nrows))
        print(f'[{name}] reproduction of published values: {"match" if ok else "MISMATCH"}', flush=True)
        repro_ok &= ok
    if not repro_ok and not a.skip_repro_gate:
        json.dump(jclean(out), open(f'{a.out}/results.json', 'w'), indent=1)
        sys.exit('Reproduction of the published values failed; no new tool is scored (use --skip-repro-gate only for plumbing tests on fake truth)')
    if a.skip_repro_gate:
        meta['repro_gate'] = 'SKIPPED (plumbing test)' if not repro_ok else 'skipped but passed'
    else:
        meta['repro_gate'] = 'passed'
    cache = f'{a.out}/task_cache' if a.resume else None
    # --- build task list
    tasks, subsets_by, flaginfo = [], {}, {}
    for name in cohorts:
        c = _CTX[name]; subsets_by[name], flaginfo[name] = flag_sets(c)
        have = [m for m in NEW_ALL if len(c.missing[m]) < len(c.srrs)]
        for sub, genes in subsets_by[name].items():
            for M in have:
                tasks.append((name, sub, 'P_' + M, EXISTING + [M], genes))
            if all(m in have for m in acc_for(name)):
                tasks.append((name, sub, 'A', EXISTING + acc_for(name), genes))
            if all(m in have for m in NEW_RED):
                tasks.append((name, sub, 'S', FRAME_S, genes))
            tasks.append((name, sub, 'B', EXISTING, genes))
    results = run_tasks(tasks, a.procs, cache)
    # --- fallback for set 'A' below the validity floor (decided on the all-core-gene rows; applied to all three subsets of that cohort-tissue)
    extra = []
    for name in cohorts:
        fa = [r for r in results if r['cohort'] == name and r['frame'] == 'A' and r['subset'] == 'all']
        if fa and not fa[0]['valid']:
            nk = [m for m in EXISTING + acc_for(name) if m != 'kmerator']
            for sub, genes in subsets_by[name].items():
                extra.append((name, sub, 'A_without_kmerator', nk, genes))
    results += run_tasks(extra, a.procs, cache)
    # --- per cohort-tissue outputs
    R = {}
    for r in results:
        R.setdefault(r['cohort'], {}).setdefault(r['frame'], {})[r['subset']] = r
    for name in cohorts:
        c = _CTX[name]; d = f'{a.out}/{name}'; os.makedirs(d, exist_ok=True)
        mr = metric_rows([r for r in results if r['cohort'] == name])
        cov = coverage_table(c, subsets_by[name]); cov_df = pd.DataFrame(cov)
        if len(mr):
            mr = mr.merge(cov_df[['subset', 'method', 'coverage', 'n_covered', 'n_core_genes_truth50_in10']], on=['subset', 'method'], how='left')
            mr['delta_ci95'] = mr.delta_ci95.map(lambda x: json.dumps(x) if x is not None else '')
        write_tsv(f'{d}/metrics.tsv', mr); write_tsv(f'{d}/coverage.tsv', cov_df)
        idf, ifl, isum = integrity(c, a.out_root)
        write_tsv(f'{d}/integrity_samples.tsv', idf); write_tsv(f'{d}/integrity_flags.tsv', ifl)
        write_tsv(f'{d}/integrity_flagged_only.tsv', ifl[ifl.flag] if len(ifl) else ifl)
        out['cohorts'][name] = dict(n_samples=len(c.srrs), n_records=len(c.recs), subsets=flaginfo[name], tools_missing_samples=c.missing,
                                    frames={f: {s: x for s, x in v.items()} for f, v in R.get(name, {}).items()},
                                    coverage=cov, integrity=isum)
        print(f'[{name}] wrote {d}', flush=True)
    out['claims'] = claim_rules(R)
    cl = []
    for rule, body in out['claims'].items():
        flat(body, [rule], cl)
    write_tsv(f'{a.out}/claims.tsv', pd.DataFrame(cl, columns=['path', 'value']))
    meta['seconds'] = round(time.time() - t0, 1)
    meta['decisions'] = DECISIONS
    json.dump(jclean(out), open(f'{a.out}/results.json', 'w'), indent=1)
    print(f'results.json written ({meta["seconds"]} s)')


def flat(o, path, out):
    if isinstance(o, dict):
        for k, v in o.items(): flat(v, path + [str(k)], out)
    else:
        out.append(('/'.join(path), json.dumps(jclean(o))))


def task_file(cache, t):
    return f'{cache}/{t[0]}__{t[1]}__{t[2]}.json' if cache else None


def run_tasks(tasks, procs, cache=None):
    """Scores the frames (2 processes by default). With a cache dir (--resume) finished tasks are stored one JSON each and skipped on a re-run
    (only valid while the inputs are unchanged). A worker that dies (out of memory) does not hang the run: the unfinished tasks are retried."""
    if not tasks: return []
    done, todo = [], []
    for t in tasks:
        f = task_file(cache, t)
        if f and os.path.exists(f):
            try:
                r = json.load(open(f))
                if r.get('methods') == t[3]:
                    done.append(r); continue
            except Exception:
                pass
        todo.append(t)
    if cache: os.makedirs(cache, exist_ok=True)
    if done: print(f'  resumed {len(done)} finished tasks from {cache}', flush=True)

    def finish(r, i, n):
        f = task_file(cache, (r['cohort'], r['subset'], r['frame']))
        if f:
            json.dump(jclean(r), open(f + '.part', 'w')); os.replace(f + '.part', f)
        done.append(r); print(f'  [{i}/{n}] {r["cohort"]} {r["subset"]} {r["frame"]} {r["seconds"]} s', flush=True)
    n_all = len(todo); k = 0
    if procs <= 1:
        for t in todo:
            k += 1; finish(score_frame(t), k, n_all)
        return done
    import concurrent.futures as cf, multiprocessing as mp
    for attempt in range(4):
        if not todo: break
        pending = list(todo); todo = []
        try:
            with cf.ProcessPoolExecutor(max_workers=procs, mp_context=mp.get_context('fork')) as ex:
                futs = {ex.submit(score_frame, t): t for t in pending}
                for fu in cf.as_completed(futs):
                    finish(fu.result(), k + 1, n_all); k += 1
        except cf.process.BrokenProcessPool:
            got = {(r['cohort'], r['subset'], r['frame']) for r in done}
            todo = [t for t in pending if (t[0], t[1], t[2]) not in got]
            print(f'  a worker died (out of memory?); retrying {len(todo)} unfinished tasks (attempt {attempt + 2})', flush=True)
    if todo: sys.exit(f'{len(todo)} tasks failed after retries')
    return done


DECISIONS = [
    'D1 Frame P(M)/A/S sample set: samples that have at least one estimate for every method of the frame (a failed tool-sample drops out of that frame, n_samples reports it).',
    'D2 Frame A validity fallback decided on the all-core-gene rows; the same frame definition (with or without kmerator) is then used for unflagged and flagged rows.',
    'D3 Noise SD / median |error| are computed on the sample-genes of the frame itself (truth>=50, B6/nonB6, estimate present for every method of the frame), with the cohort calibration (B6 rule Keele, non-hit rule Munger).',
    'D4 Coverage universe: Keele = core genes with a hybrid value and a unique symbol and truth (core_all.tsv.gz); Munger = all E5 core-gene records with truth; any stratum.',
    'D5 Subset rows: genes absent from flag_table.tsv belong to neither subset (counted in subsets.n_not_in_flag_table); they are in the all-core rows.',
    'D6 Subsample tolerance "within 1%" read as relative to nominal: |n_sub/ENA_reads / nominal - 1| <= 0.01.',
    'D7 featureCounts assigned fraction = Assigned / sum of all statuses in featurecounts.summary.txt; countTags hits per read = sum of per-gene k-mer count sums (kmerator sum column) / total_reads.',
    'D8 Cohort-tissue median for the 0.5-2x checks is taken over all samples with the value present, flagged ones included.',
    'D9 Kmerator estimate = per-gene mean count over its specific k-mers (the mean column); a gene missing from the file has no estimate.',
    'D10 Missing integrity value counts as a flag (value missing), not as pass.',
]


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    f = sub.add_parser('frames'); f.add_argument('cohort', choices=list(COHORTS)); f.add_argument('--frames', required=True)
    f.add_argument('--tabdir', help='Keele: make_tables.py output dir'); f.add_argument('--sample-table', help='Keele kidney: sample table samples.tsv')
    f.add_argument('--scorer', default=SCORER_PATH)
    s = sub.add_parser('score'); s.add_argument('--frames', required=True); s.add_argument('--out', required=True)
    s.add_argument('--out-root', default=OUT_ROOT); s.add_argument('--scorer', default=SCORER_PATH); s.add_argument('--cohorts', nargs='*', choices=list(COHORTS))
    s.add_argument('--procs', type=int, default=2); s.add_argument('--skip-repro-gate', action='store_true'); s.add_argument('--allow-scorer-mismatch', action='store_true')
    s.add_argument('--label', default=''); s.add_argument('--resume', action='store_true', help='keep one JSON per finished frame under OUT/task_cache and reuse it on a re-run (inputs must be unchanged)')
    r = sub.add_parser('repro'); r.add_argument('--frames', required=True); r.add_argument('--out', required=True); r.add_argument('--scorer', default=SCORER_PATH)
    a = ap.parse_args()
    if a.cmd == 'frames':
        cmd_frames(a)
    elif a.cmd == 'score':
        cmd_score(a)
    else:
        load_scorer(a.scorer); os.makedirs(a.out, exist_ok=True); allok = True
        for name in COHORTS:
            if not os.path.exists(f'{a.frames}/{name}/records.tsv.gz'): continue
            c = load_cohort(name, a.frames); c.V = {}; _CTX[name] = c
            ok, rows, nrows, _ = repro_one(c); allok &= ok
            write_tsv(f'{a.out}/reproduction_{name}.tsv', pd.DataFrame(rows + nrows)); print(name, 'match' if ok else 'MISMATCH')
        sys.exit(0 if allok else 1)


if __name__ == '__main__':
    main()
