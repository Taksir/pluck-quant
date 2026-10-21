#!/usr/bin/env python3
"""Build the small source-data TSVs for Fig. 5, Figs S8, S9, S12 and Tables S8, S13, S14.

Reads only existing results under revision/ in the original working folder (nothing is re-counted or re-quantified,
the scorer is not called) and writes *_u5.tsv files to figures/data/. The figure scripts
(fig5.py, figS8.py, figS9.py, figS12.py) read only those TSVs.

Usage: python3 fig5_data_u5.py CDNA.fa.gz
  CDNA.fa.gz = Ensembl release 102 Mus_musculus.GRCm38.cdna.all.fa.gz (only needed to refit the plain
  baseline classifiers B1/B2 of audit_desk/scripts/flag_baseline.py, whose fitted scores were not stored).

Recomputations, each checked against the stored value (the check is written to fig5_checks_u5.tsv):
  - ROC curves of the flag score (1 - u_g) and of the stored compact/post-count scores, from
    reliability_flag/out/scored_genes.tsv.gz; AUC must equal stage2_metrics.tsv to 4 dp.
  - ROC curves of B1 and B2, refit on Keele liver exactly as audit_desk/scripts/flag_baseline.py does;
    AUC must match audit_desk/out/flag_baseline_auc.tsv to 3 dp.
  - Per-gene error variances in the read simulation (reliability_flag/scripts/e2_check.py rule); the pooled (2% trimmed
    mean) ratios must equal reliability_flag/out/e2_check.json to 3 dp. Adds a gene-bootstrap 95%
    interval for each pooled ratio (new; 1,000 resamples, seed 20261015).
"""
import sys, gzip, json, re
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score, roc_curve

# Folder constants below point to the original working folder; set them for your copy.
REV = '/mnt/project-files/revision'
RF = f'{REV}/reliability_flag/out'
AD = f'{REV}/audit_desk/out'
OUT = '/mnt/project-files/paper/figures/data'
CDNA = sys.argv[1]
COHORTS = ['keele_kidney', 'keele_liver', 'munger', 'gse199702_liver', 'gse199702_heart', 'gse199702_kidney']
checks = []

def hdr(fn, source, notes):
    lines = [f'# file: {fn}', f'# source: {source}', '# written by: fig5_data_u5.py']
    lines += [f'# note: {n}' for n in notes]
    return '\n'.join(lines) + '\n'

def write(df, fn, source, notes, fmt=None):
    with open(f'{OUT}/{fn}', 'w') as fh:
        fh.write(hdr(fn, source, notes))
        df.to_csv(fh, sep='\t', index=False, float_format=fmt)

S = pd.read_csv(f'{RF}/scored_genes.tsv.gz', sep='\t')
S['u_score'] = 1 - S['u']
st2 = pd.read_csv(f'{RF}/stage2_metrics.tsv', sep='\t')

# ---------------------------------------------------------------- B1 / B2 refit (flag_baseline.py logic)
best = {}
name = None; seq = []
def flush():
    if name is None: return
    s = ''.join(seq); g = re.search(r'gene:(ENSMUSG\d+)', name).group(1)
    if g not in best or len(s) > best[g][0]:
        best[g] = (len(s), (s.count('G') + s.count('C')) / max(len(s), 1))
with gzip.open(CDNA, 'rt') as fh:
    for line in fh:
        if line[0] == '>': flush(); name = line; seq = []
        else: seq.append(line.strip())
flush()
GL = pd.DataFrame([(g, v[0], v[1]) for g, v in best.items()], columns=['gene_id', 'len_longest', 'gc_longest'])
sys.path.insert(0, f'{REV}/reliability_flag/scripts')
from common_rf import design_features
DF = {lay: design_features(RF, lay) for lay in ('50SE', '100SE', '101SE')}
def table(cohort):
    s = S[S.cohort == cohort]; lay = s.layout.iloc[0]
    d = s.merge(DF[lay], on='gene_id', how='left', suffixes=('', '_x')).merge(GL, on='gene_id', how='left')
    d['logT'] = np.log10(d['T'].clip(lower=1)); d['logT'] = (d['logT'] - d['logT'].mean()) / d['logT'].std()
    d['loglen'] = np.log10(d['len_longest']); d['lognt'] = d['ntx']; d['logtsl1'] = d['ntsl1']
    return d
FEATS = {'B1': ['loglen', 'lognt', 'logtsl1', 'gc_longest', 'logT'],
         'B2': ['loglen', 'lognt', 'logtsl1', 'gc_longest', 'logT', 'frac_multicopy']}
mk = lambda: make_pipeline(SimpleImputer(strategy='median'), StandardScaler(), LogisticRegression(C=1.0, max_iter=2000))
liver = table('keele_liver')
models = {k: mk().fit(liver[c].astype(float), liver['tail'].astype(int).values) for k, c in FEATS.items()}
base_auc = pd.read_csv(f'{AD}/flag_baseline_auc.tsv', sep='\t')

# ---------------------------------------------------------------- ROC curves (Fig. 5a, Fig. S8)
GRID = np.unique(np.concatenate([np.linspace(0, 0.1, 101), np.linspace(0.1, 1, 181)]))
roc_rows, op_rows = [], []
STAGE2_NAME = {'U': 'primary_u', 'compact': 'compact_logit', 'postcount': 'postcount_logit'}
for c in COHORTS:
    d = table(c); y = d['tail'].astype(int).values
    scores = {'U': d['u_score'].values, 'compact': d['s_compact'].values, 'postcount': d['s_post'].values}
    for k, cols in FEATS.items():
        scores[k] = models[k].predict_proba(d[cols].astype(float))[:, 1]
    for k, sc in scores.items():
        ok = np.isfinite(sc)
        a = roc_auc_score(y[ok], sc[ok])
        fpr, tpr, _ = roc_curve(y[ok], sc[ok])
        # linear interpolation of the ROC (ties in the score give straight segments, as roc_curve returns them)
        uf = np.unique(fpr); ut = np.array([tpr[fpr == f].max() for f in uf])
        tg = np.interp(GRID, uf, ut)
        roc_rows += [dict(cohort=c, model=k, fpr=round(g, 4), tpr=round(t, 4)) for g, t in zip(GRID, tg)]
        if k in STAGE2_NAME:
            ref = float(st2[(st2.cohort == c) & (st2.method == STAGE2_NAME[k])].auc.iloc[0]); tol = 5e-5
        else:
            ref = float(base_auc[(base_auc.cohort == c) & (base_auc.model == k)].auc.iloc[0]); tol = 1.5e-3
        checks.append(dict(item=f'AUC {k} {c}', recomputed=round(a, 4), stored=ref, ok=abs(a - ref) <= tol))
        roc_rows.append(dict(cohort=c, model=k, fpr=-1, tpr=round(a, 4)))   # fpr=-1 row carries the AUC
    for thr in (0.6, 0.7):
        fl = (d['u'] < thr).values; tail = y.astype(bool)
        op_rows.append(dict(cohort=c, threshold=thr, n_genes=len(d), n_tail=int(tail.sum()), n_flagged=int(fl.sum()),
                            flag_rate=round(fl.mean(), 4), fpr=round((fl & ~tail).sum() / (~tail).sum(), 4),
                            tpr_recall=round((fl & tail).sum() / tail.sum(), 4), precision=round((fl & tail).sum() / fl.sum(), 4),
                            tail_all=round(tail.mean(), 4), tail_unflagged=round(tail[~fl].mean(), 4), tail_flagged=round(tail[fl].mean(), 4)))
roc = pd.DataFrame(roc_rows); ops = pd.DataFrame(op_rows)
for _, r in ops[ops.threshold == 0.6].iterrows():
    s2 = st2[(st2.cohort == r.cohort) & (st2.method == 'primary_u')].iloc[0]
    checks.append(dict(item=f'recall u<0.6 {r.cohort}', recomputed=r.tpr_recall, stored=s2.recall, ok=abs(r.tpr_recall - s2.recall) < 1e-4))
    checks.append(dict(item=f'tail_unflagged u<0.6 {r.cohort}', recomputed=r.tail_unflagged, stored=s2.tail_unflagged, ok=abs(r.tail_unflagged - s2.tail_unflagged) < 1e-4))
write(roc, 'fig5a_roc_curves_u5.tsv', 'reliability_flag/out/scored_genes.tsv.gz (U = 1 - u_g, compact = s_compact, postcount = s_post); B1/B2 refit on Keele liver per audit_desk/scripts/flag_baseline.py',
      ['ROC of each score for the tail (hybrid/Salmon-pan residual-variance ratio > 2), linearly interpolated onto a fixed FPR grid.',
       'Rows with fpr = -1 carry the AUC in the tpr column. keele_liver is the fitting cohort (in-sample for every fitted model; U has no fit).',
       'Models: U = flag score; compact = 8-feature design logistic model; postcount = compact + log T + observed branch; B1 plain (length, transcripts, TSL1, GC, expression); B2 = B1 + genome multi-copy fraction.'])
write(ops, 'fig5a_operating_points_u5.tsv', 'reliability_flag/out/scored_genes.tsv.gz',
      ['Flag u_g < 0.6 (primary) and u_g < 0.7 (sensitivity) on the genes scored in each cohort-tissue; must reproduce stage2_metrics.tsv.'])

# ---------------------------------------------------------------- threshold sweep (Fig. S8b)
sw = []
for c in COHORTS:
    d = S[S.cohort == c]; tail = d['tail'].values.astype(bool); u = d['u'].values
    for thr in np.round(np.arange(0.30, 0.951, 0.025), 3):
        fl = u < thr
        sw.append(dict(cohort=c, threshold=thr, flag_rate=round(fl.mean(), 4), recall=round((fl & tail).sum() / tail.sum(), 4),
                       precision=round((fl & tail).sum() / max(fl.sum(), 1), 4), tail_unflagged=round(tail[~fl].mean(), 4)))
write(pd.DataFrame(sw), 'figS8_threshold_sweep_u5.tsv', 'reliability_flag/out/scored_genes.tsv.gz',
      ['Flag rate, recall, precision and tail rate among unflagged genes as the u_g threshold varies; 0.6 is the primary rule, 0.7 the sensitivity rule. Other thresholds are descriptive only.'])

# ---------------------------------------------------------------- per-gene ratios (Fig. 5b)
pg = S[['cohort', 'gene_id', 'ratio', 'flag', 'branch_em']].copy()
pg['flag'] = pg['flag'].astype(int); pg['branch_em'] = pg['branch_em'].astype(int)
write(pg, 'fig5b_pergene_resvar_ratio_u5.tsv', 'reliability_flag/out/scored_genes.tsv.gz',
      ['Per-gene residual-variance ratio hybrid / Salmon-pan from the eQTL scans (cis-eQTL residuals; strain-ANOVA residuals for Munger); tail = ratio > 2. flag = u_g < 0.6.'],
      fmt='%.5g')

# ---------------------------------------------------------------- simulation per-gene error variance (Fig. 5c)
D = np.load(f'{REV}/e2_sim/work/design.npz'); NG = len(D['genes']); NS = 8
T = np.stack([np.bincount(D['gene_idx'], weights=D['expected'][s], minlength=NG) for s in range(NS)])
B6 = D['hap'].astype(int) == 0
feat = {lay: pd.read_csv(f'{RF}/{fn}', sep='\t').set_index('gene_id').frac_emg.reindex(D['genes']).to_numpy()
        for lay, fn in (('se50', 'feat26.tsv.gz'), ('se101', 'feat76.tsv.gz'), ('pe', 'feat76.tsv.gz'))}
e2c = json.load(open(f'{RF}/e2_check.json'))
LAYLAB = {'se50': '50 bp SE', 'se101': '101 bp SE', 'pe': '101 bp PE'}
tr = lambda x: np.sort(x)[int(.02 * len(x)): int(.98 * len(x))].mean()
rng = np.random.default_rng(20261015)
gene_rows, pooled_rows = [], []
for lay in ('se50', 'se101', 'pe'):
    E = np.load(f'{REV}/e2_sim/out/est/clean_{lay}.npz')
    core = np.isfinite(E['hybrid']).all(0); ok = core & np.all(T >= 50, 0)
    V = {}
    for m in ('hybrid', 'salmon_pan', 'salmon_ref'):
        X = E[m]
        C = np.array([np.median(T[s, B6[s] & core & (T[s] >= 50) & (X[s] > 0)] / X[s, B6[s] & core & (T[s] >= 50) & (X[s] > 0)]) for s in range(NS)])
        V[m] = np.log2((C[:, None] * X + 1) / (T + 1)).var(0, ddof=1)
    good = ok & np.isfinite(V['hybrid']) & np.isfinite(V['salmon_pan'])
    fl = ~(feat[lay] >= 0.6)
    for i in np.where(good)[0]:
        gene_rows.append(dict(layout=LAYLAB[lay], gene_id=D['genes'][i], flag=int(fl[i]), var_hybrid=V['hybrid'][i], var_salmon_pan=V['salmon_pan'][i], var_salmon_ref=V['salmon_ref'][i]))
    for sub, msk in (('all', good), ('unflagged', good & ~fl), ('flagged', good & fl)):
        vh, vp, vr = V['hybrid'][msk], V['salmon_pan'][msk], V['salmon_ref'][msk]
        ph, pr = tr(vh) / tr(vp), tr(vr) / tr(vp)
        n = len(vh); bh, br = [], []
        for _ in range(1000):
            b = rng.integers(0, n, n); bh.append(tr(vh[b]) / tr(vp[b])); br.append(tr(vr[b]) / tr(vp[b]))
        pooled_rows.append(dict(layout=LAYLAB[lay], subset=sub, n_genes=n, ratio_hybrid=round(ph, 3), ratio_hybrid_lo=round(np.percentile(bh, 2.5), 3), ratio_hybrid_hi=round(np.percentile(bh, 97.5), 3),
                                ratio_salmon_ref=round(pr, 3), ratio_salmon_ref_lo=round(np.percentile(br, 2.5), 3), ratio_salmon_ref_hi=round(np.percentile(br, 97.5), 3),
                                frac_genes_hybrid_ratio_gt2=round(float(np.mean(vh / vp > 2)), 4)))
        ref = e2c[lay][sub]
        checks.append(dict(item=f'E2 pooled ratio hybrid {lay} {sub}', recomputed=round(ph, 3), stored=ref['pooled_var_ratio_hyb_pan'], ok=abs(round(ph, 3) - ref['pooled_var_ratio_hyb_pan']) < 1.5e-3))
        checks.append(dict(item=f'E2 pooled ratio salmon_ref {lay} {sub}', recomputed=round(pr, 3), stored=ref['salmon_ref_pooled_var_ratio'], ok=abs(round(pr, 3) - ref['salmon_ref_pooled_var_ratio']) < 1.5e-3))
        checks.append(dict(item=f'E2 n {lay} {sub}', recomputed=n, stored=ref['n'], ok=n == ref['n']))
write(pd.DataFrame(gene_rows), 'fig5c_e2_pergene_errvar_u5.tsv', 'e2_sim/work/design.npz, e2_sim/out/est/clean_{se50,se101,pe}.npz, reliability_flag/out/feat26.tsv.gz, feat76.tsv.gz',
      ['Per-gene within-gene variance (8 simulated samples) of the calibrated log2 error against the expected-count truth (reliability_flag/scripts/e2_check.py rule). Core genes with >= 50 expected reads in every sample. flag = u_g < 0.6 on the layout design.'],
      fmt='%.6g')
write(pd.DataFrame(pooled_rows), 'fig5c_e2_pooled_ratio_ci_u5.tsv', 'as fig5c_e2_pergene_errvar_u5.tsv',
      ['Pooled ratio = 2%-trimmed mean of per-gene variance (method) / same for Salmon-pan; reproduces e2_check.json. Gene-bootstrap 95% interval (1,000 resamples, seed 20261015) is added here.'])

# ---------------------------------------------------------------- SV check (Fig. S9, Table S14)
seg = pd.read_csv(f'{AD}/sv_probe_status_by_segment.tsv', sep='\t')
write(seg, 'figS9a_sv_probe_status_u5.tsv', 'audit_desk/out/sv_probe_status_by_segment.tsv (copied)',
      ['Fraction of single-copy (one B6 hit) probes absent (founder count 0) or extra (count >= 2) in each founder assembly; N = C57BL/6NJ artifact floor. Founders: A A/J, C 129S1/SvImJ, D NOD/ShiLtJ, E NZO/HlLtJ, F CAST/EiJ, G PWK/PhJ, H WSB/EiJ.'])
sve = pd.read_csv(f'{AD}/sv_vs_e5.tsv', sep='\t')
write(sve, 'figS9b_sv_vs_e5_u5.tsv', 'audit_desk/out/sv_vs_e5.tsv (copied)',
      ['Munger founder mice: SV-affected genes (>= 25% or >= 50% excess probe loss in that founder) against all other genes; gene-wise median over the founder\'s mice. method vf = hybrid (scorer code), ref = ref-unique probes.'])
svf = pd.read_csv(f'{AD}/sv_core_gene_flags.tsv', sep='\t')
ft = pd.read_csv(f'{RF}/flag_table.tsv', sep='\t')
sv = svf.merge(ft[['gene_id', 'gene_name', 'chrom', 'unique_frac_50bp', 'flag_50bp', 'unique_frac_100bp', 'flag_100bp']], on='gene_id', how='left')
sv['founders_ge50'] = sv['founders_ge50'].fillna('')
cols = ['gene_id', 'gene_name', 'chrom', 'founders_ge25', 'founders_ge50'] + [c for c in svf.columns if c.startswith('excess_absent_max_')] + ['unique_frac_50bp', 'flag_50bp', 'unique_frac_100bp', 'flag_100bp']
sv = sv[cols].sort_values(['chrom', 'gene_id'], key=lambda s: s.astype(str))
write(sv, 'tableS14_sv_risk_genes_u5.tsv', 'audit_desk/out/sv_core_gene_flags.tsv joined to reliability_flag/out/flag_table.tsv',
      ['567 core genes that lose at least 25% of their checkable probes (above the C57BL/6NJ value for the same gene) in at least one founder assembly; max over the four core estimator units. Built after scoring; reporting only, not used in any analysis.',
       'excess_absent_max_X = largest excess absent fraction over the four units in founder X (A A/J, C 129S1, D NOD, E NZO, F CAST, G PWK, H WSB). Reliability-flag columns added for cross-reference.'])
ov = []
nfl50 = int(sv.flag_50bp.astype(bool).sum()); nfl100 = int(sv.flag_100bp.astype(bool).sum())
core_rate50 = ft.flag_50bp.astype(bool).mean(); core_rate100 = ft.flag_100bp.astype(bool).mean()
ov.append(dict(set='SV-risk genes (>= 25%)', n=len(sv), flagged_50bp=nfl50, flagged_50bp_frac=round(nfl50 / len(sv), 4), flagged_100bp=nfl100, flagged_100bp_frac=round(nfl100 / len(sv), 4)))
ov.append(dict(set='all core genes', n=len(ft), flagged_50bp=int(ft.flag_50bp.astype(bool).sum()), flagged_50bp_frac=round(core_rate50, 4), flagged_100bp=int(ft.flag_100bp.astype(bool).sum()), flagged_100bp_frac=round(core_rate100, 4)))
write(pd.DataFrame(ov), 'figS9c_sv_flag_overlap_u5.tsv', 'tableS14_sv_risk_genes_u5.tsv; reliability_flag/out/flag_table.tsv',
      ['How many SV-risk genes the design-only reliability flag also marks (descriptive).'])
checks.append(dict(item='SV-risk genes (>=25%)', recomputed=len(sv), stored=567, ok=len(sv) == 567))
checks.append(dict(item='core flag rate 50bp', recomputed=round(core_rate50, 3), stored=0.149, ok=abs(core_rate50 - 0.149) < 6e-4))
checks.append(dict(item='core flag rate 100bp', recomputed=round(core_rate100, 3), stored=0.196, ok=abs(core_rate100 - 0.196) < 6e-4))

# ---------------------------------------------------------------- Table S13 (flag table)
s13 = ft[['gene_id', 'gene_name', 'chrom', 'unique_frac_50bp', 'flag_50bp', 'unique_frac_100bp', 'flag_100bp', 'frac_shared_50bp', 'frac_multicopy_50bp',
          'n_emg_50bp', 'ntx_50bp', 'gene_probes_50bp', 'gene_probes_101bp']].copy()
s13['sv_risk_ge25'] = s13.gene_id.isin(set(sv.gene_id)).astype(int)
write(s13, 'tableS13_reliability_flag_u5.tsv', 'reliability_flag/out/flag_table.tsv (sha fff978bf...d506; columns selected) + SV-risk membership from audit_desk/out/sv_core_gene_flags.tsv',
      ['All 16,885 core-tier genes. unique_frac = u_g (share of the gene\'s candidate windows, over all its annotated transcripts and free of annotated founder SNPs and small indels, that are gene-unique and single-copy in GRCm38, Section 2.3); flag = u_g < 0.6 (threshold set on Keele liver, flag_spec.json sha ffa185a7...7e3c).',
       '50 bp columns use the 26-nt-spaced design; 100 bp columns use the 76-nt design (also used for 101 bp data). n_emg_50bp = EM probes; ntx = annotated transcripts; gene_probes = gene probes at 50 / 101 bp.'])

# ---------------------------------------------------------------- Table S8 (flagged / unflagged splits, both calibrations, every pipeline)
# Rows of scoring set A come from tableS8_calib_<cohort>_u5.tsv, written by tableS8_calib_u5.py (scorer code; both constant tables are computed
# before the calibration function is replaced). Its "as scored" rows are checked against scores/acc_scores/<cohort>/metrics.tsv set A below.
# RNA-Skim rows come from its own pairwise set P(RNA-Skim) (kmer_tools/results/score/<cohort>/metrics.tsv, as scored only: its per-sample outputs
# are not in the shared folder, so no all-gene rows exist for it).
TOOL = {'vf': 'PLUCK', 'salmon_pan': 'Salmon-pan', 'salmon_ref': 'Salmon-ref', 'star_fc': 'STAR+fC', 'kallisto_pan': 'kallisto-pan', 'kallisto_ref': 'kallisto-ref',
        'kmerator': 'Kmerator+countTags', 'rnaskim': 'RNA-Skim', 'ref': 'ref-unique probes'}
LAB = {'keele_kidney': ('Keele kidney', 50), 'keele_liver': ('Keele liver', 50), 'munger': ('Munger founders', 100)}
parts = []
for coh in ('keele_kidney', 'keele_liver', 'munger'):
    t = pd.read_csv(f'{OUT}/tableS8_calib_{coh}_u5.tsv', sep='\t', comment='#')
    st = pd.read_csv(f'{REV}/final_comparison/scores/acc_scores/{coh}/metrics.tsv', sep='\t'); st = st[st.frame == 'A']
    m = t[t.calibration == 'as scored'].merge(st, on=['subset', 'method'], suffixes=('', '_s'))
    for k in ('n_delta_genes', 'median_delta', 'delta_ci95', 'spurious_rate', 'within_gene_spearman', 'noise_sd'):
        nbad = int((m[k].astype(str) != m[k + '_s'].astype(str)).sum())
        checks.append(dict(item=f'Table S8 {coh} as scored {k}: rows differing from locked frame A (of {len(m)})', recomputed=nbad, stored=0, ok=nbad == 0))
    # flagged all-gene rows of the four pipelines against final_audit/star_flag_calibration.tsv (unaffected by the monkeypatch bug)
    sc = pd.read_csv(f'{REV}/final_audit/star_flag_calibration.tsv', sep='\t'); sc = sc[(sc.cohort == coh) & (sc.calibration == 'all-gene') & (sc.subset == 'flagged')]
    for r in sc.itertuples():
        x = t[(t.subset == 'flagged') & (t.calibration == 'all-gene') & (t.method == {'hybrid': 'vf'}.get(r.method, r.method))].iloc[0]
        checks.append(dict(item=f'Table S8 {coh} flagged all-gene {r.method} delta, spurious', recomputed=f'{x.median_delta}, {x.spurious_rate}', stored=f'{r.median_delta}, {r.spurious_rate}',
                           ok=abs(x.median_delta - r.median_delta) < 1e-4 and abs(x.spurious_rate - r.spurious_rate) < 1e-4))
    t['frame'] = 'A'; parts.append(t)
    rs = pd.read_csv(f'{REV}/final_comparison/kmer_tools/results/score/{coh}/metrics.tsv', sep='\t')
    rs = rs[(rs.frame == 'P_rnaskim') & (rs.method == 'rnaskim')]
    parts.append(pd.DataFrame(dict(cohort=coh, subset=rs.subset, calibration='as scored', method='rnaskim', n_delta_genes=rs.n_delta_genes, median_delta=rs.median_delta,
                                   delta_ci95=rs.delta_ci95, spurious_rate=rs.spurious_rate, within_gene_spearman=rs.within_gene_spearman, n_genes_spearman=rs.n_genes_spearman,
                                   noise_sd=rs.noise_sd, n_genes_noise=rs.n_genes_noise, frame='P(RNA-Skim)')))
s8 = pd.concat(parts, ignore_index=True)
s8.insert(1, 'cohort_label', s8.cohort.map(lambda c: LAB[c][0])); s8.insert(2, 'read_length_bp', s8.cohort.map(lambda c: LAB[c][1]))
s8['tool'] = s8.method.map(TOOL)
ORDER = ['PLUCK', 'Salmon-pan', 'Salmon-ref', 'STAR+fC', 'kallisto-pan', 'kallisto-ref', 'Kmerator+countTags', 'RNA-Skim', 'ref-unique probes']
s8['o1'] = s8.cohort.map({'keele_kidney': 0, 'keele_liver': 1, 'munger': 2}); s8['o2'] = s8.subset.map({'all': 0, 'unflagged': 1, 'flagged': 2})
s8['o3'] = s8.tool.map({t: i for i, t in enumerate(ORDER)}); s8['o4'] = s8.calibration.map({'as scored': 0, 'all-gene': 1})
s8 = s8.sort_values(['o1', 'o2', 'o3', 'o4']).drop(columns=['o1', 'o2', 'o3', 'o4'])
s8 = s8[['cohort', 'cohort_label', 'read_length_bp', 'subset', 'tool', 'method', 'calibration', 'frame', 'n_delta_genes', 'within_gene_spearman', 'n_genes_spearman',
         'noise_sd', 'n_genes_noise', 'median_delta', 'delta_ci95', 'spurious_rate']]
write(s8, 'tableS8_flag_splits_u5.tsv', 'figures/data/tableS8_calib_<cohort>_u5.tsv (tableS8_calib_u5.py: scorer on final_comparison/frames + results_ll/acc, set A); '
      'final_comparison/kmer_tools/results/score/<cohort>/metrics.tsv (set P_rnaskim, RNA-Skim rows)',
      ['as scored = per-sample constants fitted inside the subset being scored (as in the scorer); all-gene = the constants of the all-gene score. For subset all the two are the same.',
       'Scoring set A, all core genes with a PLUCK value; RNA-Skim rows come from its own set (PLUCK and the four existing methods plus RNA-Skim, protein-coding index), as scored only.',
       'Noise SD and spurious rate are against the cohort truths, which share reads with full-read methods. Keele liver is the tuning cohort for the flag. kallisto-pan was run on Munger only.'])

# ---------------------------------------------------------------- arm 4 pilot (Fig. S12)
A4 = f'{REV}/arm4_pilot/results'
acc = pd.read_csv(f'{A4}/accuracy.tsv', sep='\t')
acc = acc[acc.method.isin(['hybrid_stored', 'arm4', 'salmon_pan_full'])].copy()
SAMPLE = {'SRR10523477': 'Keele kidney CC001', 'SRR10523526': 'Keele liver CC021', 'SRR826311': 'Munger A/J', 'SRR826340': 'Munger CAST'}
acc['sample'] = acc.srr.map(SAMPLE)
write(acc[['cohort', 'srr', 'sample', 'subset', 'method', 'n_genes', 'median_abs_lr', 'mad_sd_lr', 'median_lr', 'delta_hit_minus_nonhit']], 'figS12a_arm4_accuracy_u5.tsv',
      'arm4_pilot/results/accuracy.tsv (rows hybrid_stored, arm4, salmon_pan_full)',
      ['Per sample against the cohort truths (4 samples; no cohort noise SD). hybrid_stored = PLUCK as scored; arm4 = PLUCK with a Salmon-pan second pass on reads carrying a flagged-gene probe; salmon_pan_full = Salmon-pan on all reads.',
       'Accuracy rule S1: arm 4 has lower median |lr| and lower MAD-SD than PLUCK on flagged genes in >= 3 of 4 samples. Samples meeting it: 0 of 4.'], fmt='%.5g')
tm = pd.read_csv(f'{A4}/timing_per_rep.tsv', sep='\t')
g = tm.groupby(['cohort', 'srr'])
tt = pd.DataFrame(dict(wall_H=g.wall_H.median(), wall_A4=g.wall_A4.median(), wall_SR=g.wall_SR.median(),
                       ratio_SR_over_A4=g.ratio_SR_over_A4.median(), ratio_SR_over_A4_min=g.ratio_SR_over_A4.min(), ratio_SR_over_A4_max=g.ratio_SR_over_A4.max(),
                       ratio_SR_over_H=g.ratio_SR_over_H.median())).reset_index()
tt['sample'] = tt.srr.map(SAMPLE)
write(tt.round(3), 'figS12b_arm4_timing_u5.tsv', 'arm4_pilot/results/timing_per_rep.tsv (median and range of 3 repetitions)',
      ['Wall seconds on one 4-vCPU cloud VM, 3 shuffled repetitions per arm. H = PLUCK alone, A4 = arm 4 whole pipeline, SR = Salmon-ref. Speed rule: SR/A4 >= 10 at 4 threads. Samples meeting it: 0 of 4 (1.2 to 2.1x).'])
for _, r in tt.iterrows():
    checks.append(dict(item=f'arm4 SR/A4 median {r.srr}', recomputed=round(r.ratio_SR_over_A4, 2), stored={'SRR10523477': 2.09, 'SRR10523526': 1.73, 'SRR826311': 1.16, 'SRR826340': 1.25}[r.srr],
                       ok=abs(round(r.ratio_SR_over_A4, 2) - {'SRR10523477': 2.09, 'SRR10523526': 1.73, 'SRR826311': 1.16, 'SRR826340': 1.25}[r.srr]) < 0.006))

# ---------------------------------------------------------------- checks
ck = pd.DataFrame(checks)
write(ck, 'fig5_checks_u5.tsv', 'this script', ['Every recomputed value against its stored source value.'])
print(ck.to_string()); print('ALL OK' if ck.ok.all() else 'SOME CHECKS DO NOT MATCH')
