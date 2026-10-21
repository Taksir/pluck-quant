#!/usr/bin/env python3
"""Copy the values behind Figs S5-S7 and Tables S5-S6 from revision/ into small TSVs (values copied, not recomputed,
except columns marked "derived").

Reads (read-only):
  revision/e4b_full/results/e4b_results.tar.gz          tab/<tissue>/metrics_e4b.json (scorer compare block)
  revision/final_comparison/scores/acc_scores/<cohort>/metrics.tsv, claims.tsv          (primary scores)
  revision/final_comparison/scores/sensitivity_excl13/<cohort>/metrics.tsv, excluded13.txt
  revision/final_comparison/kmer_tools/results/score/<cohort>/metrics.tsv, claims.tsv  (RNA-Skim scoring sets, RK rules)
  revision/final_comparison/per_sample_switch/tables/*.tsv
  revision/final_comparison/r7/hits_per_read.tsv, identity_depth.tsv
Writes figures/data/figS5_e4b_tissue_metrics_u3.tsv, figS6_sensitivity_u3.tsv, figS6_hits_identity_u3.tsv,
figS7_switch_u3.tsv, figS7_switch_alltools_u3.tsv, figS7_switch_claims_u3.tsv,
figS7_switch_flagged_u3.tsv and tables/TableS5_full_metrics.tsv, TableS6_claim_rules.tsv.
Run: python3 -I figures/scripts/fig3_supp_extract.py
"""
import io, json, os, tarfile
import pandas as pd

# Folder constants below point to the original working folder; set them for your copy.
REV = '/mnt/project-files/revision'
FC = f'{REV}/final_comparison'
DATA = '/mnt/project-files/paper/figures/data'
TABLES = '/mnt/project-files/paper/tables'
COH = ['keele_kidney', 'keele_liver', 'munger']
CLAB = {'keele_kidney': 'Keele kidney, 50 bp', 'keele_liver': 'Keele liver, 50 bp', 'munger': 'Munger founders, 100 bp'}
MLAB = {'vf': 'hybrid', 'ref': 'ref-unique probes', 'salmon_pan': 'Salmon-pan', 'salmon_ref': 'Salmon-ref', 'star_fc': 'STAR+fC',
        'kallisto_ref': 'kallisto-ref', 'kallisto_pan': 'kallisto-pan', 'kmerator': 'Kmerator+countTags', 'rnaskim': 'RNA-Skim',
        'salmon_pan_p25': 'Salmon-pan 25%', 'salmon_pan_p10': 'Salmon-pan 10%', 'salmon_ref_p25': 'Salmon-ref 25%',
        'salmon_ref_p10': 'Salmon-ref 10%', 'hybrid': 'hybrid'}
STAMP = '# written by: fig3_supp_extract.py; values copied from the source unless a column says derived\n'


def write(path, df, src, notes=()):
    with open(path, 'w') as f:
        f.write(f'# file: {os.path.basename(path)}\n# source: {src}\n' + STAMP)
        for n in notes:
            f.write(f'# note: {n}\n')
        df.to_csv(f, sep='\t', index=False)
    print(path, len(df))


def e4b():
    rows = []
    with tarfile.open(f'{REV}/e4b_full/results/e4b_results.tar.gz') as t:
        for tis in ('liver', 'heart', 'kidney'):
            d = json.load(t.extractfile(f'tab/{tis}/metrics_e4b.json'))
            for m, x in d['locked_scorer']['compare'].items():
                if m.startswith('_'):
                    continue
                rows.append(dict(tissue=tis, read_length_bp=101, method=m, tool=MLAB[m], n_samples=x['n_samples'], n_delta_genes=x['n_genes'],
                                 median_delta=x['median_delta'], delta_ci95=str(x['delta_ci95']), spurious_rate=x['spurious_rate'],
                                 within_gene_spearman=x['within_gene_spearman'], n_genes_spearman=x['n_genes_spearman'],
                                 delta_wild=x['delta_wild'], delta_domestic=x['delta_domestic'], pooled_B6_med_log2=x['pooled_B6_med_log2']))
    df = pd.DataFrame(rows)
    write(f'{DATA}/figS5_e4b_tissue_metrics_u3.tsv', df, 'revision/e4b_full/results/e4b_results.tar.gz tab/<tissue>/metrics_e4b.json, key locked_scorer.compare',
          ['GSE199702, 101 bp SE, GBRS/EMASE truth (shares reads with full-read methods). Scored subsets: liver 53, heart 30, kidney 30 samples (Salmon-pan available).',
           'No noise ratio is taken from these files (GBRS truth shares reads; hybrid sits below its independent Poisson floor).'])


def metrics(root):
    out = []
    for c in COH:
        d = pd.read_csv(f'{root}/{c}/metrics.tsv', sep='\t')
        out.append(d)
    d = pd.concat(out)
    d.insert(1, 'cohort_label', d.cohort.map(CLAB))
    d.insert(d.columns.get_loc('method') + 1, 'tool', d.method.map(MLAB))
    return d


def sensitivity():
    p = metrics(f'{FC}/scores/acc_scores'); s = metrics(f'{FC}/scores/sensitivity_excl13')
    keep = ['cohort', 'cohort_label', 'frame', 'subset', 'method', 'tool', 'n_samples', 'n_delta_genes', 'median_delta', 'delta_ci95',
            'spurious_rate', 'within_gene_spearman', 'noise_sd']
    p = p[p.frame == 'A'][keep].assign(run='all libraries (primary)'); s = s[s.frame == 'A'][keep].assign(run='13 low-mapping libraries excluded')
    df = pd.concat([p, s])
    write(f'{DATA}/figS6_sensitivity_u3.tsv', df, 'final_comparison/scores/acc_scores/<cohort>/metrics.tsv and scores/sensitivity_excl13/<cohort>/metrics.tsv (set A, all methods)',
          ['Munger has no excluded library; its rows are identical in both runs. Sensitivity only; the primary result keeps all libraries.'])
    ex = pd.read_csv(f'{FC}/scores/sensitivity_excl13/excluded13.txt', sep='\t', header=None, names=['cohort', 'srr'])
    h = pd.read_csv(f'{FC}/r7/hits_per_read.tsv', sep='\t').assign(check='hits per read')
    i = pd.read_csv(f'{FC}/r7/identity_depth.tsv', sep='\t').assign(check='identity depth')
    hi = pd.concat([h[['check', 'cohort', 'srr', 'hits_per_read', 'median', 'ratio', 'flag']].rename(columns={'hits_per_read': 'value'}),
                    i[['check', 'cohort', 'srr', 'depth', 'median', 'ratio', 'flag']].rename(columns={'depth': 'value'})])
    key = set(zip(ex.cohort, ex.srr))
    hi['in_low_mapping_13'] = [(c, s) in key for c, s in zip(hi.cohort, hi.srr)]
    write(f'{DATA}/figS6_hits_identity_u3.tsv', hi, 'final_comparison/r7/hits_per_read.tsv, identity_depth.tsv; scores/sensitivity_excl13/excluded13.txt',
          ['Rule: hits per read within 0.5-2x the cohort-tissue median; identity depth >= 0.5x median. Stored Keele and Munger probe counts only, not the timed counter runs.',
           'in_low_mapping_13 (derived) marks the 13 libraries with low mapping rates in every full-read tool. The liver frame has no identity table.'])


def switch():
    t = f'{FC}/per_sample_switch/tables'
    h = pd.read_csv(f'{t}/headline_frameA_allgenes.tsv', sep='\t')
    rows = []
    for v in ('A', 'B'):
        for c in COH:
            m = pd.read_csv(f'{t}/metrics_{v}_{c}.tsv', sep='\t')
            m = m[(m.frame == 'A') & (m.subset == 'all')]
            for _, r in m.iterrows():
                rows.append(dict(cohort=c, variant=v, method=r['method'], tool=MLAB[r['method']], within_gene_spearman=r['within_gene_spearman'],
                                 noise_sd=r['noise_sd'], spurious_rate=r['spurious_rate'], median_delta=r['median_delta']))
    write(f'{DATA}/figS7_switch_u3.tsv', h, 'final_comparison/per_sample_switch/tables/headline_frameA_allgenes.tsv',
          ['primary = cohort-level switch (--switch cohort, scored); A = f per sample (--switch sample, released code); B = scaled f per sample (--switch sample-scaled, the code as timed).'])
    write(f'{DATA}/figS7_switch_alltools_u3.tsv', pd.DataFrame(rows), 'final_comparison/per_sample_switch/tables/metrics_{A,B}_<cohort>.tsv (set A, all genes)',
          ['Only the hybrid column differs between variants; other tools are listed to show they are unchanged.'])
    cl = pd.read_csv(f'{t}/claims_primary_vs_variants.tsv', sep='\t', keep_default_na=False)
    cl = cl[cl.claim.str.contains('kallisto_pan/per_cohort/munger/(?:ratio_hybrid_over_M|passes|gap_M_minus_hybrid)', regex=True)]
    write(f'{DATA}/figS7_switch_claims_u3.tsv', cl, 'final_comparison/per_sample_switch/tables/claims_primary_vs_variants.tsv (kallisto-pan Munger R1/R2 rows)',
          ['primary column = cohort switch (--switch cohort); per_sample column = --switch sample (A) or --switch sample-scaled (B). R2 limit 1.25; R1 limit 0.02.'])
    f = pd.read_csv(f'{t}/flagged_both_calibrations.tsv', sep='\t')
    f.insert(f.columns.get_loc('method') + 1, 'tool', f.method.map(MLAB))
    write(f'{DATA}/figS7_switch_flagged_u3.tsv', f, 'final_comparison/per_sample_switch/tables/flagged_both_calibrations.tsv',
          ['calibration own = constants fitted inside the subset (as scored); global = per-sample constants of the all-gene score.'])


def table_s5():
    p = metrics(f'{FC}/scores/acc_scores').assign(scorer='fc_gather_score.py 70bdc836...a9f975')
    k = metrics(f'{FC}/kmer_tools/results/score').assign(scorer='kt_score.py fbf80cee...acaaf')
    k = k[k.frame.isin(['P_rnaskim', 'P_rnaskim_z', 'K'])]
    df = pd.concat([p, k])
    os.makedirs(TABLES, exist_ok=True)
    write(f'{TABLES}/TableS5_full_metrics.tsv', df, 'final_comparison/scores/acc_scores/<cohort>/metrics.tsv; final_comparison/kmer_tools/results/score/<cohort>/metrics.tsv',
          ['Every scoring set (frame column) x subset x method of the scorer. A = all methods on common sample-genes; P_<tool> = the four existing methods plus one tool; '
           'S = reduced-read Salmon arms; B = the four existing methods (reproduction check); P_rnaskim, P_rnaskim_z (missing scored as 0) and K (with Kmerator) = RNA-Skim sets.',
           'method vf = hybrid (PLUCK); ref = ref-unique probes. Truths share reads with full-read methods. Subsets: all / unflagged / flagged (u_g < 0.6), calibrated within each subset.'])


def table_s6():
    a = pd.read_csv(f'{FC}/scores/acc_scores/claims.tsv', sep='\t', keep_default_na=False).assign(source='scores/acc_scores/claims.tsv')
    k = pd.read_csv(f'{FC}/kmer_tools/results/score/claims.tsv', sep='\t', keep_default_na=False).assign(source='kmer_tools/results/score/claims.tsv')
    df = pd.concat([a, k])
    df['value'] = df['value'].astype(str).str.strip('"')
    df.insert(0, 'rule', df.path.str.split('/').str[0])
    write(f'{TABLES}/TableS6_claim_rules.tsv', df, 'final_comparison/scores/acc_scores/claims.tsv; final_comparison/kmer_tools/results/score/claims.tsv',
          ['R1-R6 and RK1-RK3 as flattened by the scorers; decided on all-core-gene rows of set P_<tool> (R1-R3), P_rnaskim (RK1-RK2) and P_<X> (RK3). '
           'RK4 (speed) has no threshold rule; its values are in Table S7. hybrid = PLUCK.'])


if __name__ == '__main__':
    e4b(); sensitivity(); switch(); table_s5(); table_s6()
