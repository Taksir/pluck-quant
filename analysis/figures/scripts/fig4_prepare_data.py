#!/usr/bin/env python3
"""Fig. 4c source data: per-gene paired difference, PLUCK minus ref-unique probes, in the real DO animals.

Reads revision/idea2_do/out/delta_per_gene_D1.tsv and delta_per_gene_D2.tsv (per-gene class-minus-D0 shift against
GEO's individualized RSEM, written by the DO scorer) and writes the per-gene differences
delta_hybrid - delta_ref to figures/data/fig4c_do_per_gene_hybrid_minus_refunique.tsv, so that fig4.py reads only
figures/data/. The difference is truth-free in the sense that the truth term cancels within a gene.
Check: the class medians must equal fig4c_do_paired_hybrid_minus_refunique.tsv: 0.1232 (D2) and 0.0609 (D1).
"""
import time
import pandas as pd

# Folder constants below point to the original working folder; set them for your copy.
REV = '/mnt/project-files/revision/idea2_do/out'
DATA = '/mnt/project-files/paper/figures/data'
OUT = f'{DATA}/fig4c_do_per_gene_hybrid_minus_refunique.tsv'


def main():
    parts = []
    for cl in ('D2', 'D1'):
        d = pd.read_csv(f'{REV}/delta_per_gene_{cl}.tsv', sep='\t')
        parts.append(pd.DataFrame(dict(cls=cl, gene=d.gene, diff=(d.delta_hybrid - d.delta_ref).round(4))))
    df = pd.concat(parts, ignore_index=True)
    u0 = pd.read_csv(f'{DATA}/fig4c_do_paired_hybrid_minus_refunique.tsv', sep='\t', comment='#').set_index('class')
    checks = []
    for cl, g in df.groupby('cls'):
        med = round(float((g['diff']).median()), 4)
        checks.append(f'{cl} n={len(g)} median={med} (stored {u0.loc[cl, "median_diff"]}) frac_above_0={round(float((g["diff"] > 0).mean()), 3)}')
    with open(OUT, 'w') as f:
        f.write('# file: fig4c_do_per_gene_hybrid_minus_refunique.tsv\n')
        f.write('# source: idea2_do/out/delta_per_gene_D1.tsv, delta_per_gene_D2.tsv (derived: per-gene delta_hybrid - delta_ref, rounded to 4 dp)\n')
        f.write(f'# written: {time.strftime("%Y-%m-%d", time.gmtime())} by fig4_prepare_data.py\n')
        for c in checks:
            f.write(f'# check: {c}\n')
        f.write('# note: 48 DO animals, liver, 100 bp SE; D2 homozygous divergent, D1 heterozygote-like; hybrid is plotted as PLUCK\n')
        df.rename(columns={'cls': 'class'}).to_csv(f, sep='\t', index=False)
    print('\n'.join(checks))
    pooled_frame_k()
    split_frame_a()


def pooled_frame_k():
    """Fig. 4a pooled columns: scored domestic / wild-derived split per tool in scoring set K (the only set holding all six tools,
    RNA-Skim and Kmerator+countTags included), copied from kmer_tools/results/score/results.json, subset 'all'."""
    import json
    src = '/mnt/project-files/revision/final_comparison/kmer_tools/results/score/results.json'
    lab = dict(ref='ref-unique probes', vf='hybrid', salmon_ref='Salmon-ref', salmon_pan='Salmon-pan', kmerator='Kmerator+countTags', rnaskim='RNA-Skim')
    coh = dict(keele_kidney=('Keele kidney', 50), keele_liver=('Keele liver', 50), munger=('Munger founders', 100))
    R = json.load(open(src))['cohorts']
    out = f'{DATA}/fig4a_pooled_domestic_wild_frameK.tsv'
    with open(out, 'w') as f:
        f.write('# file: fig4a_pooled_domestic_wild_frameK.tsv\n')
        f.write('# source: final_comparison/kmer_tools/results/score/results.json (set K, subset all; scorer fields median_delta, delta_wild, delta_domestic)\n')
        f.write(f'# written: {time.strftime("%Y-%m-%d", time.gmtime())} by fig4_prepare_data.py; values copied, not recomputed\n')
        f.write('# note: wild-derived = CAST, PWK, WSB; domestic = A/J, 129S1, NOD, NZO (score_e9.py WILD/DOMESTIC). In Keele the founders of a CC hit gene must all be wild or all domestic to enter a split.\n')
        f.write('cohort\tread_length_bp\tframe\tmethod\ttool\tmedian_delta\tdelta_domestic\tdelta_wild\tn_delta_genes\n')
        for k, (name, rl) in coh.items():
            m = R[k]['frames']['K']['all']['metrics']
            for code in ('vf', 'salmon_pan', 'salmon_ref', 'kmerator', 'rnaskim', 'ref'):
                x = m[code]
                f.write(f'{name}\t{rl}\tK\t{code}\t{lab[code]}\t{x["median_delta"]}\t{x["delta_domestic"]}\t{x["delta_wild"]}\t{x["n_genes"]}\n')
    print('wrote', out)



def split_frame_a():
    """Fig. 4b: domestic / wild-derived split of the median founder-haplotype shift per pipeline and cohort.
    Pipelines of scoring set A from final_comparison/scores/acc_scores/<cohort>/metrics.tsv (set A, subset all; scorer fields delta_domestic,
    delta_wild, median_delta, delta_ci95); RNA-Skim from set P_rnaskim, its own gene set, as in Table 2."""
    src = '/mnt/project-files/revision/final_comparison/scores/acc_scores'
    lab = dict(vf='hybrid', salmon_pan='Salmon-pan', salmon_ref='Salmon-ref', star_fc='STAR+fC', kallisto_pan='kallisto-pan',
               kallisto_ref='kallisto-ref', kmerator='Kmerator+countTags', ref='ref-unique probes')
    coh = dict(keele_kidney=('Keele kidney', 50), keele_liver=('Keele liver', 50), munger=('Munger founders', 100))
    rows = []
    for k, (name, rl) in coh.items():
        m = pd.read_csv(f'{src}/{k}/metrics.tsv', sep='\t')
        m = m[(m.frame == 'A') & (m.subset == 'all')]
        for _, r in m.iterrows():
            rows.append(dict(cohort=name, read_length_bp=rl, frame='A', method=r.method, tool=lab[r.method], median_delta=r.median_delta,
                             delta_ci95=r.delta_ci95, delta_domestic=r.delta_domestic, delta_wild=r.delta_wild, n_delta_genes=r.n_delta_genes))
    # RNA-Skim from set P_rnaskim, the scoring set of Table 2 and Fig. 3 (set K differs by 0.0001 in Keele liver)
    import json
    R = json.load(open('/mnt/project-files/revision/final_comparison/kmer_tools/results/score/results.json'))['cohorts']
    for k_, (name, rl) in coh.items():
        fr = R[k_]['frames']['P_rnaskim']['all']; x = fr['metrics']['rnaskim']
        rows.append(dict(cohort=name, read_length_bp=rl, frame='P_rnaskim', method='rnaskim', tool='RNA-Skim', median_delta=x['median_delta'],
                         delta_ci95=json.dumps(x['delta_ci95']), delta_domestic=x['delta_domestic'], delta_wild=x['delta_wild'], n_delta_genes=fr['n_delta_genes']))
    out = f'{DATA}/fig4b_domestic_wild_split.tsv'
    with open(out, 'w') as f:
        f.write('# file: fig4b_domestic_wild_split.tsv\n')
        f.write('# source: final_comparison/scores/acc_scores/*/metrics.tsv (set A, subset all) and, for RNA-Skim, kmer_tools/results/score/results.json set P_rnaskim (the scoring set of Table 2)\n')
        f.write(f'# written: {time.strftime("%Y-%m-%d", time.gmtime())} by fig4_prepare_data.py; values copied, not recomputed\n')
        f.write('# note: delta_domestic / delta_wild = median over hit genes whose carrying founders are all domestic (A/J, 129S1, NOD, NZO) or all wild-derived (CAST, PWK, WSB), score_e9.py; kallisto-pan was run on Munger only\n')
        pd.DataFrame(rows).to_csv(f, sep='\t', index=False)
    print('wrote', out, len(rows), 'rows')


if __name__ == '__main__':
    main()
