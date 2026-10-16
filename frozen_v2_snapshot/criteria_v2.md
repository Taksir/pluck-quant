# Pre-registered pass rules for the frozen design (v2), locked 2026-10-05

The rules are copied from `strategy/plan_v3.md` §3 and `strategy/plan_v3_addendum_gse199702.md`. Where those texts and this file differ, this file governs once it is locked. The method definitions are in `design.md`. Each experiment is scored once, with `code/score_e9.py` and the noise SD of `code/score_noise.py`.

| Step | Data | Pass rule |
|---|---|---|
| **E2** simulation | T1.4 / re-spacing simulator, 6–8 mosaic samples, all 7 non-B6 founders, 10M fragments, 50 bp SE and 101 bp PE (dUTP), expected truth | Per layout, clean condition: hybrid within-gene r² ≥ 0.95; hybrid per-gene error ≤ 1.6× Salmon-pan; \|median Δ\| ≤ 0.02 on hit genes for every founder; spurious rate ≤ 2× Salmon-pan. A dirty condition (15% unspliced pre-mRNA, 101 PE) is reported, not gated. If the hybrid's r² drops by more than 0.02 while Salmon-pan's does not, the EM's sensitivity to pre-mRNA is reported. |
| **E3** cost table | One VM, one script: CC021 and CC001 at 50 bp SE, two E2 101 bp PE samples, and one GSE199702 101 bp SE sample (Liver CC001 F) | Hybrid ≥ 10× faster than Salmon-ref and ≥ 25× faster than Salmon-pan at 4 threads, on both layouts. If not, the headline multiple is whatever was measured. |
| **E4** Keele kidney | All 47 kidney samples (GSE140873), after the identity check | Within-gene Spearman gap to Salmon-pan ≤ 0.02; spurious rate ≤ 1.5× Salmon-pan's; \|median Δ\| ≤ 0.03; ref-unique probes still ≤ −0.08 (C1 form). |
| **E4b** GSE199702, 101 bp SE | 174 females (58 strains × heart/kidney/liver), after the identity check; truth = GEO GBRS/EMASE matrix | Per tissue: within-gene Spearman gap to Salmon-pan ≤ 0.02; spurious rate ≤ 1.5× Salmon-pan; \|median Δ\| ≤ 0.03; ref-unique probes still biased (C1 form, ≤ −0.08). The 174 males are an ungated replicate. |
| **E5** Munger 2014 founders | GSE45684, 8 strains × 2 mice, 100 bp SE unstranded liver; truth = GEO founder matrix | C1 form (vf \|Δ\| ≤ 0.03, ref-unique ≤ −0.08 on hit genes) and Spearman gap ≤ 0.02 to Salmon-pan across strains. Adaptations, fixed now: sum both probe orientations (unstranded); R − k + 1 = 76; join truth by Ensembl gene ID. |
| **E6** power refit | E2 + E4 (+ E4b) | Predicted vs observed within-gene SD on Keele liver and kidney (and GSE199702) within 20%. |

**Validity:**
- The scorer's own validity rules apply (≥ 30 samples, ≥ 200 genes).
- A step that fails validity is reported as invalid, not as pass or fail.
- Salmon-ref and Salmon-pan for E4 and E4b use the `e9_salmon/` recipe (Salmon 1.10, Ensembl 102 cDNA+ncRNA; pan = plus each founder's transcripts from MGP v5).

## v2 addition: gene tiers (frozen_v2 only)
- Every gate above is evaluated on **core** genes only (`gene_tiers.tsv`, tier = core). Core genes use the v1 hybrid unchanged, so v1 and v2 give identical gate results.
- **Extended** genes (tier = extended) are scored with the same scorer and metrics and reported **separately and ungated**, next to Salmon-pan where Salmon covers them. They are never pooled into a gate number.
- Joins to truth are by gene ID (E4b, E5) or symbol (Keele), as in v1. Genes outside the truth matrix are dropped in both tiers.
