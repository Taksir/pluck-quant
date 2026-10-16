# Frozen design v2: v1 hybrid (core) + coverage tier (extended), locked 2026-10-05

v2 is v1 (`revision/frozen/`, LOCK.sha256 file hash `57173f68…6031`) with one addition: a **coverage tier** for genes the v1 hybrid has no probes for. The coordinator prepared it on its recommendation, pending the user's choice. If the user keeps v1, v2 is not used. v1 is untouched. No parameter was tuned on liver for v2: the tier uses the TA_S8_P set exactly as built and scored in `probe_broaden/results.md`.

## 1. Tiers (`gene_tiers.tsv`, one row per gene)
- **Core: 16,885 genes.** These are the Ensembl 102 TSL1 protein-coding genes (`genes.json`) that have at least one v1 probe: a 26 nt variant-free gene probe (16,495) or an `emg` EM probe.
  - Estimator: the **v1 hybrid, unchanged**. Gene probes where observed capture f ≥ 0.6, `emg` EM otherwise (v1 `design.md` §2).
- **Extended: 26,258 genes.** These are every other gene on chr1–19/X with at least one probe in the `TA_S8_P` variant-free set.
  - That covers all annotated transcripts of any biotype or TSL. Probes may be in only some isoforms and may be shared by ≤8 genes with ≤8 genome copies.
  - 139 of them are `genes.json` genes with no v1 probe. The rest lie outside the TSL1 protein-coding universe.
  - Estimator: the ridge-NNLS of `probe_broaden` (`code/s6_estimate.py`). E[c_p] = Σ_t n_pt λ_t over all transcripts, ridge 1e−6 × max column norm², exact NNLS by coordinate descent. N_g = Σ_t λ_t · max(L_t − R + 1, 25).
  - The NNLS is solved on the whole TA_S8_P panel; only extended genes take their value from it.
- Every output table carries the `tier` column, and every result is reported for core and extended separately (`criteria_v2.md`).

## 2. Probe panel and counter
- Panel = the v1 panel ∪ `probe_broaden/design.npz` probes in `set_TA_S8_P__vf`, all counted antisense in one pass with the v1 counter command (`pcount scan`).
- The 50 bp extended set is fixed by the hash of `probe_broaden/design.npz` (`inputs.sha256`).

## 3. 101 bp single-end reverse-stranded (GSE199702, E4b)
- **Core:** exactly v1 §3 (77 nt gene probes `probes_s77.npz`, `emg` re-designed at SP = 77, R = 101, τ = 0.6).
- **Extended:** `code/s1_tx.py` → `s2_genome.py` → `s3_records.py` → `s4_select.py` are re-run with only these constants changed:
  - `s4_select.py`: `S, R = 77, 101`, with the baseline file `probes_s77.npz` in place of `probes_s26.npz`;
  - `s6_estimate.py`: `R = 101`.
  - The rule `TA_S8_P` is unchanged.
  - The output `design.npz` hash is recorded before any GSE199702 count is made.
- **Strand:** mRNA-antisense, as in v1.
- **Gene join and filter:** as in v1, by Ensembl gene ID to the GBRS/EMASE matrix. A tier is assigned by the same rule at 101 bp.

## 4. Sample-identity check
Unchanged from v1 §4. It runs before any scoring, the labelled strain must be the top fit with fit ≥ 0.90, and failing samples are excluded from all methods.

## 5. Scoring
- Unchanged from v1 §5 (locked `score_e9.py`, `score_noise.py` noise SD), run once per tier.
- Gates apply to core genes only.
- Liver reference for the extended tier (post hoc, Keele, not a gate):
  - 6,036 added genes with truth;
  - within-gene Spearman vs truth 0.972 (686 genes with ≥10 samples with truth ≥ 50);
  - spurious 6.6%, Δ +0.003; the same probes without the variant filter give Δ −0.06;
  - Salmon-pan 0.971 on the 143 of these genes it covers.
