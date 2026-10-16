# Deviation: streaming counter fix for gzipped FASTQs over 4 GiB

2026-10-07, thread "Kidney accuracy drop check". Approved by Taksir ("kidney: fix and retire", 2026-10-07 22:36Z). Engineering bug fix only: no design, panel, parameter or scoring change. The locked files are untouched (`frozen_v2/LOCK.sha256` and `frozen/` still verify); the fixed counter sits beside them in `frozen_v2/fix_pcount_4gib/`.

## Hashes
| File | sha256 |
|---|---|
| Old, locked: `frozen_v2/code/pcount.c` (= `frozen/code/pcount.c`, `no_msbwt/scripts/pcount.c`) | `94c8cf2aa7727ff29062834802d621b24b7edced730a0aeac77d79059a5c05de` |
| **New: `frozen_v2/fix_pcount_4gib/pcount.c`** | **`a67f0d70f394320843cbbd217f7a89aa838d13cc9e3998a110b39b0bada6c2fb`** |

Build is unchanged: `gcc -O3 -pthread -DWITH_ISAL -o pcount pcount.c -lisal` (tested with gcc 13.3, isa-l 2.31). The same command line works: `pcount gz panel.txt threads out.txt reads.fastq.gz`.

## What was wrong
In `gz` mode the counter maps the whole file and set `ist.avail_in = st.st_size` once. isa-l's `avail_in` is a 32-bit field, so a gzipped FASTQ larger than 4 GiB was inflated for only `size mod 2^32` bytes; the counter then exited 0 with the counts truncated. In E4b this hit Kidney_CC045_F and Heart_CC003_F (`e4b_full/kidney_drop.md`). The `scan` (stdin) and `store` modes never had the bug.

## What changed (`pcount_4gib.diff`, 3 hunks, `gz_read` and `main` only)
- Before every `isal_inflate` call, `avail_in` is set to the bytes left in the mapping, capped at 1 GiB (`GZ_WIN`). The mapping is contiguous, so each window simply continues from `next_in`.
- The multi-member and end-of-input checks use the bytes left in the mapping (`gze - next_in`) instead of `avail_in`.
- `main` no longer sets `avail_in` from `st_size`.
The inflated byte stream is the same, so counts are identical for any input the old code read completely.

## Test evidence (`fix_pcount_4gib/tests.sh`, output `tests.out`, count hashes `counts.sha256`)
Inputs: the real E4b panel (7,577,419 25-mers, rebuilt with `e4b_subsample/longleaf/mkpanel_ll.py`) and a synthetic 101 bp SE FASTQ chunk (`gen.py`, 4M reads, seed 20261007, two real panel k-mers planted per read, ~1% reads with N), gzipped with `igzip -1 -T 2` as in the E4b runner. No new downloads.

| Test | gz size | Old counter | Fixed counter | Result |
|---|---|---|---|---|
| A: chunk | 0.29 GB | 4,000,000 reads | 4,000,000 reads | **old = fixed, byte-identical**; also = `scan` mode on `zcat` |
| B: 12 × chunk, one gzip member | 3.54 GB (spans 1 GiB windows) | 48,000,000 reads | 48,000,000 reads | **old = fixed, byte-identical**; = 12 × A exactly |
| C: 16 × chunk, one gzip member | 4.71 GB | **5,700,828 reads** (8.9% of hits) | 64,000,000 reads | **fixed = 16 × A exactly, and = `scan` mode on `zcat`** |
| D: 16 × chunk, 16 gzip members (`cat`) | 4.71 GB | 5,700,899 reads | 64,000,000 reads | fixed = 16 × A exactly |

The old counter's 5.70M reads in C is what the mod-2^32 rule predicts (0.42 GB past 4 GiB at 73.7 bytes/read = 5.70M), the same mechanism that truncated the two E4b samples. Speed: the fixed counter was no slower (A 5.6 s vs 7.8 s; B 44.5 s vs 50.4 s, single runs, 4 threads).

## Consequences
- Every locked result stays valid: all inputs except the two E4b samples were below 4 GiB, where the outputs are byte-identical.
- The two affected E4b samples are documented in `e4b_full/kidney_drop.md` and were outside every gate subset. They were not re-counted (no reruns requested).
- From now on, runs that may meet gz FASTQs over 4 GiB (final comparison, idea 1) should build `fix_pcount_4gib/pcount.c` and record its hash. Per the coordinator, the final-comparison prep thread records the new hash in `final_comparison/prereg/ADDENDUM.md`; this thread did not edit it.
