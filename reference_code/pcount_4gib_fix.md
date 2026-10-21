# Counter fix for gzipped FASTQ files over 4 GiB

This note describes the one change made to the streaming counter `pcount.c` after the first analyses. The fixed file is `fix_pcount_4gib/pcount.c` (SHA-256 `a67f0d70f394320843cbbd217f7a89aa838d13cc9e3998a110b39b0bada6c2fb`). It is the counter in this package (`src/pluck_quant/native/pcount.c`). The earlier file is `code/pcount.c` (SHA-256 `94c8cf2aa7727ff29062834802d621b24b7edced730a0aeac77d79059a5c05de`). No design, panel, parameter or estimator changed.

Build command (unchanged): `gcc -O3 -pthread -DWITH_ISAL -o pcount pcount.c -lisal` (tested with gcc 13.3 and isa-l 2.31). Use: `pcount gz panel.txt threads out.txt reads.fastq.gz`.

## What was wrong
In `gz` mode the counter mapped the whole file and set `avail_in = st_size` once. In isa-l `avail_in` is a 32-bit field, so a gzipped FASTQ larger than 4 GiB was inflated for only `size mod 2^32` bytes. The counter then exited with code 0 and truncated counts. The `scan` (stdin) and `store` modes never had the fault.

## What changed (`fix_pcount_4gib/pcount_4gib.diff`, 3 hunks, in `gz_read` and `main`)
- Before every `isal_inflate` call, `avail_in` is set to the bytes left in the mapping, capped at 1 GiB. The mapping is contiguous, so each window continues from `next_in`.
- The multi-member and end-of-input checks use the bytes left in the mapping (`gze - next_in`) instead of `avail_in`.
- `main` no longer sets `avail_in` from `st_size`.

The inflated byte stream is the same, so counts are identical for every input the old code read completely.

## Tests (`fix_pcount_4gib/tests.sh`, output `tests.out`, count hashes `counts.sha256`)
Inputs: a panel of 7,577,419 25-mers and a synthetic 101 bp single-end FASTQ chunk made by `fix_pcount_4gib/gen.py` (4 million reads, seed 20261007, two panel k-mers planted per read, about 1% of reads with an N), gzipped with `igzip -1 -T 2`.

| Test | gz size | Old counter | Fixed counter | Result |
|---|---|---|---|---|
| A: chunk | 0.29 GB | 4,000,000 reads | 4,000,000 reads | old and fixed byte-identical; also equal to `scan` mode on `zcat` |
| B: 12 chunks, one gzip member | 3.54 GB (spans 1 GiB windows) | 48,000,000 reads | 48,000,000 reads | old and fixed byte-identical; equal to 12 times A |
| C: 16 chunks, one gzip member | 4.71 GB | 5,700,828 reads | 64,000,000 reads | fixed equals 16 times A and `scan` mode on `zcat` |
| D: 16 chunks, 16 gzip members | 4.71 GB | 5,700,899 reads | 64,000,000 reads | fixed equals 16 times A |

The 5.70 million reads of the old counter in C is what the mod-2^32 rule predicts (0.42 GB past 4 GiB at 73.7 bytes per read). The fixed counter was not slower (A 5.6 s against 7.8 s; B 44.5 s against 50.4 s; single runs, 4 threads).

## Effect on results
All inputs of the paper except two Keele-cohort tissue samples were below 4 GiB, where old and fixed outputs are byte-identical. Those two samples were not used in any reported accuracy result.
