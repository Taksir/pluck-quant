# How to verify this release candidate

The package has been tested only in one cloud container. Before anything is uploaded, run it somewhere else.

**In a fresh cloud session (or any Linux machine with Python 3.12 or newer and gcc):**

1. Put `pluck-quant-0.1.0-src.tar.gz` there and unpack it: `tar xzf pluck-quant-0.1.0-src.tar.gz && cd pluck-quant`.
2. `sudo apt-get install -y gcc libisal-dev` (libisal-dev is optional; without it the counter uses a gzip pipe).
3. `bash scripts/verify.sh`. Expected last line: `VERIFY OK`, after `10 passed, 1 skipped` (or `11 passed` on Python 3.14, where the pymsbwt test runs: `PY=python3.14 bash scripts/verify.sh`).

What it proves: the vendored sources equal the locked frozen v2 files; the compiled panel, the counts of 500,000 real reads and 300 gene estimates equal the outputs of the frozen code; the cohort switch reproduces the stored Keele kidney hybrid exactly; the wheel installs into a clean environment and runs; the pymsbwt backend gives the same counts as the streaming counter.

What it does not cover (done separately in the project cloud container, see the report): the full 50 bp design and a whole sample, and the rebuild of the design from public data.

**Optional, about 2 hours (mostly downloads), 4 GB disk, ran in a 15 GB-RAM container:** `pluck design build --layout 50 --work w50` rebuilds the design from Ensembl 102 and the Mouse Genomes Project and ends with `locked hashes: PASS` if it matches the design behind the paper.
