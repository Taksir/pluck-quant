# Changelog

## 1.0.0 (release v1.0-submission)
- Install from GitHub (`pip install git+https://github.com/Taksir/pluck-quant@v1.0-submission`); PyPI after the paper is accepted.
- `--switch sample-scaled`: the rule of the timed runs of the paper (version B, f = N_g / EM_g per sample). `--switch cohort` (accuracy results) and `--switch sample` (version A) are unchanged.
- Layout 100 (100 bp reads, probe spacing 76, the Munger design): recipe (`pluck design build --layout 100`), locked hashes, tie-order table, fetch, verify. Layout 101 and 50 as before. One `prep` run now writes the probe files of all three layouts.
- `pluck design fetch` accepts a flat folder of files (GitHub release assets, with `design_<layout>.npz`, `gc_<layout>.npy`, `genome_count_<layout>.npy`).
- The reliability flag table (`flag_table.tsv`), the gene tiers and `flag_spec.json` ship with the package: `--flags builtin:flags`, `--tiers builtin:tiers`. The flag feature script is in `flag/`.
- Tests: version B against the stored timed outputs.

## 0.1.0
- First packaging of the frozen v2 design: streaming counter (with the fix for gzip files over 4 GiB), EM, per-sample and cohort switch, panel compile for reverse, forward and unstranded libraries, design verification against the locked hashes, design rebuild and download, optional pymsbwt backend.
- Sources taken byte for byte from the frozen code: `em.c`, `em.py`, `e9rs.py`, `design.py` (two edits at run time: spacing and probe path, as in the Longleaf rebuild), `gc.py`, `gscan.py`, `gfix.py`; `pcount.c` is the 4 GiB-fixed file (sha256 a67f0d70...).
- Design rebuild from public data (`pluck design build --layout 50`) tested in a cloud container: 88 min prep (downloads) + 27 min design stages; probes_s26.npz, genes.json, design.npz arrays, gc.npy, genome_count.npy and the compiled panel all equal the locked values.
