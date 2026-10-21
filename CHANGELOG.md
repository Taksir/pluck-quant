# Changelog

## 1.0.1 (release v1.0.1-submission)
- All design files of the paper (50, 100 and 101 bp) are assets of this one release; `pluck design fetch --url .../releases/download/v1.0.1-submission` works for every layout.
- New `analysis/` folder: the scorers, the comparison and RNA-Skim scoring scripts, the comparator runners, the sample-identity check, the figure and table scripts, their source tables and `Table2.tsv`.
- `reference_code/` (renamed from the earlier snapshot folder) keeps the original probe-design, counting, EM and extended-tier scripts with a checksum list. The planning notes that sat there are removed.
- Wording: the README, the docs, the help text, the code comments and the table headers use plain descriptions. `pluck design build/verify/fetch` print "all match" or "DIFFERS" for the reference hashes.
- The package code (counter, EM, estimators, switches), the design files and every count are unchanged. `src/pluck_quant/_frozen/` is now `src/pluck_quant/_core/`. The only edit in `e9rs.py` is its first docstring line.

## 1.0.0 (release v1.0-submission)
- Install from GitHub (`pip install git+https://github.com/Taksir/pluck-quant@v1.0-submission`); PyPI after the paper is accepted.
- `--switch sample-scaled`: the rule of the timed runs of the paper (f = N_g / EM_g per sample). `--switch cohort` (accuracy results) and `--switch sample` are unchanged.
- Layout 100 (100 bp reads, probe spacing 76, the Munger design): recipe (`pluck design build --layout 100`), reference hashes, tie-order table, fetch, verify. Layouts 101 and 50 as before. One `prep` run now writes the probe files of all three layouts.
- `pluck design fetch` accepts a flat folder of files (GitHub release assets, with `design_<layout>.npz`, `gc_<layout>.npy`, `genome_count_<layout>.npy`).
- The reliability flag table (`flag_table.tsv`), the gene tiers and `flag_spec.json` ship with the package: `--flags builtin:flags`, `--tiers builtin:tiers`. The flag feature script is in `flag/`.
- The counting backend is chosen from the input: a sample with an msBWT is counted from it, otherwise the FASTQ files are streamed; without pymsbwt it falls back to the FASTQ files given with the msBWT (with a warning), or stops with a clear error. `--bwt` and `--stream` force one backend.
- Tests: `--switch sample-scaled` against the stored timed outputs.

## 0.1.0
- First packaging: streaming counter (with the fix for gzip files over 4 GiB), EM, per-sample and cohort switch, panel compile for reverse, forward and unstranded libraries, design verification against the reference hashes, design rebuild from public data, optional pymsbwt backend.
- Sources taken byte for byte from the original code: `em.c`, `em.py`, `e9rs.py`, `design.py` (two edits at run time: spacing and probe path), `gc.py`, `gscan.py`, `gfix.py`; `pcount.c` is the 4 GiB-fixed file (sha256 a67f0d70...).
- Design rebuild from public data (`pluck design build --layout 50`) tested in a cloud container: 88 min prep (downloads) + 27 min design stages; probes_s26.npz, genes.json, design.npz arrays, gc.npy, genome_count.npy and the compiled panel all equal the reference files.
