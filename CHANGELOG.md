# Changelog

## 0.1.0
- First packaging of the frozen v2 design: streaming counter (with the fix for gzip files over 4 GiB), EM, per-sample and cohort switch, panel compile for reverse, forward and unstranded libraries, design verification against the locked hashes, design rebuild and download, optional pymsbwt backend.
- Sources taken byte for byte from the frozen code: `em.c`, `em.py`, `e9rs.py`, `design.py` (two edits at run time: spacing and probe path, as in the Longleaf rebuild), `gc.py`, `gscan.py`, `gfix.py`; `pcount.c` is the 4 GiB-fixed file (sha256 a67f0d70...).
- Design rebuild from public data (`pluck design build --layout 50`) tested in a cloud container: 88 min prep (downloads) + 27 min design stages; probes_s26.npz, genes.json, design.npz arrays, gc.npy, genome_count.npy and the compiled panel all equal the locked values.
