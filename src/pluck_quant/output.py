"""Write the gene table."""
import gzip
from pathlib import Path

import numpy as np


def _read_table(path, key="gene_id"):
    """Tab-separated table by gene id. path may be 'builtin:flags' / 'builtin:tiers' (the tables shipped with the package:
    flag_table.tsv, sha256 fff978bf...d506, and gene_tiers.tsv of the frozen design) or a file (plain or .gz)."""
    import csv
    from importlib import resources
    if str(path).startswith("builtin:"):
        name = {"builtin:flags": "flag_table.tsv.gz", "builtin:tiers": "gene_tiers.tsv.gz"}[str(path)]
        path = resources.files("pluck_quant") / "_data" / name
        f = gzip.open(path.open("rb"), "rt", newline="")
    else:
        f = gzip.open(path, "rt", newline="") if str(path).endswith(".gz") else open(path, newline="")
    with f:
        rd = csv.DictReader(f, delimiter="\t")
        return {r[key]: r for r in rd}


def write_table(out, panel, names, est, probe, tiers=None, flags=None, read_length=50, branch=False, digits=None):
    """est, probe: (n_samples, n_genes). Columns: gene_id, gene_name, tier, flag (if tables given), then one column per sample
    (and '<sample>.branch' = probe / em, if branch)."""
    tier_map = _read_table(tiers) if tiers else None
    flag_map = _read_table(flags) if flags else None
    fcol = "flag_50bp" if read_length <= 60 else "flag_100bp"
    ucol = "unique_frac_50bp" if read_length <= 60 else "unique_frac_100bp"
    head = ["gene_id", "gene_name"]
    if tier_map is not None:
        head.append("tier")
    if flag_map is not None:
        head += ["unique_frac", "flag"]
    for n in names:
        head.append(n)
        if branch:
            head.append(n + ".branch")
    op = gzip.open if str(out).endswith(".gz") else open
    with op(out, "wt") as f:
        f.write("\t".join(head) + "\n")
        for i, g in enumerate(panel.GID):
            row = [g, panel.names[i]]
            if tier_map is not None:
                row.append(tier_map[g]["tier"] if g in tier_map else "NA")
            if flag_map is not None:
                r = flag_map.get(g)
                row += [r[ucol], r[fcol]] if r else ["NA", "NA"]
            for j in range(len(names)):
                e = est[j, i]
                row.append("NA" if np.isnan(e) else (repr(float(e)) if digits is None else f"{e:.{digits}g}"))
                if branch:
                    row.append("NA" if np.isnan(e) else ("probe" if probe[j, i] else "em"))
            f.write("\t".join(row) + "\n")
