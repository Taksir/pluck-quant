"""End-to-end tests against outputs produced by the frozen v2 code (see tests/data/expected_mini)."""
import gzip
import hashlib
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from pluck_quant import count, estimate, native, output, panel as panel_mod

DATA = Path(__file__).parent / "data"
EXP = DATA / "expected_mini"
GENES = [l.strip() for l in open(DATA / "mini_design" / "genes_selected.txt") if l.strip()]


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def load_tsv(p):
    d = {}
    op = gzip.open if str(p).endswith(".gz") else open
    for l in list(op(p, "rt"))[1:]:
        t = l.rstrip("\n").split("\t")
        d[t[0]] = (None if t[1] == "NA" else float(t[1]), t[2])
    return d


def test_vendored_sources_are_the_locked_ones():
    lock = {}
    for l in (Path(__file__).parent.parent / "frozen_v2_snapshot" / "LOCK.sha256").read_text().splitlines():
        if l.startswith("#") or not l.strip():
            continue
        h, f = l.split()
        lock[Path(f).name] = h
    here = Path(panel_mod.__file__).parent
    assert sha(here / "_frozen" / "em.py") == lock["em.py"]
    assert sha(here / "_frozen" / "e9rs.py") == lock["e9rs.py"]
    assert sha(here / "native" / "em.c") == lock["em.c"]
    assert sha(here / "native" / "pcount.c") == "a67f0d70f394320843cbbd217f7a89aa838d13cc9e3998a110b39b0bada6c2fb"  # 4 GiB fix


def test_compile_matches_frozen(panel_dir):
    assert sha(panel_dir / "panel.txt") == (EXP / "panel.sha256").read_text().strip()
    got = panel_mod.array_hashes(panel_dir / "static.npz")
    want = {}
    for l in (EXP / "static_arrays.tsv").read_text().splitlines():
        f, k, dt, sh, h = l.split("\t")
        want[k] = h
    # arrhash.py wrote 0-d arrays as shape (1,); compare by name and hash
    assert {k[0]: h for k, h in got.items()} == want


@pytest.mark.parametrize("threads", [1, 4])
def test_counts_identical_to_frozen(panel_dir, tmp_path, threads):
    out = tmp_path / "c.txt"
    s = count.count_stream(panel_dir / "panel.txt", DATA / "reads" / "CC021_first500k.fastq.gz", out, threads=threads)
    assert sha(out) == (EXP / "counts.sha256").read_text().strip()
    assert s["reads"] == 500000


def test_counts_via_stdin_path(panel_dir, tmp_path, monkeypatch):
    real = native.build()
    monkeypatch.setattr(native, "build", lambda force=False: {**real, "isal": False})
    out = tmp_path / "c.txt"
    s = count.count_stream(panel_dir / "panel.txt", DATA / "reads" / "CC021_first500k.fastq.gz", out, threads=2)
    assert "stdin" in s["counter"]
    assert sha(out) == (EXP / "counts.sha256").read_text().strip()


@pytest.mark.parametrize("version", ["A"])
def test_estimates_match_frozen(panel_dir, counts_file, tmp_path, version):
    pan = estimate.Panel(panel_dir)
    res = estimate.sample_estimates(pan, estimate.read_counts(pan, counts_file))
    est, probe = estimate.hybrid_sample(res)
    out = tmp_path / "t.tsv"
    output.write_table(out, pan, ["CC021"], est[None, :], probe[None, :], branch=True)
    got = {}
    for l in list(open(out))[1:]:
        t = l.rstrip("\n").split("\t")
        got[t[0]] = (None if t[2] == "NA" else float(t[2]), t[3])
    want = load_tsv(EXP / "estimates_A.tsv")
    assert set(want) <= set(got)
    for g, (e, b) in want.items():
        ge, gb = got[g]
        if e:  # zero counts: both estimates are 0 and the old reader labels such genes 'probe'
            assert gb == b, g
        assert (e is None) == (ge is None)
        if e is not None:
            assert abs(ge - e) <= 1e-9 * abs(e) + 1e-12, (g, e, ge)


def test_cohort_rule_reproduces_stored_hybrid():
    z = np.load(DATA / "cohort_fixture_keele_kidney.npz")
    res = [{"Sg": z["S"][i], "em": z["em"][i], "Ng": z["Ng"][i]} for i in range(z["S"].shape[0])]
    h, probe, f = estimate.hybrid_cohort(res)
    ok = ~np.isnan(z["hyb"])
    assert np.array_equal(h[ok], z["hyb"][ok])
    assert np.allclose(f, z["fo"], equal_nan=True, rtol=0, atol=0)


def test_unstranded_panel_counts(panel_dir, tmp_path):
    from pluck_quant import panel as pm
    out = tmp_path / "pan_u"
    pm.compile_panel(DATA / "mini_design" / "inputs", DATA / "mini_design" / "inputs" / "probes_s26.npz", DATA / "genes.json", 50, out, "unstranded")
    info = __import__("json").loads((out / "panel.json").read_text())
    n = info["n_probes"]
    assert info["panel_lines"] == 2 * n
    c = tmp_path / "c.txt"
    count.count_stream(out / "panel.txt", DATA / "reads" / "CC021_first500k.fastq.gz", c, threads=2)
    both = np.loadtxt(c, dtype=np.int64)
    # antisense half must equal the reverse-stranded counts of the same probes
    c_rev = tmp_path / "rev.txt"
    count.count_stream(panel_dir / "panel.txt", DATA / "reads" / "CC021_first500k.fastq.gz", c_rev, threads=2)
    assert np.array_equal(both[n:], np.loadtxt(c_rev, dtype=np.int64))
    pan = estimate.Panel(out)
    assert np.array_equal(estimate.read_counts(pan, c), both[:n] + both[n:])


def test_cli_run_end_to_end(panel_dir, tmp_path):
    out = tmp_path / "genes.tsv"
    subprocess.run([sys.executable, "-m", "pluck_quant", "run", "--panel", str(panel_dir), "--sample",
                    f"CC021={DATA / 'reads' / 'CC021_first500k.fastq.gz'}", "--branch", "-o", str(out)], check=True)
    rows = [l.split("\t") for l in out.read_text().splitlines()]
    assert rows[0][:2] == ["gene_id", "gene_name"] and rows[0][-2:] == ["CC021", "CC021.branch"]
    d = {r[0]: r for r in rows[1:]}
    want = load_tsv(EXP / "estimates_A.tsv")
    g = GENES[0]
    assert abs(float(d[g][-2]) - want[g][0]) <= 1e-9 * abs(want[g][0]) + 1e-12


def test_msbwt_backend_equals_stream(panel_dir, counts_file, tmp_path):
    pytest.importorskip("MUSCython.MultiStringBWTCython")
    import shutil
    exe = shutil.which("msbwt") or str(Path(sys.executable).parent / "msbwt")
    fq = tmp_path / "r.fq"
    with gzip.open(DATA / "reads" / "CC021_first500k.fastq.gz", "rb") as f:
        fq.write_bytes(f.read())
    bwt = tmp_path / "bwt"
    bwt.mkdir()
    subprocess.run([exe, "cffq", "-p", "2", "-u", "-c", str(bwt), str(fq)], check=True, capture_output=True)
    out = tmp_path / "bwt.counts.txt"
    count.count_msbwt(panel_dir / "panel.txt", bwt, out)
    assert sha(out) == (EXP / "counts.sha256").read_text().strip()


def test_several_files_add_up(panel_dir, tmp_path):
    gz = DATA / "reads" / "CC021_first500k.fastq.gz"
    out = tmp_path / "two.txt"
    s = count.count_stream(panel_dir / "panel.txt", [gz, gz], out, threads=2)
    one = np.loadtxt(DATA.parent.parent / "tests" / "data" / "expected_mini" / "counts.txt.gz", dtype=np.int64)
    assert s["reads"] == 1000000
    assert np.array_equal(np.loadtxt(out, dtype=np.int64), 2 * one)
