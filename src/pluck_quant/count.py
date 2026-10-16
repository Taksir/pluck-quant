"""Count the panel 25-mers in reads: streaming counter (default) or an existing msBWT (optional)."""
import gzip
import json
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from . import native


def _decompressor():
    for exe in ("igzip", "pigz", "gzip"):
        p = shutil.which(exe)
        if p:
            return [p, "-dc"]
    return None


def _feed(files, stdin):
    """Write the decompressed reads of all files to stdin (used when the counter cannot inflate gzip itself), then close it."""
    dec = _decompressor()
    try:
        for f in files:
            if str(f).endswith(".gz"):
                if dec:
                    r = subprocess.run(dec + [str(f)], stdout=stdin)
                    if r.returncode:
                        raise RuntimeError(f"{dec[0]} failed on {f}")
                else:
                    with gzip.open(f, "rb") as g:
                        shutil.copyfileobj(g, stdin, 1 << 20)
            else:
                with open(f, "rb") as g:
                    shutil.copyfileobj(g, stdin, 1 << 20)
    finally:
        stdin.close()


def count_stream(panel_txt, reads, out, threads=4):
    """Run pcount on one or more FASTQ files (plain or gzip). Counts of several files add up.
    Returns {'reads': n, 'hits': total, 'hits_per_read': x, 'counter': ...} and writes OUT.json with the same."""
    info = native.build()
    reads = [str(r) for r in (reads if isinstance(reads, (list, tuple)) else [reads])]
    if len(reads) == 1 and reads[0].endswith(".gz") and info["isal"]:
        cmd = [info["pcount"], "gz", str(panel_txt), str(threads), str(out), reads[0]]
        r = subprocess.run(cmd, capture_output=True, text=True)
    else:
        cmd = [info["pcount"], "scan", str(panel_txt), str(threads), str(out)]
        with tempfile.TemporaryFile("w+") as so, tempfile.TemporaryFile("w+") as se:
            p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=so, stderr=se)
            _feed(reads, p.stdin)  # closes stdin
            rc = p.wait()
            so.seek(0)
            se.seek(0)
            r = subprocess.CompletedProcess(cmd, rc, so.read(), se.read())
    if r.returncode:
        raise RuntimeError("pcount failed:\n" + r.stdout + r.stderr)
    m = re.search(r"reads\s+(\d+)", r.stdout + r.stderr)
    n_reads = int(m.group(1)) if m else None
    hits = int(sum(int(x) for x in Path(out).read_text().split()))
    s = {"reads": n_reads, "hits": hits, "hits_per_read": (hits / n_reads) if n_reads else None,
         "counter": "pcount gz (isa-l)" if cmd[1] == "gz" else "pcount scan (stdin)"}
    Path(str(out) + ".json").write_text(json.dumps(s, indent=1))
    return s


def count_msbwt(panel_txt, bwt_dir, out):
    """Count every panel line as an exact substring of the reads in an existing msBWT (pymsbwt, Python >= 3.14).
    The result equals the streaming counter's (checked on a test set)."""
    try:
        import MUSCython.MultiStringBWTCython as MSBWT
    except ImportError as e:  # pragma: no cover
        raise RuntimeError("the msBWT backend needs pymsbwt: pip install 'pluck-quant[msbwt]' (Python >= 3.14)") from e
    bwt = MSBWT.loadBWT(str(bwt_dir))
    hits = 0
    with open(panel_txt) as f, open(out, "w") as o:
        for line in f:
            s = line.strip()
            if not s:
                continue
            c = int(bwt.countOccurrencesOfSeq(s))
            o.write(f"{c}\n")
            hits += c
    s = {"reads": None, "hits": hits, "hits_per_read": None, "counter": "pymsbwt"}
    Path(str(out) + ".json").write_text(json.dumps(s, indent=1))
    return s
