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


BWT_FILES = ("msbwt.npy", "comp_msbwt.npy", "comp_msbwt.dat")


def is_msbwt_dir(path):
    """True if PATH is a directory holding an msBWT (the files pymsbwt's loadBWT looks for)."""
    p = Path(path)
    return p.is_dir() and any((p / f).exists() for f in BWT_FILES)


def msbwt_available():
    try:
        import MUSCython.MultiStringBWTCython  # noqa: F401
        return True
    except ImportError:
        return False


def split_inputs(paths):
    """Sort the input paths of one sample into msBWT directories and read files. Anything else is an error."""
    bwts, reads = [], []
    for p in paths:
        if is_msbwt_dir(p):
            bwts.append(p)
        elif Path(p).is_file():
            reads.append(p)
        elif Path(p).is_dir():
            raise ValueError(f"{p} is a directory but not an msBWT (none of {', '.join(BWT_FILES)} in it)")
        else:
            raise ValueError(f"{p} does not exist")
    return bwts, reads


def count_msbwt(panel_txt, bwt_dirs, out):
    """Count every panel line as an exact substring of the reads in existing msBWTs (pymsbwt, Python >= 3.14); the counts of
    several msBWTs of one sample add up. The result equals the streaming counter's (checked on a test set)."""
    try:
        import MUSCython.MultiStringBWTCython as MSBWT
    except ImportError as e:  # pragma: no cover
        raise RuntimeError("an msBWT was given but pymsbwt is not installed: pip install 'pluck-quant[msbwt]' (Python >= 3.14)") from e
    if isinstance(bwt_dirs, (str, Path)):
        bwt_dirs = [bwt_dirs]
    bwts = [MSBWT.loadBWT(str(d)) for d in bwt_dirs]
    hits = 0
    with open(panel_txt) as f, open(out, "w") as o:
        for line in f:
            s = line.strip()
            if not s:
                continue
            c = sum(int(b.countOccurrencesOfSeq(s)) for b in bwts)
            o.write(f"{c}\n")
            hits += c
    s = {"reads": None, "hits": hits, "hits_per_read": None, "counter": "pymsbwt"}
    Path(str(out) + ".json").write_text(json.dumps(s, indent=1))
    return s


def count_sample(panel_txt, paths, out, threads=4, mode="auto", warn=print):
    """Count the panel in one sample. PATHS are read files (FASTQ, gzip or plain) and/or msBWT directories.
    mode 'auto' (default): use the msBWT if the sample has one, otherwise stream the reads. If an msBWT is given but pymsbwt is not
        installed, warn and stream the reads when read files are given too, else stop with an error.
    mode 'bwt': msBWT only (error if a read file is given or pymsbwt is missing).
    mode 'stream': read files only (an msBWT directory in PATHS is ignored with a warning; error if there are no read files)."""
    bwts, reads = split_inputs(paths)
    if mode == "bwt":
        if reads or not bwts:
            raise ValueError("--bwt: every input must be an msBWT directory")
        return count_msbwt(panel_txt, bwts, out)
    if mode == "stream":
        if bwts:
            warn(f"[pluck] --stream: ignoring msBWT {', '.join(map(str, bwts))}")
        if not reads:
            raise ValueError("--stream: no read file given")
        return count_stream(panel_txt, reads, out, threads)
    if bwts:
        if msbwt_available():
            if reads:
                warn("[pluck] msBWT found; counting from it (the read files given with it are not read)")
            return count_msbwt(panel_txt, bwts, out)
        if reads:
            warn(f"[pluck] WARNING: msBWT {', '.join(map(str, bwts))} found but pymsbwt is not installed "
                 f"(pip install 'pluck-quant[msbwt]', Python >= 3.14); streaming the read files instead")
            return count_stream(panel_txt, reads, out, threads)
        raise RuntimeError(f"{bwts[0]} is an msBWT but pymsbwt is not installed: pip install 'pluck-quant[msbwt]' (Python >= 3.14), "
                           f"or give the FASTQ file(s) of the sample")
    return count_stream(panel_txt, reads, out, threads)
