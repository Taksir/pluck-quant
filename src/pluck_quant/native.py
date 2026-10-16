"""Build and locate the two C programs: the streaming probe counter (pcount) and the EM (em).

The sources ship inside the package. They are compiled once with the system C compiler into a cache
directory (PLUCK_CACHE, default ~/.cache/pluck-quant) keyed by the hash of the sources."""
import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path

SRC = Path(__file__).parent / "native"


class NativeBuildError(RuntimeError):
    pass


def cache_root() -> Path:
    return Path(os.environ.get("PLUCK_CACHE") or Path.home() / ".cache" / "pluck-quant")


def _run(cmd):
    r = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
    return r


def build(force: bool = False) -> dict:
    """Return {'pcount': path, 'em': path, 'isal': bool}. Compiles on first use."""
    key = hashlib.sha256(b"".join((SRC / n).read_bytes() for n in ("em.c", "pcount.c"))).hexdigest()[:12]
    d = cache_root() / key
    meta = d / "built.json"
    if meta.exists() and not force:
        info = json.loads(meta.read_text())
        if Path(info["pcount"]).exists() and Path(info["em"]).exists():
            return info
    cc = os.environ.get("CC") or shutil.which("gcc") or shutil.which("cc") or shutil.which("clang")
    if not cc:
        raise NativeBuildError("no C compiler found (install gcc, or set CC)")
    d.mkdir(parents=True, exist_ok=True)
    r = _run([cc, "-O3", "-o", d / "em", SRC / "em.c", "-lm"])
    if r.returncode:
        raise NativeBuildError("compiling em.c failed:\n" + r.stderr)
    # with isa-l the counter inflates gzip in-process (fast); without it, gzip is piped in (see count.py)
    r = _run([cc, "-O3", "-pthread", "-DWITH_ISAL", "-o", d / "pcount", SRC / "pcount.c", "-lisal"])
    isal = r.returncode == 0
    if not isal:
        r = _run([cc, "-O3", "-pthread", "-o", d / "pcount", SRC / "pcount.c"])
        if r.returncode:
            raise NativeBuildError("compiling pcount.c failed:\n" + r.stderr)
    info = {"pcount": str(d / "pcount"), "em": str(d / "em"), "isal": isal, "compiler": cc,
            "source_sha256": {n: hashlib.sha256((SRC / n).read_bytes()).hexdigest() for n in ("em.c", "pcount.c")}}
    meta.write_text(json.dumps(info, indent=1))
    return info
