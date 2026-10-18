"""Probe design -> counter panel + EM arrays ("panel directory"); design hash checks; design rebuild recipe."""
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from importlib import resources
from pathlib import Path

import numpy as np

K = 25
STRANDS = ("reverse", "forward", "unstranded")


def array_hashes(npz_path):
    """Content hash of every array of an .npz (zip timestamps and compression do not matter)."""
    d = np.load(npz_path, allow_pickle=False)
    out = {}
    for k in sorted(d.files):
        a = np.ascontiguousarray(d[k])
        out[(k, str(a.dtype), str(a.shape))] = hashlib.sha256(a.tobytes()).hexdigest()
    return out


def sha256_file(p, bufsize=1 << 24):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        while True:
            b = f.read(bufsize)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def verify_design(layout, design_dir):
    """Compare design.npz (array contents), gc.npy and genome_count.npy of DESIGN_DIR with the locked values
    for layout 50 (50 bp reads, spacing 26), 100 (100 bp reads, spacing 76) or 101 (101 bp reads, spacing 77). Returns (ok, lines)."""
    ref = resources.files("pluck_quant") / "_design"
    lines, bad = [], 0
    want = {}
    for l in (ref / f"expected_arrays_design{layout}.tsv").read_text().splitlines():
        f, k, dt, sh, h = l.split("\t")
        want[(k, dt, sh)] = h
    got = array_hashes(Path(design_dir) / "design.npz")
    for key, h in want.items():
        ok = got.get(key) == h
        bad += not ok
        if not ok:
            lines.append(f"FAIL design.npz {key} got {(got.get(key) or 'missing')[:12]} want {h[:12]}")
    lines.append(f"design.npz arrays: {len(want)} compared, " + ("all PASS" if not bad else f"{bad} FAIL"))
    for l in (ref / f"expected_files_design{layout}.sha256").read_text().splitlines():
        if not l.strip():
            continue
        h, name = l.split()
        g = sha256_file(Path(design_dir) / name)
        ok = g == h
        bad += not ok
        lines.append(f"{'PASS' if ok else 'FAIL'} {name} {g[:16]} want {h[:16]}")
    return bad == 0, lines


def _dec(codes):
    c = np.asarray(codes, dtype=np.uint64).copy()
    res = np.zeros((len(c), K), np.uint8)
    for i in range(K - 1, -1, -1):
        res[:, i] = c & np.uint64(3)
        c >>= np.uint64(2)
    return res


def compile_panel(design_dir, probes_npz, genes_json, read_length, out_dir, strand="reverse"):
    """Turn a design (design.npz, genome_count.npy, gc.npy + the gene-probe file) into the panel directory:
    panel.txt (what the counter reads), static.npz (EM arrays, gene-probe index), GID.npy, panel.json.

    strand: 'reverse' (read 1 antisense to the mRNA: dUTP libraries; the panel holds the reverse complement of each probe),
            'forward' (read 1 is the mRNA strand; sense strings), 'unstranded' (sense strings then antisense strings,
            the two counts of a probe are added). Same arithmetic as the locked compile step for 'reverse'."""
    if strand not in STRANDS:
        raise ValueError(f"strand must be one of {STRANDS}")
    os.environ.setdefault("NUMBA_CACHE_DIR", str(Path.home() / ".cache" / "pluck-quant" / "numba"))
    from ._frozen import em, e9rs  # noqa: E402 (lazy: numba)
    t0 = time.time()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    ddir = Path(design_dir)
    em.R = int(read_length)
    D = dict(np.load(ddir / "design.npz"))
    gc = np.load(ddir / "gc.npy")
    gcount = np.load(ddir / "genome_count.npy")
    tg_u, tg_i = np.unique(D["tgene"], return_inverse=True)
    tg_i = tg_i.ravel()
    pairs = np.unique(np.stack([D["cp"], tg_i[D["ct"]]], 1), axis=0)
    ngp = np.bincount(pairs[:, 0], minlength=len(D["U"]))
    use = D["inB"] & (gcount <= np.maximum(ngp, 1))
    use_t = use & (ngp == 1)
    M = em.Model(D, use_t, gc)
    h5 = np.ones(em.NB)
    h3 = np.ones(em.NB)
    hg = np.ones(8)
    a = M.rates(h5, h3, hg)
    A = np.bincount(M.kt, weights=a, minlength=M.ntx)
    m = M.ec_of_probe >= 0
    nprb = np.bincount(M.ec_of_probe[m], minlength=M.nec).astype(np.float64)
    ea = a[M.e_ix] * nprb[M.e_ec]
    et = M.kt[M.e_ix].astype(np.int32)
    Pg = np.load(probes_npz)
    vfm = Pg["in_vf"]
    Gcodes = np.array([e9rs.seqcodes(e9rs.enc(k))[0][0] for k in Pg["kmer"][vfm].tolist()], np.uint64)
    ag = Pg["gene"][vfm].astype(np.int64)
    G = json.load(open(genes_json))
    GID = np.array([g["gene_id"] for g in G])
    ng = len(G)
    LT = np.array([np.mean([max(L - read_length + 1, K) for L in g["tx_lens"]]) for g in G])
    mg = np.bincount(ag, minlength=ng)
    gpos = {g: i for i, g in enumerate(GID)}
    gmap = np.array([gpos.get(g, -1) for g in tg_u])
    em_idx = np.nonzero(m)[0]
    U = D["U"]
    pcodes = np.unique(np.concatenate([U[em_idx], Gcodes]))
    b = _dec(pcodes)
    acgt = np.frombuffer(b"ACGT", np.uint8)
    sense = acgt[b]
    anti = acgt[(3 - b)[:, ::-1]]
    lines = {"reverse": [anti], "forward": [sense], "unstranded": [sense, anti]}[strand]
    with open(out / "panel.txt", "wb") as f:
        f.write(b"\n".join(bytes(r) for x in lines for r in x) + b"\n")
    pos = np.searchsorted(U, pcodes)
    pos[pos == len(U)] = 0
    inU = U[pos] == pcodes
    uidx = pos[inU]
    uline = np.nonzero(inU)[0]
    gline = np.searchsorted(pcodes, Gcodes)
    np.savez(out / "static.npz", ec_of_probe=M.ec_of_probe.astype(np.int32), e_ec=M.e_ec.astype(np.int32), et=et, ea=ea, A=A,
             gof=tg_i.astype(np.int64), gline=gline, ag=ag, LT=LT, mg=mg, gmap=gmap, uidx=uidx, uline=uline, pcodes=pcodes,
             npanel=np.int64(len(pcodes)), nU=np.int64(len(U)), ng=np.int64(ng), nec=np.int64(M.nec), ntx=np.int64(M.ntx),
             R=np.int64(read_length), ntg=np.int64(len(tg_u)))
    np.save(out / "GID.npy", GID)
    # gene names, for the output table
    json.dump({"gene_id": [g["gene_id"] for g in G], "gene_name": [g.get("gene_name", "") for g in G]}, open(out / "genes.names.json", "w"))
    info = {"read_length": int(read_length), "strand": strand, "n_probes": int(len(pcodes)), "n_gene_probes": int(len(Gcodes)),
            "n_em_probes": int(len(em_idx)), "n_genes": int(ng), "panel_lines": int(len(pcodes) * len(lines)),
            "unstranded_collisions": None, "compile_seconds": round(time.time() - t0, 1),
            "design_hashes": {k: sha256_file(ddir / k) for k in ("gc.npy", "genome_count.npy")}}
    if strand == "unstranded":
        vs = sense.view(np.dtype((np.void, K))).ravel()
        va = anti.view(np.dtype((np.void, K))).ravel()
        info["unstranded_collisions"] = int(np.isin(vs, va).sum())  # probes whose sense string is also some antisense string
    info["panel_sha256"] = sha256_file(out / "panel.txt")
    json.dump(info, open(out / "panel.json", "w"), indent=1)
    return info


# ----------------------------------------------------------------------------- download a prebuilt design
_GENES_SHA = "25e9297094ce681d42405e8d9a108c728a6d830d5f3c21c98b08b95f0b760148"
_PROBES = {50: ("probes_s26.npz", "f6c171db3987207df6dc543e9ffbf211de4bea22159d2771034e967d7e26df20"),
           100: ("probes_s76.npz", "533d70b66640ece3616c237e8d1acfebc33bdfd01ac0ff0a1c7960974e993cea"),
           101: ("probes_s77.npz", "053bfc6f3a3ae1b28f8cb460203f559150e95f69b65b7ec2fd163e9886087c94")}
# layout -> (probe spacing, design folder, gene-probe file); the spacing is the read length minus 24
_LAYOUTS = {50: (26, "dsn26", "probes_s26.npz"), 100: (76, "dsn76", "probes_s76.npz"), 101: (77, "dsn77", "probes_s77.npz")}


def fetch_design(layout, base_url, dest):
    """Download design.npz, gc.npy, genome_count.npy, the gene-probe file and genes.json from BASE_URL (files in BASE_URL/<layout>/, in BASE_URL, or named design_<layout>.npz, gc_<layout>.npy, genome_count_<layout>.npy) and check them
    against the locked hashes. BASE_URL is where the release assets were put."""
    import urllib.error
    import urllib.request
    layout = int(layout)
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)
    ref = resources.files("pluck_quant") / "_design"
    want = {}
    for l in (ref / f"expected_files_design{layout}.sha256").read_text().splitlines():
        if l.strip():
            h, n = l.split()
            want[n] = h
    pn, ph = _PROBES[layout]
    want[pn] = ph
    want["genes.json"] = _GENES_SHA
    want["design.npz"] = None  # checked by array content in verify_design
    for n, h in want.items():
        target = dest / n
        if not target.exists():
            b = base_url.rstrip("/")
            stem, dot, ext = n.rpartition(".")
            # a folder per layout, one flat folder, or a flat folder whose files carry the layout (GitHub release assets are flat)
            cands = [f"{b}/{layout}/{n}", f"{b}/{n}", f"{b}/{stem}_{layout}{dot}{ext}"]
            for i, u in enumerate(cands):
                try:
                    urllib.request.urlretrieve(u, target)
                    break
                except urllib.error.HTTPError as e:
                    if e.code != 404 or i == len(cands) - 1:
                        raise
        if h and sha256_file(target) != h:
            raise RuntimeError(f"{n}: sha256 does not match the locked value")
    ok, lines = verify_design(layout, dest)
    if not ok:
        raise RuntimeError("design.npz does not match the locked arrays:\n" + "\n".join(lines))
    return dest


# ----------------------------------------------------------------------------- rebuild the design from public data
_STAGES = ("prep", "design", "gc", "gscan", "gfix", "verify")


def build_design(layout, work, threads=None, log=print, stages=_STAGES):
    """Rebuild a design from Ensembl 102 and the Mouse Genomes Project REL-1505 (needs internet, ~8 GB disk, a few hours,
    about 15 GB RAM) with the same stages and edits as the Longleaf rebuild, then check the locked hashes.
    layout 50 -> SP 26 (50 bp reads); 100 -> SP 76 (100 bp reads); 101 -> SP 77 (101 bp reads). Writes
    WORK/design/dsn{26,76,77}/{design.npz,gc.npy,genome_count.npy} and WORK/ref/probes_s{26,76,77}.npz, WORK/ref/genes.json.
    The reference download and the probe selection (stage prep) are shared: running a second layout in the same WORK folder
    with stages without 'prep' reuses them."""
    sp, dname, pname = _LAYOUTS[int(layout)]
    w = Path(work).resolve()
    src = Path(__file__).parent
    py = sys.executable
    (w / "ref").mkdir(parents=True, exist_ok=True)

    def sh(cmd, cwd):
        log("+ " + " ".join(map(str, cmd)))
        t = time.time()
        subprocess.run([str(c) for c in cmd], cwd=cwd, check=True)
        log(f"  done in {time.time() - t:.0f} s")

    if "prep" in stages:
        shutil.copy(src / "_frozen" / "e9rs.py", w / "e9rs.py")
        # the frozen prep, unchanged, with the list of probe spacings set to the three the paper uses (76 = 100 bp reads);
        # probes_s26 and probes_s77 come out identical to the files of the default list
        (w / "prep_run.py").write_text("import sys\nsys.path.insert(0, '.')\nimport e9rs\ne9rs.SPACINGS = [26, 76, 77]\n"
                                       "sys.argv = ['e9rs.py', 'prep'] + sys.argv[1:]\ne9rs.main()\n")
        sh([py, "prep_run.py", "--work", w] + (["--threads", threads] if threads else []), w)
        for need in ("probes_s26.npz", "probes_s76.npz", "probes_s77.npz", "genes.json"):
            if not (w / "ref" / need).exists():
                raise RuntimeError(f"prep did not write ref/{need}")
    d = w / "design" / dname
    if "design" in stages:
        (d / "ref").mkdir(parents=True, exist_ok=True)
        (d / "w").mkdir(parents=True, exist_ok=True)
        for link, tgt in (("ref/Mus_musculus.GRCm38.102.gtf.gz", "genes.gtf.gz"), ("ref/Mus_musculus.GRCm38.cdna.all.fa.gz", "cdna.fa.gz"),
                          ("ref/Mus_musculus.GRCm38.ncrna.fa.gz", "ncrna.fa.gz"), ("ref/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz", "dna.fa.gz")):
            p = d / link
            if p.is_symlink() or p.exists():
                p.unlink()
            p.symlink_to(w / "ref" / tgt)
        if (d / "w" / "ref").is_symlink():
            (d / "w" / "ref").unlink()
        (d / "w" / "ref").symlink_to(w / "ref")
        shutil.copy(src / "_frozen" / "e9rs.py", d / "e9rs.py")
        text = (src / "_design" / "design.py").read_text()
        text = text.replace("SP = 26", f"SP = {sp}", 1)
        text = text.replace("/mnt/project-files/revision/respacing/probes/probes_s26.npz", str(w / "ref" / "probes_s26.npz"))
        text = text.split("# antisense k-mers for KMC")[0]
        (d / "design.py").write_text(text)
        shutil.copy(src / "_design" / "gc.py", d / "gc.py")
        shutil.copy(src / "_design" / "gscan.py", d / "gscan.py")
        sh([py, "design.py"], d)
    if "gc" in stages:
        sh([py, "gc.py"], d)
    if "gscan" in stages:
        sh([py, "gscan.py"], d)
    if "gfix" in stages:
        # gscan sorts with an unstable argsort; gfix puts each tie group's total on the member the locked design recorded
        for lay in (50, 100, 101):
            (d / f"carriers{lay}.b64").write_bytes((src / "_design" / f"carriers{lay}.b64").read_bytes())
        shutil.copy(src / "_design" / "gfix.py", d / "gfix.py")
        sh([py, "gfix.py", layout, d], d)
    ok = True
    if "verify" in stages:
        ok, lines = verify_design(int(layout), d)
        for l in lines:
            log(l)
    return ok, d
