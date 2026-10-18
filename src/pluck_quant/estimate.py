"""Per-sample estimation (gene-probe estimate N_g, EM estimate) and the hybrid switch.

Arithmetic is that of the frozen v2 code (run_sample2_frozenf.py for the per-sample rule, build_variants.py variant C0 for the
cohort rule), so that outputs are numerically identical to the paper's pipeline."""
import json
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from . import native

TAU = 0.6
K = 25


class Panel:
    def __init__(self, d):
        d = Path(d)
        z = np.load(d / "static.npz")
        self.S = {k: z[k] for k in z.files}
        self.GID = np.load(d / "GID.npy", allow_pickle=True)
        self.info = json.loads((d / "panel.json").read_text()) if (d / "panel.json").exists() else {}
        nm = d / "genes.names.json"
        self.names = json.loads(nm.read_text())["gene_name"] if nm.exists() else [""] * len(self.GID)
        self.dir = d
        self.strand = self.info.get("strand", "reverse")

    @property
    def n_probes(self):
        return int(self.S["npanel"])

    @property
    def panel_txt(self):
        return self.dir / "panel.txt"


def read_counts(panel: Panel, path):
    """Counts file of the counter (one integer per panel line). Unstranded panels list every probe twice (sense, antisense):
    the two counts of a probe are added."""
    c = np.loadtxt(path, dtype=np.int64, ndmin=1)
    n = panel.n_probes
    if panel.strand == "unstranded":
        if len(c) != 2 * n:
            raise ValueError(f"{path}: {len(c)} counts, panel has {2 * n} lines")
        return (c[:n] + c[n:2 * n]).astype(np.float64)
    if len(c) < n:
        raise ValueError(f"{path}: {len(c)} counts, panel has {n} lines")
    return c[:n].astype(np.float64)  # extra lines would belong to other panel segments


def sample_estimates(panel: Panel, c_sel, em_exe=None):
    """Returns dict: Sg (summed gene-probe counts), Ng (gene-probe estimate, Eq. 1), em (EM estimate, NaN if none), mg."""
    S = panel.S
    em_exe = em_exe or native.build()["em"]
    c = np.zeros(int(S["nU"]))
    c[S["uidx"]] = c_sel[S["uline"]]
    ec = S["ec_of_probe"]
    m = ec >= 0
    nec = int(S["nec"])
    cec = np.bincount(ec[m], weights=c[m], minlength=nec)
    e_ec = S["e_ec"]
    keep = cec[e_ec] > 0
    nz = np.nonzero(cec > 0)[0]
    remap = -np.ones(nec, np.int64)
    remap[nz] = np.arange(len(nz))
    e_ec2 = remap[e_ec[keep]]
    ea = S["ea"][keep]
    et = S["et"][keep].astype(np.int64)
    cec2 = cec[nz]
    eoff = np.searchsorted(e_ec2, np.arange(len(nz) + 1)).astype(np.int64)
    with tempfile.TemporaryDirectory(prefix="pluck_em_") as wd:
        for n, x in (("eoff", eoff), ("et", et), ("ea", ea.astype(np.float64)), ("cec", cec2), ("A", S["A"]), ("gof", S["gof"])):
            x.tofile(f"{wd}/{n}.bin")
        r = subprocess.run([em_exe, wd, "1e-3"], capture_output=True, text=True)
        if r.returncode:
            raise RuntimeError("em failed: " + r.stdout + r.stderr)
        X = np.fromfile(f"{wd}/X.bin", dtype=np.float64)
    ntg, ng = int(S["ntg"]), int(S["ng"])
    gsum = np.bincount(S["gof"], X, ntg)
    gmap = S["gmap"]
    em_est = np.full(ng, np.nan)
    ok = gmap >= 0
    em_est[gmap[ok]] = gsum[ok]
    Sg = np.bincount(S["ag"], weights=c_sel[S["gline"]], minlength=ng)
    mg = S["mg"]
    Ng = np.where(mg > 0, Sg * S["LT"] / (np.maximum(mg, 1) * (int(S["R"]) - K + 1)), np.nan)
    return {"Sg": Sg, "Ng": Ng, "em": em_est, "mg": mg, "em_iterations": r.stdout.strip()}


def hybrid_sample(res, tau=TAU, scaled=False):
    """Per-sample switch. Default (version A, the released rule): f = S_g / EM_g within the sample, gene-probe estimate if f >= tau.
    scaled=True (version B, the rule of the timed pipeline, run_sample2.py of the frozen runs): f = N_g / EM_g with the
    length-scaled gene-probe estimate N_g in the numerator."""
    num = res["Ng"] if scaled else res["Sg"]
    with np.errstate(invalid="ignore", divide="ignore"):
        f = np.where(res["em"] > 0, num / res["em"], np.nan)
        est = np.where((res["mg"] > 0) & (f >= tau), res["Ng"], res["em"])
    est = np.where(np.isnan(est), res["em"], est)
    probe = (res["mg"] > 0) & (f >= tau) & ~np.isnan(res["Ng"])
    return est, probe


def hybrid_cohort(results, tau=TAU):
    """Cohort switch, the rule behind every accuracy number in the paper: f_g = sum_i S_gi / sum_i EM_gi over the samples of
    one cohort (same tissue and layout), gene-probe estimate for all samples if f_g >= tau."""
    S = np.stack([r["Sg"] for r in results])
    EM = np.stack([r["em"] for r in results])
    Ng = np.stack([r["Ng"] for r in results])
    with np.errstate(invalid="ignore", divide="ignore"):
        f = np.nansum(S, 0) / np.maximum(np.nansum(EM, 0), 1e-9)
        sw = f[None, :] >= tau
        h = np.where(sw, Ng, EM)
    h = np.where(np.isnan(h), EM, h)
    probe = np.broadcast_to(sw, Ng.shape) & ~np.isnan(Ng)
    return h, probe, f
