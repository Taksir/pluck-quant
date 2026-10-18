"""Command line: pluck {run,count,quant,design,native,version}."""
import argparse
import json
import sys
import tempfile
from pathlib import Path

import numpy as np

from . import __version__, count as count_mod, estimate, native, output, panel as panel_mod


def _samples(items):
    out = []
    for it in items:
        if "=" not in it:
            raise SystemExit(f"--sample needs NAME=reads1[,reads2,...], got {it!r}")
        n, v = it.split("=", 1)
        out.append((n, v.split(",")))
    names = [n for n, _ in out]
    if len(set(names)) != len(names):
        raise SystemExit("sample names must be unique")
    return out


def _quant(args, names, count_files):
    pan = estimate.Panel(args.panel)
    res = []
    for n, cf in zip(names, count_files):
        c = estimate.read_counts(pan, cf)
        res.append(estimate.sample_estimates(pan, c))
        print(f"[pluck] {n}: EM {res[-1]['em_iterations'][:60]}", file=sys.stderr)
    if args.switch == "cohort":
        est, probe, f = estimate.hybrid_cohort(res, args.tau)
    else:
        pr = [estimate.hybrid_sample(r, args.tau, scaled=args.switch == "sample-scaled") for r in res]
        est = np.stack([p[0] for p in pr])
        probe = np.stack([p[1] for p in pr])
    output.write_table(args.out, pan, names, est, probe, tiers=args.tiers, flags=args.flags,
                       read_length=pan.info.get("read_length", int(pan.S["R"])), branch=args.branch)
    print(f"[pluck] wrote {args.out}", file=sys.stderr)


def _count(pan, paths, out, threads, mode):
    try:
        return count_mod.count_sample(pan.panel_txt, paths, out, threads, mode=mode, warn=lambda m: print(m, file=sys.stderr))
    except (ValueError, RuntimeError) as e:
        raise SystemExit(f"pluck: {e}")


def _mode(args):
    return "bwt" if args.bwt else "stream" if args.stream else "auto"


def cmd_run(args):
    pan = estimate.Panel(args.panel)
    samples = _samples(args.sample)
    wd = Path(args.counts_dir) if args.counts_dir else Path(tempfile.mkdtemp(prefix="pluck_counts_"))
    wd.mkdir(parents=True, exist_ok=True)
    cfs = []
    for n, reads in samples:
        cf = wd / f"{n}.counts.txt"
        s = _count(pan, reads, cf, args.threads, _mode(args))
        print(f"[pluck] {n}: {s}", file=sys.stderr)
        if s.get("hits_per_read") is not None and s["hits_per_read"] < args.min_hits_per_read:
            print(f"[pluck] WARNING {n}: only {s['hits_per_read']:.4f} panel hits per read; "
                  f"check strandedness, read length and that the whole file was read", file=sys.stderr)
        cfs.append(cf)
    _quant(args, [n for n, _ in samples], cfs)


def cmd_count(args):
    pan = estimate.Panel(args.panel)
    paths = list(args.reads or []) + ([args.bwt] if args.bwt else [])
    if not paths:
        raise SystemExit("give --reads (read files and/or msBWT directories) or --bwt DIR")
    mode = "bwt" if args.bwt and not args.reads else "stream" if args.stream else "auto"
    s = _count(pan, paths, args.out, args.threads, mode)
    print(json.dumps(s))


def cmd_quant(args):
    names = args.name or [Path(c).name.split(".")[0] for c in args.counts]
    if len(names) != len(args.counts):
        raise SystemExit("--name must be given once per counts file")
    _quant(args, names, args.counts)


def cmd_design(args):
    if args.dcmd == "verify":
        ok, lines = panel_mod.verify_design(args.layout, args.design_dir)
        print("\n".join(lines))
        sys.exit(0 if ok else 1)
    if args.dcmd == "compile":
        info = panel_mod.compile_panel(args.design_dir, args.probes, args.genes, args.read_length, args.out, args.strand)
        print(json.dumps(info, indent=1))
    if args.dcmd == "fetch":
        d = panel_mod.fetch_design(args.layout, args.url, args.dest)
        print("design verified in", d)
    if args.dcmd == "build":
        ok, d = panel_mod.build_design(args.layout, args.work, threads=args.threads)
        print("design dir:", d, "| locked hashes:", "PASS" if ok else "FAIL")
        sys.exit(0 if ok else 1)


def cmd_native(args):
    print(json.dumps(native.build(force=args.force), indent=1))


def _add_quant_opts(p):
    p.add_argument("--panel", required=True, help="panel directory made by 'pluck design compile'")
    p.add_argument("--switch", choices=["sample", "cohort", "sample-scaled"], default="sample",
                   help="hybrid switch: 'sample' = per sample (version A, the released rule); 'cohort' = summed over all samples given "
                        "(the rule behind the accuracy results of the paper; give the samples of one tissue and layout); "
                        "'sample-scaled' = per sample with the length-scaled gene-probe estimate in the ratio (version B, the rule of "
                        "the timed runs)")
    p.add_argument("--tau", type=float, default=estimate.TAU)
    p.add_argument("--tiers", help="gene_tiers.tsv, or builtin:tiers for the table shipped with the package (adds a tier column)")
    p.add_argument("--flags", help="flag_table.tsv, or builtin:flags for the table shipped with the package (adds unique_frac and flag columns)")
    p.add_argument("--branch", action="store_true", help="add a column per sample telling which estimate was used")
    p.add_argument("-o", "--out", required=True, help="output table (TSV; .gz allowed)")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="pluck", description=__doc__)
    ap.add_argument("--version", action="version", version=f"pluck-quant {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("run", help="count the panel in reads and write the gene table")
    p.add_argument("--sample", action="append", required=True, metavar="NAME=path[,path...]",
                   help="one sample: FASTQ file(s) (gzip or plain) and/or an existing msBWT directory; repeat for more samples. "
                        "If the sample has an msBWT, PLUCK counts from it (needs pymsbwt, Python >= 3.14); otherwise it streams the "
                        "FASTQ file(s). An msBWT without pymsbwt falls back to the FASTQ file(s) given with it, with a warning")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--bwt", action="store_true", help="force msBWT: every --sample path must be an msBWT directory")
    g.add_argument("--stream", action="store_true", help="force streaming: ignore msBWT directories, read FASTQ files only")
    p.add_argument("-t", "--threads", type=int, default=4)
    p.add_argument("--counts-dir", help="keep the counts files here (default: temporary)")
    p.add_argument("--min-hits-per-read", type=float, default=0.001)
    _add_quant_opts(p)
    p.set_defaults(fn=cmd_run)

    p = sub.add_parser("count", help="only count the panel (one sample)")
    p.add_argument("--panel", required=True)
    p.add_argument("--reads", nargs="+", help="read files and/or msBWT directories (an msBWT is used when present, see 'run')")
    p.add_argument("--bwt", help="existing msBWT directory (forces msBWT)")
    p.add_argument("--stream", action="store_true", help="force streaming of the read files")
    p.add_argument("-t", "--threads", type=int, default=4)
    p.add_argument("-o", "--out", required=True)
    p.set_defaults(fn=cmd_count)

    p = sub.add_parser("quant", help="gene table from counts files")
    p.add_argument("--counts", nargs="+", required=True)
    p.add_argument("--name", action="append", help="sample name for each counts file")
    _add_quant_opts(p)
    p.set_defaults(fn=cmd_quant)

    p = sub.add_parser("design", help="design files: verify, compile, build")
    ds = p.add_subparsers(dest="dcmd", required=True)
    q = ds.add_parser("verify", help="check design.npz, gc.npy, genome_count.npy against the locked hashes")
    q.add_argument("--layout", type=int, choices=[50, 100, 101], required=True)
    q.add_argument("--design-dir", required=True)
    q = ds.add_parser("compile", help="design -> panel directory")
    q.add_argument("--design-dir", required=True, help="folder with design.npz, gc.npy, genome_count.npy")
    q.add_argument("--probes", required=True, help="gene-probe file (probes_s26.npz for 50 bp, probes_s76.npz for 100 bp, probes_s77.npz for 101 bp)")
    q.add_argument("--genes", required=True, help="genes.json")
    q.add_argument("--read-length", type=int, required=True)
    q.add_argument("--strand", choices=list(panel_mod.STRANDS), default="reverse")
    q.add_argument("-o", "--out", required=True)
    q = ds.add_parser("fetch", help="download a prebuilt design (release assets) and check the locked hashes")
    q.add_argument("--layout", type=int, choices=[50, 100, 101], required=True)
    q.add_argument("--url", required=True, help="base URL holding [<layout>/]design.npz, gc.npy, genome_count.npy, probes_*.npz, genes.json")
    q.add_argument("--dest", required=True)
    q = ds.add_parser("build", help="rebuild a design from Ensembl 102 + MGP REL-1505 (internet, hours) and verify the hashes")
    q.add_argument("--layout", type=int, choices=[50, 100, 101], required=True)
    q.add_argument("--work", required=True)
    q.add_argument("--threads", type=int)
    p.set_defaults(fn=cmd_design)

    p = sub.add_parser("native", help="compile (or show) the C programs")
    p.add_argument("--force", action="store_true")
    p.set_defaults(fn=cmd_native)

    args = ap.parse_args(argv)
    args.fn(args)


if __name__ == "__main__":
    main()
