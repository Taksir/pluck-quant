# pluck-quant

PLUCK estimates gene-level expression in Collaborative Cross (CC) and Diversity Outbred (DO) mice by counting exact occurrences of a fixed set of 25-mers (probes) in the reads, without alignment and without a founder genome or founder transcriptome index. The probes are designed once from Ensembl release 102 and the Mouse Genomes Project (REL-1505) founder variants, and are free of annotated SNPs and small indels of the seven non-B6 CC founders (structural variants are not screened). One streaming pass over the FASTQ file gives the counts; an equivalence-class EM and a gene-probe estimate are combined by a switch.

PLUCK = Probe Levels from Unaligned reads by Counting K-mers. The distribution name is `pluck-quant` (the name `pluck` on PyPI belongs to another project); the command is `pluck` and the import name is `pluck_quant`.

Release `v1.0-submission` is the code that accompanies the PLUCK paper. It does not change any estimator, parameter or probe of the frozen design.

## Install

PLUCK is installed from GitHub. It will be put on PyPI when the paper is accepted.

```
python3 -m venv pluck-env && . pluck-env/bin/activate
pip install "git+https://github.com/Taksir/pluck-quant@v1.0-submission"        # Python >= 3.12; numpy, pandas, numba
pip install "pluck-quant[msbwt] @ git+https://github.com/Taksir/pluck-quant@v1.0-submission"   # optional: count from an existing msBWT (pymsbwt, Python >= 3.14)
```

To work on the source, `git clone https://github.com/Taksir/pluck-quant && pip install -e "pluck-quant[dev]"` and run `pytest`. A C compiler is needed once (gcc, cc or clang): the counter and the EM are compiled at first use into `~/.cache/pluck-quant`. With `libisal-dev` (isa-l 2.31) installed the counter inflates gzip itself (fastest); without it, gzip is piped in.

## What you need besides the package: a design

The design is the set of probes and the EM structure for one read length. It is large (design.npz is 359 MB for 50 bp) and is not in the package. There are three layouts:

| `--layout` | reads | probe spacing | gene-probe file | where it is used in the paper |
|---|---|---|---|---|
| 50 | 50 bp | 26 | `probes_s26.npz` | Keele CC, 50 bp |
| 100 | 100 bp | 76 | `probes_s76.npz` | Munger founder cohort |
| 101 | 101 bp | 77 | `probes_s77.npz` | GSE199702 |

Get a design in one of two ways.

```
# 1. download the prebuilt files and check them against the locked hashes
pluck design fetch --layout 50  --url https://github.com/Taksir/pluck-quant/releases/download/design-v1       --dest design50
pluck design fetch --layout 100 --url https://github.com/Taksir/pluck-quant/releases/download/v1.0-submission --dest design100
pluck design fetch --layout 101 --url https://github.com/Taksir/pluck-quant/releases/download/v1.0-submission --dest design101
# 2. or rebuild from public data (Ensembl 102, Mouse Genomes REL-1505; internet, a few hours, about 10 GB disk, 15 GB RAM); the locked hashes are checked at the end
pluck design build --layout 50 --work work
```

`pluck design build` can be run for a second and third layout in the same `--work` folder; the reference download and probe selection are shared and done once. The builds of all three layouts were run from scratch for this release and reproduce the locked design arrays, `gc.npy`, `genome_count.npy` and the gene-probe files (see `VERIFY.md`).

Then compile the panel (a few minutes, about 3 GB RAM). Choose the strand of your library:

```
pluck design compile --design-dir design50 --probes design50/probes_s26.npz --genes design50/genes.json \
      --read-length 50 --strand reverse -o panel50
```

`--strand reverse`: read 1 is antisense to the mRNA (dUTP libraries; Keele CC, GSE199702). `forward`: read 1 is the mRNA strand. `unstranded`: both orientations are counted and added (Munger). For the 100 bp Munger design use `--read-length 100 --probes design100/probes_s76.npz --strand unstranded`; for GSE199702 use `--read-length 101 --probes design101/probes_s77.npz --strand reverse`.

## Run

```
pluck run --panel panel50 -t 4 -o genes.tsv --branch \
      --sample CC021=CC021.fastq.gz --sample CC001=CC001_a.fastq.gz,CC001_b.fastq.gz
```

Output: one row per gene of `genes.json` (17,267 genes), columns `gene_id`, `gene_name`, then one column per sample (fragment units up to one constant per sample, as in the paper; `NA` if the gene has no probe), and with `--branch` a `<sample>.branch` column telling whether the gene-probe estimate (`probe`) or the EM estimate (`em`) was used. `--tiers builtin:tiers` adds the tier (core: 16,885 genes with a gene or EM probe) and `--flags builtin:flags` adds the design-only reliability flag (`unique_frac`, `flag`; flagged genes have unique fraction below 0.6). Both tables ship with the package (`flag_table.tsv`, sha256 `fff978bf...d506`; the flag rule is frozen in `flag_spec.json`, sha256 `ffa185a7...7e3c`, and its feature script is in `flag/`); a file of your own can be given instead.

### The switch

`--switch sample` (default) decides per sample and per gene: the gene-probe estimate if the observed capture f = S_g / EM_g is at least 0.6. This is the rule of the released code (version A in the paper). `--switch cohort` sums S_g and EM_g over all samples given in the call and decides once per gene; this is the rule behind every accuracy number of the paper. Pass the samples of one tissue and one layout in one call. `--switch sample-scaled` is the rule of the timed runs of the paper (version B): per sample, with the length-scaled gene-probe estimate N_g in the ratio f = N_g / EM_g. The accuracy results and the timed runs can therefore both be reproduced with this code: `--switch cohort` for the first, `--switch sample-scaled` for the second.

The EM stops when no gene total above 10 fragments changes by more than 1e-3 between rounds, over all genes of the run; an EM estimate can therefore depend slightly (about 1e-4 relative) on which other genes are in the design. Counts and gene-probe estimates are exact.

### Other commands

```
pluck count --panel panel50 --reads a.fastq.gz -o a.counts.txt        # only the counts
pluck quant --panel panel50 --counts a.counts.txt b.counts.txt -o genes.tsv
pluck run --panel panel50 --bwt --sample S1=path/to/msbwt_dir -o genes.tsv   # existing msBWT, needs pymsbwt
pluck design verify --layout 50 --design-dir design50
pluck native                                                           # compile / show the C programs
```

Each `pluck run` prints panel hits per read for every sample and warns below 0.001, which usually means wrong strandedness or a truncated file.

## Not in this release

The extended tier (26,258 further genes, non-negative least squares), the identity check against founder variants, the comparator scorers and the benchmark pipelines are not part of the package. Their frozen scripts are kept unchanged in `frozen_v2_snapshot/` of the source repository.

## Tests

`pytest` in the source repository runs the frozen-code comparison: a 300-gene cut of the real 50 bp design and 500,000 public reads (ENA SRR10523526) must give panel, counts and gene estimates identical to the frozen v2 code (counts byte-identical, estimates to 1e-9). On a whole sample (CC021, 18.1 M reads) all 5,025,384 probe counts and the estimates equal those stored for the paper.

## Licence and citation

Code: MIT, copyright MD TAKSIR HASAN MAJUMDER. Paper, figures and data tables: CC BY 4.0. Please cite the PLUCK paper (reference to be added at publication) and, for the optional msBWT backend, pymsbwt (github.com/Taksir/msbwt; original msBWT by J. M. Holt).
