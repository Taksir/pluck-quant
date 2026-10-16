# pluck-quant

PLUCK estimates gene-level expression in Collaborative Cross (CC) and Diversity Outbred (DO) mice by counting exact occurrences of a fixed set of 25-mers (probes) in the reads, without alignment and without a founder genome or founder transcriptome index. The probes are designed once from Ensembl release 102 and the Mouse Genomes Project (REL-1505) founder variants, and are free of annotated SNPs and small indels of the seven non-B6 CC founders (structural variants are not screened). One streaming pass over the FASTQ file gives the counts; an equivalence-class EM and a gene-probe estimate are combined by a switch.

PLUCK = Probe Levels from Unaligned reads by Counting K-mers. The PyPI name is `pluck-quant` (the name `pluck` belongs to another project); the command is `pluck`, the import name is `pluck_quant`.

Version 0.1 packages the frozen v2 design of the PLUCK paper. It does not change any estimator, parameter or probe.

## Install

```
pip install pluck-quant                 # Python >= 3.12; numpy, pandas, numba
pip install "pluck-quant[msbwt]"        # optional: count from an existing msBWT (pymsbwt, Python >= 3.14)
```

A C compiler is needed once (gcc, cc or clang): the counter and the EM are compiled at first use into `~/.cache/pluck-quant`. With `libisal-dev` (isa-l 2.31) installed the counter inflates gzip itself (fastest); without it, gzip is piped in.

## What you need besides the package: a design

The design is the set of probes and the EM structure for one read length. It is large (design.npz 359 MB for 50 bp) and is not in the package. Get it in one of two ways.

```
# 1. download the prebuilt files from the release and check them against the locked hashes
pluck design fetch --layout 50 --url https://.../releases/download/v0.1.0 --dest design50
# 2. or rebuild from public data (internet; 2 hours in the test container, mostly download time; about 4 GB disk; ran in a 15 GB-RAM container); the locked hashes are checked at the end
pluck design build --layout 50 --work work50
```

Then compile the panel (a few minutes, about 3 GB RAM). Choose the strand of your library:

```
pluck design compile --design-dir design50 --probes design50/probes_s26.npz --genes design50/genes.json \
      --read-length 50 --strand reverse -o panel50
```

`--strand reverse`: read 1 is antisense to the mRNA (dUTP libraries; Keele CC, GSE199702). `forward`: read 1 is the mRNA strand. `unstranded`: both orientations are counted and added (written as in the paper's Methods; checked here on the structure of the counts only, not against the Munger cohort). Layouts: 50 bp uses `probes_s26.npz`; 101 bp uses `probes_s77.npz` and `--layout 101`.

## Run

```
pluck run --panel panel50 -t 4 -o genes.tsv --branch \
      --sample CC021=CC021.fastq.gz --sample CC001=CC001_a.fastq.gz,CC001_b.fastq.gz
```

Output: one row per gene of `genes.json` (17,267 genes), columns `gene_id`, `gene_name`, then one column per sample (fragment units up to one constant per sample, as in the paper; `NA` if the gene has no probe), and with `--branch` a `<sample>.branch` column telling whether the gene-probe estimate (`probe`) or the EM estimate (`em`) was used. `--tiers gene_tiers.tsv` adds the tier (core: 16,885 genes with a gene or EM probe) and `--flags flag_table.tsv` adds the design-only reliability flag (`unique_frac`, `flag`; flagged genes have unique fraction below 0.6). Both tables are in the release.

### The switch

`--switch sample` (default) decides per sample and per gene: the gene-probe estimate if the observed capture f = S_g / EM_g is at least 0.6. This is the rule of the released code (version A in the paper). `--switch cohort` sums S_g and EM_g over all samples given in the call and decides once per gene; this is the rule behind every accuracy number of the paper. Pass the samples of one tissue and one layout in one call. The timed pipeline of the paper used a third, scaled per-sample rule (version B) that is not offered here.

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

## Not in 0.1

The extended tier (26,258 further genes, non-negative least squares), the identity check against founder variants, the scorers and the benchmark pipelines are not part of the package. Their frozen scripts are kept unchanged in `frozen_v2_snapshot/` of the source repository. Design recipes ship for 50 bp (spacing 26) and 101 bp (spacing 77); the 100 bp Munger design (spacing 76) has no recipe here.

## Tests

`pytest` in the source repository runs the frozen-code comparison: a 300-gene cut of the real 50 bp design and 500,000 public reads (ENA SRR10523526) must give panel, counts and gene estimates identical to the frozen v2 code (counts byte-identical, estimates to 1e-9). On a whole sample (CC021, 18.1 M reads) all 5,025,384 probe counts and the estimates equal those stored for the paper.

## Licence and citation

Code: MIT, copyright MD TAKSIR HASAN MAJUMDER. Paper, figures and data tables: CC BY 4.0. Please cite the PLUCK paper (reference to be added at publication) and, for the optional msBWT backend, pymsbwt (github.com/Taksir/msbwt; original msBWT by J. M. Holt).
