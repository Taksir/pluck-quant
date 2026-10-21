# Checking this release

The checks below were first run on release `v1.0-submission` and repeated on `v1.0.1-submission` (same package code, same design files) in a cloud container (Ubuntu 24.04, 4 cores, 15 GB RAM, Python 3.12, 3.13 and 3.14, gcc 13.3, isa-l 2.31). The commands let you repeat each check on your own machine.

## 1. Tests and a clean install (about 3 minutes)

```
git clone https://github.com/Taksir/pluck-quant && cd pluck-quant
sudo apt-get install -y gcc libisal-dev      # libisal-dev is optional (faster gzip reading)
bash scripts/verify.sh                       # last line: VERIFY OK
```

`verify.sh` runs the test suite, builds the wheel and source archive, installs the wheel into a clean environment, and runs `pluck` on 500,000 public reads (ENA SRR10523526, Keele CC021 liver). It compares with the outputs of the original analysis code (`reference_code/`):

- the compiled panel, the counts (byte for byte, any thread count) and the 300 test-gene estimates are identical;
- the three switches (`sample`, `sample-scaled`, `cohort`) reproduce the stored outputs of the paper's pipeline;
- counting from an msBWT and streaming the FASTQ give identical counts (this test runs with Python 3.14: `PY=python3.14 bash scripts/verify.sh`).

## 2. Whole-sample and whole-cohort checks against the stored results

| Check | Result |
|---|---|
| One whole sample (CC021, 18.1 M reads), 50 bp design | all 5,025,384 probe counts and all 17,267 gene estimates equal the stored ones |
| Keele kidney cohort (47 samples), `--switch cohort`, from the stored counts | 811,549 gene-sample values: probe-branch values equal to 5e-16, EM-branch values within 1.3e-6 relative, the same probe/EM decision for every gene |
| Munger founder cohort (47 samples, 100 bp, unstranded panel), `--switch cohort`, from the stored counts | 811,549 values: probe-branch equal to 5e-16, EM-branch within 1.5e-4 relative, no value missing or extra |
| One Munger FASTQ (SRR826311, 9.8 M reads) counted with the 100 bp unstranded panel | all 2,074,435 probe counts equal the stored counts |

The EM stops when no gene total above 10 fragments changes by more than 1e-3. That rule looks at every gene of the run, so an EM estimate shifts slightly with the set of genes in the design. Counts and gene-probe estimates are exact.

## 3. Rebuilding the designs from public data

`pluck design build --layout L --work W` downloads Ensembl 102 and the Mouse Genomes Project REL-1505 and rebuilds the probes and the design, then compares with the design behind the paper: the content of every array of `design.npz`, and the hashes of `gc.npy`, `genome_count.npy` and the gene-probe file.

| Layout | Stages | Time in the test container | Result |
|---|---|---|---|
| 50 bp (probe spacing 26) | design, GC, genome copies, tie order, check | about 29 min | all arrays and hashes match |
| 100 bp (spacing 76, Munger) | same | about 27 min | all arrays and hashes match |
| 101 bp (spacing 77, GSE199702) | same | about 29 min | all arrays and hashes match |

The reference download and probe selection (`prep`) are shared by the three layouts and took 111 minutes (once, for all three layouts) (mostly the genome-uniqueness scan and variant processing). Disk use is about 10 GB, peak memory about 7 GB.

The genome-copy scan sorts k-mer codes with an unstable sort, so the copy count of a probe and its reverse complement can land on either member depending on the machine. `gfix.py` puts the count on the member recorded in the paper's design. On the test machine the step moved 1,804 members of the 1,865 tie groups on the 100 bp layout and none on the 50 and 101 bp layouts; another machine can give a different answer, so the step stays in.

## 4. Downloads

The prebuilt design files for all three layouts are assets of release `v1.0.1-submission`. `pluck design fetch` checks them against the hashes above, so a download can be compared with a rebuild.
