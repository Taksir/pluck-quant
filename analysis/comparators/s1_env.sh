#!/bin/bash
# Stage 1: micromamba env at the versions in tool_versions.txt + Salmon 1.10.0 release binary (as in e9_salmon). Resumable.
set -uo pipefail
W=${W:-/home/user/fc}; cd $W; export MAMBA_ROOT_PREFIX=$W/mm
if [ ! -x bin/micromamba ]; then
  curl -fsSL -o mm.tar.bz2 https://conda.anaconda.org/conda-forge/linux-64/micromamba-2.3.0-0.tar.bz2 && tar xjf mm.tar.bz2 bin/micromamba && rm mm.tar.bz2
fi
if [ ! -x $W/mm/envs/tools/bin/kallisto ]; then
  bin/micromamba create -y -q -n tools -c conda-forge -c bioconda \
    kallisto=0.52.0 star=2.7.11b subread=2.1.1 kmer-jellyfish=2.3.1 bowtie=1.3.1 samtools=1.24 seqtk isa-l aria2 pigz \
    numpy scipy pandas numba pysam python=3.12 make
fi
T=$W/mm/envs/tools/bin
if [ ! -x $W/salmon-1.10.0/bin/salmon ]; then
  mkdir -p $W/salmon-1.10.0
  curl -fsSL https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz | tar xz -C $W/salmon-1.10.0 --strip-components=1
fi
$W/salmon-1.10.0/bin/salmon --version
for t in kallisto STAR featureCounts jellyfish bowtie samtools seqtk python; do printf "%s: " $t; ($T/$t --version 2>&1 || $T/$t -v 2>&1) | head -2 | tr '\n' ' '; echo; done
[ -x $T/kallisto ] && echo S1_ENV_DONE || echo S1_ENV_FAILED
