#!/bin/bash
# Rebuild the Keele kidney, Keele liver and Munger scorer inputs from stored files + GEO, with the original scripts unchanged, then export the
# plain input files that fc_gather_score.py reads and check that the four existing methods reproduce the published compare blocks.
# Needs ~0.6 GB disk, 2 cores (use 2 jobs).
# usage: bash rebuild_frames.sh [FRAMES_WORKDIR, default under the original working folder]       Takes about 45 min on 2 cores (estimates 13 min x2 in parallel, make_tables 2+5 min, scorer 1.5 min, gates.py 2 min).
set -euo pipefail
# Folder constants below point to the original working folder; set them for your copy.
FR=${1:-/home/user/fc/frames}; R=/mnt/project-files/revision; GA=$(cd "$(dirname "$0")" && pwd)
PY=${PY:-/home/user/fc/mm/envs/tools/bin/python}; export PYTHONDONTWRITEBYTECODE=1 NUMBA_CACHE_DIR=$FR/numba_cache NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p $FR/{e4/ref,e4/truth,e4/salmon/fa,e4/dk,e5/w/ref,e5/dl,e5/panel,e5/out,out,frames}; cd $FR
# the original scripts hard-code two folders (e4, e5) in the original working folder: point them at the work trees (they did not exist before)
for t in e4 e5; do [ -e /home/user/$t ] && [ ! -L /home/user/$t ] && { echo "/home/user/$t exists and is not a symlink; stop"; exit 1; }; ln -sfn $FR/$t /home/user/$t; done
ln -sfn $R/frozen_v2/code e4/code; ln -sfn $R/frozen_v2/code e5/code
C=$R/fast_competitors/e9_prep_cache/ref
for f in genes.json runinfo.tsv geo_filelist.txt cc.zip; do ln -sf $C/$f e4/ref/$f; ln -sf $C/$f e5/w/ref/$f; done
zcat $R/e2_sim/work/ref/tx2gene.tsv.gz > e4/salmon/fa/tx2gene.tsv
cp $R/e4_kidney/dk_sites.npz e4/dk/sites.npz
# truth_files.txt = the 94 GEO gene_quant files of the 47 kidney + 47 liver samples (list written once from runinfo.tsv + geo_filelist.txt)
[ -s truth_files.txt ] || $PY - <<'PYE'
import re, os
R='/mnt/project-files/revision'; C=f'{R}/fast_competitors/e9_prep_cache/ref'
geo={}
for l in open(f'{C}/geo_filelist.txt'):
    m=re.search(r'(GSM\d+)_(CC\d+)_(\w+?)_gene_quant\.txt\.gz',l)
    if m: geo[(m.group(2),m.group(3))]=m.group(0)
title={l.split('\t')[0]:l.split('\t')[1] for l in open(f'{C}/runinfo.tsv') if l.split('\t')[1].endswith('_rna')}
kid=[l.split('\t')[0] for l in open(f'{R}/e4_kidney/prereg/samples.tsv').read().splitlines()[1:]]
liv=sorted(f[:-4] for f in os.listdir(f'{R}/probe_em/counts') if f.endswith('.npz'))
need={geo[(title[s].split('_')[0],title[s].split('_')[1])] for s in kid+liv}
open('truth_files.txt','w').write('\n'.join(sorted(need))+'\n')
PYE
bash $GA/dl_inputs.sh      # GTF (Ensembl 102, 33 MB), 94 truth files (42 MB), Munger matrix (7 MB); run from $FR
ln -sf $FR/e4/ref/genes.gtf.gz e5/w/ref/genes.gtf.gz
$PY $GA/mkpanel_nowrite.py                                       # map.npz only (copy of mkpanel.py without writing panel.txt)
S=$R/e4_kidney/scripts
nice -n 10 $PY $S/estimates.py liver  out/core_liver.npz  > out/est_liver.log 2>&1 &
nice -n 10 $PY $S/estimates.py kidney out/core_kidney.npz > out/est_kidney.log 2>&1 &
wait
nice -n 10 $PY $S/make_tables.py liver  out/core_liver.npz  liver_tab  > out/mt_liver.log 2>&1 &
nice -n 10 $PY $S/make_tables.py kidney out/core_kidney.npz kidney_tab > out/mt_kidney.log 2>&1 &
wait
for t in kidney liver; do nice -n 10 $PY $R/frozen_v2/code/score_e9.py ${t}_tab --extra ${t}_tab/extra.tsv --out ${t}_tab/prereg_out > out/score_$t.log 2>&1 & done; wait
for t in kidney liver; do nice -n 10 $PY $S/gates.py ${t}_tab $t > out/gates_$t.log 2>&1 & done; wait
for t in kidney liver; do cmp ${t}_tab/prereg_out/metrics.json $R/e4_kidney/${t}_frame/prereg_out/metrics.json && echo "$t metrics.json byte-identical to published"; cmp ${t}_tab/e4_gates.json $R/e4_kidney/${t}_frame/e4_gates.json && echo "$t e4_gates.json byte-identical"; done
# Munger: stored est.npz (hybrid etc.) + salmon.pkl rebuilt from the stored Salmon tables by the script's own salmon(); script unchanged
cd $FR/e5; ln -sf $R/e5_munger/samples.tsv samples.tsv
ln -sfn $R/e5_munger/design/dsn76 dsn76; ln -sfn $R/e5_munger/design/v2e v2e; ln -sfn $R/e5_munger/salmon sal
ln -sf $R/e5_munger/scripts/panel.json panel/panel.json; ln -sf $R/e5_munger/probes_s76.npz w/ref/probes_s76.npz; ln -sf $R/e5_munger/est.npz out/est.npz
$PY $GA/mk_salmon_pkl.py
nice -n 10 $PY $R/e5_munger/scripts/score_e5_post_lock_fix.py score > $FR/out/score_e5.log 2>&1
cmp out/metrics.json $R/e5_munger/metrics.json && echo "munger metrics.json identical to published"
# plain input files + reproduction of the four existing methods through the gather script's own path
G=$GA/fc_gather_score.py; cd $FR
$PY $G frames keele_kidney --frames frames --tabdir kidney_tab --sample-table $R/e4_kidney/prereg/samples.tsv
$PY $G frames keele_liver  --frames frames --tabdir liver_tab
(cd e5 && $PY $G frames munger --frames ../frames)
$PY $G repro --frames frames --out repro_out
echo REBUILD_DONE
