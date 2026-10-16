#!/bin/bash
# counter fix tests; every count file is the panel-order text pcount writes
set -euo pipefail; cd /home/user/fix; P=w/panel.txt; T=4
run(){ local a=$(date +%s.%N); ./pcount_$1 gz $P $T $2 $3 2> $2.log; echo "$1 $3 $(awk -v a=$a -v b=$(date +%s.%N) "BEGIN{printf \"%.1f s\", b-a}") $(grep ^reads $2.log)"; }
echo "== A chunk $(stat -c%s chunk.fq.gz) bytes"
run old A_old.txt chunk.fq.gz; run fix A_fix.txt chunk.fq.gz
cmp A_old.txt A_fix.txt && echo "A: old == fix (byte-identical)"
zcat chunk.fq.gz | ./pcount_fix scan $P $T A_scan.txt 2> A_scan.log; cmp A_fix.txt A_scan.txt && echo "A: gz == scan(zcat)"
echo "== B 12x single member"
for i in $(seq 12); do zcat chunk.fq.gz; done | igzip -1 -T 2 > mid.fq.gz; echo "B bytes $(stat -c%s mid.fq.gz)"
run old B_old.txt mid.fq.gz; run fix B_fix.txt mid.fq.gz
cmp B_old.txt B_fix.txt && echo "B: old == fix (byte-identical)"
python3 -I -c "import numpy as np;a=np.loadtxt('A_fix.txt',dtype=np.int64);b=np.loadtxt('B_fix.txt',dtype=np.int64);print('B: fix == 12 x chunk', np.array_equal(b,12*a))"
rm -f mid.fq.gz
echo "== C 16x single member"
for i in $(seq 16); do zcat chunk.fq.gz; done | igzip -1 -T 2 > big.fq.gz; echo "C bytes $(stat -c%s big.fq.gz)"
run old C_old.txt big.fq.gz; run fix C_fix.txt big.fq.gz
zcat big.fq.gz | ./pcount_fix scan $P $T C_scan.txt 2> C_scan.log; cmp C_fix.txt C_scan.txt && echo "C: fix gz == scan(zcat)"
python3 -I -c "import numpy as np;a=np.loadtxt('A_fix.txt',dtype=np.int64);f=np.loadtxt('C_fix.txt',dtype=np.int64);o=np.loadtxt('C_old.txt',dtype=np.int64);print('C: fix == 16 x chunk', np.array_equal(f,16*a), 'fix sum', f.sum(), 'old sum', o.sum(), 'old/fix', round(o.sum()/f.sum(),4))"
rm -f big.fq.gz
echo "== D 16x multi member (cat)"
for i in $(seq 16); do cat chunk.fq.gz; done > multi.fq.gz; echo "D bytes $(stat -c%s multi.fq.gz)"
run old D_old.txt multi.fq.gz; run fix D_fix.txt multi.fq.gz
python3 -I -c "import numpy as np;a=np.loadtxt('A_fix.txt',dtype=np.int64);f=np.loadtxt('D_fix.txt',dtype=np.int64);o=np.loadtxt('D_old.txt',dtype=np.int64);print('D: fix == 16 x chunk', np.array_equal(f,16*a), 'old sum', o.sum(), 'fix sum', f.sum())"
rm -f multi.fq.gz
grep -h "^reads" *_old.txt.log *_fix.txt.log *_scan.log | head -20
echo ALL_TESTS_DONE
