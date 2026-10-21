#!/usr/bin/env bash
# One-command check of the source tree (about 3 minutes). Run from the repository root:  bash scripts/verify.sh
# Needs: Python >= 3.12, gcc, internet for pip. Optional: libisal-dev (apt) for the fast gzip path; Python 3.14 to include the msBWT test.
set -euo pipefail
cd "$(dirname "$0")/.."
PY=${PY:-python3}
$PY -m venv .venv_verify
. .venv_verify/bin/activate
pip install -q -e ".[dev]" build twine
echo "== 1. unit and end-to-end tests (compare with outputs of the original analysis code)"
python -m pytest -q
echo "== 2. build sdist and wheel, check metadata"
rm -rf dist && python -m build -q && twine check dist/*
echo "== 3. install the wheel into a clean environment and run it on the test data"
$PY -m venv .venv_clean && ./.venv_clean/bin/pip install -q dist/*.whl
D=tests/data; C=./.venv_clean/bin/pluck; O=$(mktemp -d)
$C design compile --design-dir $D/mini_design/inputs --probes $D/mini_design/inputs/probes_s26.npz --genes $D/genes.json --read-length 50 -o $O/panel > $O/panel.json
test "$(sha256sum $O/panel/panel.txt | cut -d' ' -f1)" = "$(cat $D/expected_mini/panel.sha256)"
$C run --panel $O/panel --sample CC021=$D/reads/CC021_first500k.fastq.gz --counts-dir $O -t 4 --branch -o $O/genes.tsv
test "$(sha256sum $O/CC021.counts.txt | cut -d' ' -f1)" = "$(cat $D/expected_mini/counts.sha256)"
python - "$O/genes.tsv" "$D/expected_mini/estimates_A.tsv" <<'PY'
import sys
got = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in list(open(sys.argv[1]))[1:]}
n = bad = 0
for l in list(open(sys.argv[2]))[1:]:
    g, e, b = l.rstrip("\n").split("\t")
    r = got[g]; n += 1
    if e == "NA" or r[2] == "NA":
        bad += (e == "NA") != (r[2] == "NA"); continue
    bad += abs(float(r[2]) - float(e)) > 1e-9 * abs(float(e)) + 1e-12
print(f"{n} genes compared, {bad} differ")
sys.exit(1 if bad else 0)
PY
echo "VERIFY OK"
