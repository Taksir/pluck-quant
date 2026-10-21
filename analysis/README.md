# analysis/

Scripts and source tables behind the numbers, tables and figures of the PLUCK paper. The installed `pluck` command does not use this folder.

| Folder | Content |
|---|---|
| `scoring/` | `score_e9.py` and `score_noise.py` (scorers), `fc_gather_score.py` (gathers and scores all methods), `kt_score.py` (RNA-Skim scoring), scripts that rebuild the truth tables, `SHA256SUMS` of the four main scripts |
| `comparators/` | scripts and sample lists used to run the other methods and to time them |
| `identity/` | sample-identity check with founder-diagnostic 25-mers |
| `figures/scripts/` | one script per figure, and the scripts for the tables made from the source tables |
| `figures/data/` | source tables of the figures |
| `tables/` | the data files of the tables of the paper (`Table2.tsv`, Supplementary Data files `TableS*.tsv`) |

How to use them:
- Folder constants at the top of a script (`REV`, `D`, `OUT`, `/home/user/...`) point to the original working folders. Set them for your copy.
- The raw reads, the truth tables from GEO and the outputs of the compared methods are not stored here. `scoring/dl_inputs.sh` and the header of `scoring/fc_gather_score.py` say what is read.
- `scoring/fc_gather_score.py` loads `scoring/score_e9.py` and checks its SHA-256 against `SCORER_SHA`.
- Table `.tsv` files: lines that start with `#` are notes; read with `pandas.read_csv(path, sep="\t", comment="#")`.
- The extended-tier estimator is `reference_code/code/s6_estimate.py` and the identity scoring is `reference_code/code/dk.py`.
