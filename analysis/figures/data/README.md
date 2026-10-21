# Source tables of the figures and tables

Each `.tsv` file holds the values behind one figure panel or table. Lines that start with `#` at the top of a file say what the file holds and give caveats for reading it. Read a file with

```
pandas.read_csv(path, sep="\t", comment="#")
```

- The file name starts with the figure or table it feeds (`fig2a_...` is Figure 2 panel a, `figS10_...` is Supplementary Figure S10).
- Values were copied from the output tables of the analyses. A column or note says "derived" where a value was computed afterwards.
- `vf` means gene probes alone in the read-simulation files and the hybrid (PLUCK) in the scored real-data files. The `tool` column names the method.
- Scripts that draw the figures are in `../scripts/`. Folder constants at the top of each script point to the original working folders; set them for your copy.
