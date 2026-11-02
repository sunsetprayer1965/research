# Chromosome linking analysis

## Workflow

1. **Map interactions to genomic bins.** `locus.py` uses the BED intervals to map intra- and inter-chromosome interaction records to 10 kb bins, then writes the contact matrix used by Pastis.
2. **Reconstruct 3D chromosome structures.** `run_all.ps1` runs Pastis with PM2, PM1, MDS, and NMDS for each dataset under `dataset_linking_number/`.
3. **Prepare chromosome coordinates.** `coordinates_blank.py` formats the coordinates; `java/chromosome/split.py` writes each chromosome to a separate numbered file.
4. **Estimate linking frequency.** `java/computewritec/ComputeWriteC3.java` compares chromosome pairs using random closure sampling and writes one percentage per pair under `runs/<dataset>/linking/`.
5. **Run optional statistical follow-up.** `EB.py` performs empirical Bayes estimation. `semi-threshold.py` and `semi2.py` compare two sets of linking proportions using semi-soft thresholding. These are separate analyses and are not automatically run by `run_all.ps1`.

## Panel 1 vs Panel 2 example

Both semi-soft scripts read `panel_1_max.txt` and `panel_2_max.txt` from the same directory as the scripts. The files contain percentages; the scripts divide values by 100 before analysis. Run either script from PowerShell in this directory:

```powershell
python .\semi-threshold.py
python .\semi2.py
```

The output below was provided for **`semi-threshold.py`** using those two files:

| Scale | Rule | c-hat | Minimum F-hat | Near-zero after shrinkage | Baseline near-zero |
|---|---:|---:|---:|---:|---:|
| Untransformed | G1 | -0.001 | -5.6763 | 3.33% | 2.50% |
| Untransformed | G2 | 0.001 | -5.6763 | 3.33% | 2.50% |
| Transformed | G1 | -0.001 | -7.5970 | 1.67% | Not reported |
| Transformed | G2 | 0.001 | -7.5970 | 1.67% | Not reported |

The untransformed results report 3.33% near-zero differences after shrinkage, compared with 2.50% before shrinkage. The transformed results report 1.67% near-zero differences after shrinkage. G1 and G2 have the same displayed values at the shown precision for each scale.

**Parameter-range note:** the current `semi-threshold.py` call searches `c_range=(-10, 10)`. Since `c` is a positive scale parameter, this range includes invalid negative values and values very close to zero. The reported `c-hat` values of ±0.001 are at the grid's nearest points to zero, so treat them as a parameter-range artifact, not a reliable optimum. Constrain the search to positive `c` values before interpreting or reporting `c-hat` as a scientific result.
