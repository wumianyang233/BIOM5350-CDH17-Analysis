# CDH17 statistical analysis notes

## Methods

The script selects samples from the preserved `_study`, `_primary_site`, `_sample_type`, and `detailed_category` fields, then checks the resulting rows against `analysis_group`. Expression values are used unchanged on the reported `log2(TPM + 0.001)` scale. All comparisons are unpaired; no patient matching or pairing is used.

Each comparison uses a two-sided Mann–Whitney U test. The reported median difference is group 1 minus group 2. Rank-biserial correlation is defined as `2 × U1 / (n1 × n2) − 1`, where U1 is SciPy’s Mann–Whitney statistic for group 1; positive values indicate that group 1 tends to rank higher than group 2, with ties receiving half credit. Holm adjustment is applied across the two prespecified primary TCGA tests only. The GTEx comparison is exploratory and its p-value is unadjusted.

## Sample selection and results

| Role | Comparison (group 1 vs group 2) | n1 | n2 | Median 1 | Median 2 | Median difference | Rank-biserial r | Raw two-sided p | Holm-adjusted p |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| primary | TCGA colon primary tumors vs TCGA colon solid-tissue normals | 288 | 41 | 7.4535 | 8.2430 | -0.7895 | -0.5812 | 1.73368e-09 | 3.46737e-09 |
| primary | TCGA rectum primary tumors vs TCGA rectum solid-tissue normals | 92 | 10 | 7.6750 | 7.9195 | -0.2445 | -0.2293 | 0.237362 | 0.237362 |
| secondary_exploratory | GTEx normal sigmoid colon vs GTEx normal transverse colon | 141 | 167 | -1.0560 | 6.9080 | -7.9640 | -0.7317 | 1.8915e-28 | — |

## Interpretation and limitations

The two TCGA comparisons are the prespecified primary analyses; the GTEx sigmoid-versus-transverse comparison is secondary and exploratory. A nonsignificant p-value should not be interpreted as evidence that groups are equivalent. Median differences describe the observed sample medians on the transformed expression scale; the rank-biserial statistic describes relative rank tendency and is not a fold change.

The TCGA solid-tissue-normal groups are modest in size, especially rectum normals, so estimates and p-values may be imprecise. Mann–Whitney tests compare distributions/ranks and do not automatically test equality of medians if distribution shapes differ. TCGA and GTEx samples come from different cohorts; the GTEx result is not a substitute for a TCGA matched control and is not adjusted for source effects. No batch correction, imputation, or sample exclusion was applied.

This is a focused course-project analysis of RNA expression. It does not establish protein abundance, ADC efficacy, or clinical safety.

## Reproduction

From the project root run `python scripts/03_statistical_analysis.py`. The script reads the saved annotated CSV and rewrites only `results/statistical_summary.csv`, this notes file, and `figures/figure_D_tcga_statistical_comparisons.png`.
