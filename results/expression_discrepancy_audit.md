# CDH17 expression discrepancy audit

**Project:** BIOM5350 Fall 2026 Final Project  
**Audit date:** 2026-11-25  
**Scope:** Independent descriptive audit of the saved Phase 3 annotated data and Phase 4 median findings. No hypothesis tests, batch correction, figures, or source-data edits were performed.

## Summary

The six reported medians reproduce exactly from `data/cdh17_expression_annotated.csv` at the saved precision. The 7,860 rows are classified into the same six groups and have no duplicate sample IDs, missing expression values, or failed joins. The Phase 3 provenance and acquisition script specify one Xena expression dataset, one feature (`ENSG00000079112.9`), and one reported scale (`log2(TPM + 0.001)`). The two expression columns in the CSV are identical row by row; `CDH17_expression_api_value` is the returned transformed value, not an untransformed TPM measurement.

The large TCGA-normal versus GTEx-colon median gap is real in this file: 8.243 versus 1.420 (6.823 units on the reported log2 scale). However, GTEx colon is internally heterogeneous: transverse colon has median 6.908, while sigmoid colon has median -1.056. The current TCGA normal annotations do not distinguish these colon subsites. Therefore, the pooled cross-source contrast cannot isolate a source effect from subsite or other cohort/collection differences. Treat it as exploratory pending review.

## Data and provenance checks

- Saved analysis file: 7,860 rows × 14 columns.
- Expression column: `CDH17_expression_log2_TPM_plus_0.001` (`float64`); duplicate API-return column: `CDH17_expression_api_value` (`float64`). The two columns match exactly for all rows.
- Feature recorded in the Phase 3 acquisition script and provenance: `ENSG00000079112.9` (Ensembl gene identifier with version suffix).
- Dataset recorded in both: expression `TcgaTargetGtex_rsem_gene_tpm`; phenotype `TcgaTargetGTEX_phenotype.txt`; UCSC Xena Toil hub `https://toil.xenahubs.net`.
- Documented expression scale: Toil RSEM TPM transformed as `log2(TPM + 0.001)`. The values were not converted or re-normalized in this audit.
- `sampleID`: 0 missing, 0 duplicates. `_join_status`: `both` for all 7,860 rows. Expression: 0 missing and all values parse as numeric.
- Source annotations were checked in `_study`, `_primary_site`, `_sample_type`, and `detailed_category`. All TCGA rows in the four TCGA groups have the expected site and sample type; tumor detailed labels are `Colon Adenocarcinoma` or `Rectum Adenocarcinoma` as appropriate. GTEx colon consists of 167 `Colon - Transverse` and 141 `Colon - Sigmoid` samples.
- Four retained GTEx non-colon samples have missing `_primary_site` but usable detailed tissue labels (`Stomach`, two `Esophagus - Mucosa`, and `Skin - Sun Exposed (Lower Leg)`). In addition, `_gender` is missing for 3 TCGA rows; `detailed_category` and `primary disease or tissue` are complete in this saved analysis file. They remain in Other GTEx normal tissues by the documented detailed-label fallback. No row lacks both usable tissue labels. Reconstructing groups while applying that fallback reproduces all six saved labels; there are no group mismatches.
- Group counts exactly match the Phase 3 feasibility counts; no count discrepancies require reconciliation.

## Independent six-group summaries

Expression units in this table are the dataset's reported `log2(TPM + 0.001)` values. The floor column counts values saved as `-9.966`, the rounded value expected for the zero-TPM floor. It is a floor-equivalent proxy, not proof that every such rounded value came from mathematically exact zero TPM.

| Analysis group | n | Median | Q1 | Q3 | Minimum | Maximum | Missing expression | Floor-equivalent n (%) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| TCGA colon primary tumors | 288 | 7.4535 | 6.9165 | 7.9783 | -0.394 | 9.469 | 0 | 0 (0.00%) |
| TCGA rectum primary tumors | 92 | 7.6750 | 7.2170 | 7.9890 | 5.171 | 9.059 | 0 | 0 (0.00%) |
| TCGA colon solid-tissue normals | 41 | 8.2430 | 7.9350 | 8.4110 | 5.390 | 8.955 | 0 | 0 (0.00%) |
| TCGA rectum solid-tissue normals | 10 | 7.9195 | 7.5003 | 8.2770 | 6.870 | 8.958 | 0 | 0 (0.00%) |
| GTEx normal colon | 308 | 1.4200 | -1.2830 | 7.1483 | -6.506 | 9.302 | 0 | 0 (0.00%) |
| Other GTEx normal tissues | 7,121 | -4.2930 | -6.5060 | -2.9320 | -9.966 | 8.047 | 0 | 1,447 (20.32%) |

The reported EDA medians all match. A useful descriptive contrast is the 6.823-unit difference between the TCGA colon-normal median (8.243) and pooled GTEx-colon median (1.420); this is not a source-adjusted or causal estimate.

## GTEx colon subsite audit

| Detailed source label | n | Median | Q1 | Q3 | Minimum | Maximum | Floor-equivalent n (%) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Colon - Sigmoid | 141 | -1.056 | -2.053 | 0.5568 | -6.506 | 7.781 | 0 (0.00%) |
| Colon - Transverse | 167 | 6.908 | 1.947 | 7.395 | -2.315 | 9.302 | 0 (0.00%) |

The median separation between the two GTEx colon labels is 7.964 expression units. These labels are retained separately because pooling them conceals substantial heterogeneity. No floor-equivalent GTEx colon observations were found. GTEx normal samples span 51 distinct detailed tissue categories in this file.

## TCGA distribution and outlier checks

Group-specific Tukey fences below are a descriptive flagging rule (1.5 × IQR), not a data-removal rule or a test. No flagged values were removed. Floor-equivalent values were absent in every TCGA group.

| TCGA group | n | Median | Q1 | Q3 | Range | Tukey fences | Flagged values |
|---|---:|---:|---:|---:|---:|---:|---:|
| Colon primary tumors | 288 | 7.4535 | 6.9165 | 7.9783 | -0.394 to 9.469 | 5.324 to 9.571 | 13 below lower fence; 0 above |
| Colon solid-tissue normals | 41 | 8.2430 | 7.9350 | 8.4110 | 5.390 to 8.955 | 7.221 to 9.125 | 3 below lower fence; 0 above |
| Rectum primary tumors | 92 | 7.6750 | 7.2170 | 7.9890 | 5.171 to 9.059 | 6.059 to 9.147 | 4 below lower fence; 0 above |
| Rectum solid-tissue normals | 10 | 7.9195 | 7.5003 | 8.2770 | 6.870 to 8.958 | 6.335 to 9.442 | 0 |

All Tukey flags are on the low side. The 13 colon-tumor flagged observations span -0.394 to 5.269; the lowest two are `TCGA-CK-6746-01` (-0.394) and `TCGA-CA-6718-01` (0.030). The three colon-normal flags span 5.390 to 7.002 (`TCGA-AZ-6605-11`, `TCGA-AA-3489-11`, `TCGA-AZ-6603-11`). The four rectum-tumor flags span 5.171 to 6.038 (`TCGA-DC-6155-01`, `TCGA-EI-7004-01`, `TCGA-DY-A1H8-01`, `TCGA-AF-2687-01`). These are unusual relative to their respective group IQRs, but the saved phenotype fields do not show an obvious group-label inconsistency for them. They should remain in the data unless later sample-level review identifies a documented exclusion reason.

## Detailed GTEx tissue categories with the ten highest medians

GTEx normal tissue categories are ranked by their own median; labels remain at their supplied anatomical resolution. These are not combined-organ summaries.

| Detailed GTEx category | n | Median | Q1 | Q3 |
|---|---:|---:|---:|---:|
| Colon - Transverse | 167 | 6.9080 | 1.9470 | 7.3950 |
| Small Intestine - Terminal Ileum | 92 | 6.2965 | 4.0318 | 7.1990 |
| Colon - Sigmoid | 141 | -1.0560 | -2.0530 | 0.5568 |
| Bladder | 9 | -1.7810 | -3.8160 | -1.5520 |
| Esophagus - Muscularis | 245 | -1.8310 | -2.6350 | -1.1810 |
| Spleen | 100 | -1.9380 | -2.8525 | -1.1170 |
| Fallopian Tube | 5 | -2.1780 | -2.8260 | -1.8310 |
| Artery - Tibial | 281 | -2.2450 | -2.9320 | -1.5100 |
| Testis | 165 | -2.3150 | -2.7270 | -1.8310 |
| Prostate | 100 | -2.3150 | -3.1710 | -1.6850 |

Some categories have small sample counts (for example, Fallopian Tube n=5), so their rank and median are imprecise descriptive summaries. This ranking is not a formal tissue-selectivity result.

## Interpreting the TCGA–GTEx normal-colon gap

### Documented in the available files and source documentation

- Every observation uses the same Xena Toil expression dataset and reported feature/scale recorded above. The Xena dataset metadata identifies a shared processing resource; that supports a common computational pipeline, not equivalence of specimen collection or a guarantee that all cohort/source effects have been removed.
- The source annotations separate TCGA from GTEx and distinguish TCGA `Solid Tissue Normal` from GTEx `Normal Tissue`. GTEx is documented as a postmortem donor tissue resource; TCGA COAD/READ belongs to a separate cancer biospecimen program. These are documented cohort-framework differences, not evidence that a particular collection mechanism caused the expression gap. TCGA colon normals in this file are 41 samples; GTEx colon normals are 308.
- GTEx colon is split into transverse and sigmoid labels, with markedly different observed distributions. The TCGA normal records have broad `_primary_site=Colon`; their detailed category is `Colon Adenocarcinoma`, so these saved fields do not allow equivalent sigmoid/transverse stratification.
- The available sample-level CSV does not contain donor ischemic interval, procurement details, collection site within each specimen, histologic cell composition, tumor purity, or a harmonized technical batch variable. It therefore cannot quantify those mechanisms.

### Plausible explanations that remain hypotheses

- **Subsite mixture:** the GTEx colon pool contains both sigmoid and transverse tissue, whose medians differ by 7.964 units. Different subsite composition could move the pooled GTEx median. The data do not provide an equivalent subsite split for TCGA normals, so this alone cannot explain or correct the cross-source contrast.
- **Collection and specimen handling:** TCGA cancer specimens and GTEx donor tissues arise from distinct cohort collection frameworks. Differences in procurement, ischemic interval, tissue handling, or sampling could affect bulk RNA measurements. Such variables are not present here, so no specific handling difference can be assigned as the cause.
- **Tissue composition and pathology:** bulk expression can vary with epithelial, stromal, immune, and other cell proportions; TCGA tumors also vary in tumor content. The current table lacks cell-composition and tumor-purity measures suitable for evaluating this explanation.
- **Residual source/batch effects:** uniform Toil processing reduces some avoidable computational differences, but source remains aligned with cohort and specimen type. Without documented batch covariates or a design that crosses source with tissue/status, the observed TCGA–GTEx contrast cannot separate technical source effects from biological and collection differences.

The sigmoid/transverse result is the clearest immediate audit finding: the single pooled GTEx-colon summary hides strong internal variation. It raises a subsite/composition question; it does not establish a cause for that variation.

## Defensible next analysis plan

1. **Primary comparisons within TCGA:** analyze colon primary tumors versus colon solid-tissue normals, and rectum primary tumors versus rectum solid-tissue normals, as two separate contrasts on the unchanged `log2(TPM + 0.001)` scale. Do not pair samples unless patient barcodes independently verify matching individuals. Report group distributions and effect sizes with confidence intervals. A rank-based two-sample procedure is a reasonable primary inferential option given low-side outliers and unequal group sizes; a Welch mean comparison can be a sensitivity analysis. Interpret rank tests as distributional comparisons, not automatically as tests of median equality. Correct the two planned primary comparisons for multiplicity (for example, Holm adjustment). The rectum-normal group is small (n=10), so emphasize uncertainty.
2. **Retain all observations:** do not discard the flagged values solely because they meet a Tukey rule. Inspect sample identifiers and available clinical/pathology covariates if they become available. Keep the floor-equivalent count visible; do not impute floor values or convert the scale.
3. **Cross-source analyses:** keep TCGA and GTEx contrasts exploratory. Report sigmoid and transverse GTEx colon separately, and do not present pooled GTEx colon as a matched normal control. If a later model includes source, tissue subsite, or sample class, first check whether the design can identify those effects; the current TCGA broad-site labels and missing collection/batch covariates limit that possibility. Seek additional matched or better annotated data if a source-adjusted claim is important.
4. **Reproducibility and interpretation:** preserve raw phenotype columns and detailed tissue labels, pre-specify contrasts, report counts/missingness/floor proportions, and avoid interpreting RNA abundance as protein abundance or clinical efficacy/safety.

## Sources and audit trail

- Local Phase 3 retrieval and dataset details: [`data/data_provenance.md`](../data/data_provenance.md) and [`scripts/01_download_cdh17.py`](../scripts/01_download_cdh17.py). The script specifies `TcgaTargetGtex_rsem_gene_tpm`, `ENSG00000079112.9`, and targeted retrieval of the single feature.
- UCSC Xena Toil hub: https://toil.xenahubs.net/ . Dataset metadata records the Toil release as 2016-09-03; the local provenance document records the expression and phenotype dataset IDs and reported scale.
- GTEx Portal documentation: https://gtexportal.org/home/documentationPage . GTEx Consortium atlas paper: https://pubmed.ncbi.nlm.nih.gov/32913098/ (Science, 2020, doi:10.1126/science.aaz1776).
- TCGA colon and rectal cancer cohort paper: https://pubmed.ncbi.nlm.nih.gov/22810696/ (Nature, 2012, doi:10.1038/nature11252).
- Toil workflow paper: Vivian et al., “Toil enables reproducible, open source, big biomedical data analyses,” https://doi.org/10.1371/journal.pcbi.1005397 .

This audit used the saved annotated CSV and the local Phase 3 provenance/script for data checks. It did not re-download expression data or alter the CSV, acquisition script, Phase 4 notebook, or figures.
