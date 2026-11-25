# CDH17 data provenance and retrieval

## Source and datasets

- UCSC Xena hub: https://toil.xenahubs.net
- Xena data API: https://toil.xenahubs.net/data/
- Expression dataset ID/name: `TcgaTargetGtex_rsem_gene_tpm`
- Matching phenotype dataset ID/name: `TcgaTargetGTEX_phenotype.txt`
- Expression file endpoint (reference only; the complete file was not downloaded): https://toil.xenahubs.net/download/TcgaTargetGtex_rsem_gene_tpm.gz
- Phenotype file endpoint (reference only; annotations were fetched through the API): https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz
- Gene feature: `ENSG00000079112.9`, Ensembl gene ID for CDH17 (version suffix retained as present in the matrix).
- Expression scale: Toil RSEM TPM, documented as `log2(TPM + 0.001)`. Values are kept exactly as returned by Xena; no conversion to `log2(TPM + 1)`, additional normalization, or batch correction is applied.
- Toil release metadata: `2016-09-03` (dataset metadata).
- Access date: 2026-11-25 (UTC date at the most recent script run; the QC report records the same run date).

## Retrieved annotations

The script fetches the matching `sampleID` plus these original phenotype field names from `TcgaTargetGTEX_phenotype.txt`: `_study, _primary_site, _sample_type, detailed_category, primary disease or tissue, _gender`. Categorical API codes are decoded using the same dataset's Xena code tables, and the source fields are retained in the CSV. `analysis_group`, `anatomical_tissue_label`, `exclusion_reason`, and join/QC helper fields are derived fields; they do not replace the source phenotype columns.

## Retrieval and join procedure

1. Query Xena's `/data/` API for sample IDs advertised by the expression and phenotype datasets.
2. Calculate and record the shared sample IDs. Independently request only feature `ENSG00000079112.9` for every expression sample and the annotation fields listed above for every phenotype sample, in batches of 500 samples, so unmatched IDs are preserved for QC.
3. The `fetch` endpoint returns requested values in the order of the supplied sample IDs. Each response vector is immediately labeled with those exact IDs. The expression and annotation tables are then outer-joined explicitly on `sampleID` with a one-to-one validation.
4. Save all fetched rows, including excluded groups, to `data/cdh17_expression_raw_all_samples.csv`; save only the six defined groups to `data/cdh17_expression_annotated.csv`.

Reproduce with `python scripts/01_download_cdh17.py` from this project directory. The script uses `requests` and `pandas`, applies request timeouts and up to four attempts with backoff, validates response dimensions, checks duplicate IDs and joins, and atomically replaces its own output files only after the retrieval and validation steps succeed. It does not download the full expression matrix.

## Sample group rules

- TCGA colon primary tumors: `_study=TCGA`, `_primary_site=Colon`, `_sample_type=Primary Tumor`, `detailed_category=Colon Adenocarcinoma`.
- TCGA rectum primary tumors: `_study=TCGA`, `_primary_site=Rectum`, `_sample_type=Primary Tumor`, `detailed_category=Rectum Adenocarcinoma`.
- TCGA colon/rectum solid-tissue normals: matching `_primary_site` and `_sample_type=Solid Tissue Normal`.
- GTEx normal colon: `_study=GTEX`, `_sample_type=Normal Tissue`, and `_primary_site=Colon`; if broad site is missing, a source detailed tissue label beginning `Colon` is used. Source `detailed_category` values (including sigmoid/transverse where supplied) are retained.
- Other GTEx normal tissues: `_study=GTEX`, `_sample_type=Normal Tissue`, and a non-Colon `_primary_site`; if broad site is missing, a nonblank detailed tissue label identifies the tissue. Detailed labels are retained.
- TARGET, non-normal GTEx, metastatic/recurrent/other colorectal TCGA types, non-colorectal TCGA, and ambiguous or unclassified records are excluded from the main file but retained in the raw intermediate when returned by the datasets.
- TCGA normal samples are not treated as paired to tumors based on these annotations.

## Documentation and limitations

Xena hub: https://toil.xenahubs.net

Xena dataset metadata and downloadable references can be inspected via the hub UI/API. The full expression file is roughly 1.32 GB compressed according to the prior range-header feasibility check; this workflow makes targeted `/data/` API requests instead. The phenotype file's prior range-header size was 135,753 bytes. These are reference sizes, not downloaded files. Verify current endpoint behavior if reproducing at a later date.

Toil's shared processing improves technical consistency but does not guarantee removal of TCGA/GTEx cohort and collection effects. The source fields should be retained for downstream design decisions. No statistical comparisons are made by this acquisition script.
