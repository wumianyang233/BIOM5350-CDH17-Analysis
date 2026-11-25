# Reproducibility guide

Run commands from the repository root. The workflow expects Python 3.14.4 and the packages pinned in `requirements.txt`.

## 1. Create the environment

```bash
python3 -m venv .venv
source .venv/bin/activate        # Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

## 2. Retrieve the CDH17 feature and phenotype annotations

```bash
python scripts/01_download_cdh17.py
```

This makes targeted requests to the UCSC Xena Toil `/data/` API for `ENSG00000079112.9` and the matching phenotype fields. It does not fetch the full expression matrix. It writes:

- `data/cdh17_expression_annotated.csv`
- `data/cdh17_expression_raw_all_samples.csv` (local intermediate; excluded from GitHub)
- `data/data_provenance.md`
- `results/data_qc_summary.csv`
- `results/data_qc_report.md`

The acquisition script atomically replaces these generated files after a successful retrieval and validation. Run it only when refreshing the public source extraction; it requires internet access.

## 3. Recreate descriptive statistics and Figures A-C

```bash
jupyter nbconvert --to notebook --execute --inplace notebooks/02_CDH17_EDA.ipynb
```

The notebook checks the schema, source annotations, group labels, missingness, expression floor, and tissue labels before summarizing the six analysis groups. It writes `results/eda_summary.csv` and:

- `figures/figure_A_tcga_colorectal_distributions.png`
- `figures/figure_B_tcga_tumors_vs_gtex_colon.png`
- `figures/figure_C_gtex_normal_tissues.png`

To run interactively instead, start Jupyter from the project root and open `notebooks/02_CDH17_EDA.ipynb`, then run all cells in order. If the notebook server is started elsewhere, set `CDH17_PROJECT_ROOT` to the project root before launching it.

## 4. Recreate the prespecified statistical analysis and Figure D

```bash
python scripts/03_statistical_analysis.py
```

The script independently selects samples from the preserved phenotype columns, verifies those selections against `analysis_group`, and runs two primary within-TCGA comparisons plus one secondary exploratory GTEx colon-subsite comparison. It writes:

- `results/statistical_summary.csv`
- `results/statistical_analysis_notes.md`
- `figures/figure_D_tcga_statistical_comparisons.png`

## Execution order and inputs

Run acquisition first, then the EDA notebook, then the statistical script. Both analysis steps consume `data/cdh17_expression_annotated.csv`; the notebook also verifies `data/data_provenance.md`. The statistical script needs no network connection. The scripts locate their project root relative to their own paths, and the notebook uses `CDH17_PROJECT_ROOT` or searches the current directory and its parents. No project inputs outside this repository are required; public UCSC Xena is the only external source.

The statistical script and notebook rewrite their listed outputs on execution. The acquisition script rewrites the listed data/QC outputs on successful completion. To preserve the validated versions, copy the project to a temporary directory before testing a full rerun.

## Reproducibility check note

A temporary end-to-end rerun reproduced the annotated data, descriptive-statistics CSV, statistical summary, statistical notes, and Figures C-D exactly. Figures A-B use randomly jittered sample points, so their PNG pixels can differ slightly across reruns even when the underlying data and summaries are unchanged. The validated Figures A-D were left untouched during the check.
