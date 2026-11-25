#!/usr/bin/env python3
"""Retrieve CDH17 expression and matching Toil Xena phenotype annotations.

This script requests only the CDH17 feature and selected phenotype fields from
the Xena ``/data/`` API. It never downloads the complete expression matrix.
Run from the project root with:

    python scripts/01_download_cdh17.py

The script also supports an absolute path to ``scripts/01_download_cdh17.py``
when invoked from another working directory.

The output files are written relative to this script's project root. Existing
outputs from this script are replaced atomically after a successful retrieval.
"""

from __future__ import annotations

import json
import math
import os
import sys
import tempfile
import time
from datetime import date
from pathlib import Path
from typing import Any

import pandas as pd
import requests


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"

HUB_URL = "https://toil.xenahubs.net"
DATA_API = f"{HUB_URL}/data/"
EXPRESSION_DATASET = "TcgaTargetGtex_rsem_gene_tpm"
PHENOTYPE_DATASET = "TcgaTargetGTEX_phenotype.txt"
GENE_ID = "ENSG00000079112.9"
EXPRESSION_COLUMN = "CDH17_expression_log2_TPM_plus_0.001"
EXPRESSION_RAW_COLUMN = "CDH17_expression_api_value"

PHENOTYPE_FIELDS = [
    "_study",
    "_primary_site",
    "_sample_type",
    "detailed_category",
    "primary disease or tissue",
    "_gender",
]
SAMPLE_FIELD = "sampleID"
FETCH_BATCH_SIZE = 500
MAX_ATTEMPTS = 4
TIMEOUT = (20, 120)

GROUPS = [
    "TCGA colon primary tumors",
    "TCGA rectum primary tumors",
    "TCGA colon solid-tissue normals",
    "TCGA rectum solid-tissue normals",
    "GTEx normal colon",
    "Other GTEx normal tissues",
]

# Counts reported in the feasibility check, used only as a comparison point.
FEASIBILITY_COUNTS = {
    "TCGA colon primary tumors": 288,
    "TCGA rectum primary tumors": 92,
    "TCGA colon solid-tissue normals": 41,
    "TCGA rectum solid-tissue normals": 10,
    "GTEx normal colon": 308,
    "Other GTEx normal tissues": 7121,
}


class XenaError(RuntimeError):
    """Raised when Xena returns invalid or incomplete data."""


def request_xena(query: str, description: str) -> Any:
    """POST a read-only Xena query with timeout and bounded retry handling."""
    last_error: Exception | None = None
    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            response = requests.post(
                DATA_API,
                data=query.encode("utf-8"),
                headers={"Content-Type": "text/plain; charset=utf-8"},
                timeout=TIMEOUT,
            )
            if response.status_code == 429 or response.status_code >= 500:
                raise requests.HTTPError(
                    f"HTTP {response.status_code}: {response.text[:300]}",
                    response=response,
                )
            response.raise_for_status()
            try:
                return response.json()
            except ValueError as exc:
                raise XenaError(
                    f"Xena returned invalid JSON for {description}: {response.text[:300]}"
                ) from exc
        except (requests.RequestException, XenaError) as exc:
            last_error = exc
            if attempt == MAX_ATTEMPTS:
                break
            time.sleep(min(2 ** (attempt - 1), 8))
    raise XenaError(
        f"Could not retrieve {description} from {DATA_API} after {MAX_ATTEMPTS} attempts: "
        f"{last_error}"
    )


def qstr(value: str) -> str:
    """Return a safely quoted string literal for the Xena query language."""
    return json.dumps(value, ensure_ascii=False)


def get_dataset_sample_ids(dataset: str) -> list[str]:
    query = (
        "(query {:select [:code.value] :from [:dataset] "
        ":join [:field [:= :dataset.id :dataset_id] :code [:= :field.id :field_id]] "
        f":where [:and [:= :dataset.name {qstr(dataset)}] [:= :field.name {qstr(SAMPLE_FIELD)}]]}})"
    )
    response = request_xena(query, f"sample IDs for {dataset}")
    if not isinstance(response, list):
        raise XenaError(f"Unexpected sample-ID response for {dataset}: {type(response).__name__}")
    ids = [item.get("value") if isinstance(item, dict) else None for item in response]
    if any(not isinstance(sample_id, str) or not sample_id for sample_id in ids):
        raise XenaError(f"Dataset {dataset} returned missing or malformed sample IDs.")
    return ids


def get_codebook() -> dict[str, list[str]]:
    field_list = " ".join(qstr(name) for name in PHENOTYPE_FIELDS)
    query = (
        "(query {:select [:field.name :code.ordering :code.value] :from [:field] "
        ":join [:dataset [:= :dataset.id :dataset_id] :code [:= :field.id :field_id]] "
        f":where [:and [:= :dataset.name {qstr(PHENOTYPE_DATASET)}] "
        f"[:in :field.name [{field_list}]]] :order-by [:field.name :code.ordering]}})"
    )
    response = request_xena(query, "phenotype code tables")
    if not isinstance(response, list):
        raise XenaError("Unexpected phenotype code-table response.")

    grouped: dict[str, list[tuple[int, str]]] = {field: [] for field in PHENOTYPE_FIELDS}
    for item in response:
        if not isinstance(item, dict):
            raise XenaError(f"Malformed code-table row: {item!r}")
        field = item.get("name")
        if field not in grouped:
            continue
        try:
            ordering = int(item["ordering"])
        except (KeyError, TypeError, ValueError) as exc:
            raise XenaError(f"Invalid code ordering for phenotype field {field!r}.") from exc
        value = item.get("value")
        grouped[field].append((ordering, "" if value is None else str(value)))

    codebook: dict[str, list[str]] = {}
    for field, entries in grouped.items():
        entries.sort(key=lambda pair: pair[0])
        if not entries:
            raise XenaError(f"Xena returned no code table for phenotype field {field!r}.")
        # The fetch API returns zero-based code indices; ordering must be contiguous.
        orderings = [ordering for ordering, _ in entries]
        if orderings != list(range(len(entries))):
            raise XenaError(f"Non-contiguous Xena code ordering for phenotype field {field!r}.")
        codebook[field] = [value for _, value in entries]
    return codebook


def decode_code(value: Any, labels: list[str], field: str, sample_id: str) -> str:
    if value is None or value == "":
        return ""
    try:
        numeric = float(value)
        if math.isnan(numeric):
            return ""
        if not numeric.is_integer():
            raise ValueError
        index = int(numeric)
    except (TypeError, ValueError) as exc:
        raise XenaError(
            f"Non-integer code {value!r} for phenotype field {field!r}, sample {sample_id}."
        ) from exc
    if index < 0 or index >= len(labels):
        raise XenaError(
            f"Code {index} is outside the code table for phenotype field {field!r}, "
            f"sample {sample_id}."
        )
    return labels[index]


def fetch_gene_values(sample_ids: list[str]) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    for start in range(0, len(sample_ids), FETCH_BATCH_SIZE):
        batch = sample_ids[start : start + FETCH_BATCH_SIZE]
        query = (
            f"(fetch [{{:table {qstr(EXPRESSION_DATASET)} :columns [{qstr(GENE_ID)}] "
            f":samples [{ ' '.join(qstr(sample_id) for sample_id in batch) }]}}])"
        )
        result = request_xena(query, f"{GENE_ID} expression values, samples {start + 1}-{start + len(batch)}")
        if not isinstance(result, list) or len(result) != 1 or not isinstance(result[0], list):
            raise XenaError(f"Unexpected expression response for sample batch starting at {start}.")
        values = result[0]
        if len(values) != len(batch):
            raise XenaError(
                f"Expression response has {len(values)} values for {len(batch)} requested samples "
                f"(batch starting at {start})."
            )
        records.extend(
            {SAMPLE_FIELD: sample_id, EXPRESSION_RAW_COLUMN: value}
            for sample_id, value in zip(batch, values, strict=True)
        )
    return pd.DataFrame.from_records(records)


def fetch_phenotypes(sample_ids: list[str], codebook: dict[str, list[str]]) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    fields = " ".join(qstr(field) for field in PHENOTYPE_FIELDS)
    for start in range(0, len(sample_ids), FETCH_BATCH_SIZE):
        batch = sample_ids[start : start + FETCH_BATCH_SIZE]
        sample_expr = " ".join(qstr(sample_id) for sample_id in batch)
        query = (
            f"(fetch [{{:table {qstr(PHENOTYPE_DATASET)} :columns [{fields}] "
            f":samples [{sample_expr}]}}])"
        )
        result = request_xena(
            query,
            f"phenotype annotations, samples {start + 1}-{start + len(batch)}",
        )
        if not isinstance(result, list) or len(result) != len(PHENOTYPE_FIELDS):
            raise XenaError(f"Unexpected phenotype response for sample batch starting at {start}.")
        if any(not isinstance(row, list) or len(row) != len(batch) for row in result):
            raise XenaError(
                f"Phenotype response dimensions do not match {len(PHENOTYPE_FIELDS)} fields "
                f"and {len(batch)} requested samples (batch starting at {start})."
            )
        for sample_index, sample_id in enumerate(batch):
            record: dict[str, Any] = {SAMPLE_FIELD: sample_id}
            for field_index, field in enumerate(PHENOTYPE_FIELDS):
                record[field] = decode_code(
                    result[field_index][sample_index], codebook[field], field, sample_id
                )
            records.append(record)
    return pd.DataFrame.from_records(records, columns=[SAMPLE_FIELD, *PHENOTYPE_FIELDS])


def classify_row(row: pd.Series) -> tuple[str, str]:
    """Return analysis group and a mutually exclusive reason when excluded."""
    def field_text(name: str) -> str:
        value = row.get(name, "")
        return "" if pd.isna(value) else str(value).strip()

    study = field_text("_study")
    site = field_text("_primary_site")
    sample_type = field_text("_sample_type")
    detail = field_text("detailed_category")

    if study == "TCGA":
        if site == "Colon" and sample_type == "Primary Tumor" and detail == "Colon Adenocarcinoma":
            return "TCGA colon primary tumors", ""
        if site == "Rectum" and sample_type == "Primary Tumor" and detail == "Rectum Adenocarcinoma":
            return "TCGA rectum primary tumors", ""
        if site == "Colon" and sample_type == "Solid Tissue Normal":
            return "TCGA colon solid-tissue normals", ""
        if site == "Rectum" and sample_type == "Solid Tissue Normal":
            return "TCGA rectum solid-tissue normals", ""
        if site in {"Colon", "Rectum"}:
            return "", "TCGA colorectal metastatic, recurrent, other, or ambiguous sample"
        return "", "TCGA non-colorectal sample"

    if study == "GTEX":
        if sample_type != "Normal Tissue":
            return "", "GTEx non-normal or unclassified sample"
        if site == "Colon":
            return "GTEx normal colon", ""
        if site:
            return "Other GTEx normal tissues", ""
        # The phenotype table has a small number of normal GTEx rows with a
        # blank broad site field but a specific detailed tissue label. Retain
        # these identifiable tissues; leave them blank only if both source
        # labels fail to identify anatomy.
        detailed_site = detail or field_text("primary disease or tissue")
        if detailed_site:
            if detailed_site.casefold().startswith("colon"):
                return "GTEx normal colon", ""
            return "Other GTEx normal tissues", ""
        return "", "GTEx normal tissue with unclassified primary site"

    if study == "TARGET":
        return "", "TARGET sample"
    if not study:
        return "", "Missing study/source annotation"
    return "", "Other or unclassified study/source"


def atomic_write_csv(frame: pd.DataFrame, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(prefix=f".{target.name}.", suffix=".tmp", dir=target.parent)
    os.close(fd)
    temp_path = Path(temp_name)
    try:
        frame.to_csv(temp_path, index=False, lineterminator="\n")
        os.replace(temp_path, target)
    finally:
        temp_path.unlink(missing_ok=True)


def atomic_write_text(content: str, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(prefix=f".{target.name}.", suffix=".tmp", dir=target.parent)
    os.close(fd)
    temp_path = Path(temp_name)
    try:
        temp_path.write_text(content, encoding="utf-8")
        os.replace(temp_path, target)
    finally:
        temp_path.unlink(missing_ok=True)


def fmt_bytes(path: Path) -> str:
    return str(path.stat().st_size) if path.exists() else "not created"


def main() -> int:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Querying Xena hub {HUB_URL}; no full matrix download will be made.", flush=True)
    expression_ids = get_dataset_sample_ids(EXPRESSION_DATASET)
    phenotype_ids = get_dataset_sample_ids(PHENOTYPE_DATASET)
    expression_duplicates = pd.Series(expression_ids).duplicated(keep=False)
    phenotype_duplicates = pd.Series(phenotype_ids).duplicated(keep=False)
    expression_duplicate_ids = sorted(set(pd.Series(expression_ids)[expression_duplicates]))
    phenotype_duplicate_ids = sorted(set(pd.Series(phenotype_ids)[phenotype_duplicates]))

    expression_set, phenotype_set = set(expression_ids), set(phenotype_ids)
    shared_ids = sorted(expression_set & phenotype_set)
    if not shared_ids:
        raise XenaError("The expression and phenotype datasets have no matching sample IDs.")

    print(
        f"Dataset sample IDs: expression={len(expression_ids):,}, phenotype={len(phenotype_ids):,}, "
        f"shared={len(shared_ids):,}; fetching each dataset's IDs independently in batches "
        f"of {FETCH_BATCH_SIZE}.",
        flush=True,
    )
    codebook = get_codebook()
    expression = fetch_gene_values(sorted(expression_ids))
    phenotypes = fetch_phenotypes(sorted(phenotype_ids), codebook)

    # Explicit keyed merge; Xena fetch vectors are first labeled with the exact
    # requested IDs, then joined by sampleID with one-to-one validation.
    if expression[SAMPLE_FIELD].duplicated().any() or phenotypes[SAMPLE_FIELD].duplicated().any():
        raise XenaError("Duplicate sample IDs appeared in fetched rows; refusing an ambiguous join.")
    joined = expression.merge(
        phenotypes,
        on=SAMPLE_FIELD,
        how="outer",
        validate="one_to_one",
        indicator="_join_status",
        sort=True,
    )
    joined["_join_status"] = joined["_join_status"].astype(str)
    joined[EXPRESSION_COLUMN] = joined[EXPRESSION_RAW_COLUMN]

    # All selected phenotype fields are retained with their original field names.
    classifications = joined.apply(classify_row, axis=1, result_type="expand")
    classifications.columns = ["analysis_group", "exclusion_reason"]
    joined = pd.concat([joined, classifications], axis=1)

    numeric_expression = pd.to_numeric(joined[EXPRESSION_RAW_COLUMN], errors="coerce")
    raw_not_missing = joined[EXPRESSION_RAW_COLUMN].notna() & joined[EXPRESSION_RAW_COLUMN].astype(str).ne("")
    nonnumeric_mask = raw_not_missing & numeric_expression.isna()
    numeric_missing_mask = ~raw_not_missing
    joined["_expression_numeric_check"] = numeric_expression

    main = joined[(joined["_join_status"] == "both") & joined["analysis_group"].isin(GROUPS)].copy()
    # Add group-specific anatomical labels while retaining every source field.
    main["anatomical_tissue_label"] = main["detailed_category"].fillna("")
    no_detail = main["anatomical_tissue_label"].eq("")
    main.loc[no_detail, "anatomical_tissue_label"] = main.loc[no_detail, "primary disease or tissue"].fillna("")

    group_counts = main["analysis_group"].value_counts().to_dict()
    main_count_by_group = {group: int(group_counts.get(group, 0)) for group in GROUPS}
    all_group_rows = joined[joined["analysis_group"].isin(GROUPS)]
    group_missing = {
        group: int(
            ((all_group_rows["analysis_group"] == group) & numeric_missing_mask).sum()
        )
        for group in GROUPS
    }
    group_nonnumeric = {
        group: int(
            ((all_group_rows["analysis_group"] == group) & nonnumeric_mask).sum()
        )
        for group in GROUPS
    }

    join_counts = joined["_join_status"].value_counts().to_dict()
    exclusion_counts = joined.loc[~joined["analysis_group"].isin(GROUPS), "exclusion_reason"].value_counts().to_dict()
    unexpected = []
    if expression_duplicate_ids:
        unexpected.append(f"Expression sample IDs contain duplicates ({len(expression_duplicate_ids)} distinct IDs).")
    if phenotype_duplicate_ids:
        unexpected.append(f"Phenotype sample IDs contain duplicates ({len(phenotype_duplicate_ids)} distinct IDs).")
    if join_counts.get("left_only", 0) or join_counts.get("right_only", 0):
        unexpected.append(
            f"Join is incomplete: left_only={join_counts.get('left_only', 0)}, "
            f"right_only={join_counts.get('right_only', 0)}."
        )
    if numeric_missing_mask.sum():
        unexpected.append(f"Expression is missing for {int(numeric_missing_mask.sum())} joined rows.")
    if nonnumeric_mask.sum():
        unexpected.append(f"Expression has {int(nonnumeric_mask.sum())} nonnumeric non-missing values.")
    if len(expression_ids) != len(set(expression_ids)):
        unexpected.append("Expression dataset sample count includes duplicate IDs.")
    if len(phenotype_ids) != len(set(phenotype_ids)):
        unexpected.append("Phenotype dataset sample count includes duplicate IDs.")
    for group, expected in FEASIBILITY_COUNTS.items():
        if main_count_by_group[group] != expected:
            unexpected.append(
                f"{group}: live count {main_count_by_group[group]} differs from feasibility count {expected}."
            )

    summary_rows: list[dict[str, Any]] = []

    def add(metric: str, value: Any, details: str = "") -> None:
        summary_rows.append({"metric": metric, "value": value, "details": details})

    add("expression_dataset_sample_ids", len(expression_ids), EXPRESSION_DATASET)
    add("phenotype_dataset_sample_ids", len(phenotype_ids), PHENOTYPE_DATASET)
    add("shared_sample_ids", len(shared_ids), "Unique intersection used for targeted fetch")
    add("expression_dataset_duplicate_ids", len(expression_duplicate_ids), ", ".join(expression_duplicate_ids[:20]))
    add("phenotype_dataset_duplicate_ids", len(phenotype_duplicate_ids), ", ".join(phenotype_duplicate_ids[:20]))
    add("expression_values_fetched", len(expression))
    add("phenotype_rows_fetched", len(phenotypes))
    add("joined_rows", len(joined), "Outer merge on sampleID, validate=one_to_one")
    for status in ["both", "left_only", "right_only"]:
        add(f"join_{status}_rows", int(join_counts.get(status, 0)))
    add("join_success_rate_percent", round(100 * join_counts.get("both", 0) / max(len(joined), 1), 6))
    add("expression_missing_rows", int(numeric_missing_mask.sum()), "Missing values are retained, not imputed")
    add("expression_nonnumeric_rows", int(nonnumeric_mask.sum()), "Nonnumeric values are retained in raw expression column")
    for field in PHENOTYPE_FIELDS:
        add(
            f"phenotype_missing::{field}",
            int(joined[field].isna().sum() + joined[field].astype(str).eq("").sum()),
            "Missing or blank source annotation values; rows are retained",
        )
    add("annotated_main_dataset_rows", len(main), "Only six requested groups and successful joins")
    add("annotated_main_dataset_columns", len(main.columns), "Includes raw expression, annotations, group and QC helpers")
    add("raw_intermediate_rows", len(joined), "All expression/phenotype rows represented by outer join")
    for group in GROUPS:
        add(f"group_rows::{group}", main_count_by_group[group], f"feasibility reference={FEASIBILITY_COUNTS[group]}")
        add(f"group_expression_missing::{group}", group_missing[group])
        add(f"group_expression_nonnumeric::{group}", group_nonnumeric[group])
    for group in GROUPS:
        fallback_count = int(
            (
                (main["analysis_group"] == group)
                & main["_primary_site"].fillna("").astype(str).eq("")
            ).sum()
        )
        add(f"group_primary_site_missing_but_detailed_label_used::{group}", fallback_count)
    for reason, count in sorted(exclusion_counts.items()):
        add(f"excluded::{reason}", int(count))
    add("access_date_utc", date.today().isoformat(), "Date at script execution; environment-local date")

    main_columns = [
        SAMPLE_FIELD,
        EXPRESSION_COLUMN,
        EXPRESSION_RAW_COLUMN,
        *PHENOTYPE_FIELDS,
        "analysis_group",
        "anatomical_tissue_label",
        "exclusion_reason",
        "_join_status",
        "_expression_numeric_check",
    ]
    raw_columns = [
        SAMPLE_FIELD,
        EXPRESSION_COLUMN,
        EXPRESSION_RAW_COLUMN,
        *PHENOTYPE_FIELDS,
        "analysis_group",
        "exclusion_reason",
        "_join_status",
        "_expression_numeric_check",
    ]
    main = main.reindex(columns=main_columns).sort_values(SAMPLE_FIELD, kind="stable")
    raw = joined.reindex(columns=raw_columns).sort_values(SAMPLE_FIELD, kind="stable")
    summary = pd.DataFrame(summary_rows)

    report = build_report(
        expression_ids,
        phenotype_ids,
        shared_ids,
        expression_duplicate_ids,
        phenotype_duplicate_ids,
        joined,
        main,
        join_counts,
        main_count_by_group,
        group_missing,
        group_nonnumeric,
        exclusion_counts,
        unexpected,
    )

    # Stage each output to a temporary sibling and replace only after all
    # retrieval, validation, and report rendering succeeded.
    atomic_write_csv(raw, DATA_DIR / "cdh17_expression_raw_all_samples.csv")
    atomic_write_csv(main, DATA_DIR / "cdh17_expression_annotated.csv")
    atomic_write_csv(summary, RESULTS_DIR / "data_qc_summary.csv")
    atomic_write_text(report, RESULTS_DIR / "data_qc_report.md")
    atomic_write_text(build_provenance(), DATA_DIR / "data_provenance.md")

    print("\nRetrieval and QC complete.")
    print(f"Main dataset: {len(main):,} rows × {len(main.columns)} columns")
    print(f"Raw intermediate: {len(raw):,} rows × {len(raw.columns)} columns")
    print("Main group counts:")
    for group in GROUPS:
        print(f"  {group}: {main_count_by_group[group]:,}")
    for path in [
        DATA_DIR / "cdh17_expression_annotated.csv",
        DATA_DIR / "cdh17_expression_raw_all_samples.csv",
        RESULTS_DIR / "data_qc_summary.csv",
        RESULTS_DIR / "data_qc_report.md",
        DATA_DIR / "data_provenance.md",
    ]:
        print(f"{path.relative_to(PROJECT_ROOT)}: {fmt_bytes(path)} bytes")
    if unexpected:
        print("Unexpected results: " + " | ".join(unexpected))
    else:
        print("Unexpected results: none detected.")
    return 0


def build_report(
    expression_ids: list[str],
    phenotype_ids: list[str],
    shared_ids: list[str],
    expression_duplicate_ids: list[str],
    phenotype_duplicate_ids: list[str],
    joined: pd.DataFrame,
    main: pd.DataFrame,
    join_counts: dict[str, int],
    group_counts: dict[str, int],
    group_missing: dict[str, int],
    group_nonnumeric: dict[str, int],
    exclusion_counts: dict[str, int],
    unexpected: list[str],
) -> str:
    rows = []
    for group in GROUPS:
        expected = FEASIBILITY_COUNTS[group]
        observed = group_counts[group]
        difference = observed - expected
        rows.append(
            f"| {group} | {observed:,} | {expected:,} | {difference:+,} | "
            f"{group_missing[group]} | {group_nonnumeric[group]} |"
        )
    exclusions = "\n".join(
        f"| {reason} | {count:,} |" for reason, count in sorted(exclusion_counts.items())
    ) or "| None | 0 |"
    if unexpected:
        unexpected_text = "\n".join(f"- {item}" for item in unexpected)
    else:
        unexpected_text = "- None detected."
    both = int(join_counts.get("both", 0))
    total = len(joined)
    rate = 100 * both / max(total, 1)
    blank_site = joined["_primary_site"].isna() | joined["_primary_site"].astype(str).eq("")
    fallback_rows = main[main["_primary_site"].isna() | main["_primary_site"].astype(str).eq("")]
    fallback_labels = fallback_rows["anatomical_tissue_label"].value_counts().to_dict()
    fallback_summary = (
        ", ".join(f"{label} ({count})" for label, count in sorted(fallback_labels.items()))
        if fallback_labels
        else "none"
    )
    return f"""# CDH17 targeted data acquisition and QC

Generated by `scripts/01_download_cdh17.py` on {date.today().isoformat()} (UTC date).
This report covers targeted retrieval and data integrity only; no statistical comparisons or figures were produced.

## Retrieval and join

- Expression dataset: `{EXPRESSION_DATASET}`; feature `{GENE_ID}`.
- Phenotype dataset: `{PHENOTYPE_DATASET}`.
- Expression sample IDs advertised by Xena: {len(expression_ids):,}.
- Phenotype sample IDs advertised by Xena: {len(phenotype_ids):,}.
- Shared IDs requested: {len(shared_ids):,}.
- Outer-joined rows: {total:,}; one-to-one merge on `sampleID` succeeded for {both:,} rows ({rate:.4f}%).
- Join statuses: both={join_counts.get('both', 0):,}, left_only={join_counts.get('left_only', 0):,}, right_only={join_counts.get('right_only', 0):,}.
- Duplicate IDs: expression={len(expression_duplicate_ids):,} distinct duplicate IDs; phenotype={len(phenotype_duplicate_ids):,} distinct duplicate IDs.
- Successful joins are not inferred from ordering: the script labels every fetched API value with the corresponding requested `sampleID`, then merges using that key and `validate="one_to_one"`. Expression and phenotype IDs are fetched independently so future unmatched IDs remain visible in the outer-joined raw file.

## Main analysis dataset groups

| Group | Observed rows | Feasibility reference | Difference | Missing expression | Nonnumeric expression |
|---|---:|---:|---:|---:|---:|
{chr(10).join(rows)}

The main file contains {len(main):,} rows and {len(main.columns)} columns. Tumor groups require the corresponding TCGA site, `Primary Tumor`, and exact detailed category (`Colon Adenocarcinoma` or `Rectum Adenocarcinoma`). TCGA normals require the matching site and `Solid Tissue Normal`. GTEx normal colon requires `_study=GTEX` and `_sample_type=Normal Tissue`; `_primary_site=Colon` is used when present, with a detailed source tissue label beginning `Colon` as fallback when the broad site is blank. Other GTEx normals use a non-Colon broad site or, when it is blank, a nonblank detailed tissue label. Detailed source labels remain in the phenotype columns and `anatomical_tissue_label`; only records lacking both usable tissue labels are excluded as unclassified.

## Exclusions from the main dataset

| Mutually exclusive exclusion reason | Rows |
|---|---:|
{exclusions}

The raw intermediate retains rows from excluded categories and any unmatched rows. It does not infer that TCGA normal samples are paired with tumors.

## Missingness and checks

- Missing expression values: {int(joined[EXPRESSION_RAW_COLUMN].isna().sum()):,}.
- Nonnumeric, nonmissing expression values: {int((joined[EXPRESSION_RAW_COLUMN].notna() & pd.to_numeric(joined[EXPRESSION_RAW_COLUMN], errors='coerce').isna()).sum()):,}.
- Missing or blank values in each source phenotype field are counted individually in `results/data_qc_summary.csv`; missing annotations remain blank in the saved raw rows.
- `_primary_site` is blank for {int(blank_site.sum()):,} joined rows. {len(fallback_rows):,} records are retained using detailed tissue labels when available ({fallback_summary}); records without usable tissue labels remain excluded as unclassified.
- No expression values were imputed, transformed, or batch-corrected. The Xena-reported `log2(TPM + 0.001)` values are preserved.
- Duplicate checks, sample-count consistency, join checks, and per-group missingness are recorded in `results/data_qc_summary.csv`.

## Unexpected results and count comparison

{unexpected_text}

Group differences from the earlier feasibility check, if any, are reported above with observed counts and signed differences; the script does not force the reference counts.

## Batch and interpretation limitations

TCGA and GTEx were uniformly reprocessed through Toil, but that does not establish removal of all source, collection, tissue handling, or cohort effects. The source label remains explicit and no batch correction was applied. TCGA and GTEx sample definitions and collection procedures differ. The Toil dataset is an older release (2016-09-03). Expression is already transformed as `log2(TPM + 0.001)`; do not reinterpret it as `log2(TPM + 1)`.
"""


def build_provenance() -> str:
    return f"""# CDH17 data provenance and retrieval

## Source and datasets

- UCSC Xena hub: {HUB_URL}
- Xena data API: {DATA_API}
- Expression dataset ID/name: `{EXPRESSION_DATASET}`
- Matching phenotype dataset ID/name: `{PHENOTYPE_DATASET}`
- Expression file endpoint (reference only; the complete file was not downloaded): {HUB_URL}/download/{EXPRESSION_DATASET}.gz
- Phenotype file endpoint (reference only; annotations were fetched through the API): {HUB_URL}/download/{PHENOTYPE_DATASET}.gz
- Gene feature: `{GENE_ID}`, Ensembl gene ID for CDH17 (version suffix retained as present in the matrix).
- Expression scale: Toil RSEM TPM, documented as `log2(TPM + 0.001)`. Values are kept exactly as returned by Xena; no conversion to `log2(TPM + 1)`, additional normalization, or batch correction is applied.
- Toil release metadata: `2016-09-03` (dataset metadata).
- Access date: {date.today().isoformat()} (UTC date at the most recent script run; the QC report records the same run date).

## Retrieved annotations

The script fetches the matching `sampleID` plus these original phenotype field names from `{PHENOTYPE_DATASET}`: `{', '.join(PHENOTYPE_FIELDS)}`. Categorical API codes are decoded using the same dataset's Xena code tables, and the source fields are retained in the CSV. `analysis_group`, `anatomical_tissue_label`, `exclusion_reason`, and join/QC helper fields are derived fields; they do not replace the source phenotype columns.

## Retrieval and join procedure

1. Query Xena's `/data/` API for sample IDs advertised by the expression and phenotype datasets.
2. Calculate and record the shared sample IDs. Independently request only feature `{GENE_ID}` for every expression sample and the annotation fields listed above for every phenotype sample, in batches of 500 samples, so unmatched IDs are preserved for QC.
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

Xena hub: {HUB_URL}

Xena dataset metadata and downloadable references can be inspected via the hub UI/API. The full expression file is roughly 1.32 GB compressed according to the prior range-header feasibility check; this workflow makes targeted `/data/` API requests instead. The phenotype file's prior range-header size was 135,753 bytes. These are reference sizes, not downloaded files. Verify current endpoint behavior if reproducing at a later date.

Toil's shared processing improves technical consistency but does not guarantee removal of TCGA/GTEx cohort and collection effects. The source fields should be retained for downstream design decisions. No statistical comparisons are made by this acquisition script.
"""


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (XenaError, requests.RequestException, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
