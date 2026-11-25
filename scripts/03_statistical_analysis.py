#!/usr/bin/env python3
"""Run the prespecified CDH17 group comparisons on the saved Xena extract.

Run from the project root with:
    python scripts/03_statistical_analysis.py

The script also supports an absolute path to ``scripts/03_statistical_analysis.py``
when invoked from another working directory.

The three comparisons are unpaired, two-sided Mann–Whitney U tests. Expression
values remain on the reported log2(TPM + 0.001) scale. This script does not
modify input data or perform batch correction.
"""

from __future__ import annotations

from pathlib import Path
import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_CSV = PROJECT_ROOT / "data" / "cdh17_expression_annotated.csv"
RESULTS_DIR = PROJECT_ROOT / "results"
FIGURES_DIR = PROJECT_ROOT / "figures"
SUMMARY_CSV = RESULTS_DIR / "statistical_summary.csv"
NOTES_MD = RESULTS_DIR / "statistical_analysis_notes.md"
FIGURE_PNG = FIGURES_DIR / "figure_D_tcga_statistical_comparisons.png"

EXPRESSION_COLUMN = "CDH17_expression_log2_TPM_plus_0.001"
EXPRESSION_ALIAS_COLUMN = "CDH17_expression_api_value"
EXPRESSION_SCALE = "log2(TPM + 0.001)"


def _require_columns(frame: pd.DataFrame, columns: list[str]) -> None:
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise ValueError(f"Input is missing required columns: {missing}")


def _selection(frame: pd.DataFrame, *, study: str, site: str, sample_type: str,
               detail: str, expected_group: str) -> pd.DataFrame:
    """Select samples only from preserved source phenotype fields."""
    mask = (
        frame["_study"].eq(study)
        & frame["_primary_site"].eq(site)
        & frame["_sample_type"].eq(sample_type)
        & frame["detailed_category"].eq(detail)
    )
    selected = frame.loc[mask].copy()
    if selected.empty:
        raise ValueError(
            "No rows matched original annotations: "
            f"study={study!r}, site={site!r}, sample_type={sample_type!r}, detail={detail!r}"
        )
    labels = selected["analysis_group"].dropna().unique().tolist()
    if labels != [expected_group]:
        raise ValueError(
            f"Original-field selection expected analysis_group={expected_group!r}; "
            f"found {labels!r}"
        )
    return selected


def _median_difference_ci_free_summary(
    comparison_id: str,
    role: str,
    group1: str,
    group2: str,
    first: pd.DataFrame,
    second: pd.DataFrame,
) -> dict[str, object]:
    x = first[EXPRESSION_COLUMN].to_numpy(dtype=float)
    y = second[EXPRESSION_COLUMN].to_numpy(dtype=float)
    # scipy returns U for its first input. U/(n1*n2) is the probability of a
    # first-group value exceeding a second-group value, with half credit for ties.
    result = mannwhitneyu(x, y, alternative="two-sided", method="auto")
    n1, n2 = len(x), len(y)
    rank_biserial = 2.0 * float(result.statistic) / (n1 * n2) - 1.0
    median1 = float(np.median(x))
    median2 = float(np.median(y))
    return {
        "comparison_id": comparison_id,
        "analysis_role": role,
        "group_1": group1,
        "group_2": group2,
        "n_1": n1,
        "n_2": n2,
        "median_1_log2_TPM_plus_0.001": median1,
        "median_2_log2_TPM_plus_0.001": median2,
        "median_difference_1_minus_2": median1 - median2,
        "mann_whitney_U_1": float(result.statistic),
        "rank_biserial_correlation_1_vs_2": rank_biserial,
        "p_value_two_sided": float(result.pvalue),
        "p_value_holm_primary_only": np.nan,
        "expression_scale": EXPRESSION_SCALE,
        "test": "Mann-Whitney U, two-sided, unpaired",
    }


def _holm_adjust_primary(rows: list[dict[str, object]]) -> None:
    """Apply Holm step-down adjustment only to rows marked primary."""
    indices = [i for i, row in enumerate(rows) if row["analysis_role"] == "primary"]
    ordered = sorted(indices, key=lambda i: float(rows[i]["p_value_two_sided"]))
    running_max = 0.0
    total = len(ordered)
    for rank, index in enumerate(ordered):
        raw_p = float(rows[index]["p_value_two_sided"])
        adjusted = min(1.0, (total - rank) * raw_p)
        running_max = max(running_max, adjusted)
        rows[index]["p_value_holm_primary_only"] = running_max


def _format_p(value: float) -> str:
    return f"{value:.3g}" if value >= 0.001 else f"{value:.2e}"


def _make_figure(
    comparisons: list[tuple[str, pd.DataFrame, pd.DataFrame]],
    rows: list[dict[str, object]],
) -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    colors = ["#376996", "#E07A5F"]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 6.6), sharey=True)

    primary_rows = {str(row["comparison_id"]): row for row in rows}
    for ax, (comparison_id, first, second) in zip(axes, comparisons, strict=True):
        groups = [first[EXPRESSION_COLUMN].to_numpy(float), second[EXPRESSION_COLUMN].to_numpy(float)]
        bp = ax.boxplot(
            groups,
            positions=[0, 1],
            widths=0.48,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "#17202A", "linewidth": 2},
            whiskerprops={"color": "#495057", "linewidth": 1.2},
            capprops={"color": "#495057", "linewidth": 1.2},
        )
        for box, color in zip(bp["boxes"], colors, strict=True):
            box.set_facecolor(color)
            box.set_alpha(0.32)
            box.set_edgecolor(color)
            box.set_linewidth(1.6)

        rng = np.random.default_rng(5350 + (0 if comparison_id == "tcga_colon" else 1))
        for position, values, color in zip([0, 1], groups, colors, strict=True):
            jitter = rng.uniform(-0.17, 0.17, size=len(values))
            ax.scatter(
                position + jitter,
                values,
                s=18,
                color=color,
                alpha=0.48,
                linewidths=0,
                zorder=3,
            )

        row = primary_rows[comparison_id]
        p_raw = float(row["p_value_two_sided"])
        p_holm = float(row["p_value_holm_primary_only"])
        site = "Colon" if comparison_id == "tcga_colon" else "Rectum"
        ax.set_title(f"TCGA {site}: primary tumor vs solid-tissue normal", fontsize=12, pad=12)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(
            [f"Primary tumor\n(n={len(first)})", f"Solid-tissue normal\n(n={len(second)})"],
            fontsize=10,
        )
        ax.set_xlabel("")
        ax.grid(axis="y", color="#D9DEE3", linewidth=0.8, alpha=0.8)
        ax.set_axisbelow(True)
        ax.text(
            0.03,
            0.97,
            f"Two-sided Mann–Whitney U\nraw p = {_format_p(p_raw)}\nHolm-adjusted p = {_format_p(p_holm)}",
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=9,
            bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "#CED4DA", "alpha": 0.9},
        )

    axes[0].set_ylabel(f"CDH17 expression ({EXPRESSION_SCALE})", fontsize=11)
    axes[1].set_ylabel("")
    axes[0].set_ylim(-1.25, 10.1)
    fig.suptitle("CDH17 expression in TCGA colorectal tumors and solid-tissue normals", fontsize=14, y=0.99)
    fig.text(
        0.5,
        0.015,
        "Points are individual samples; boxplots show median and interquartile range. Samples are treated as unpaired.",
        ha="center",
        fontsize=9,
        color="#495057",
    )
    fig.tight_layout(rect=(0, 0.045, 1, 0.94))
    fig.savefig(FIGURE_PNG, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _write_notes(summary: pd.DataFrame) -> None:
    lines = [
        "# CDH17 statistical analysis notes",
        "",
        "## Methods",
        "",
        f"The script selects samples from the preserved `_study`, `_primary_site`, `_sample_type`, and `detailed_category` fields, then checks the resulting rows against `analysis_group`. Expression values are used unchanged on the reported `{EXPRESSION_SCALE}` scale. All comparisons are unpaired; no patient matching or pairing is used.",
        "",
        "Each comparison uses a two-sided Mann–Whitney U test. The reported median difference is group 1 minus group 2. Rank-biserial correlation is defined as `2 × U1 / (n1 × n2) − 1`, where U1 is SciPy’s Mann–Whitney statistic for group 1; positive values indicate that group 1 tends to rank higher than group 2, with ties receiving half credit. Holm adjustment is applied across the two prespecified primary TCGA tests only. The GTEx comparison is exploratory and its p-value is unadjusted.",
        "",
        "## Sample selection and results",
        "",
        "| Role | Comparison (group 1 vs group 2) | n1 | n2 | Median 1 | Median 2 | Median difference | Rank-biserial r | Raw two-sided p | Holm-adjusted p |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in summary.to_dict(orient="records"):
        adjusted = row["p_value_holm_primary_only"]
        adjusted_text = "—" if pd.isna(adjusted) else f"{adjusted:.6g}"
        lines.append(
            "| {role} | {g1} vs {g2} | {n1} | {n2} | {m1:.4f} | {m2:.4f} | {diff:.4f} | {r:.4f} | {p:.6g} | {adj} |".format(
                role=row["analysis_role"],
                g1=row["group_1"],
                g2=row["group_2"],
                n1=int(row["n_1"]),
                n2=int(row["n_2"]),
                m1=float(row["median_1_log2_TPM_plus_0.001"]),
                m2=float(row["median_2_log2_TPM_plus_0.001"]),
                diff=float(row["median_difference_1_minus_2"]),
                r=float(row["rank_biserial_correlation_1_vs_2"]),
                p=float(row["p_value_two_sided"]),
                adj=adjusted_text,
            )
        )

    primary = summary[summary["analysis_role"].eq("primary")]
    lines += [
        "",
        "## Interpretation and limitations",
        "",
        "The two TCGA comparisons are the prespecified primary analyses; the GTEx sigmoid-versus-transverse comparison is secondary and exploratory. A nonsignificant p-value should not be interpreted as evidence that groups are equivalent. Median differences describe the observed sample medians on the transformed expression scale; the rank-biserial statistic describes relative rank tendency and is not a fold change.",
        "",
        "The TCGA solid-tissue-normal groups are modest in size, especially rectum normals, so estimates and p-values may be imprecise. Mann–Whitney tests compare distributions/ranks and do not automatically test equality of medians if distribution shapes differ. TCGA and GTEx samples come from different cohorts; the GTEx result is not a substitute for a TCGA matched control and is not adjusted for source effects. No batch correction, imputation, or sample exclusion was applied.",
        "",
        "This is a focused course-project analysis of RNA expression. It does not establish protein abundance, ADC efficacy, or clinical safety.",
        "",
        "## Reproduction",
        "",
        "From the project root run `python scripts/03_statistical_analysis.py`. The script reads the saved annotated CSV and rewrites only `results/statistical_summary.csv`, this notes file, and `figures/figure_D_tcga_statistical_comparisons.png`.",
        "",
    ]
    NOTES_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    if not INPUT_CSV.is_file():
        raise FileNotFoundError(f"Annotated input CSV not found: {INPUT_CSV}")
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    frame = pd.read_csv(INPUT_CSV)
    _require_columns(
        frame,
        [
            "sampleID",
            EXPRESSION_COLUMN,
            EXPRESSION_ALIAS_COLUMN,
            "_study",
            "_primary_site",
            "_sample_type",
            "detailed_category",
            "analysis_group",
        ],
    )
    if frame["sampleID"].isna().any() or frame["sampleID"].duplicated().any():
        raise ValueError("Input sampleID values must be present and unique.")
    expression = pd.to_numeric(frame[EXPRESSION_COLUMN], errors="coerce")
    if expression.isna().any() or not np.isfinite(expression.to_numpy()).all():
        raise ValueError("Expression values contain missing, nonnumeric, or nonfinite entries.")
    alias = pd.to_numeric(frame[EXPRESSION_ALIAS_COLUMN], errors="coerce")
    if alias.isna().any() or not np.array_equal(expression.to_numpy(), alias.to_numpy()):
        raise ValueError("The stored expression value and API-return column do not match exactly.")

    frame[EXPRESSION_COLUMN] = expression
    colon_tumor = _selection(
        frame, study="TCGA", site="Colon", sample_type="Primary Tumor",
        detail="Colon Adenocarcinoma", expected_group="TCGA colon primary tumors",
    )
    colon_normal = _selection(
        frame, study="TCGA", site="Colon", sample_type="Solid Tissue Normal",
        detail="Colon Adenocarcinoma", expected_group="TCGA colon solid-tissue normals",
    )
    rectum_tumor = _selection(
        frame, study="TCGA", site="Rectum", sample_type="Primary Tumor",
        detail="Rectum Adenocarcinoma", expected_group="TCGA rectum primary tumors",
    )
    rectum_normal = _selection(
        frame, study="TCGA", site="Rectum", sample_type="Solid Tissue Normal",
        detail="Rectum Adenocarcinoma", expected_group="TCGA rectum solid-tissue normals",
    )
    gtex_sigmoid = _selection(
        frame, study="GTEX", site="Colon", sample_type="Normal Tissue",
        detail="Colon - Sigmoid", expected_group="GTEx normal colon",
    )
    gtex_transverse = _selection(
        frame, study="GTEX", site="Colon", sample_type="Normal Tissue",
        detail="Colon - Transverse", expected_group="GTEx normal colon",
    )

    selections = [colon_tumor, colon_normal, rectum_tumor, rectum_normal, gtex_sigmoid, gtex_transverse]
    selected_ids = pd.concat([part["sampleID"] for part in selections], ignore_index=True)
    if selected_ids.duplicated().any():
        raise ValueError("One or more samples were selected into multiple comparison groups.")

    rows: list[dict[str, object]] = [
        _median_difference_ci_free_summary(
            "tcga_colon", "primary", "TCGA colon primary tumors", "TCGA colon solid-tissue normals",
            colon_tumor, colon_normal,
        ),
        _median_difference_ci_free_summary(
            "tcga_rectum", "primary", "TCGA rectum primary tumors", "TCGA rectum solid-tissue normals",
            rectum_tumor, rectum_normal,
        ),
        _median_difference_ci_free_summary(
            "gtex_colon_subsite", "secondary_exploratory", "GTEx normal sigmoid colon", "GTEx normal transverse colon",
            gtex_sigmoid, gtex_transverse,
        ),
    ]
    _holm_adjust_primary(rows)
    summary = pd.DataFrame(rows)
    summary.to_csv(SUMMARY_CSV, index=False, float_format="%.10g")

    _make_figure(
        [
            ("tcga_colon", colon_tumor, colon_normal),
            ("tcga_rectum", rectum_tumor, rectum_normal),
        ],
        rows,
    )
    _write_notes(summary)

    print(f"Input: {INPUT_CSV}")
    print(f"Validated rows: {len(frame)}; columns: {len(frame.columns)}")
    print("Verified source-field selections:")
    for group, selected in [
        ("TCGA colon primary tumors", colon_tumor),
        ("TCGA colon solid-tissue normals", colon_normal),
        ("TCGA rectum primary tumors", rectum_tumor),
        ("TCGA rectum solid-tissue normals", rectum_normal),
        ("GTEx sigmoid colon", gtex_sigmoid),
        ("GTEx transverse colon", gtex_transverse),
    ]:
        print(f"  {group}: n={len(selected)}")
    print("\nResults:")
    print(summary.to_string(index=False))
    print(f"\nWrote {SUMMARY_CSV}")
    print(f"Wrote {NOTES_MD}")
    print(f"Wrote {FIGURE_PNG}")


if __name__ == "__main__":
    main()
