# AI-use implementation outline

This document summarizes the project workflow and review points from the existing project files and recorded project progression. It is not a transcript and does not reconstruct exact prompts or chat messages.

## Workflow completed

1. **Environment and source feasibility:** established the project layout and checked that UCSC Xena's Toil TCGA/GTEx expression source and phenotype annotations supported targeted retrieval of CDH17.
2. **Targeted data retrieval and QC:** implemented a script that requests the CDH17 feature and selected annotations, joins by sample ID, preserves raw annotations and excluded sample records, and writes acquisition QC and provenance documentation.
3. **Exploratory analysis:** built and ran the EDA notebook to validate the saved schema and labels, calculate descriptive summaries, review expression floors/outliers and tissue annotations, and create Figures A-C.
4. **Discrepancy review:** independently recalculated group summaries and examined the TCGA versus GTEx normal-colon discrepancy, including colon subsites and detailed GTEx tissue labels.
5. **Prespecified statistical comparisons:** implemented three unpaired two-sided Mann–Whitney comparisons, with Holm adjustment applied only to the two primary TCGA tests, plus notes and Figure D.
6. **Submission preparation:** audited paths, dependencies, data files, and documentation for a proposed public repository. Repository creation or upload is a separate human review decision and has not been performed in this preparation step.

## Human review checkpoints

- The project owner defined the scope and prohibited downloading the full matrix, changing the expression transformation, starting inferential analysis before EDA review, and making clinical efficacy or safety claims from RNA expression.
- Data-source feasibility, acquisition/QC findings, EDA outputs, discrepancy audit, and statistical results were presented for review at their respective stages.
- The project owner specified the three statistical comparisons and their unpaired design, test, effect-size direction, and multiple-testing plan before analysis.
- The proposed repository contents and reproducibility check are to be reviewed before any GitHub repository is created or files are uploaded.

## Documentation still subject to course policy

The local `ai_documentation/` directory was empty before this outline was added, and no instructor-provided AI policy or disclosure template was present in the project directory. This outline records the workflow at a factual level only; it does not claim or invent exact prompts, a prompt transcript, model citations, or checks not recorded in the project. If the course requires those items, add them from the authentic project chat or instructor template before submission.
