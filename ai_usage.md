# AI Usage Documentation

## Project

**Title:** Investigating the Effects of Vancomycin on Gene Expression in
Two *Streptococcus anginosus* Strains

**Dataset:** NCBI Gene Expression Omnibus (GEO), accession GSE163717

## 1. AI Tools Used

ChatGPT (OpenAI) was used as a supporting tool during project planning,
implementation, troubleshooting, and documentation.

AI assistance was used for: - Planning and organizing the computational
workflow. - Explaining gene-expression analysis concepts, including
FPKM, log2 fold change, PCA, and hierarchical clustering. - Drafting and
troubleshooting Python scripts for expression-matrix construction,
fold-change calculation, visualization, and validation. - Helping
interpret code output and organize preliminary findings. - Improving
project documentation and the draft report.

AI was used as an aid. It was not treated as an experimental data source
or as a substitute for running and evaluating the analysis.

## 2. Project Components Where AI Assistance Was Used

  -----------------------------------------------------------------------
  Project component                   Description of assistance
  ----------------------------------- -----------------------------------
  Workflow planning                   Helped organize processing of the
                                      eight GEO expression files and
                                      comparison of the two strains.

  Data preprocessing                  Assisted with a Python script to
                                      combine sample-level FPKM files
                                      into an expression matrix and
                                      create sample metadata.

  Fold-change analysis                Assisted with mean-expression
                                      calculations, a 0.1 pseudocount,
                                      log2 fold changes, and expression
                                      filters.

  Visualization                       Assisted with code for PCA, an
                                      expression heatmap, hierarchical
                                      clustering, and a fold-change
                                      comparison scatter plot.

  Validation                          Assisted with checks for missing
                                      values, duplicate gene identifiers,
                                      invalid expression values, metadata
                                      consistency, and recalculated fold
                                      changes.

  Documentation                       Assisted with drafting the README,
                                      this AI-use statement, and sections
                                      of the project report.
  -----------------------------------------------------------------------

## 3. Verification of AI-Generated Outputs

AI-generated code and explanations were reviewed and tested rather than
accepted without verification.

The following checks were performed: - Python scripts were run against
the project dataset. - The expression matrix was checked for its
expected dimensions, missing values, duplicate gene identifiers, and
invalid expression values. - Sample metadata was checked against the
sample columns in the expression matrix. - Fold changes were
independently recalculated by the validation script and compared with
the saved results. The recalculated values agreed with the saved values
within floating-point precision. - The four visualizations were
generated from the analysis outputs and reviewed for inclusion in the
report. - Numerical statements in the report were checked against
generated result files and summary outputs. - Scientific background and
dataset details were associated with the original publication and NCBI
GEO record, rather than relying on AI-generated claims alone.

AI-generated suggestions that were not supported by the data were not
presented as findings.

## 4. Limitations and Human Review

The analysis uses processed FPKM expression files from GEO, not raw
sequencing reads. Therefore, raw-read quality control and read alignment
were not independently repeated in this project.

The current comparison is descriptive. It uses two biological replicates
per strain and treatment condition and does not include a formal
statistical differential-expression test or adjusted p-values.
Accordingly, fold-change directions and visualization patterns are not
described as statistically significant, and they do not by themselves
establish a mechanism of vancomycin tolerance or resistance.

The project team remains responsible for reviewing the code, checking
the results, interpreting the biological findings, and approving the
final report. AI assistance does not replace the team's responsibility
for the accuracy and integrity of the work.

## 5. Ethical and Data-Privacy Considerations

The project analyzes publicly available processed bacterial
gene-expression data from NCBI GEO. It does not involve recruiting human
participants, collecting new human samples, or processing identifiable
personal information.

The dataset source and original research publication are acknowledged in
the report. The analysis is intended for academic purposes, and
limitations are stated to avoid overstating the conclusions.

No confidential participant information or private personal data was
intentionally provided to the AI tool as part of this workflow.

## 6. Reproducibility

The project repository should include: - Python scripts for
preprocessing, fold-change analysis, visualization, and validation. -
`README.md`, describing the project, dataset, repository structure,
software requirements, and instructions for running the workflow. - This
`ai_usage.md` file, documenting AI assistance and verification. -
Relevant output tables and figures, where appropriate and permitted by
the dataset's terms.

The software versions recorded for the analysis are: - Python
3.10.0rc2 - pandas 2.0.3 - NumPy 1.26.4 - Matplotlib 3.8.0 - Seaborn
0.13.0 - scikit-learn 1.3.0

## 7. Declaration

AI tools were used as supporting resources for coding assistance,
troubleshooting, explanation, and documentation. The project team
reviewed and tested the generated code and checked the reported results
against the actual analysis outputs. The final interpretation and
submission remain the responsibility of the project team.
