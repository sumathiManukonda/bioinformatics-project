# Vancomycin Effects on Gene Expression in *Streptococcus anginosus*

## Project overview

This project investigates changes in gene expression following vancomycin exposure in two *Streptococcus anginosus* strains: KHUD_S1 and KHUD_S3.

The analysis uses processed FPKM expression data from NCBI GEO accession GSE163717.

## Dataset

Eight samples were analyzed:

- KHUD_S1: two untreated controls and two vancomycin-treated samples.
- KHUD_S3: two untreated controls and two vancomycin-treated samples.

Vancomycin treatment: 8 µg/mL for one hour.

Each sample contains gene identifiers and FPKM expression values.

## Project structure

- `data/`: Original compressed expression files.
- `results/`: Expression matrices, fold changes, filtered gene lists, validation information, and software versions.
- `figures/`: PCA plot, expression heatmap, and additional visualizations.
- `notebooks/`: Reserved for analysis notebooks.
- `build_expression_matrix.py`: Combines the eight expression files.
- `calculate_fold_changes.py`: Calculates descriptive log2 fold changes.
- `create_visualizations.py`: Generates PCA and heatmap figures.
- `create_additional_visualizations.py`: Generates additional figures.
- `validate_results.py`: Checks data integrity and fold-change calculations.

## Analysis workflow

1. Read the eight compressed expression files.
2. Match genes across all samples and construct an expression matrix.
3. Calculate mean control and treated FPKM for each strain.
4. Calculate log2 fold change using a pseudocount of 0.1.
5. Compare expression-change directions between strains.
6. Apply expression filters:
   - FPKM >= 1 in at least two samples.
   - FPKM >= 1 in at least four samples.
7. Generate PCA and heatmap visualizations.
8. Generate additional visualizations.
9. Validate the expression matrix and fold-change calculations.

## Main results

- 1,757 genes included in the expression matrix.
- 8 samples analyzed.
- 1,746 genes passed the two-sample expression filter.
- 1,729 genes passed the four-sample expression filter.
- The first two PCA components explained 80.21% of the variation.

## Running the analysis

Run these commands from the project root:

```bash
python3 build_expression_matrix.py
python3 calculate_fold_changes.py
python3 create_visualizations.py
python3 create_additional_visualizations.py
python3 validate_results.py
```
