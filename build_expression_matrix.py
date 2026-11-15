
import pandas as pd
from pathlib import Path

# Project directories
DATA_DIR = Path("data")
RESULTS_DIR = Path("results")
RESULTS_DIR.mkdir(exist_ok=True)

# Find all eight compressed expression files
files = sorted(DATA_DIR.glob("*.txt.gz"))

if len(files) != 8:
    raise ValueError(f"Expected 8 files, found {len(files)}")

expression_matrix = None
reference_genes = None
sample_records = []

print("Reading expression files...\n")

for file in files:

    # Read compressed tab-separated file
    df = pd.read_csv(
        file,
        sep="\t",
        compression="gzip",
        skip_blank_lines=True
    )

    # Validate structure
    if df.shape[1] != 2 or df.columns[0] != "gene":
        raise ValueError(f"Unexpected file structure: {file.name}")

    sample_name = df.columns[1]

    # Check for duplicate gene identifiers
    if df["gene"].duplicated().any():
        raise ValueError(f"Duplicate gene IDs in {file.name}")

    # Check gene identifiers against first sample
    current_genes = set(df["gene"])

    if reference_genes is None:
        reference_genes = current_genes
    elif current_genes != reference_genes:
        raise ValueError(f"Gene identifiers differ in {file.name}")

    # Rename expression column using sample name
    df = df.rename(columns={sample_name: sample_name})

    # Merge samples using gene identifiers
    if expression_matrix is None:
        expression_matrix = df
    else:
        expression_matrix = expression_matrix.merge(
            df,
            on="gene",
            how="inner",
            validate="one_to_one"
        )

    # Extract sample metadata
    parts = sample_name.split(".")

    sample_records.append({
        "sample": sample_name,
        "strain": parts[0],
        "treatment": (
            "Control" if parts[1] == "CON" else "Vancomycin"
        ),
        "replicate": parts[2],
        "filename": file.name
    })

    print(f"Processed: {file.name}")

# Save combined expression matrix
matrix_path = RESULTS_DIR / "expression_matrix.csv"
expression_matrix.to_csv(matrix_path, index=False)

# Save sample metadata
metadata = pd.DataFrame(sample_records)
metadata_path = RESULTS_DIR / "sample_metadata.csv"
metadata.to_csv(metadata_path, index=False)

# Final validation
print("\n========== ANALYSIS SUMMARY ==========")
print(f"Samples processed: {len(files)}")
print(f"Genes matched: {len(expression_matrix)}")
print(f"Expression columns: {expression_matrix.shape[1] - 1}")
print(f"Missing values: {expression_matrix.isna().sum().sum()}")

print(f"\nExpression matrix saved: {matrix_path}")
print(f"Sample metadata saved: {metadata_path}")

print("\nExpression matrix preview:")
print(expression_matrix.head())

print("\nSample metadata:")
print(metadata)
