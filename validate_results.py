
import pandas as pd
import numpy as np
from pathlib import Path

RESULTS = Path("results")

# Load analysis files
matrix = pd.read_csv(RESULTS / "expression_matrix.csv")
metadata = pd.read_csv(RESULTS / "sample_metadata.csv")
fold_changes = pd.read_csv(RESULTS / "gene_log2fc_comparison.csv")

print("========== DATA VALIDATION ==========")

# 1. Check matrix dimensions
print("\n1. Expression matrix")
print("Genes:", len(matrix))
print("Columns:", len(matrix.columns))
print("Unique gene IDs:", matrix["gene"].nunique())
print("Missing values:", matrix.isna().sum().sum())
print("Duplicate genes:", matrix["gene"].duplicated().sum())

# 2. Check expression values
sample_columns = [c for c in matrix.columns if c != "gene"]
values = matrix[sample_columns].to_numpy()

print("\n2. Expression values")
print("All values finite:", np.isfinite(values).all())
print("Negative values:", (values < 0).sum())
print("Minimum FPKM:", values.min())
print("Maximum FPKM:", values.max())

# 3. Check sample metadata
print("\n3. Sample metadata")
print("Number of samples:", len(metadata))
print("\nSamples per strain and treatment:")
print(pd.crosstab(metadata["strain"], metadata["treatment"]))

print("\nSample names match matrix:",
      set(metadata["sample"]) == set(sample_columns))

# 4. Recalculate fold changes independently
groups = {
    "KHUD_S1": {
        "control": ["KHUD_S1.CON.R1", "KHUD_S1.CON.R2"],
        "treated": ["KHUD_S1.VAN.R1", "KHUD_S1.VAN.R2"]
    },
    "KHUD_S3": {
        "control": ["KHUD_S3.CON.R1", "KHUD_S3.CON.R2"],
        "treated": ["KHUD_S3.VAN.R1", "KHUD_S3.VAN.R2"]
    }
}

for strain, conditions in groups.items():
    control = matrix[conditions["control"]].mean(axis=1)
    treated = matrix[conditions["treated"]].mean(axis=1)
    recalculated = np.log2((treated + 0.1) / (control + 0.1))

    saved = fold_changes[f"{strain}_log2FC"]

    print(f"\n4. {strain} fold-change check")
    print("Maximum absolute difference:",
          np.max(np.abs(recalculated - saved)))
    print("Calculations match:",
          np.allclose(recalculated, saved))

# 5. Check output files
print("\n5. Required output files")
for filename in [
    "expression_matrix.csv",
    "sample_metadata.csv",
    "gene_log2fc_comparison.csv",
    "filtered_genes_2_samples.csv",
    "filtered_genes_4_samples.csv"
]:
    path = RESULTS / filename
    print(filename, "EXISTS" if path.exists() else "MISSING")

print("\n========== VALIDATION COMPLETE ==========")
