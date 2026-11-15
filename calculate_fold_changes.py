
import pandas as pd
import numpy as np
from pathlib import Path

# Directories
RESULTS = Path("results")
RESULTS.mkdir(exist_ok=True)

# Load expression matrix
df = pd.read_csv(RESULTS / "expression_matrix.csv")

# Sample groups
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

# Pseudocount from the project proposal
PSEUDOCOUNT = 0.1

# Calculate mean expression and log2 fold changes
results = pd.DataFrame()
results["gene"] = df["gene"]

for strain, conditions in groups.items():

    control_mean = df[conditions["control"]].mean(axis=1)
    treated_mean = df[conditions["treated"]].mean(axis=1)

    log2fc = np.log2(
        (treated_mean + PSEUDOCOUNT) /
        (control_mean + PSEUDOCOUNT)
    )

    results[f"{strain}_control_mean"] = control_mean
    results[f"{strain}_treated_mean"] = treated_mean
    results[f"{strain}_log2FC"] = log2fc

# Classify response direction in each strain
s1 = results["KHUD_S1_log2FC"]
s3 = results["KHUD_S3_log2FC"]

results["response_pattern"] = np.select(
    [
        (s1 > 0) & (s3 > 0),
        (s1 < 0) & (s3 < 0),
        (s1 > 0) & (s3 < 0),
        (s1 < 0) & (s3 > 0)
    ],
    [
        "Increased in both strains",
        "Decreased in both strains",
        "Increased in S1, decreased in S3",
        "Decreased in S1, increased in S3"
    ],
    default="No directional change in at least one strain"
)

# Save complete comparison
results.to_csv(
    RESULTS / "gene_log2fc_comparison.csv",
    index=False
)

# Apply the two expression filters
all_sample_columns = [
    col for col in df.columns if col != "gene"
]

expressed_count = (df[all_sample_columns] >= 1).sum(axis=1)

for threshold in [2, 4]:

    keep = expressed_count >= threshold

    filtered = results.loc[keep].copy()

    filtered.to_csv(
        RESULTS / f"filtered_genes_{threshold}_samples.csv",
        index=False
    )

    print(
        f"Genes retained by filter "
        f"(FPKM >= 1 in at least {threshold} samples): "
        f"{keep.sum()}"
    )

# Summarise response patterns
print("\n========== RESPONSE SUMMARY ==========")

print(results["response_pattern"].value_counts())

# Display strongest changes
print("\nTop 10 increases in KHUD_S1:")
print(
    results.nlargest(10, "KHUD_S1_log2FC")[
        ["gene", "KHUD_S1_log2FC", "KHUD_S3_log2FC"]
    ].to_string(index=False)
)

print("\nTop 10 decreases in KHUD_S1:")
print(
    results.nsmallest(10, "KHUD_S1_log2FC")[
        ["gene", "KHUD_S1_log2FC", "KHUD_S3_log2FC"]
    ].to_string(index=False)
)

print("\nTop 10 increases in KHUD_S3:")
print(
    results.nlargest(10, "KHUD_S3_log2FC")[
        ["gene", "KHUD_S1_log2FC", "KHUD_S3_log2FC"]
    ].to_string(index=False)
)

print("\nTop 10 decreases in KHUD_S3:")
print(
    results.nsmallest(10, "KHUD_S3_log2FC")[
        ["gene", "KHUD_S1_log2FC", "KHUD_S3_log2FC"]
    ].to_string(index=False)
)

print("\nFull gene comparison saved successfully.")
