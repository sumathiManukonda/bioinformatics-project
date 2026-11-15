
import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# --------------------------------------------------
# Configuration
# --------------------------------------------------

RESULTS_DIR = "results"
FIGURES_DIR = "figures"

EXPRESSION_FILE = os.path.join(
    RESULTS_DIR, "expression_matrix.csv"
)

os.makedirs(FIGURES_DIR, exist_ok=True)

sns.set_theme(style="white", font_scale=0.9)


# --------------------------------------------------
# Load expression matrix
# --------------------------------------------------

print("Loading expression matrix...")

df = pd.read_csv(EXPRESSION_FILE)

gene_column = df.columns[0]
df = df.set_index(gene_column)

# Convert expression values to numeric
df = df.apply(pd.to_numeric, errors="coerce")

if df.isna().any().any():
    raise ValueError("Missing or non-numeric expression values found.")

print(f"Loaded {df.shape[0]} genes and {df.shape[1]} samples.")


# --------------------------------------------------
# Identify samples by strain and treatment
# --------------------------------------------------

def identify_sample(sample_name):
    name = str(sample_name).upper()

    if "KS1" in name or "S1" in name:
        strain = "KHUD_S1"
    elif "KS3" in name or "S3" in name:
        strain = "KHUD_S3"
    else:
        raise ValueError(
            f"Cannot identify strain from sample name: {sample_name}"
        )

    if "VAN" in name:
        treatment = "VAN"
    elif "CON" in name:
        treatment = "CON"
    else:
        raise ValueError(
            f"Cannot identify treatment from sample name: {sample_name}"
        )

    return strain, treatment


sample_info = {}

for sample in df.columns:
    sample_info[sample] = identify_sample(sample)

print("Sample groups identified successfully.")


# --------------------------------------------------
# 1. Hierarchical clustering heatmap
# --------------------------------------------------

print("\nGenerating hierarchical clustering heatmap...")

# Select the 50 most variable genes
top_genes = (
    df.var(axis=1)
    .sort_values(ascending=False)
    .head(50)
    .index
)

selected = df.loc[top_genes]

# Log transformation reduces the influence of very high values
log_expression = np.log2(selected + 1)

# Standardize expression values for each gene
standardized = log_expression.sub(
    log_expression.mean(axis=1), axis=0
)

standardized = standardized.div(
    log_expression.std(axis=1).replace(0, 1), axis=0
)

# Sample annotation colours
strain_colors = {
    "KHUD_S1": "#4C78A8",
    "KHUD_S3": "#F58518"
}

treatment_colors = {
    "CON": "#72B7B2",
    "VAN": "#E45756"
}

column_colors = pd.DataFrame({
    "Strain": [
        strain_colors[sample_info[s][0]] for s in df.columns
    ],
    "Treatment": [
        treatment_colors[sample_info[s][1]] for s in df.columns
    ]
}, index=df.columns)

# Generate clustered heatmap
cluster = sns.clustermap(
    standardized,
    method="average",
    metric="euclidean",
    cmap="vlag",
    center=0,
    row_cluster=True,
    col_cluster=True,
    col_colors=column_colors,
    figsize=(13, 12),
    xticklabels=True,
    yticklabels=True,
    dendrogram_ratio=(0.15, 0.15),
    cbar_kws={
        "label": "Standardized expression"
    }
)

cluster.fig.suptitle(
    "Hierarchical Clustering of Gene Expression",
    y=1.02,
    fontsize=15
)

cluster.savefig(
    os.path.join(FIGURES_DIR, "hierarchical_clustering.png"),
    dpi=300,
    bbox_inches="tight"
)

plt.close(cluster.fig)

print("Saved: figures/hierarchical_clustering.png")


# --------------------------------------------------
# 2. Fold-change comparison scatter plot
# --------------------------------------------------

print("\nCalculating fold changes for both strains...")

fold_changes = pd.DataFrame(index=df.index)

for strain in ["KHUD_S1", "KHUD_S3"]:

    control_samples = [
        sample for sample in df.columns
        if sample_info[sample] == (strain, "CON")
    ]

    treated_samples = [
        sample for sample in df.columns
        if sample_info[sample] == (strain, "VAN")
    ]

    if len(control_samples) != 2 or len(treated_samples) != 2:
        raise ValueError(
            f"{strain} must have two control and two treated samples. "
            f"Found {len(control_samples)} controls and "
            f"{len(treated_samples)} treated samples."
        )

    control_mean = df[control_samples].mean(axis=1)
    treated_mean = df[treated_samples].mean(axis=1)

    fold_changes[f"log2FC_{strain}"] = np.log2(
        (treated_mean + 0.1) / (control_mean + 0.1)
    )


# Classify genes by direction of change
x = fold_changes["log2FC_KHUD_S1"]
y = fold_changes["log2FC_KHUD_S3"]

categories = np.select(
    [
        (x > 0) & (y > 0),
        (x < 0) & (y < 0),
        (x > 0) & (y < 0),
        (x < 0) & (y > 0)
    ],
    [
        "Increased in both strains",
        "Decreased in both strains",
        "Increased in S1, decreased in S3",
        "Decreased in S1, increased in S3"
    ],
    default="No directional change in at least one strain"
)

plot_data = pd.DataFrame({
    "KHUD_S1": x,
    "KHUD_S3": y,
    "Response": categories
})

palette = {
    "Increased in both strains": "#2A9D8F",
    "Decreased in both strains": "#457B9D",
    "Increased in S1, decreased in S3": "#E76F51",
    "Decreased in S1, increased in S3": "#E9C46A",
    "No directional change in at least one strain": "#999999"
}

plt.figure(figsize=(10, 8))

sns.scatterplot(
    data=plot_data,
    x="KHUD_S1",
    y="KHUD_S3",
    hue="Response",
    palette=palette,
    s=35,
    alpha=0.75,
    edgecolor=None
)

# Reference lines
limits = [
    min(x.min(), y.min()),
    max(x.max(), y.max())
]

plt.plot(
    limits,
    limits,
    linestyle="--",
    color="black",
    linewidth=1,
    label="Equal fold change"
)

plt.axhline(0, color="grey", linewidth=0.7)
plt.axvline(0, color="grey", linewidth=0.7)

plt.xlabel("KHUD_S1 log2 fold change")
plt.ylabel("KHUD_S3 log2 fold change")

plt.title(
    "Comparison of Vancomycin-Induced Gene Expression Changes",
    fontsize=14
)

plt.legend(
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    title="Gene response",
    fontsize=8
)

plt.tight_layout()

plt.savefig(
    os.path.join(FIGURES_DIR, "fold_change_comparison.png"),
    dpi=300,
    bbox_inches="tight"
)

plt.close()

# Save calculated values for future reference
fold_changes["Response"] = categories

fold_changes.to_csv(
    os.path.join(RESULTS_DIR, "fold_change_comparison.csv")
)

print("Saved: figures/fold_change_comparison.png")
print("Saved: results/fold_change_comparison.csv")

print("\nGene response summary:")
print(pd.Series(categories).value_counts())

print("\nBoth visualizations generated successfully.")
