
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from pathlib import Path

RESULTS = Path("results")
FIGURES = Path("figures")
FIGURES.mkdir(exist_ok=True)

# Load expression data
df = pd.read_csv(RESULTS / "expression_matrix.csv")
metadata = pd.read_csv(RESULTS / "sample_metadata.csv")

# Prepare expression matrix: genes as rows, samples as columns
expression = df.set_index("gene")

# Log-transform FPKM values
log_expression = np.log2(expression + 0.1)

# PCA: samples are observations, genes are features
sample_data = log_expression.T

scaled_data = StandardScaler().fit_transform(sample_data)

pca = PCA(n_components=2)
coordinates = pca.fit_transform(scaled_data)

pca_df = pd.DataFrame({
    "PC1": coordinates[:, 0],
    "PC2": coordinates[:, 1],
    "sample": sample_data.index
})

pca_df = pca_df.merge(
    metadata[["sample", "strain", "treatment"]],
    on="sample"
)

# PCA plot
plt.figure(figsize=(9, 7))
sns.scatterplot(
    data=pca_df,
    x="PC1",
    y="PC2",
    hue="strain",
    style="treatment",
    s=130
)

label_names = {
    "KHUD_S1.CON.R1": "S1 Control R1",
    "KHUD_S1.CON.R2": "S1 Control R2",
    "KHUD_S1.VAN.R1": "S1 Vancomycin R1",
    "KHUD_S1.VAN.R2": "S1 Vancomycin R2",
    "KHUD_S3.CON.R1": "S3 Control R1",
    "KHUD_S3.CON.R2": "S3 Control R2",
    "KHUD_S3.VAN.R1": "S3 Vancomycin R1",
    "KHUD_S3.VAN.R2": "S3 Vancomycin R2"
}

label_offsets = {
    "KHUD_S1.CON.R1": (-100, 10),
    "KHUD_S1.CON.R2": (10, 10),
    "KHUD_S1.VAN.R1": (-100, -20),
    "KHUD_S1.VAN.R2": (10, 10),
    "KHUD_S3.CON.R1": (10, 10),
    "KHUD_S3.CON.R2": (10, -20),
    "KHUD_S3.VAN.R1": (-110, -35),
    "KHUD_S3.VAN.R2": (10, 20)
}

for _, row in pca_df.iterrows():
    plt.annotate(
        label_names[row["sample"]],
        (row["PC1"], row["PC2"]),
        xytext=label_offsets[row["sample"]],
        textcoords="offset points",
        fontsize=9,
        arrowprops=dict(arrowstyle="-", color="gray", lw=0.5)
    )

plt.title("PCA of Streptococcus anginosus Gene Expression")
plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)")
plt.tight_layout()
plt.savefig(FIGURES / "pca_plot.png", dpi=300)
plt.close()

# Heatmap using the 50 genes with greatest variation
top_genes = log_expression.var(axis=1).nlargest(50).index
heatmap_data = log_expression.loc[top_genes]

# Standardize each gene for visualization
heatmap_data = heatmap_data.sub(
    heatmap_data.mean(axis=1), axis=0
).div(
    heatmap_data.std(axis=1).replace(0, 1), axis=0
)

plt.figure(figsize=(11, 13))
sns.heatmap(
    heatmap_data,
    cmap="vlag",
    center=0,
    yticklabels=False
)

plt.title("Expression Heatmap of 50 Most Variable Genes")
plt.xlabel("Samples")
plt.ylabel("Genes")
plt.tight_layout()
plt.savefig(FIGURES / "expression_heatmap.png", dpi=300)
plt.close()

print("Visualizations created successfully.")
print(f"PCA variance explained: {pca.explained_variance_ratio_.sum()*100:.2f}%")
print("Saved: figures/pca_plot.png")
print("Saved: figures/expression_heatmap.png")
