import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from skbio.stats.distance import DistanceMatrix
from skbio.stats.ordination import pcoa
from skbio.stats.distance import permanova
from scipy.spatial.distance import pdist, squareform

try:
    arg_data = pd.read_excel("arg_data.xlsx")
    metadata = pd.read_excel("metadata.xlsx")
except FileNotFoundError:
    print("Data files not found")
    exit()



def draw_conf_ellipse(x, y, ax, n_std=2, facecolor='none', **kwargs):
    if len(x) < 2:   #not enough points
        return
    cov = np.cov(x, y)
    mean_x, mean_y = np.mean(x), np.mean(y)
    eigvals, eigvecs = np.linalg.eigh(cov)

    #SORT eigenvalues
    order = eigvals.argsort()[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]

    angle = np.degrees(np.arctan2(eigvecs[1, 0], eigvecs[0, 0]))
    width, height = 2 * n_std * np.sqrt(eigvals)

    ellipse = Ellipse((mean_x, mean_y), width, height, angle=angle,
                      facecolor=facecolor, **kwargs)
    ax.add_patch(ellipse)

arg_matrix = arg_data.pivot_table(
    index="sample_id",
    columns="ARG",
    values="RNum_Gi",
    aggfunc="sum",
    fill_value=0
)

meta_matrix = metadata.set_index("sample_id").loc[arg_matrix.index]

brayArr = pdist(arg_matrix.values, metric="braycurtis")
brayMat = DistanceMatrix(squareform(brayArr), ids=arg_matrix.index)
xx = brayMat 

pcoa_out = pcoa(brayMat)

#take first 2 coords
coords = pcoa_out.samples.iloc[:, :2].copy()
coords.columns = ["PCoA1", "PCoA2"]
coords["group"] = meta_matrix["group"].values
coords["Site"] = meta_matrix["Site"].values

#permanova
perm = permanova(distance_matrix=brayMat, grouping=coords["group"], permutations=999)
varExpl = pcoa_out.proportion_explained * 100

merged = arg_data.merge(metadata[["sample_id", "Site", "group"]], on="sample_id")

pivot_site = merged.pivot_table(
    index="sample_id",
    columns="ARG",
    values="RNum_Gi",
    aggfunc="sum",
    fill_value=0
)
pivot_site = pivot_site.merge(metadata[["sample_id", "Site", "group"]],
                              left_index=True, right_on="sample_id")

site_means = pivot_site.groupby(["Site", "group"]).mean(numeric_only=True).reset_index()
site_means["sample_id"] = site_means["Site"] + " | " + site_means["group"]

#bray again
site_mat = site_means.drop(columns=["Site", "group"]).set_index("sample_id")
braySite = pdist(site_mat.values, metric="braycurtis")
braySiteMat = DistanceMatrix(squareform(braySite), ids=site_mat.index)

#pcoa again
pcoa_site = pcoa(braySiteMat)
coordsSite = pcoa_site.samples.iloc[:, :2].copy()
coordsSite.columns = ["PCoA1", "PCoA2"]
coordsSite = coordsSite.reset_index().rename(columns={"index": "sample_id"})
coordsSite[["Site", "group"]] = coordsSite["sample_id"].str.split(" \| ", expand=True)
varExpl_site = pcoa_site.proportion_explained * 100


#plot
plt.rcParams.update({
    "font.size": 12,
    "axes.linewidth": 1.2,
    "xtick.major.width": 1.2,
    "ytick.major.width": 1.2
})

colors = {"influent": "#1b7837", "final effluent": "#762a83"}

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# plot 1
ax = axes[0]
ax.set_facecolor("white")
ax.grid(True, color="lightgrey", linestyle="--", alpha=0.7)
ax.set_axisbelow(True)

for g, d in coords.groupby("group"):
    c = colors.get(g, "#777777")
    ax.scatter(d["PCoA1"], d["PCoA2"], color=c, edgecolor="black", s=90,
               alpha=0.9, label=g.capitalize(), zorder=3)
    draw_conf_ellipse(d["PCoA1"], d["PCoA2"], ax, n_std=2,
                      edgecolor=c, linewidth=2, facecolor=c, alpha=0.2, zorder=2)

ax.set_xlabel("PCoA1 (" + str(round(varExpl[0]*1.0,1)) + "%)", fontweight="bold")
ax.set_ylabel(f"PCoA2 ({varExpl[1]:.1f}%)", fontweight="bold")
ax.legend(title="Sample Type", fontsize=11, title_fontsize=12,
          frameon=True, loc="lower right")
ax.text(0.02, 0.98,
        f"PERMANOVA\nPseudo-F = {perm['test statistic']:.2f}\np = {perm['p-value']:.3f}",
        transform=ax.transAxes, verticalalignment="top", fontsize=11,
        bbox=dict(boxstyle="round,pad=0.5", facecolor="white", alpha=0.9,
                  edgecolor="gray", linewidth=1))
ax.text(-0.1, 1.05, "a", transform=ax.transAxes, fontsize=16, fontweight="bold")

# plot 2
ax = axes[1]
ax.set_facecolor("white")
ax.grid(True, color="lightgrey", linestyle="--", alpha=0.7)
ax.set_axisbelow(True)

for g, d in coordsSite.groupby("group"):
    c = colors.get(g, "#777777")
    ax.scatter(d["PCoA1"], d["PCoA2"], color=c, edgecolor="black", s=90,
               alpha=0.9, label=g.capitalize(), zorder=3)
    draw_conf_ellipse(d["PCoA1"], d["PCoA2"], ax, n_std=2,
                      edgecolor=c, linewidth=2, facecolor=c, alpha=0.2, zorder=2)

ax.set_xlabel(f"PCoA1 ({varExpl_site[0]:.1f}%)", fontweight="bold")
ax.set_ylabel("PCoA2 (" + str(round(varExpl_site[1]*100/100, 1)) + "%)", fontweight="bold")
ax.text(-0.1, 1.05, "b", transform=ax.transAxes, fontsize=16, fontweight="bold")

plt.tight_layout()
plt.show()

print("PERMANOVA results:")
print("Pseudo-F = " + str(round(perm["test statistic"], 2)))
print(f"p-value = {perm['p-value']:.3f}")
