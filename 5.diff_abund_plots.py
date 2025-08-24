import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests
from matplotlib.patches import Patch

try:
    arg_data = pd.read_excel("arg_data.xlsx")   #replace these
    metadata = pd.read_excel("metadata.xlsx")
except FileNotFoundError:
    print("Data files not found")
    exit()



arg_matrix = arg_data.pivot_table(
    index="sample_id",
    columns="ARG",
    values="RNum_Gi",
    aggfunc="sum",
    fill_value=0
)


arg_matrix = arg_matrix.merge(
    metadata[["sample_id", "group"]],
    left_index=True,
    right_on="sample_id"
).set_index("sample_id")

#wilcoxon rank-sum test
stats_results = []

for gene in arg_matrix.columns.drop("group"):
    infl_vals = arg_matrix[arg_matrix["group"] == "influent"][gene]
    effl_vals = arg_matrix[arg_matrix["group"] == "final effluent"][gene]

    #skipping genes with 0 counts
    if infl_vals.sum() == 0 and effl_vals.sum() == 0:
        continue

    stat, pval = ranksums(infl_vals, effl_vals)
    log2_fc = np.log2((effl_vals.mean() + 1e-6) / (infl_vals.mean() + 1e-6))

    stats_results.append({
        "ARG": gene,
        "pval": pval,
        "log2FC": log2_fc
    })

results_df = pd.DataFrame(stats_results)

#FDR
results_df["FDR"] = multipletests(results_df["pval"], method="fdr_bh")[1]
results_df["significant"] = results_df["FDR"] < 0.05
results_df["direction"] = np.where(results_df["log2FC"] > 0, "gained", "lost")



colors_map = {"gained": "#B22222", "lost": "#1F4E79", "ns": "grey"}

def pick_color(row):
    if row["significant"] and row["direction"] == "gained":
        return colors_map["gained"]
    elif row["significant"] and row["direction"] == "lost":
        return colors_map["lost"]
    else:
        return colors_map["ns"]

results_df["color"] = results_df.apply(pick_color, axis=1)


highlight_df = pd.read_csv("129_differential_ARGs.csv")
highlight_df["log2FC"] = np.log2(
    (highlight_df["mean_effluent"] + 1e-6) / (highlight_df["mean_influent"] + 1e-6)
)

def sig_to_stars(fdr):
    if fdr < 0.001: return "***"
    elif fdr < 0.01: return "**"
    elif fdr < 0.05: return "*"
    else: return ""

highlight_df["sig_star"] = highlight_df["FDR"].apply(sig_to_stars)

#choose top gained & lost
top_gained = highlight_df[highlight_df["direction"] == "gained"].nlargest(20, "log2FC")
top_lost = highlight_df[highlight_df["direction"] == "lost"].nsmallest(20, "log2FC")
plot_subset = pd.concat([top_gained, top_lost])

name_map = {
    "ESCHERICHIA COLI EF-TU MUTANTS CONFERRING RESISTANCE TO PULVOMYCIN": "Ecol_EFTu_PLV",
    "VANS": "vanS", "ERMB": "ermB", "MEFA": "mefA", "VANU": "vanU", "MEL": "mel",
    "EFRA": "efrA", "EFRB": "efrB", "TOLC": "tolC", "BCRA": "bcrA", "TET34": "tet34",
    "MACB": "macB", "TRUNCATED_PUTATIVE_RESPONSE_REGULATOR_ARLR": "arlR", "MPHE": "mphE",
    "LSAE": "lsaE", "LNUB": "lnuB", "DNA-BINDING_PROTEIN_H-NS": "H-NS",
    "LLMA_23S_RIBOSOMAL_RNA_METHYLTRANSFERASE": "LlmA_23S_CLI",
    "ACRB": "acrB", "OXA-333": "OXA-333", "BRUCELLA SUIS MPRF": "mprF",
    "STAPHYLOCOCCUS MUPA CONFERRING RESISTANCE TO MUPIROCIN": "mupA",
    "CCRA": "ccrA", "TRIC": "triC"
}
plot_subset["ARG"] = plot_subset["ARG"].map(name_map).fillna(plot_subset["ARG"])
plot_subset = plot_subset.sort_values("log2FC")


plt.rcParams.update({
    "font.size": 10,
    "axes.linewidth": 1.2,
    "xtick.major.width": 1.2,
    "ytick.major.width": 1.2,
    "font.family": "sans-serif"
})

palette = {"gained": "#B22222", "lost": "#1F4E79"}

fig, (ax_volcano, ax_bar) = plt.subplots(
    1, 2,
    figsize=(14, 7),
    gridspec_kw={"width_ratios": [1.1, 1], "wspace": 0.3}
)

for ax in [ax_volcano, ax_bar]:
    ax.set_facecolor("white")
    ax.grid(axis="y", linestyle="--", linewidth=1, alpha=0.7)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.2)

# Volcano
ax_volcano.scatter(
    results_df["log2FC"],
    -np.log10(results_df["FDR"]),
    c=results_df["color"],
    edgecolor="black",
    linewidth=0.7,
    s=60,
    alpha=1,
    zorder=3
)
ax_volcano.axhline(-np.log10(0.05), color="black", linestyle="--", linewidth=1)
ax_volcano.axvline(0, color="black", linestyle="-", linewidth=0.8)
ax_volcano.set_xlabel("Log₂ Fold Change", fontsize=10, fontweight="bold", labelpad=8)
ax_volcano.set_ylabel("-Log₁₀ Adjusted p-value (FDR)", fontsize=10, fontweight="bold", labelpad=8)
ax_volcano.tick_params(axis="both", which="major", labelsize=9)
ax_volcano.text(-0.07, 1.02, "a", transform=ax_volcano.transAxes,
    fontsize=12, fontweight="bold", va="top", ha="left")

top_gained_sig = results_df[(results_df["significant"]) & (results_df["direction"] == "gained")]
top_gained_sig = top_gained_sig.nlargest(3, "log2FC")
top_gained_sig["short_name"] = top_gained_sig["ARG"].map(name_map).fillna(top_gained_sig["ARG"])

for _, row in top_gained_sig.iterrows():
    ax_volcano.text(
        row["log2FC"],
        -np.log10(row["FDR"]),
        row["short_name"],
        fontsize=8,
        ha="left",
        va="bottom",
        rotation=25,
        color="black"
    )

#Bars
bars = ax_bar.barh(
    y=plot_subset["ARG"],
    width=plot_subset["log2FC"],
    color=plot_subset["direction"].map(palette),
    edgecolor="black",
    linewidth=1.0
)

for bar, (_, row) in zip(bars, plot_subset.iterrows()):
    if row.sig_star == "":
        continue
    x_pos = bar.get_width() + (0.1 if bar.get_width() > 0 else -0.1)
    ha = "left" if bar.get_width() > 0 else "right"
    ax_bar.text(
        x_pos,
        bar.get_y() + bar.get_height() / 2,
        row.sig_star,
        color="black",
        va="center",
        ha=ha,
        fontsize=9,
        fontweight="bold"
    )

ax_bar.axvline(0, color="black", linestyle="--", linewidth=1)
ax_bar.set_xlabel("Log₂ Fold Change", fontsize=10, fontweight="bold", labelpad=8)
ax_bar.set_ylabel("")
ax_bar.tick_params(axis="both", which="major", labelsize=9)
ax_bar.text(-0.07, 1.02, "b", transform=ax_bar.transAxes,
    fontsize=12, fontweight="bold", va="top", ha="left")

legend_items = [
    Patch(facecolor=palette["gained"], edgecolor="black", label="gained"),
    Patch(facecolor=palette["lost"], edgecolor="black", label="lost"),
    Patch(facecolor=colors_map["ns"], edgecolor="black", label="not significant")
]
fig.legend(
    handles=legend_items,
    loc="lower center",
    ncol=3,
    frameon=False,
    fontsize=9,
    bbox_to_anchor=(0.5, -0.03)
)

plt.tight_layout(rect=[0, 0.05, 1, 1])
plt.show()
