import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

try:
    arg_data = pd.read_excel("arg_data.xlsx")   #replace these
    metadata = pd.read_excel("metadata.xlsx")
except FileNotFoundError:
    print("Data files not found")
    exit()


agg_table = arg_data.groupby(["sample_id", "AMR_category"], as_index=False)["RNum_Gi"].sum()

wide_table = agg_table.pivot(index="sample_id", columns="AMR_category", values="RNum_Gi").fillna(0)

tidy_meta = metadata.set_index("sample_id").loc[wide_table.index]
wide_table = wide_table.merge(tidy_meta[["Site", "group"]], left_index=True, right_index=True)

influent_tbl = wide_table[wide_table["group"] == "influent"]
effluent_tbl = wide_table[wide_table["group"] == "final effluent"]

influent_sum = influent_tbl.groupby("Site").sum().drop(columns=["group"], errors="ignore")
effluent_sum = effluent_tbl.groupby("Site").sum().drop(columns=["group"], errors="ignore")

#15 top classes
combined_total = influent_sum.add(effluent_sum, fill_value=0)
top_classes = combined_total.sum(axis=0).sort_values(ascending=False).head(15).index.tolist()

def squash_others(df):
    df = df.copy()
    df["Other"] = df.drop(columns=top_classes, errors="ignore").sum(axis=1)
    return df[top_classes + ["Other"]]

influent_sum = squash_others(influent_sum)
effluent_sum = squash_others(effluent_sum)

site_order = sorted(set(influent_sum.index) | set(effluent_sum.index))
influent_sum = influent_sum.reindex(site_order).fillna(0)
effluent_sum = effluent_sum.reindex(site_order).fillna(0)

#for percent reduction
influent_totals = influent_sum.sum(axis=1)
effluent_totals = effluent_sum.sum(axis=1)

drop_pct = ((influent_totals - effluent_totals) / influent_totals) * 100

print("\nRelative Burden Reduction by Site:")
for s, val in drop_pct.items():
    print(s + ": " + str(round(val, 2)) + "% reduction")

wide_rel = agg_table.pivot(index="sample_id", columns="AMR_category", values="RNum_Gi").fillna(0)
meta_for_rel = metadata.set_index("sample_id").loc[wide_rel.index]

rel_abund = wide_rel.div(wide_rel.sum(axis=1), axis=0) * 100
rel_abund["Other"] = rel_abund.drop(columns=top_classes, errors="ignore").sum(axis=1)

grand_total = rel_abund.sum().sum()
top15_total = rel_abund[top_classes].sum().sum()
prop_top15 = (top15_total / grand_total) * 100
print(f"\ntop 15 classes account for {prop_top15:.2f}% of burden\n")

rel_abund_top = rel_abund[top_classes + ["Other"]].merge(
    meta_for_rel[["Site", "group"]], left_index=True, right_index=True
)

influent_rel = rel_abund_top[rel_abund_top["group"] == "influent"]
effluent_rel = rel_abund_top[rel_abund_top["group"] == "final effluent"]

#mean by site (relative)
influent_mean = influent_rel.groupby("Site")[top_classes + ["Other"]].mean()
effluent_mean = effluent_rel.groupby("Site")[top_classes + ["Other"]].mean()

print("Influent mean (rounded):")
print(influent_mean.round(2))

print("\nEffluent mean (rounded):")
print(np.round(effluent_mean, 2))

#plot
plt.rcParams.update({
    "font.size": 14,
    "axes.linewidth": 1.2,
    "xtick.major.width": 1.2,
    "ytick.major.width": 1.2,
    "font.family": "sans-serif"
})

bar_colors = sns.color_palette("tab20", n_colors=20)
grey_color = [(0.7, 0.7, 0.7)]
palette = bar_colors[:15] + grey_color

fig, axes = plt.subplots(2, 2, figsize=(20, 12), sharey="row")
panel_ids = ["a", "b", "c", "d"]

#plots a and b
for idx, (ax, df) in enumerate(zip(axes[0], [influent_sum, effluent_sum])):
    base = np.zeros(len(df))
    for j, cat in enumerate(df.columns):
        ax.bar(df.index, df[cat], bottom=base, color=palette[j],
               label=cat, width=0.7, edgecolor="white", linewidth=1)
        base += df[cat].values
    ax.set_facecolor("white")
    ax.grid(axis="y", linestyle="--", linewidth=1, alpha=0.7)
    ax.set_axisbelow(True)
    for sp in ax.spines.values():
        sp.set_visible(True)
        sp.set_color("black")
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    ax.text(0.01, 0.95, panel_ids[idx], transform=ax.transAxes,
            fontsize=18, fontweight="bold", va="top", ha="left")

axes[0][0].set_ylabel("Total Relative Burden", fontsize=16, fontweight="bold", labelpad=15)

# plots c and d
for idx, (ax, df) in enumerate(zip(axes[1], [influent_mean, effluent_mean]), start=2):
    base = np.zeros(len(df))
    for j, cat in enumerate(df.columns):
        ax.bar(df.index, df[cat], bottom=base, color=palette[j],
               label=cat, width=0.7, edgecolor="white", linewidth=1)
        base += df[cat].values
    ax.set_facecolor("white")
    ax.grid(axis="y", linestyle="--", linewidth=1, alpha=0.7)
    ax.set_axisbelow(True)
    for sp in ax.spines.values():
        sp.set_visible(True)
        sp.set_color("black")
    ax.set_xticklabels(df.index, rotation=45, ha="right", fontsize=14)
    ax.set_xlabel("Site", fontsize=16, fontweight="bold", labelpad=15)
    ax.text(0.01, 0.95, panel_ids[idx], transform=ax.transAxes,
            fontsize=18, fontweight="bold", va="top", ha="left")

axes[1][0].set_ylabel("Mean Relative Abundance (%)", fontsize=16, fontweight="bold", labelpad=15)

handles, labels = axes[1][0].get_legend_handles_labels()
axes[0][1].legend(handles, labels, loc="upper left",
                  bbox_to_anchor=(1.05, 1), fontsize=14, frameon=False)

plt.tight_layout(rect=[0, 0, 0.92, 1])
plt.show()
