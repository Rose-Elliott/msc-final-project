import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
import seaborn as sns
from skbio.diversity.alpha import shannon

try:
    metadata = pd.read_excel('metadata.xlsx')    #replace these
    arg_data = pd.read_excel('arg_data.xlsx')
except FileNotFoundError:
    print("Data files not found")
    exit()


args_merged = pd.merge(arg_data, metadata[['sample_id', 'group']], on='sample_id', how='left')
args_filtered = args_merged[args_merged['group'].isin(['influent', 'final effluent'])]
pivot = args_filtered.pivot_table(index='sample_id', columns='ARG', values='RNum_Gi', fill_value=0)

argonly_abund = pivot.copy()

def richness(counts):
    return np.sum(counts > 0)

pivot['richness'] = argonly_abund.apply(richness, axis=1)
pivot['shannon'] = argonly_abund.apply(lambda row: shannon(row.values), axis=1)
alphadf = pivot[['richness', 'shannon']].merge(metadata[['sample_id', 'group']], left_index=True, right_on='sample_id')


plt.rcParams.update({
    'font.size': 20,
    'axes.linewidth': 1.0,
    'xtick.major.width': 1.0,
    'ytick.major.width': 1.0,
    'font.family': 'sans-serif'
})

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

palette = {
    'influent': '#1b7837',
    'final effluent': '#762a83'
}

for ax in axes:
    ax.set_facecolor('white')
    ax.grid(axis='y', linestyle='--', linewidth=1, alpha=0.7)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color('black')

#plot a
sns.boxplot(
    data=alphadf, x='group', y='richness', ax=axes[0],
    palette=palette,
    boxprops=dict(edgecolor='black', linewidth=1),
    medianprops=dict(color='black'),
    whiskerprops=dict(color='black'),
    capprops=dict(color='black')
)
axes[0].set_ylabel('Richness', fontsize=20, fontweight='bold', labelpad=15)
axes[0].set_xlabel('')
axes[0].set_xticklabels(['Influent', 'Effluent'], fontsize=18)
axes[0].text(0.01, 0.95, 'a', transform=axes[0].transAxes,
             fontsize=20, fontweight='bold', va='top', ha='left')

#plot b
sns.boxplot(
    data=alphadf, x='group', y='shannon', ax=axes[1],
    palette=palette,
    boxprops=dict(edgecolor='black', linewidth=1),
    medianprops=dict(color='black'),
    whiskerprops=dict(color='black'),
    capprops=dict(color='black')
)
axes[1].set_ylabel('Shannon Diversity Index', fontsize=20, fontweight='bold', labelpad=15)
axes[1].set_xlabel('')
axes[1].set_xticklabels(['Influent', 'Effluent'], fontsize=18)
axes[1].text(0.01, 0.95, 'b', transform=axes[1].transAxes,
             fontsize=20, fontweight='bold', va='top', ha='left')

for metric in ['richness', 'shannon']:
    group1 = alphadf[alphadf['group'] == 'influent'][metric]
    group2 = alphadf[alphadf['group'] == 'final effluent'][metric]
    stat, p = mannwhitneyu(group1, group2, alternative='two-sided')
    print(f"{metric.capitalize()} comparison between in and eff")
    print(f"U statistic = {stat:.2f}, p-value = {p:.4f}")
    print('sig difference' if p < 0.05 else 'no sig difference')
    print()

plt.tight_layout()
plt.show()
