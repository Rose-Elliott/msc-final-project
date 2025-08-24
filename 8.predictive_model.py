import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import OneHotEncoder
import matplotlib.pyplot as plt
import matplotlib.patches as patches

try:
    arg_df = pd.read_excel("arg_data.xlsx")    #replace these
    meta_df = pd.read_excel("metadata.xlsx")
except FileNotFoundError:
    print("Data files not found")
    exit()



full_arg = pd.merge(
    arg_df,
    meta_df[['sample_id', 'group', 'Site']],
    on='sample_id',
    how='left'
)

agg_counts = full_arg.groupby(['Site', 'ARG', 'group'])['RNum_Gi'].sum().reset_index()

pivot_table = agg_counts.pivot_table(
    index=['Site', 'ARG'],
    columns='group',
    values='RNum_Gi',
    fill_value=0
).reset_index()

pivot_table = pivot_table.rename(
    columns={'influent': 'influent_abundance', 'final effluent': 'effluent_abundance'}
)


pivot_table['persist'] = np.where(
    (pivot_table['influent_abundance'] > 0) & (pivot_table['effluent_abundance'] > 0), 1, 0
)
pivot_table = pivot_table[pivot_table['influent_abundance'] > 0]  # ignore zero-influent

print("Pivot table preview:")
print(pivot_table.head())

influent_only = full_arg[full_arg['group'] == 'influent']

mge_mode = influent_only.groupby(['Site', 'ARG'])['MGE_prediction'].agg(
    lambda x: x.mode().iloc[0] if not x.mode().empty else 'unknown'
).reset_index()

taxa_mode = influent_only.groupby(['Site', 'ARG'])['Taxa'].agg(
    lambda x: x.mode().iloc[0] if not x.mode().empty else 'unknown'
).reset_index()

#Merge features into main table 
features_df = pd.merge(mge_mode, taxa_mode, on=['Site', 'ARG'], how='left')
model_df = pd.merge(pivot_table, features_df, on=['Site', 'ARG'], how='left')


model_df['MGE_bin'] = model_df['MGE_prediction'].apply(lambda x: 1 if x in ['plasmid', 'phage'] else 0)
model_df['Phylum'] = model_df['Taxa'].apply(
    lambda x: x.split(';')[2] if isinstance(x, str) and len(x.split(';')) > 2 else 'Unknown'
)

#one-hot encode Phylum 
ohe = OneHotEncoder(handle_unknown='ignore', sparse_output=False)
phylum_encoded = ohe.fit_transform(model_df[['Phylum']])
phylum_cols = ohe.get_feature_names_out(['Phylum'])
phylum_df = pd.DataFrame(phylum_encoded, columns=phylum_cols, index=model_df.index)
model_df = pd.concat([model_df, phylum_df], axis=1)

feat_cols = ['influent_abundance', 'MGE_bin'] + list(phylum_cols)
X = model_df[feat_cols].copy()
y = model_df['persist']
X['influent_abundance'] = np.log1p(X['influent_abundance'])

print(f"Feature matrix: {X.shape}, Target: {y.shape}")

#train/test split
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, stratify=y, random_state=42
)
print(f"Training set: {X_train.shape}, Test set: {X_test.shape}")

# Random Forest 
rf_model = RandomForestClassifier(n_estimators=100, random_state=42)
rf_model.fit(X_train, y_train)
y_pred = rf_model.predict(X_test)
y_prob = rf_model.predict_proba(X_test)[:, 1]


feat_importances = pd.Series(rf_model.feature_importances_, index=feat_cols).sort_values(ascending=False)
top_feats = feat_importances.head(15)

def tidy_label(name):
    if name == "influent_abundance":
        return "Influent relative abundance"
    elif name == "MGE_bin":
        return "MGE association"
    elif name.startswith("Phylum_"):
        return name.replace("Phylum_", "")
    else:
        return name

clean_labels = [tidy_label(f) for f in top_feats.index]

#ROC
fpr, tpr, _ = roc_curve(y_test, y_prob)
roc_auc = auc(fpr, tpr)
print(f"ROC AUC: {roc_auc:.3f}")

plt.rcParams.update({
    'font.size': 14,
    'axes.linewidth': 1.0,
    'xtick.major.width': 1.0,
    'ytick.major.width': 1.0,
    'font.family': 'sans-serif'
})

fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.0), gridspec_kw={'wspace': 0.3})

axs[0].barh(clean_labels, top_feats.values, color='#8B0000', edgecolor='black')
axs[0].set_xlabel('Feature Importance (log scale)', fontsize=16, fontweight='bold')
axs[0].set_xscale('log')
axs[0].invert_yaxis()
axs[0].grid(axis='x', linestyle='--', alpha=0.7)
axs[0].set_facecolor('white')
axs[0].tick_params(labelsize=13)
axs[0].text(0.01, 0.98, 'a', transform=axs[0].transAxes, fontsize=18, fontweight='bold', va='top', ha='left')
axs[0].add_patch(patches.Rectangle((0, 0), 1, 1, transform=axs[0].transAxes,
                                    fill=False, color='black', linewidth=1.5))

# ROC curve
axs[1].plot(fpr, tpr, color='darkorange', lw=2.5, label=f'ROC curve (AUC = {roc_auc:.2f})')
axs[1].plot([0, 1], [0, 1], color='grey', lw=1, linestyle='--')
axs[1].set_xlim([0, 1])
axs[1].set_ylim([0, 1.05])
axs[1].set_xlabel('False Positive Rate', fontsize=16, fontweight='bold')
axs[1].set_ylabel('True Positive Rate', fontsize=16, fontweight='bold')
axs[1].legend(loc='lower right', fontsize=13)
axs[1].grid(True, linestyle='--', alpha=0.7)
axs[1].set_facecolor('white')
axs[1].tick_params(labelsize=13)
axs[1].text(0.01, 0.98, 'b', transform=axs[1].transAxes, fontsize=18, fontweight='bold', va='top', ha='left')
axs[1].add_patch(patches.Rectangle((0, 0), 1, 1, transform=axs[1].transAxes,
                                    fill=False, color='black', linewidth=1.5))

plt.tight_layout()
plt.show()

# check 
print("\nTop 5 predicted persistence probabilities:")
print(y_prob[:5])
