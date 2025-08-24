import pandas as pd
import numpy as np
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests

try:
    arg_data = pd.read_excel("arg_data.xlsx") #replace these
    metadata = pd.read_excel("metadata.xlsx")
except FileNotFoundError:
    print("...")
    exit()

pivot = arg_data.pivot_table(
    index='sample_id', 
    columns='ARG', 
    values='RNum_Gi', 
    aggfunc='sum', 
    fill_value=0
)


pivot = pivot.merge(
    metadata[['sample_id', 'group']], 
    left_index=True, right_on='sample_id'
).set_index('sample_id')

#rank sum tests
results_list = []  

for arg_name in pivot.columns[:-1]:  
    # grab values for each group
    influent_vals = pivot[pivot['group'] == 'influent'][arg_name]
    effluent_vals = pivot[pivot['group'] == 'final effluent'][arg_name]

    stat, p_value = ranksums(influent_vals, effluent_vals)

    # mean differences
    mean_influent = influent_vals.mean()
    mean_effluent = effluent_vals.mean()
    mean_diff = mean_effluent - mean_influent

    direction = 'gained' if mean_diff > 0 else 'lost'

    results_list.append({
        'ARG': arg_name,
        'p-value': p_value,
        'mean_influent': mean_influent,
        'mean_effluent': mean_effluent,
        'mean_diff': mean_diff,
        'direction': direction
    })


results_df = pd.DataFrame(results_list)

fdr_values = multipletests(results_df['p-value'], method='fdr_bh')[1]
results_df['FDR'] = fdr_values

sig_df = results_df[results_df['FDR'] < 0.05].sort_values('FDR')

sig_df.to_csv("129_differential_ARGs.csv", index=False)
print(sig_df)
