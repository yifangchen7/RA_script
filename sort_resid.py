#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

file_path = '/Users/chenyifang/Desktop/NEW_AQP_LIP/R1/summed_contact_frequencies.tsv'
df = pd.read_csv(file_path, sep='\t')
df.columns = ['residue_1', 'residue_2', 'summed_contact_frequency']

# Convert the summed_contact_frequency column to numeric data type
df['summed_contact_frequency'] = pd.to_numeric(df['summed_contact_frequency'], errors='coerce')

# Divide the summed_contact_frequency values by 4 to get the average
df['average_contact_frequency'] = df['summed_contact_frequency'] / 4

# Handle duplicates by averaging the contact frequencies for the same residue pairs
df = df.groupby(['residue_1', 'residue_2'], as_index=False).agg({'average_contact_frequency': 'mean'})

# Extract numeric values from residue IDs for sorting, while keeping the original names intact
df['residue_1_num'] = df['residue_1'].str.extract('(\d+)').astype(int)
df['residue_2_num'] = df['residue_2'].str.extract('(\d+)').astype(int)

# Sort the dataframe by the numeric part of residue_1 and residue_2
sorted_df = df.sort_values(by=['residue_1_num', 'residue_2_num'])

# Save the sorted results to a new TSV file with original residue names intact
sorted_df[['residue_1', 'residue_2', 'average_contact_frequency']].to_csv('/Users/chenyifang/Desktop/NEW_AQP_LIP/R1/sorted_averaged_contact_frequencies_by_residue_numeric.tsv', sep='\t', index=False)

# Draw a heatmap using the original residue names, with residues sorted by numeric order
pivot_df = sorted_df.pivot('residue_1', 'residue_2', 'average_contact_frequency')

# Customize heatmap for publication
plt.figure(figsize=(30, 20))  # Adjust size for clarity
heatmap = sns.heatmap(
    pivot_df, 
    cmap='YlOrBr', 
    annot=False, 
    fmt=".3f", 
    cbar_kws={
        'label': 'Contact Frequency',
        'shrink': 0.5,  
        'aspect': 30,
        'pad': 0.02
    }, 
    linewidths=0.5,  
    linecolor='gray',  
    square=True  
)

# Adjust the label size of the color bar
cbar = heatmap.collections[0].colorbar
cbar.ax.yaxis.label.set_size(14) 
cbar.ax.yaxis.label.set_weight('bold')

# Improve titles and labels
plt.title('Repeat1', fontsize=18, weight='bold', pad=20)
plt.xlabel('Residue 2 (LIP)', fontsize=14, weight='bold')
plt.ylabel('Residue 1 (AQP)', fontsize=14, weight='bold')
plt.xticks(fontsize=12, rotation=45, ha='right', weight='bold')
plt.yticks(fontsize=12, rotation=0, weight='bold')
plt.tight_layout()
#plt.savefig('/Users/chenyifang/Desktop/AQP_LIP/contact_frequency_R1_by_residue_numeric.png', dpi=600, bbox_inches='tight')
plt.show()

