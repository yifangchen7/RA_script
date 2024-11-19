#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the summed contact frequency file
file_path = '/Users/chenyifang/Desktop/NEW_AQP_LIP/R1/summed_contact_frequencies.tsv'
df = pd.read_csv(file_path, sep='\t')
df.columns = ['residue_1', 'residue_2', 'summed_contact_frequency']

# Convert the summed_contact_frequency column to numeric data type
df['summed_contact_frequency'] = pd.to_numeric(df['summed_contact_frequency'], errors='coerce')

# Divide the summed_contact_frequency values by 4 to get the average
df['average_contact_frequency'] = df['summed_contact_frequency'] / 4

# Handle duplicates by averaging the contact frequencies for the same residue pairs
df = df.groupby(['residue_1', 'residue_2'], as_index=False).agg({'average_contact_frequency': 'mean'})

# Filter out contact frequencies less than 0.3
filtered_df = df[df['average_contact_frequency'] >= 0.3]

# Function to extract residue numbers from a string (ignoring the residue name)
def extract_residue_number(string):
    return int(''.join(filter(str.isdigit, string)))

# Extract numeric values from residue IDs for sorting
filtered_df['residue_1_num'] = filtered_df['residue_1'].apply(extract_residue_number)
filtered_df['residue_2_num'] = filtered_df['residue_2'].apply(extract_residue_number)
sorted_filtered_df = filtered_df.sort_values(by=['residue_1_num', 'residue_2_num'])

# Set the order of categories based on the sorted residue numbers
residue_1_categories = sorted_filtered_df['residue_1'].unique()
residue_2_categories = sorted_filtered_df['residue_2'].unique()

# Ensure residue_2_categories are also sorted by their residue numbers
residue_2_categories_sorted = sorted(residue_2_categories, key=lambda x: extract_residue_number(x))

# Pivot the DataFrame using residue names
pivot_df = sorted_filtered_df.pivot(index='residue_1', columns='residue_2', values='average_contact_frequency')

# Reindex to enforce the correct order
pivot_df = pivot_df.reindex(index=residue_1_categories, columns=residue_2_categories_sorted)

# Draw a heatmap using the residue names sorted by their numbers
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
cbar.ax.yaxis.label.set_size(14)  # Set the font size of the color bar label
cbar.ax.yaxis.label.set_weight('bold')  # Make the label bold

# Improve titles and labels
plt.title('Repeat1', fontsize=18, weight='bold', pad=20)
plt.xlabel('Residue 2 (LIP)', fontsize=14, weight='bold')
plt.ylabel('Residue 1 (AQP)', fontsize=14, weight='bold')
plt.xticks(fontsize=12, rotation=45, ha='right', weight='bold')
plt.yticks(fontsize=12, rotation=0, weight='bold')
plt.tight_layout()
plt.savefig('/Users/chenyifang/Desktop/NEW_AQP_LIP/contact_frequency_R1_sorted_by_residue_number.png', dpi=600, bbox_inches='tight')
plt.show()

