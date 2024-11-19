#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

file_path = '/Users/chenyifang/Desktop/output.tsv'

# Read the TSV file into a DataFrame
df = pd.read_csv(file_path, sep='\t', header=None, names=['Residue', 'Ligand', 'Frequency'])

# Sort the DataFrame by the 'Frequency' column in descending order
df_sorted = df.sort_values(by='Frequency', ascending=False)

# To ensure the correct order of residues, set the order of columns manually based on sorted DataFrame
residues_order = df_sorted['Residue'].tolist()

# Pivot the DataFrame for heatmap
heatmap_data = df_sorted.pivot(index='Ligand', columns='Residue', values='Frequency')
heatmap_data = heatmap_data[residues_order]

# Plotting the heatmap
plt.figure(figsize=(15, 2))
sns.heatmap(heatmap_data, cmap='Blues', annot=False)
plt.xticks(rotation=45)
plt.title('Heatmap of Contact Frequency')
plt.show()


# In[ ]:




