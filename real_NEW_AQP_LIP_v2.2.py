#!/usr/bin/env python
# coding: utf-8

# In[ ]:


###Step 1###
#use GetContacts to get corresponding tsv files#
#..//getcontacts/get_dynamic_contacts.py --topology ref.pdb --trajectory MD.xtc --itypes all --output contacts_D_H.tsv --sele "chain D" --sele2 "chain H"
#..//data/public/software/getcontacts/get_contact_frequencies.py --input_files contacts_D_H.tsv --output_file resfrequencies_D_H.tsv"

###Step 2###
#mapping three-letter residue codes to one-letter codes#
import csv
def process_tsv(input_file, output_file):
    #with open(input_file, 'r') as tsv_in, open(output_file, 'w', newline='') as tsv_out:
        reader = csv.reader(tsv_in, delimiter='\t')
        writer = csv.writer(tsv_out, delimiter='\t')
        for row in reader:
            row = [item.replace('D:', '') if item.startswith('D:') else item for item in row]
            row = [item.replace('H:', '') if item.startswith('H:') else item for item in row]
            row = [convert_residue_name(item.split(":")[0]) + item.split(":")[1] if ':' in item else item for item in row]
            writer.writerow(row)

def convert_residue_name(residue):
    residue_map = {
        'ASN': 'N',
        'GLY': 'G',
        'TYR': 'Y',
        'ASP': 'D',
        'SER': 'S',
        'ALA': 'A',
        'ARG': 'R',
        'CYS': 'C',
        'GLN': 'Q',
        'GLU': 'E',
        'HIS': 'H',
        'ILE': 'I',
        'LEU': 'L',
        'LYS': 'K',
        'MET': 'M',
        'PHE': 'F',
        'PRO': 'P',
        'THR': 'T',
        'TRP': 'W',
        'VAL': 'V'
    }
    return residue_map.get(residue, residue)
#input_file = '..//resfrequencies_D_H.tsv'
#output_file = '..//converted_D_H.tsv'
#process_tsv(input_file, output_file)

###Step3###
#sum up and average the caculations for different chains and plot for one repeat#
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from collections import defaultdict

input_files = [
    './R1/converted_A_E_1.tsv',
    './R1/converted_B_F_1.tsv',
    './R1/converted_C_G_1.tsv',
    './R1/converted_D_H_1.tsv'
]

# Dictionary to store summed contact frequencies
contact_freq = defaultdict(float)

# Process each file to sum contact frequencies
for file in input_files:
    df = pd.read_csv(file, sep='\t', skiprows=2, header=None, names=['residue_1', 'residue_2', 'contact_frequency'])
    
    for _, row in df.iterrows():
        pair = (row['residue_1'], row['residue_2'])
        contact_freq[pair] += row['contact_frequency']

# Convert summed frequencies to DataFrame and calculate the average
output_df = pd.DataFrame(
    [(res1, res2, freq/4) for (res1, res2), freq in contact_freq.items()],
    columns=['residue_1', 'residue_2', 'average_contact_frequency']
)

# Handle duplicates by averaging the contact frequencies for the same residue pairs
output_df = output_df.groupby(['residue_1', 'residue_2'], as_index=False).agg({'average_contact_frequency': 'mean'})

# Sort the dataframe by average_contact_frequency in descending order
#sorted_df = output_df.sort_values(by='average_contact_frequency', ascending=False)
#sorted_df.to_csv('sorted_averaged_contact_frequencies.tsv', sep='\t', index=False)
#pivot_df = sorted_df.pivot('residue_1', 'residue_2', 'average_contact_frequency')

#---------------------------------------------------------------
# Filter out contact frequencies less than 0.3
filtered_df = output_df[output_df['average_contact_frequency'] >= 0.2]

# Function to extract residue numbers from a string (ignoring the residue name)
def extract_residue_number(string):
    return int(''.join(filter(str.isdigit, string)))

# Extract numeric values from residue IDs for sorting
filtered_df['residue_1_num'] = filtered_df['residue_1'].apply(extract_residue_number)
filtered_df['residue_2_num'] = filtered_df['residue_2'].apply(extract_residue_number)

# Sort the dataframe by the numeric part of residue_1 and residue_2
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


# Save the filtered and sorted results to a new TSV file with original residue names intact
sorted_filtered_df[['residue_1_num', 'residue_2_num', 'average_contact_frequency']].to_csv(
    './R1/sorted_filtered_averaged_contact_frequencies_by_residue_numeric.tsv',
    sep='\t', index=False)

#---------------------------------------------------------------


# Generate a heatmap
plt.figure(figsize=(15, 10))  # Adjust size for clarity
heatmap = sns.heatmap(
    pivot_df, 
    #cmap='YlOrBr', 
    cmap='Blues',
    vmin=0.0, vmax=1,
    #cmap='PiYG',
    #center=0.6,
    annot=True,
    fmt=".2f", 
    cbar_kws={
        'label': 'Contact Frequency',
        'shrink': 0.5,  
        'aspect': 10,
        'pad': 0.02
    }, 
    linewidths=0.5,  
    linecolor='gray',  
    square=True  
)

# Adjust the label size of the color bar
cbar = heatmap.collections[0].colorbar
cbar.ax.yaxis.label.set_size(20) 
cbar.ax.yaxis.label.set_weight('bold')
plt.title('RUN1', fontsize=20, weight='bold', pad=20)
plt.xlabel('Residue Index (LIP5)', fontsize=20, weight='bold')
plt.ylabel('Residue Index (AQP2)', fontsize=20, weight='bold')
plt.xticks(fontsize=20, rotation=45, ha='right') #, weight='bold')
plt.yticks(fontsize=20, rotation=45) #, weight='bold')
plt.tight_layout()
plt.savefig('./contact_frequency_R1_cutoff.pdf', dpi=600, bbox_inches='tight')
#plt.show()













input_files = [
    './R2/converted_A_E_2.tsv',
    './R2/converted_B_F_2.tsv',
    './R2/converted_C_G_2.tsv',
    './R2/converted_D_H_2.tsv'
]

# Dictionary to store summed contact frequencies
contact_freq = defaultdict(float)

# Process each file to sum contact frequencies
for file in input_files:
    df = pd.read_csv(file, sep='\t', skiprows=2, header=None, names=['residue_1', 'residue_2', 'contact_frequency'])
    
    for _, row in df.iterrows():
        pair = (row['residue_1'], row['residue_2'])
        contact_freq[pair] += row['contact_frequency']

# Convert summed frequencies to DataFrame and calculate the average
output_df = pd.DataFrame(
    [(res1, res2, freq/4) for (res1, res2), freq in contact_freq.items()],
    columns=['residue_1', 'residue_2', 'average_contact_frequency']
)

# Handle duplicates by averaging the contact frequencies for the same residue pairs
output_df = output_df.groupby(['residue_1', 'residue_2'], as_index=False).agg({'average_contact_frequency': 'mean'})

# Sort the dataframe by average_contact_frequency in descending order
#sorted_df = output_df.sort_values(by='average_contact_frequency', ascending=False)
#sorted_df.to_csv('sorted_averaged_contact_frequencies.tsv', sep='\t', index=False)
#pivot_df = sorted_df.pivot('residue_1', 'residue_2', 'average_contact_frequency')

#---------------------------------------------------------------
# Filter out contact frequencies less than 0.3
filtered_df = output_df[output_df['average_contact_frequency'] >= 0.2]

# Function to extract residue numbers from a string (ignoring the residue name)
def extract_residue_number(string):
    return int(''.join(filter(str.isdigit, string)))

# Extract numeric values from residue IDs for sorting
filtered_df['residue_1_num'] = filtered_df['residue_1'].apply(extract_residue_number)
filtered_df['residue_2_num'] = filtered_df['residue_2'].apply(extract_residue_number)

# Sort the dataframe by the numeric part of residue_1 and residue_2
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


# Save the filtered and sorted results to a new TSV file with original residue names intact
sorted_filtered_df[['residue_1_num', 'residue_2_num', 'average_contact_frequency']].to_csv(
    './R2/sorted_filtered_averaged_contact_frequencies_by_residue_numeric.tsv',
    sep='\t', index=False)

#---------------------------------------------------------------

#---------------------------------------------------------------
# Generate a heatmap
plt.figure(figsize=(15, 10))  # Adjust size for clarity
heatmap = sns.heatmap(
    pivot_df,
    #cmap='YlOrBr', 
    cmap='Blues',
    vmin=0.0, vmax=1,
    #cmap='PiYG',
    #center=0.6,
    annot=True,
    fmt=".2f",
    cbar_kws={
        'label': 'Contact Frequency',
        'shrink': 0.5,
        'aspect': 10,
        'pad': 0.02
    },
    linewidths=0.5,
    linecolor='gray',
    square=True
)

# Adjust the label size of the color bar
cbar = heatmap.collections[0].colorbar
cbar.ax.yaxis.label.set_size(20)
cbar.ax.yaxis.label.set_weight('bold')
plt.title('RUN2', fontsize=20, weight='bold', pad=20)
plt.xlabel('Residue Index (LIP5)', fontsize=20, weight='bold')
plt.ylabel('Residue Index (AQP2)', fontsize=20, weight='bold')
plt.xticks(fontsize=20, rotation=45, ha='right') #, weight='bold')
plt.yticks(fontsize=20, rotation=45) #, weight='bold')
plt.tight_layout()

plt.savefig('./contact_frequency_R2_cutoff.pdf', dpi=600, bbox_inches='tight')
#plt.show()
















###Step 4###
#averaged value for different trajectories#
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

file1 = 'R1/sorted_averaged_contact_frequencies.tsv'
file2 = 'R2/sorted_averaged_contact_frequencies.tsv'
df1 = pd.read_csv(file1, sep='\t')
df2 = pd.read_csv(file2, sep='\t')
df1.columns = ['residue_1', 'residue_2', 'summed_contact_frequency']
df2.columns = ['residue_1', 'residue_2', 'summed_contact_frequency']

# Convert the summed_contact_frequency columns to numeric data types
df1['summed_contact_frequency'] = pd.to_numeric(df1['summed_contact_frequency'], errors='coerce')
df2['summed_contact_frequency'] = pd.to_numeric(df2['summed_contact_frequency'], errors='coerce')

# Merge the dataframes on residue_1 and residue_2
merged_df = pd.merge(df1, df2, on=['residue_1', 'residue_2'], suffixes=('_repeat1', '_repeat2'))

# Calculate the average contact frequency
merged_df['average_contact_frequency'] = merged_df[['summed_contact_frequency_repeat1', 'summed_contact_frequency_repeat2']].mean(axis=1)

# Handle any duplicates after averaging (if necessary)
merged_df = merged_df.groupby(['residue_1', 'residue_2']).agg({'average_contact_frequency': 'mean'}).reset_index()

# Sort by the average contact frequency in descending order
sorted_df = merged_df.sort_values(by='average_contact_frequency', ascending=False)

# Save the sorted and averaged results to a new TSV file
sorted_df.to_csv('LIP_AQP_averaged_contact_frequencies_cutoff.tsv', sep='\t', index=False)

# Filter out contact frequencies less than 0.3
filtered_df = sorted_df[sorted_df['average_contact_frequency'] >= 0.2]

# Function to extract residue numbers from a string (ignoring the residue name)
def extract_residue_number(string):
    return int(''.join(filter(str.isdigit, string)))

# Extract numeric values from residue IDs for sorting
filtered_df['residue_1_num'] = filtered_df['residue_1'].apply(extract_residue_number)
filtered_df['residue_2_num'] = filtered_df['residue_2'].apply(extract_residue_number)

# Sort the dataframe by the numeric part of residue_1 and residue_2
sorted_filtered_df = filtered_df.sort_values(by=['residue_1_num', 'residue_2_num'])

# Save the filtered and sorted results to a new TSV file
sorted_filtered_df.to_csv(
    './LIP_AQP_sorted_filtered_averaged_contact_frequencies_by_residue_numeric.tsv',
    sep='\t', index=False
)

# Set the order of categories based on the sorted residue numbers
residue_1_categories = sorted_filtered_df['residue_1'].unique()
residue_2_categories = sorted_filtered_df['residue_2'].unique()

# Ensure residue_2_categories are also sorted by their residue numbers
residue_2_categories_sorted = sorted(residue_2_categories, key=lambda x: extract_residue_number(x))

# Pivot the DataFrame using residue names
pivot_df = sorted_filtered_df.pivot(index='residue_1', columns='residue_2', values='average_contact_frequency')

# Reindex to enforce the correct order
pivot_df = pivot_df.reindex(index=residue_1_categories, columns=residue_2_categories_sorted)

# Generate a heatmap
plt.figure(figsize=(15, 10))  # Adjust size for clarity
heatmap = sns.heatmap(
    pivot_df,
    cmap='Blues',
    vmin=0.0, vmax=1,
    annot=True,
    fmt=".2f",
    cbar_kws={
        'label': 'Contact Frequency',
        'shrink': 0.5,
        'aspect': 10,
        'pad': 0.02
    },
    linewidths=0.5,
    linecolor='gray',
    square=True
)

# Adjust the label size of the color bar
cbar = heatmap.collections[0].colorbar
cbar.ax.yaxis.label.set_size(20)
cbar.ax.yaxis.label.set_weight('bold')

plt.xlabel('Residue Index (LIP5)', fontsize=20, weight='bold')
plt.ylabel('Residue Index (AQP2)', fontsize=20, weight='bold')
plt.xticks(fontsize=20, rotation=45, ha='right')
plt.yticks(fontsize=20, rotation=45)
plt.tight_layout()

# Save the heatmap
plt.savefig('./APQ2-LIP5_contact_frequency_avg_cutoff.pdf', dpi=600, bbox_inches='tight')
# plt.show()

