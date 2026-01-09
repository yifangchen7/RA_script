import csv
def convert_residue_name(residue):
    # Dictionary mapping three-letter residue codes to one-letter codes
    residue_map = {
        'ASN': 'N',
        'UNL': 'Lig',
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

def process_tsv(input_file, output_file):
    with open(input_file, 'r') as tsv_in, open(output_file, 'w', newline='') as tsv_out:
        reader = csv.reader(tsv_in, delimiter='\t')
        writer = csv.writer(tsv_out, delimiter='\t')
        
        # Skip the first two rows
        next(reader)  # Skip the first row
        next(reader)  # Skip the second row
        
        for row in reader:
            # Apply find and replace operations
            row = [item.replace('X:', '') if item.startswith('X:') else item for item in row]
            row = [item.replace('UNL:', 'Lig:') if item.startswith('UNL:') else item for item in row]
            row = [convert_residue_name(item.split(":")[0]) + item.split(":")[1] if ':' in item else item for item in row]
            writer.writerow(row)

input_file = './frequency_1.tsv'
output_file = './frequency_1_converted.tsv'
process_tsv(input_file, output_file)

#plotting
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Read the TSV file into a DataFrame with explicit column names
file_path = './frequency_1_converted.tsv'
column_names = ['residue', 'ligand', 'contact_frequency']
df = pd.read_csv(file_path, delimiter='\t', names=column_names)

# Pivot the DataFrame for heatmap plotting
pivot_data = df.pivot(index='ligand', columns='residue', values='contact_frequency')

# Plot heatmap
plt.figure(figsize=(15, 1))
sns.heatmap(pivot_data, cmap="Blues", annot=False)
plt.xticks(rotation=45)
plt.xlabel('Residue')
plt.ylabel('Ligand')
plt.title('Contact Frequency Heatmap')
plt.show()

######reorder the residues based on the id######
import pandas as pd
# Read the TSV file into a DataFrame
df = pd.read_csv("./frequency_1_converted.tsv", delimiter="\t", skiprows=2)
df.columns = ['Lig', 'residue', 'frequency']

# extract residue numbers from a string
def extract_residue_number(string):
    numbers = ''.join(filter(str.isdigit, string))
    return int(numbers) if numbers else 0

# Extract residue numbers from the 'residue' column
df['residue_numbers'] = df['residue'].apply(extract_residue_number)

# Sort the DataFrame by residue_numbers
df_sorted = df.sort_values(by='residue_numbers')
df_sorted.drop(columns=['residue_numbers'], inplace=True)

# Write the sorted DataFrame back to a CSV file
df_sorted.to_csv("./new.csv", index=False)

