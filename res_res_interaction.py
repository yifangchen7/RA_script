import numpy as np
import matplotlib.pyplot as plt

cutoff = 0.475        # two objects will be considered as contacting when their distance is below this cutoff

name_1 = 'VMP1'       # the name of the protein you want on the x axis
n_res_1 = 385         # number of residues for the protein you want on the x axis of the plot
res_offset_1 = 22     # protein residue will start from this number

name_2 = 'ATG2'        # the name of the protein you want on the y axis
n_res_2 = 307         # number of residues for the protein you want on the y axis of the plot
res_offset_2 = 166

file = '/Users/chenyifang/Desktop/RA/VMP1_ATG2/R1/min_res_VMP1_ATG2.xvg'   #gmx pairdist
save_dir = '/Users/chenyifang/Desktop/RA/VMP1_ATG2/R1/PPI_VMP1_ATG2.png'

n = 10                # top n contacting indices for output

#####################################
###### no changes needed below ######
#####################################

interaction_matrix = np.zeros((n_res_2+res_offset_2, n_res_1+res_offset_1))

n_frames_total = 0
n_contacts_total = 0

data = np.genfromtxt(file, skip_header=24, unpack=True)
n_frames = len(data[0, :])
dist = data[1:, :]

for i in range(n_res_1):
    for j in range(n_res_2):
        contact_frames = np.sum(dist[(i*n_res_2)+j,:] < cutoff)
        occupancy = contact_frames / n_frames * 100
        interaction_matrix[j+res_offset_2, i+res_offset_1] = occupancy

print(interaction_matrix.max())
flattend = interaction_matrix.flatten()
sorted_indices = np.argsort(flattend)
top_n_indices = sorted_indices[-n:]
top_n_indices = top_n_indices[::-1]
with open(f'top_{n}_indices_and_values.txt', 'w') as file:
    for index in top_n_indices:
        value = flattend[index]
        original_indices = np.unravel_index(index, interaction_matrix.shape)
        print(f'index: {original_indices}, Value: {value}')
        file.write(f'index: {original_indices}, Value: {value}/n')

plt.figure(figsize=(12, 10))
plt.imshow(interaction_matrix, aspect='auto', origin='lower', cmap='Reds', vmin=0)

cbar = plt.colorbar(aspect=15)
cbar.ax.yaxis.set_tick_params(labelsize=12)
cbar.set_label('Occupancy (%)', fontsize=12)

plt.xlabel(name_1, fontsize=14)
plt.ylabel(name_2, fontsize=14)

plt.xlim(res_offset_1)
plt.ylim(res_offset_2)

#plt.xticks(ticks=np.arange(0, n_res, 100), labels=np.arange(0, n_res, 100), fontsize=12)
#plt.yticks(ticks=range(len(lip_list)), labels=lip_list, fontsize=14)

plt.tight_layout()
plt.savefig(save_dir)
plt.show()