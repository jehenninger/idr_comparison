import os
from Bio import SeqIO
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
from bresenham import bresenham

# to do: new graph type, plot amino acids and link interactions
# to do: just give uniprot numbers and look them up
# to do: ability to do multiple IDRs
# to do: add PTMs
# to do: check loops to make sure you reach the end of the protein
# to do: figure out bug where the smaller IDR is being used to plot. Also need to align the amino
# acid graphs to line up with the centered IDR

## test script to compare two IDRs to each

idr_uniprot_A = ""
idr_uniprot_B = ""

reverse_A = False
reverse_B = False

fasta_file_path = "/Users/jon/Google Drive (WIBR)/analysis/181020 IDR interaction in silico analysis/med1_polIICTD_IDRs.txt"

## read in fasta file
identifier = []
seq = []
for seq_record in SeqIO.parse(fasta_file_path, "fasta"):
    identifier.append(seq_record.id)
    seq.append(str(seq_record.seq))

## retrieve individual IDR sequences
seq_A_name = identifier[0]
seq_B_name = identifier[1]
seq_A = seq[0]
if reverse_A:
    seq_A = seq_A[::-1]
    seq_A_name = seq_A_name + "_reverse"
seq_B = seq[1]
if reverse_B:
    seq_B = seq_B[::-1]
    seq_B_name = seq_B_name + "_reverse"

seq_A_length = len(seq_A)
seq_B_length = len(seq_B)

print("IDR A: " + seq_A_name)
print(seq_A[0:20] + " ...")
print("Length: ", len(seq_A))

print("IDR B: " + seq_B_name)
print(seq_B[0:20] + " ...")
print("Length: ", len(seq_B))

## generate amino acid interaction matrix
charge_weight = 1
cation_pi_weight = 1
pi_stack_weight = 1
h_bond_weight = 1
same_charge_weight = -1

amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

# define regular expressions
charge_pos_regex = re.compile('[RHK]')
charge_neg_regex = re.compile('[DE]')

cat_pi_pos_regex = re.compile('[RK]')
cat_pi_ring_regex = re.compile('[YFW]')

pi_stack_regex = re.compile('[YFWH]')

h_bond_donor_regex = re.compile('[RNQHKSTWY]')
h_bond_acceptor_regex = re.compile('[NDQEHSTY]')

interaction_lut = pd.DataFrame(np.zeros(shape=(len(amino_acids), len(amino_acids))),
                               index=amino_acids, columns=amino_acids)

for index, row in interaction_lut.iterrows():
    # index is the amino acid

    # charge-charge
    if re.match(charge_pos_regex, index):
        interaction_lut.loc[index, 'D'] = interaction_lut.loc[index, 'D'] + charge_weight
        interaction_lut.loc[index, 'E'] = interaction_lut.loc[index, 'E'] + charge_weight

        interaction_lut.loc[index, 'R'] = interaction_lut.loc[index, 'R'] + same_charge_weight
        interaction_lut.loc[index, 'H'] = interaction_lut.loc[index, 'H'] + same_charge_weight
        interaction_lut.loc[index, 'K'] = interaction_lut.loc[index, 'K'] + same_charge_weight

    if re.match(charge_neg_regex, index):
        interaction_lut.loc[index, 'R'] = interaction_lut.loc[index, 'R'] + charge_weight
        interaction_lut.loc[index, 'H'] = interaction_lut.loc[index, 'H'] + charge_weight
        interaction_lut.loc[index, 'K'] = interaction_lut.loc[index, 'K'] + charge_weight

        interaction_lut.loc[index, 'D'] = interaction_lut.loc[index, 'D'] + same_charge_weight
        interaction_lut.loc[index, 'E'] = interaction_lut.loc[index, 'E'] + same_charge_weight

    # pi-cation
    if re.match(cat_pi_pos_regex, index):
        interaction_lut.loc[index, 'Y'] = interaction_lut.loc[index, 'Y'] + cation_pi_weight
        interaction_lut.loc[index, 'F'] = interaction_lut.loc[index, 'F'] + cation_pi_weight
        interaction_lut.loc[index, 'W'] = interaction_lut.loc[index, 'W'] + cation_pi_weight

    if re.match(cat_pi_ring_regex, index):
        interaction_lut.loc[index, 'R'] = interaction_lut.loc[index, 'R'] + cation_pi_weight
        interaction_lut.loc[index, 'K'] = interaction_lut.loc[index, 'K'] + cation_pi_weight

    # pi-stacking
    if re.match(pi_stack_regex, index):
        interaction_lut.loc[index, 'Y'] = interaction_lut.loc[index, 'Y'] + pi_stack_weight
        interaction_lut.loc[index, 'F'] = interaction_lut.loc[index, 'F'] + pi_stack_weight
        interaction_lut.loc[index, 'W'] = interaction_lut.loc[index, 'W'] + pi_stack_weight
        interaction_lut.loc[index, 'H'] = interaction_lut.loc[index, 'H'] + pi_stack_weight

    # h-bonding
    if re.match(h_bond_donor_regex, index):
        interaction_lut.loc[index, 'N'] = interaction_lut.loc[index, 'N'] + h_bond_weight
        interaction_lut.loc[index, 'D'] = interaction_lut.loc[index, 'D'] + h_bond_weight
        interaction_lut.loc[index, 'Q'] = interaction_lut.loc[index, 'Q'] + h_bond_weight
        interaction_lut.loc[index, 'E'] = interaction_lut.loc[index, 'E'] + h_bond_weight
        interaction_lut.loc[index, 'H'] = interaction_lut.loc[index, 'H'] + h_bond_weight
        interaction_lut.loc[index, 'S'] = interaction_lut.loc[index, 'S'] + h_bond_weight
        interaction_lut.loc[index, 'T'] = interaction_lut.loc[index, 'T'] + h_bond_weight
        interaction_lut.loc[index, 'Y'] = interaction_lut.loc[index, 'Y'] + h_bond_weight

    if re.match(h_bond_acceptor_regex, index):
        interaction_lut.loc[index, 'R'] = interaction_lut.loc[index, 'R'] + h_bond_weight
        interaction_lut.loc[index, 'N'] = interaction_lut.loc[index, 'N'] + h_bond_weight
        interaction_lut.loc[index, 'Q'] = interaction_lut.loc[index, 'Q'] + h_bond_weight
        interaction_lut.loc[index, 'H'] = interaction_lut.loc[index, 'H'] + h_bond_weight
        interaction_lut.loc[index, 'K'] = interaction_lut.loc[index, 'K'] + h_bond_weight
        interaction_lut.loc[index, 'S'] = interaction_lut.loc[index, 'S'] + h_bond_weight
        interaction_lut.loc[index, 'T'] = interaction_lut.loc[index, 'T'] + h_bond_weight
        interaction_lut.loc[index, 'W'] = interaction_lut.loc[index, 'W'] + h_bond_weight
        interaction_lut.loc[index, 'Y'] = interaction_lut.loc[index, 'Y'] + h_bond_weight

# find amino acids that have all 0 scores (we can ignore these later to make things faster)
aa_to_ignore = interaction_lut[(interaction_lut == 0).all(axis=1)].index.tolist()

## sliding window of k-mers and get score of interaction

kmer_len = 9

num_of_windows_A = len(seq_A) - kmer_len + 1
num_of_windows_B = len(seq_B) - kmer_len + 1

score = np.zeros(shape=(num_of_windows_A, num_of_windows_B))
print(score.shape)
for i in range(num_of_windows_A):
    for j in range(num_of_windows_B):
        kmer_A = seq_A[i:(i + kmer_len)]
        kmer_B = seq_B[j:(j + kmer_len)]

        for n in kmer_B:
            if n not in aa_to_ignore:  # might speed things up to not calculate when there aren't interactions
                for k in kmer_A:
                    if k not in aa_to_ignore:  # might speed things up to not calculate when there aren't interactions
                        score[i, j] = score[i, j] + interaction_lut.loc[n, k]

print(score[0:10, 0:10])

heatmap_fig, heatmap_ax = plt.subplots(figsize=(10, 10))
im = heatmap_ax.imshow(score, cmap='coolwarm_r', vmin=-100, vmax=100)
heatmap_ax.invert_yaxis()
plt.colorbar(im, ax=heatmap_ax)

plt.xlabel(seq_B_name)
plt.ylabel(seq_A_name)


## plot interaction strength between two IDRs

# default y-grid space
y_grid_space = 100

# test which IDR is longer, will go on bottom and will center the top IDR
seq_A_on_bottom = False
if seq_A_length > seq_B_length:
    seq_A_on_bottom = True
    longest_idr_size = num_of_windows_A
else:
    longest_idr_size = num_of_windows_B

# initial x coordinates of lines
idr_A_x = list(range(num_of_windows_A))
idr_B_x = list(range(num_of_windows_B))

# modify x coordinates and generate y coordinates of lines
center_correction = round(abs((num_of_windows_A / 2) - (num_of_windows_B / 2)))  # centers top line to bottom line
if seq_A_on_bottom:
    idr_B_x = [x + center_correction for x in idr_B_x]
    idr_A_y = np.ones(num_of_windows_A)
    idr_A_y.astype(int)
    idr_B_y = y_grid_space * np.ones(num_of_windows_B)
    idr_B_y.astype(int)
else:
    idr_A_x = [x + center_correction for x in idr_A_x]
    idr_A_y = y_grid_space * np.ones(num_of_windows_A)
    idr_A_y.astype(int)
    idr_B_y = np.ones(num_of_windows_B)
    idr_B_y.astype(int)

idr_A_y.tolist()
idr_B_y.tolist()

# generate heatmap of interacting pairs
idr_grid = np.zeros(shape=(y_grid_space + 1, longest_idr_size))
print(idr_grid.shape)
print(score.shape)
for i in range(score.shape[0]):
    for j in range(score.shape[1]):

        # get all grid points that intersect the line between amino acids in IDR A and IDR B
        # using Bresenham algorithm
        interacting_grids = list(bresenham(int(idr_A_x[i]), int(idr_A_y[i]),
                                           int(idr_B_x[j]), int(idr_B_y[j])))

        for k in interacting_grids:
            idr_grid[int(k[1]), int(k[0])] = idr_grid[int(k[1]), int(k[0])] + score[i, j]

idr_grid = np.divide(idr_grid, score.size)  # normalize grid score to the total number of elements

# make axes
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 7.5))

# plot top sequence grid
if seq_A_on_bottom:
    grid_top_idr_length = seq_B_length
    seq_top = seq_B
    grid_bottom_idr_length = seq_A_length
    seq_bottom = seq_A
else:
    grid_top_idr_length = seq_A_length
    seq_top = seq_A
    grid_bottom_idr_length = seq_B_length
    seq_bottom = seq_B

amino_acid_grid_top = np.zeros(shape=(len(amino_acids), grid_top_idr_length))
count = 0
for a in seq_top:
    aa_index = amino_acids.index(a)
    amino_acid_grid_top[aa_index, count] = 1

    count = count + 1

im_top = ax1.imshow(amino_acid_grid_top, cmap='binary', vmin=0, vmax=1)
plt.colorbar(im_top, ax=ax1)
ax1.set_xlim(left=-5, right=grid_top_idr_length + 10)
ax1.set_yticks(range(len(amino_acids)))
ax1.set_yticklabels(amino_acids, fontsize=6)
ax1.set_aspect('auto')

# plot bottom sequence grid
amino_acid_grid_bottom = np.zeros(shape=(len(amino_acids), grid_bottom_idr_length))
count = 0
for a in seq_bottom:
    aa_index = amino_acids.index(a)
    amino_acid_grid_bottom[aa_index, count] = 1

    count = count + 1

im_bottom = ax3.imshow(amino_acid_grid_bottom, cmap='binary', vmin=0, vmax=1)
plt.colorbar(im_bottom, ax=ax3)
ax3.set_xlim(left=-5, right=grid_bottom_idr_length + 10)
ax3.set_yticks(range(len(amino_acids)))
ax3.set_yticklabels(amino_acids, fontsize=6)
ax3.set_aspect('auto')

# plot IDR structure lines
ax2.plot(idr_A_x, idr_A_y, 'g')
ax2.plot(idr_B_x, idr_B_y, 'b')

# set max and min for alpha (or line width)
plot_weight_min = 0
plot_weight_max = 100

im2 = ax2.imshow(idr_grid, cmap='coolwarm_r', vmin=-0.5, vmax=0.5)
ax2.invert_yaxis()
ax2.set_ylim(bottom=-30, top=y_grid_space + 30)
ax2.set_xlim(left=-5, right=grid_top_idr_length + 10)
plt.colorbar(im2, ax=ax2)
ax2.set_aspect('auto')

if seq_A_on_bottom:
    ax2.text(num_of_windows_A / 2, -15, seq_A_name, horizontalalignment='center')
    ax2.text(num_of_windows_A / 2, y_grid_space + 15, seq_B_name, horizontalalignment='center')
else:
    ax2.text(num_of_windows_B / 2, -15, seq_B_name, horizontalalignment='center')
    ax2.text(num_of_windows_B / 2, y_grid_space + 15, seq_A_name, horizontalalignment='center')

output_dir = os.path.dirname(fasta_file_path)

plt.savefig(os.path.join(output_dir, seq_A_name + "_" + seq_B_name + "_" + str(kmer_len) + "mer" + ".eps"))
plt.savefig(os.path.join(output_dir, seq_A_name + "_" + seq_B_name + "_" + str(kmer_len) + "mer" + ".png"))

plt.show()
