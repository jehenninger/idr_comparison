import os
from Bio import SeqIO
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt

## script to find mean charge of protein region by sliding window

fasta_file_path = "/Users/jon/data_analysis/young/Ozgur_O/Chromatin_IDR_selected.fasta"
kmer_len = 5

## read in fasta file
identifier = []
seq = []
for seq_record in SeqIO.parse(fasta_file_path, "fasta"):
    identifier.append(seq_record.id)
    seq.append(str(seq_record.seq))

## generate amino acid interaction matrix
pos_charge_weight = 1.0
neg_charge_weight = -1.0
histidine_charge_weight = 0.5

# amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

amino_acids = ['G', 'A', 'S', 'T', 'C', 'V', 'L', 'I', 'M', 'P', 'F', 'Y', 'W', 'D', 'E', 'N', 'Q', 'H', 'K', 'R']

# define regular expressions
charge_pos_regex = re.compile('[RK]')
charge_neg_regex = re.compile('[DE]')
histidine_regex = re.compile('[H]')

## retrieve individual IDR sequences
for idx, s in enumerate(seq):
    seq_name = identifier[idx]
    seq_length = len(s)

    print("Protein: " + seq_name)

    print("Length: ", seq_length)

    ## sliding window of k-mers and get score of interaction

    num_of_windows = seq_length - kmer_len + 1

    score = np.zeros(shape=(num_of_windows, 1))

    for i in range(num_of_windows):
        kmer = s[i:(i + kmer_len)]
        kmer_window_score = np.zeros(shape=(kmer_len,1))
        for count, k in enumerate(kmer):
            if re.match(charge_pos_regex, k):
                kmer_window_score[count] = kmer_window_score[count] + pos_charge_weight
            elif re.match(charge_neg_regex, k):
                kmer_window_score[count] = kmer_window_score[count] + neg_charge_weight
            elif re.match(histidine_regex, k):
                kmer_window_score[count] = kmer_window_score[count] + histidine_charge_weight

        score[i] = score[i] + np.mean(kmer_window_score)  # calculate mean charge of window

    # make axes
    fig, (ax1, ax2) = plt.subplots(2, 1)

    # plot top sequence grid
    grid_length = seq_length
    sequence = s

    amino_acid_grid = np.zeros(shape=(len(amino_acids), grid_length))
    count = 0
    for a in sequence:
        aa_index = amino_acids.index(a)
        amino_acid_grid[aa_index, count] = 1

        count = count + 1

    im = ax1.imshow(amino_acid_grid, cmap='binary', vmin=0, vmax=1)

    ax1.set_xlim(left=-5, right=grid_length + 10)
    ax1.set_yticks(range(len(amino_acids)))
    ax1.set_yticklabels(amino_acids, fontsize=6)
    ax1.set_aspect('auto')

    # plot IDR structure lines

    upper_charge = 0.25
    lower_charge = -0.25
    mask_upper = np.ma.masked_where(score < upper_charge, score)
    mask_lower = np.ma.masked_where(score > lower_charge, score)
    # mask_middle = np.ma.masked_where(np.logical_or(score < -0.15, score > 0.15), score)
    x = range(seq_length - kmer_len + 1)
    ax2.plot(x, score, '-', color='0.8')
    ax2.plot(x, mask_upper, '-b')
    ax2.plot(x, mask_lower, '-r')

    ax2.lines[0].set_linewidth(0.5)
    ax2.lines[1].set_linewidth(1.0)
    ax2.lines[2].set_linewidth(1.0)
    ax2.lines[1].set_alpha(0.3)
    ax2.lines[2].set_alpha(0.3)

    ax2.set_ylim(bottom=-1, top=1)
    ax2.set_xlim(left=-5, right=grid_length + 10)
    ax2.set_aspect('auto')

    fig.suptitle(seq_name)

    output_dir = os.path.dirname(fasta_file_path)

    plt.savefig(os.path.join(output_dir, seq_name + "_" + str(kmer_len) + "mer" + ".eps"))
    plt.savefig(os.path.join(output_dir, seq_name + "_" + str(kmer_len) + "mer" + ".png"), dpi=300)

    plt.close()

