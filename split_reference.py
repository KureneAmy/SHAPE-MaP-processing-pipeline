#!/usr/bin/env python3

import argparse
import os
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description='Split a FASTA file into DNA and RNA version per sequence.')
    parser.add_argument('-i', '--input_fasta', required=True, help='Path to the input multi-FASTA file.')
    parser.add_argument('-d', '--dna_dir', required=True, help='Directory to output DNA version FASTA files.')
    parser.add_argument('-r', '--rna_dir', required=True, help='Directory to output RNA version FASTA files.')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite output files if they exist.')

    args = parser.parse_args()

    os.makedirs(args.dna_dir, exist_ok=True)
    os.makedirs(args.rna_dir, exist_ok=True)

    all_sequences = list(SeqIO.parse(args.input_fasta, "fasta"))
    for record in all_sequences:
        seq_name = record.id
        # DNA
        dna_seq = str(record.seq).replace('U', 'T').replace('u', 't')
        dna_file = os.path.join(args.dna_dir, f"{seq_name}.fa")
        with open(dna_file, 'w') as f:
            f.write(f">{seq_name}\n{dna_seq}\n")

        # RNA
        rna_seq = str(record.seq).replace('T', 'U').replace('t', 'u')
        rna_file = os.path.join(args.rna_dir, f"{seq_name}.fa")
        with open(rna_file, 'w') as f:
            f.write(f">{seq_name}\n{rna_seq}\n")

if __name__ == "__main__":
    main()