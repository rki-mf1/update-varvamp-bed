# -*- coding: utf-8 -*-
"""correct_primer_positions.py
Namuun Battur
naominamuun@gmail.com
required packages: pandas, biopython, fuzzysearch, sys, os.path
"""

from fuzzysearch import find_near_matches
from Bio import SeqIO
import pandas as pd
from sys import argv
import os


def get_absolute_path(relative_path):
    '''Convert the given relative path to an absolute path'''
    return os.path.abspath(relative_path)


def verify(sequence):
    '''This code verifies if a sequence is a DNA or contains ambiguous bases according to the IUPAC code.'''
    seq = set(sequence)

    valid_bases = {"A", "T", "C", "G", "Y", "R", "S", "W", "K", "M", "B", "D", "H", "V", "N"}

    if seq.issubset(valid_bases):
        return "DNA"
    else:
        return "Invalid DNA sequence"


def rev_comp_if(seq):
    '''Reverse complement the DNA sequence if valid, with ambiguity handling'''
    comp = []
    if verify(seq) == "DNA":
        for base in seq:
            if base == "A": comp.append("T")
            elif base == "G": comp.append("C")
            elif base == "T": comp.append("A")
            elif base == "C": comp.append("G")
            elif base in ["Y", "R", "S", "W", "K", "M", "B", "D", "H", "V", "N"]:
                comp.append("N")  # Mask ambiguous bases with 'N'
    else:
        return "Invalid DNA Sequence"

    return "".join(comp[::-1])  # Return the reverse complement


def update_positions(primer_df, record_dict):
    '''Update the primer positions in the provided reference sequence'''
    try:
        start_pos, end_pos = [], []

        for _, row in primer_df.iterrows():
            primer = str(row['seq']).upper()
            revprimer = rev_comp_if(primer)
            reference = str(record_dict[list(record_dict)[0]].seq).upper()

            primer_match = find_near_matches(primer, reference, max_l_dist=2)
            rev_primer_match = find_near_matches(revprimer, reference, max_l_dist=2)

            if len(primer_match) == 1:
                start_pos.append(primer_match[0].start)
                end_pos.append(primer_match[0].end)
            elif len(rev_primer_match) == 1:
                start_pos.append(rev_primer_match[0].start)
                end_pos.append(rev_primer_match[0].end)
            elif len(primer_match) > 1 or len(rev_primer_match) > 1:
                raise Exception("Primer sequence found multiple times in reference sequence.")
            else:
                raise Exception("Primer sequence not found in reference sequence.")

        seq_size = [a - b for a, b in zip(end_pos, start_pos)]
        if seq_size != list(primer_df['size']):
            raise Exception("The primer size before and after the update are not identical")

        return start_pos, end_pos
    except Exception as e:
        raise Exception(f"Error in updating positions: {e}")


def main(argv):
    if len(argv) != 4:
        print('Please add required input files \n e.g.: python3 correct_primer_positions.py [path/to/Reference/.fasta] [path/to/primer/.tsv] [path/to/primer/.bed]')
        return
    else:
        try:
            # Ensure paths are absolute
            fasta_path = get_absolute_path(argv[1])
            primer_tsv_path = get_absolute_path(argv[2])
            primer_bed_path = get_absolute_path(argv[3])

            # Check if files exist
            if not os.path.exists(fasta_path):
                print(f"FASTA file not available: {fasta_path}")
                return
            if not os.path.exists(primer_tsv_path):
                print(f"TSV file not available: {primer_tsv_path}")
                return
            if not os.path.exists(primer_bed_path):
                print(f"BED file not available: {primer_bed_path}")
                return

            # Read files
            record_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
            primer_df = pd.read_csv(primer_tsv_path, sep="\t")
            primer_bed_df = pd.read_csv(primer_bed_path, sep="\t", header=None)

            # Update primer positions
            start_pos, end_pos = update_positions(primer_df, record_dict)

            # Update the BED file with new positions
            primer_bed_df.iloc[:, 1] = start_pos
            primer_bed_df.iloc[:, 2] = end_pos

            # Ensure the output directory exists
            output_dir = os.path.dirname(primer_bed_path)
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)

            output_file = f"{output_dir}/{os.path.basename(primer_bed_path)}.position.corrected.bed"
            primer_bed_df.to_csv(output_file, index=False, sep='\t', header=False)

            print(f"Primer.bed file has been updated. See: {output_file}")

        except Exception as e:
            print(f"An error has occurred: {e}")


if __name__ == "__main__":
    main(argv)