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
import os.path

def verify(sequence):
	'''This code verifies if a sequence is a DNA'''

	# set the input sequence
	seq = set(sequence)

	# confirm if its elements is equal to the set of valid DNA bases
	# Use a union method to ensure the sequence is verified if does not
	# contain all the bases
	if seq == {"A", "T", "C", "G"}.union(seq):
		return "DNA"
	else:
		return "Invalid DNA sequence"

def rev_comp_if(seq):
# adjusted code from https://www.geeksforgeeks.org/reverse-complement-of-dna-strand-using-python/
	comp = []
	if verify(seq) == "DNA":
		for base in seq:
			if base == "A": comp.append("T")
			elif base == "G": comp.append("C")
			elif base == "T": comp.append("A")
			elif base == "C": comp.append("G")
			### mask ambigious bases with N
			elif base == "Y": comp.append("N")
			elif base == "R": comp.append("N")
			elif base == "S": comp.append("N")
			elif base == "W": comp.append("N")
			elif base == "K": comp.append("N")
			elif base == "M": comp.append("N")
			elif base == "B": comp.append("N")
			elif base == "D": comp.append("N")
			elif base == "H": comp.append("N")
			elif base == "V": comp.append("N")
			elif base == "N": comp.append("N")
	else:
		return "Invalid DNA Sequence"

	# reverse the sequence
	comp_rev = comp[::-1]

	# convert list to string
	comp_rev = "".join(comp_rev)
	return comp_rev

"""considered ambigious bases (https://www.bioinformatics.org/sms/iupac.html)
"""

def update_positions(primer_df, record_dict):
    try:
        start_pos = []
        end_pos = []

        for index, row in primer_df.iterrows():
            primer = str(row['seq']).upper()

            # reverse compliment the primer sequence
            revprimer = rev_comp_if(str(row['seq']).upper())

            # get reference sequence of the 1st entry from fasta file
            reference=str(record_dict[list(record_dict)[0]].seq).upper()

            # fuzzy search primer sequence in reference sequence
            # allowing max. 2 mismatches (Levenshtein distance of 2)
            primer_match = find_near_matches(primer, reference, max_l_dist=2)
            rev_primer_match = find_near_matches(revprimer, reference, max_l_dist=2)

            if len(primer_match) == 1:
                start_pos.append(primer_match[0].start)
                end_pos.append(primer_match[0].end)
            elif len(rev_primer_match) == 1:
                start_pos.append(rev_primer_match[0].start)
                end_pos.append(rev_primer_match[0].end)
            elif (len(primer_match) > 1 or len(rev_primer_match) > 1):
                ## Maybe instead of raising an exception here choose the match nearest to initial primer position in varvamp bed file
                raise Exception("Primer sequence found multiple times in reference sequence.")
            else:
                raise Exception("Primer sequence not found in reference sequence.")

        # check if primer size is still the same
        seq_size = list(map(lambda a, b: a-b, end_pos, start_pos))
        if seq_size == list(primer_df['size']):
            pass
        else:
            raise Exception("The primer size before and after the update are not identical")
    except:
        raise Exception("update_position not possible")
    return start_pos, end_pos


def main(argv):
    
    if len(argv) != 4:
        print('Please add required input files \n e.g.: python3 correct_primer_positions.py [path/to/Reference/.fasta] [path/to/primer/.tsv] [path/to/primer/.bed]')
        return
    else:  
        try:
			#  read reference fasta file
            if os.path.exists(argv[1]):
                record_dict = SeqIO.to_dict(SeqIO.parse(argv[1],"fasta"))
            else:
                print("FASTA file not available.")

			# read the primer tsv and bed files
            if os.path.exists(argv[2]):
                primer_df = pd.read_csv(argv[2], sep="\t")
            else:
                print("TSV file not available.")
            if os.path.exists(argv[3]):
                primer_bed_df = pd.read_csv(argv[3], sep="\t", header=None)
            else:
                print("BED file not available.")

			# update positions in bed file
            s_pos, e_pos = update_positions(primer_df, record_dict)

            primer_bed_df.iloc[:,1] = s_pos
            primer_bed_df.iloc[:,2] = e_pos

            primer_bed_df.to_csv(os.path.dirname(argv[3])+'/'+os.path.basename(argv[3])+'.position.corrected.bed', index=False, sep='\t', header=False)
            print("Primer.bed file has been updated. See: \n"+os.path.dirname(argv[3])+'/'+os.path.basename(argv[3])+'.position.corrected.bed')
        except Exception as e:  # pragma: nocover
            return "An error has occured: " + str(e)
if __name__ == "__main__":  # pragma: nocover
    main(argv)
