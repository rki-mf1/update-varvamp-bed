# Update primer BED file from VarVAMP

A small Python script to adjust the primer positions in a BED amplicon design file produced by the [VarVAMP pipeline](https://github.com/jonas-fuchs/varVAMP).

![](https://img.shields.io/badge/python-3.12.1-brightgreen)
![](https://img.shields.io/badge/uses-conda-yellow.svg)

## Objective

[VarVAMP](https://github.com/jonas-fuchs/varVAMP) is a pipeline to design amplicon primer schemes based on multiple sequence alingments (MSA) of virus genomes. However, building an MSA can introduce gaps in the sequences and the resulting primer positions will be based on the alginment. 

If the primer scheme is now used for sequencing and bioinformatics analysis, it is important to remove the artificial primer sequences. To do that, most tools need a reference sequence and the primer BED file with the start and stop coordinates of the designed primers as input. 

Now, these primer positions need to be adjusted to match the reference sequenced used in the analysis.

## Install the script

```bash
git clone https://github.com/rki-mf1/update-varvamp-bed.git
conda env create --file=update-varvamp-bed/py_env_bed_pos_correction.yml
# or use mamba instead of conda
conda activate py_bed
```

## Run the script

The reference file in FASTA format and the primer TSV and BED files from [VarVAMP](https://github.com/jonas-fuchs/varVAMP) are taken as input and the updated BED file is obtained as output. Briefly summarized, the primer sequences from the TSV file are searched for in the reference, with max. mismatch of 2, and then the match positions are taken as new primer positions.

```bash
python3 correct_primer_positions.py [path/to/Reference/.fasta] [path/to/primer/.tsv] [path/to/primer/.bed]
```
