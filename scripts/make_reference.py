"""
Use PyBedTools to generate a fasta file containing the sequence for a single region
"""
from Bio import SeqIO
import pybedtools

# create a bed tool for the required region
bed_tool = pybedtools.BedTool([(snakemake.params.chrom, snakemake.params.start, snakemake.params.end)])

# associate the bedtool with the reference genome fasta
bed_tool = bed_tool.sequence(fi=snakemake.input[0])

# get the sequence and save it
with open(snakemake.output[0], "w") as outfile:
    sequence = SeqIO.read(bed_tool.seqfn, "fasta")
    SeqIO.write(sequence, outfile, "fasta")
