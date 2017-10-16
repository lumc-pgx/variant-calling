"""
Use PyBedTools to generate a fasta file containing the sequence for a single region
"""
from Bio import SeqIO
import pybedtools
import locus_processing
import yaml


locus = locus_processing.load_locus_yaml(snakemake.input.locus)

try:
    with open(snakemake.config["EXPERIMENT"], "r") as infile:
        experiment = yaml.safe_load(infile)
        start_pos = experiment["targets"][0]["primers"][0]["forward"]["start"]
        end_pos = experiment["targets"][0]["primers"][0]["reverse"]["end"]
except (KeyError, IOError):
    start_pos = locus.coordinates.start
    end_pos = locus.coordinates.end


# create a bed tool for the required region
bed_tool = pybedtools.BedTool([(locus.chromosome.name, start_pos - 1, end_pos)])

# associate the bedtool with the reference genome fasta
bed_tool = bed_tool.sequence(fi=snakemake.input.genome)

# get the sequence and save it
with open(snakemake.output[0], "w") as outfile:
    sequence = SeqIO.read(bed_tool.seqfn, "fasta")
    SeqIO.write(sequence, outfile, "fasta")
