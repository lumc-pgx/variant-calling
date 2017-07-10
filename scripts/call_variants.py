"""
Call variants for allele sequences
"""
from Bio import SeqIO
from variant_tools import alignment, variant_calling
import edlib

ref_seq = str(SeqIO.read(snakemake.input.reference, "fasta").seq)

with open(snakemake.input.alleles, "r") as alleles, \
     open(snakemake.output.variants, "w") as varfile, \
     open(snakemake.output.alignments, "w") as alnfile:

    for allele in SeqIO.parse(alleles, "fasta"):
        allele_seq = str(allele.seq)
        aln = alignment.get_alignment(ref_seq, allele_seq, semi_global=True, adjust_strand=True)
        variants = variant_calling.call_variants(aln, offset=snakemake.params.offset)
        print(allele.id, file=varfile)
        print("\n".join([str(v) for v in variants]), file=varfile)
        print(allele.id, file=alnfile)
        print(alignment.pretty_alignment(aln[0], aln[1], width=100, match="="), file=alnfile)

