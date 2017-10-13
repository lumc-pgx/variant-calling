"""
Call variants for allele sequences
"""
from Bio import SeqIO
from variant_tools import alignment, variant_calling
from collections import OrderedDict
import json
import locus_processing
import yaml


locus = locus_processing.load_locus_yaml(snakemake.input.gene)

try:
    with open(snakemake.config["EXPERIMENT"], "r") as infile:
        experiment = yaml.safe_load(infile)
        start_pos = experiment["targets"][0]["primers"][0]["forward"]["start"]
except (KeyError, IOError):
    start_pos = locus.coordinates.start

# load the reference sequence
ref_seq = str(SeqIO.read(snakemake.input.reference, "fasta").seq)

# process the alleles
with open(snakemake.input.alleles, "r") as alleles, \
     open(snakemake.output.variants, "w") as varfile, \
     open(snakemake.output.alignments, "w") as alnfile:

    # somewhere to put the results
    result_list = []
    
    for allele in SeqIO.parse(alleles, "fasta"):
        # align the allele to the ref
        allele_seq = str(allele.seq)
        aln = alignment.align(ref_seq, allele_seq)
        
        # dump the alignment to file
        print(allele.id, file=alnfile)
        print(alignment.pretty_alignment(aln[0], aln[1], width=100, match="="), file=alnfile)
        print("", file=alnfile)

        # calculate similarity
        distance = min(alignment.get_strand_dists(allele_seq, ref_seq))
        identity = 1 - distance / ((len(allele_seq) + len(ref_seq)) / 2)
       
        # create result dict and add it to result list
        # use an ordered dict so that field order is consistent
        result = OrderedDict()
        result["sequence_id"] = allele.id
        result["identity"] = identity
        result["variants"] = [str(v) for v in variant_calling.call_variants(aln, start_pos - 1)]
        result_list.append(result)

    # dump results to json
    print(json.dumps(result_list, indent=4, separators=(',', ': ')), file=varfile)
