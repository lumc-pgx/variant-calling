"""
Call variants for allele sequences
"""
from Bio import SeqIO
from variant_tools import alignment, variant_calling, str_search
from collections import OrderedDict
import json
import locus_processing
import yaml
from snakemake import WorkflowError

locus = locus_processing.load_locus_yaml(snakemake.input.locus)

try:
    with open(snakemake.config["EXPERIMENT"], "r") as infile:
        experiment = yaml.safe_load(infile)
        start_pos = experiment["targets"][0]["primers"][0]["forward"]["start"]
except (KeyError, IOError):
    start_pos = locus.coordinates.start


def find_and_mask_repeats(query_seq, ref_seq):
    locus_strs = [str_search.parse_hgvs(s.g_notation) for s in locus.snps if s.tags is not None and 'short_tandem_repeat' in s.tags]
    strs = set()
    masked_query = query_seq
    for repeat in set((s["start"], s["unit"]) for s in locus_strs):
        
        # currently, we can only handle str nomenclature that INCLUDES the repeat motif
        if repeat[1] == "":
            raise WorkflowError("Cannot parse str notation without reference string")
        
        # find all occurrences of the specified repeat
        found = str_search.str_search_aligned(repeat[1], query_seq)
        
        # identify which of the found repeats correspond to the query repeat
        for f in found:
            ref_start = len([_ for _ in ref_seq[:f["start"]] if _ != "-"]) + start_pos
            if ref_start == repeat[0]:
                strs.add(
                    "{}_{}{}[{}]".format(
                        ref_start, 
                        ref_start + len(f["canonical_unit"]) - 1,
                        f["canonical_unit"],
                        f["num_units"]
                    )
                )
                
                # mask this repeat region in the query sequence
                masked_query = masked_query[:f["start"]] + "$" * (f["end"] + 1 - f["start"]) + masked_query[f["end"] + 1:]
    
    return list(strs), masked_query


def variant_dict(variant):
    """
    Create a dict from a variant
    """
    return {
        "g_notation": str(variant),
        "tags": []
    }
    
  
# process the alleles
with open(snakemake.input.alleles, "r") as alleles, \
     open(snakemake.output.variants, "w") as varfile, \
     open(snakemake.output.alignments, "w") as alnfile:

    # somewhere to put the results
    result_list = []
    
    # load the reference sequence
    ref_seq = str(SeqIO.read(snakemake.input.reference, "fasta").seq)
    
    for allele in SeqIO.parse(alleles, "fasta"):
        # align the allele to the ref
        allele_seq = str(allele.seq)
        aln = alignment.align(ref_seq, allele_seq)
        
        # dump the alignment to file
        print(allele.id, file=alnfile)
        print(alignment.pretty_alignment(aln[0], aln[1], width=100, match="="), file=alnfile)
        print("", file=alnfile)
        
        # search and mask repeats in the aligned allele sequence
        strs, masked_allele = find_and_mask_repeats(aln[1], aln[0])
        aln = (aln[0], masked_allele)

        # calculate similarity
        distance = min(alignment.get_strand_dists(allele_seq, ref_seq))
        identity = 1 - distance / ((len(allele_seq) + len(ref_seq)) / 2)
       
        # create result dict and add it to result list
        # use an ordered dict so that field order is consistent
        result = OrderedDict()
        result["sequence_id"] = allele.id
        result["identity"] = identity
        result["variants"] = [
            variant_dict(v) for v in variant_calling.call_variants(aln, start_pos - 1) + strs
        ]
        result_list.append(result)

    # dump results to json
    print(json.dumps(result_list, indent=4, separators=(',', ': ')), file=varfile)
