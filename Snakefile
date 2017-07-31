# imports
import os
import glob
import yaml
import datetime

# yaml representer for dumping config
from yaml.representer import Representer
import collections
yaml.add_representer(collections.defaultdict, Representer.represent_dict)


with open(config["BARCODES"], "r") as bc_file:
    BARCODE_IDS = [line.strip()[1:] for line in bc_file if line.startswith(">")]

with open(config["GENE"], "r") as infile:
    GENE = yaml.safe_load(infile)

GENE_NAME = GENE["name"]
CHROMOSOME = GENE["chromosome"]

try:
    with open(config["EXPERIMENT"], "r") as infile:
        EXPERIMENT = yaml.safe_load(infile)

    START = EXPERIMENT["targets"][0]["primers"][0]["forward"]["start"]
    END = EXPERIMENT["targets"][0]["primers"][0]["reverse"]["end"]
except (KeyError, IOError):
    START = GENE["coordinates"]["start"]
    END = GENE["coordinates"]["end"]


# handlers for workflow exit status
onsuccess:
    print("Variant calling workflow completed successfully")
    config_file = "config.{}.yaml".format("{:%Y-%m-%d_%H:%M:%S}".format(datetime.datetime.now()))
    with open(config_file, "w") as outfile:
        print(yaml.dump(config, default_flow_style=False), file=outfile)

onerror:
    print("Error encountered while executing workflow")
    shell("cat {log}")


# main workflow
localrules:
    all

rule all:
    input:
        expand("variants/{barcodes}.json", barcodes=BARCODE_IDS)


rule reference:
    input:
        reference = config["GENOME"]
    output:
        "reference/{}.fasta".format(GENE_NAME)
    params:
        chrom = CHROMOSOME["name"],
        start = START,
        end =  END
    conda:
        "envs/make_reference.yaml"
    script:
        "scripts/make_reference.py"


rule variants:
    input:
        reference = "reference/{}.fasta".format(GENE_NAME),
        alleles = config["ALLELE_FASTA_FOLDER"] + "/{barcode}.fasta"
    output:
        variants = "variants/{barcode}.json",
        alignments = "variants/{barcode}.aln"
    params:
        offset = START
    conda:
        "envs/call_variants.yaml"
    script:
        "scripts/call_variants.py"

