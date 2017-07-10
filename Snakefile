# load config file
configfile: srcdir("config.yaml")

# imports
import os
import glob
import yaml

INPUT_FILES = glob.glob(config["LAA_DATA_PATH"] + "/*.fa") + glob.glob(config["LAA_DATA_PATH"] + "/*.fasta")
BARCODE_IDS = [".".join(os.path.basename(f).split(".")[:-1]) for f in INPUT_FILES]

with open(config["GENE"], "r") as infile:
    GENE = yaml.safe_load(infile)

with open(config["EXPERIMENT"], "r") as infile:
    EXPERIMENT = yaml.safe_load(infile)

GENE_NAME = GENE["name"]
CHROMOSOME = GENE["chromosome"]["name"]
START = EXPERIMENT["targets"][0]["primers"][0]["forward"]["start"]
END = EXPERIMENT["targets"][0]["primers"][0]["reverse"]["end"]

# handlers for workflow exit status
onsuccess:
    print("Variant calling workflow completed successfully")
onerror:
    print("Error encountered while executing workflow")
    shell("cat {log}")


# main workflow
localrules:
    all

rule all:
    input:
        #expand("variants/{barcodes}.txt", barcodes=BARCODE_IDS)
        "reference/CYP2D6.fasta"


rule reference:
    input:
        reference = config["GENOME"]
    output:
        "reference/{}.fasta".format(GENE_NAME)
    params:
        chrom = CHROMOSOME,
        start = START,
        end =  END
    script:
        "scripts/make_reference.py"

