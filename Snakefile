include: "helper.snake"

PARAMS = VariantCalling(config, "Variant Calling")

onsuccess: PARAMS.onsuccess()
onerror: PARAMS.onerror()

# main workflow
localrules:
    all

rule all:
    input:
        expand("variants/{barcodes}.json", barcodes=PARAMS.BARCODE_IDS)


rule reference:
    input:
        reference = config["GENOME"]
    output:
        "reference/{}.fasta".format(PARAMS.GENE_NAME)
    params:
        chrom = PARAMS.CHROMOSOME["name"],
        start = PARAMS.START,
        end = PARAMS.END
    conda:
        "envs/make_reference.yaml"
    script:
        "scripts/make_reference.py"


rule variants:
    input:
        reference = "reference/{}.fasta".format(PARAMS.GENE_NAME),
        alleles = config["ALLELE_FASTA_FOLDER"] + "/{barcode}.fasta"
    output:
        variants = "variants/{barcode}.json",
        alignments = "variants/{barcode}.aln"
    params:
        offset = PARAMS.START
    conda:
        "envs/call_variants.yaml"
    script:
        "scripts/call_variants.py"

