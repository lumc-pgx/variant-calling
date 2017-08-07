include: "rules/gene_reference.snake"
include: "rules/call_variants.snake"

include: "helper.snake"

PARAMS = VariantCalling(config, "Variant Calling")

onsuccess: PARAMS.onsuccess()
onerror: PARAMS.onerror()


# main workflow
localrules:
    all


rule all:
    input:
        PARAMS.outputs