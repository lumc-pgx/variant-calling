include: "helper.snake"

PARAMS = VariantCallingHelper(config, "Variant Calling")

onsuccess: PARAMS.onsuccess()
onerror: PARAMS.onerror()


# main workflow
localrules:
    all, link_sources


rule all:
    input:
        PARAMS.outputs


include: "rules/gene_reference.snake"
include: "rules/call_variants.snake"
include: "rules/link_sources.snake"


