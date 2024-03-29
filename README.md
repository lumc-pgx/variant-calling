# Variant calling workflow

**Calling of normalized hgvs variants for amplicon sequences**

The variant calling workflow module performs the following operations:  
- Creation of a reference fasta sequence for the targetted region.
- Alignment of amplicon sequences to the reference sequence.
- 'Calling' of normalized hgvs-like variants.

The pipeline outputs two files per barcode:
- A {barcode}.json file which contains the found variants per allele sequence
- A {barcode}.aln file  which contains the alignments used for the variant calling

![rule graph](static/rulegraph.png)
   
## Requirements
- [Conda/Miniconda](https://conda.io/miniconda.html)  

## Installation
- Clone the repository
  - `git clone https://github.com/lumc-pgx/variant-calling.git`

- Change to the variant_calling directory
  - `cd variant_calling`

- Create a conda environment for running the pipeline
  - `conda env create -n variant_calling -f environment.yaml`

## Configuration
The pipeline configuration settings are specified in [config.yaml](config.yaml).  
Edit the configfile with run-specific paths and settings.  

## Execution
- Activate the conda environment
  - `source activate variant_calling`
- Run the pipeline using Snakemake
         
