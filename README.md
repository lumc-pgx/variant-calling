# Variant calling workflow

**Calling of normalized hgvs variants for amplicon sequences**

The variant calling workflow module performs the following operations:  
- Creation of a reference fasta sequence for the targetted region.
- Alignment of amplicon sequences to the reference sequence.
- 'Calling' of normalized hgvs-like variants.

The pipeline outputs two files per barcode:
- A {barcode}.json file which contains the found variants per allele sequence
- A {barcode}.aln file  which contains the alignments used for the variant calling

## Requirements
- [Conda/Miniconda](https://conda.io/miniconda.html)  

## Installation
- Clone the repository
  - `git clone https://wgallard@git.lumc.nl/PharmacogenomicsPipe/variant_calling.git`

- Change to the variant_calling directory
  - `cd variant_calling`

- Create a conda environment for running the pipeline
  - `conda env create -n variant_calling -f environment.yaml`

- In order to use the pipeline on the cluster, update your .profile to use the drmaa library:
  - `echo "export DRMAA_LIBRARY_PATH=libdrmaa.so.1.0" >> ~/.profile`
  - `source ~/.profile`

## Configuration
Pipeline configuration settings can be altered by editing [config.yaml](config.yaml).  

## Execution
- Activate the conda environment
  - `source activate variant_calling`
- For parallel execution on the cluster
  - `./run_cluster.sh`
- To specify that the pipeline should write output to a location other than the default:
  - `./run_cluster.sh -d path/to/output/directory`

