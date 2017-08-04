# Variant calling workflow

**Calling of normalized hgvs variants for amplicon sequences**

The variant calling workflow module performs the following operations:  
- Creation of a reference fasta sequence for the targetted region.
- Alignment of amplicon sequences to the reference sequence.
- 'Calling' of normalized hgvs-like variants.

The pipeline outputs two files per barcode:
- A {barcode}.json file which contains the found variants per allele sequence
- A {barcode}.aln file  which contains the alignments used for the variant calling

```plantuml
digraph snakemake_dag {
    graph[bgcolor=white, margin=0];
    node[shape=box, style=rounded, fontname=sans,                 fontsize=10, penwidth=2];
    edge[penwidth=2, color=grey];
	0[label = "variants", color = "0.22 0.6 0.85", style="rounded"];
	1[label = "all", color = "0.44 0.6 0.85", style="rounded"];
	2[label = "reference", color = "0.00 0.6 0.85", style="rounded"];
	2 -> 0
	0 -> 1
}
```
   
## Requirements
- [Conda/Miniconda](https://conda.io/miniconda.html)  

## Installation
- Clone the repository
  - `git clone https://git.lumc.nl/PharmacogenomicsPipe/variant_calling.git`

- Change to the variant_calling directory
  - `cd variant_calling`

- Create a conda environment for running the pipeline
  - `conda env create -n variant_calling -f environment.yaml`

- In order to use the pipeline on the cluster, update your .profile to use the drmaa library:
  - `echo "export DRMAA_LIBRARY_PATH=libdrmaa.so.1.0" >> ~/.profile`
  - `source ~/.profile`

## Configuration
The pipeline configuration settings are specified in [config.yaml](config.yaml).  
Edit the configfile with run-specific paths and settings.  

## Execution
- Activate the conda environment
  - `source activate variant_calling`
- For parallel execution on the cluster
  - `pipe-runner --configfile config.yaml`
- To specify that the pipeline should write output to a location other than the default:
  - `pipe-runner --configfile config.yaml --directory path/to/output/directory`
         
