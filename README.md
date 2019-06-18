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
	graph [bb="0,0,180,180",
		bgcolor=white,
		margin=0
	];
	node [fontname=sans,
		fontsize=10,
		label="\N",
		penwidth=2,
		shape=box,
		style=rounded
	];
	edge [color=grey,
		penwidth=2
	];
	0	 [color="0.50 0.6 0.85",
		height=0.5,
		label=variants,
		pos="95,90",
		width=0.75694];
	3	 [color="0.33 0.6 0.85",
		height=0.5,
		label=all,
		pos="95,18",
		width=0.75];
	0 -> 3	 [pos="e,95,36.104 95,71.697 95,63.983 95,54.712 95,46.112"];
	1	 [color="0.00 0.6 0.85",
		height=0.5,
		label=gene_reference,
		pos="46,162",
		width=1.2778];
	1 -> 0	 [pos="e,83.027,108.1 58.112,143.7 63.868,135.47 70.862,125.48 77.206,116.42"];
	2	 [color="0.17 0.6 0.85",
		height=0.5,
		label=link_source,
		pos="145,162",
		width=0.97222];
	2 -> 0	 [pos="e,107.22,108.1 132.64,143.7 126.77,135.47 119.63,125.48 113.16,116.42"];
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
         
