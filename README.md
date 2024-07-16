# BradyrhizobiumPDS
Identifying phage defense systems in Bradyrhizobium. Simms Lab at UC Berkeley. This pipeline was written by Margot Lavitt, designed alongside Dr. Ellen Simms for internal use in the Simms Lab.

## Introduction
Identification of phage defense systems encoded in genomes of Bradyrhizobia. 

## Contact
For general questions or assistance installing or running, contact margotl@berkeley.edu
## Citations
#### Pipeline: Snakemake
Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., & Köster, J. (2021). Sustainable data analysis with Snakemake. F1000Research, 10, 33. https://doi.org/10.12688/f1000research.29032.2
#### Tool: Defense Finder
Payne, L. J., Todeschini, T. C., Wu, Y., Perry, B. J., Ronson, C. W., Fineran, P. C., Nobrega, F. L., & Jackson, S. A. (2021). Identification and classification of antiviral defence systems in bacteria and archaea with PADLOC reveals new system types. Nucleic acids research, 49(19), 10868–10878. https://doi.org/10.1093/nar/gkab883
#### Tool: PADLOC
J., Ronson, C. W., Fineran, P. C., Nobrega, F. L., & Jackson, S. A. (2021). Identification and classification of antiviral defence systems in bacteria and archaea with PADLOC reveals new system types. Nucleic acids research, 49(19), 10868–10878. https://doi.org/10.1093/nar/gkab883

## Set Up
Requires use of computing cluster and optionally GitHub.
### 1. Set Up on Cluster
1. Set up [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#full-installation) environment using [mamba](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html#mamba) (You only need to run this once)

	```mamba create -y -n simmslab_BPDS -c conda-forge -c bioconda mamba snakemake snakefmt snakedeploy git=2.34.1```
2. Make a new directory where you want to store and run this pipeline and its output
	
	```mkdir <insert directory name>```
3. Make a new directory where you want to store and run this pipeline and its output

	```cd <insert directory name>```
4. To populate directory with the BradyrhizobiumPDS pipeline

	```snakedeploy deploy-workflow https://github.com/margotl9/SimmsLab_BradyrhizobiumPDS.git . --branch main```
	
	_If you run the command `ls`, you should now see 'config' and 'workflow' directories_

### 2. Modify Config
Either edit it directly on cluster by opening an interactive session on https://ood.brc.berkeley.edu [instructions here](https://docs-research-it.berkeley.edu/services/high-performance-computing/user-guide/ood/#code-server-vs-code) or by creating a Git repository for current directory to edit locally. I personally used the latter because the interactive session constantly prompted me to log in.

**Fill in config** - _none of the following fields can be left blank!_
* ```config/config.yaml```: full paths for resources and results directories -- make sure these directories have already been created before running.
* ```config/simple/config.yaml```: `sbatch` and other job submission fields as necessary
* ```config/samples.tsv```: sample information 
	* location: RT
	* site_plant: 7b
	* fasta: path to file **fasta files must be stored in `resources/input`**



### 3. Run on cluster
1. Activate `simmslab_BPDS`
	
	```mamba activate simmslab_BPDS```
1. Dry-run to ensure proper rule chaining and input paths
	
	```snakemake --profile config/simple -p --dry-run```
3. Run pipeline, specifying an integer (Ex 4) of threads
	
	```snakemake --profile config/simple --max-threads <insert number of threads> -p```

## Output
For every input sequence, with its unique `location` + `site_plant` fields from `samples.tsv`, there should be:
#### DefenseFinder 
_See [DefenseFinder GitHub](https://github.com/mdmparis/defense-finder?tab=readme-ov-file#outputs) for in-depth descriptions of files and column headers_
| File        | Description |
| ----------- | ----------- |
| _defense_finder_analysis.csv	|	summary of systems found with descriptions from [DefenseFinder's predicted structures database](https://defensefinder.mdmlab.fr/wiki/structure)	|
| .defense_finder_systems.tsv	|	summary of systems found	|
| .defense_finder_genes.tsv		|	summary of genes found	|
| .defense_finder_hmmer.tsv		|	summary of HMM hits	|
| .prt	and .prt.idx			|	implicitly created by DefenseFinder	|

#### PADLOC
_Taken from [PADLOC's GitHub](https://github.com/padlocbio/padloc#output) output descriptions_
| File        | Description |
| ----------- | ----------- |
|	.domtblout		|	Domain table file generated by HMMER	|
|	_padloc.csv		|	PADLOC output file for identified defense systems	|
|	_padloc.gff		|	GFF annotation file for identified defense systems	|
