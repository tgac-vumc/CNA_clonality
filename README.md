# CNA_clonality

<p align="center">
  <img width="50%" height="50%" src="https://github.com/tgac-vumc/CNA_clonality/blob/main/dag.svg">
</p>

Code associated with manuscript "The added value of genome-wide copy nubers to objectively resolve clonality of multiple lung tumors with ambiguous or inconclusive mutational diagnosis"

To run this pipeline Snakemake is required.

for easy installation you need (mini)conda.

Miniconda installation from folder where you want to install miniconda:

```
cd </path/to/files/dir/>
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

follow the instructions of the installation process, give the location where you want miniconda to be installed and answer YES to add miniconda to your path.

go to the directory where the analysis need to be performed

```
cd </path/to/analysis/dir>
git clone https://github.com/tgac-vumc/CNA_clonality.git
cd CNA_clonality
```

install the snakemake environment,

```

conda env create --name snakemake --file envs/snakemake.yaml

```
activate the snakemake environment
```
conda activate snakemake

```

Setup your files correctly:
	- Specify directories of input and output files in config.yaml   

navigate to CNA_clonality directory and start snakemake.

```
snakemake --use-conda

```
Useful snakemake options

-j , --cores, --jobs : Use at most N cores in parallel (default: 1). If N is omitted, the limit is set to the number of available cores.   
-n , --dryrun : Do not execute anything. but show rules which are planned to be performed.    
-k , --keep-going : Go on with independent jobs if a job fails.    
-f , --force : Force the execution of the selected target or the first rule regardless of already created output.   
-U , --until : Runs the pipeline until it reaches the specified rules or files. Only runs jobs that are dependencies of the specified rule or files, does not run sibling DAGs.   
-T , --timestamp : Add a timestamp to all logging output   
 
for all options go to http://snakemake.readthedocs.io/en/stable/executable.html#all-options
