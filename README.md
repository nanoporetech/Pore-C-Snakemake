# Snakemake workflow for pore-C

A snakemake pipeline to drive the analysis of pore-C data. This workflow follows the format of those in the [Snakemake-Workflows](https://github.com/snakemake-workflows) project.


#### Step 1: Obtain a copy of this workflow

1. Create a new github repository using this workflow [as a template](https://help.github.com/en/articles/creating-a-repository-from-a-template).
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the newly created repository to your local system, into the place where you want to perform the data analysis.

#### Step 2: Create the execution environment

Create the conda environment with snakemake and other packages required to run the pipeline.
     
     conda env create
     source activate pore_c_snakemake

#### Step 3: Configure workflow

Configure the workflow according to your needs via editing the file `config.yaml`.

#### Step 3: Execute workflow

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --use-conda --cluster qsub --jobs 100

or

    snakemake --use-conda --drmaa --jobs 100

If you not only want to fix the software stack but also the underlying OS, use

    snakemake --use-conda --use-singularity

in combination with any of the modes above.
See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

## Testing

Test cases are in the subfolder `.test`. They are automtically executed via continuous integration with Travis CI.



