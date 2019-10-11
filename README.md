![ONT](./Images/ONT_logo.png "Oxford Nanopore Technologies")

**************************

# 1. Introduction 

### Overview:

This pipeline manages a pore-c workflow starting from raw fastq files and converting
them to standard file formats for use by downstream tools. The steps involved are:

* Pre-processing a reference genome or draft assembly to generate auxiliary files used in downstream analyses
* Creating virtual digests of the genome
* Filtering the raw reads to remove any that might break downstream tools
* Align against a reference genome
* Processing results to filter spurious alignments, detect ligation junctions and assign fragments. The results are stored in a [parquet](http://parquet.apache.org/) table for downstream processing.
* Converting the results to the following formats:
  - [pairs format](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md)
  - [cooler format](https://mirnylab.github.io/cooler/) 
  - [hic medium format](https://github.com/aidenlab/juicer/wiki/Pre#medium-format-most-common)
  - [salsa2 bed format](https://github.com/marbl/SALSA)


# 2. Getting started

In most cases, it is best to pre-install conda before starting. All other dependencies will be installed automatically when running the pipeline for the first time. 

### Requirements:

This pipleine requires a computer running Linux (Ubuntu 16). >64Gb of memory would be recommended. The pipeline has been tested on minimal server installs of these operating systems.

Most software dependecies are managed using *conda*. To install conda, please install [miniconda3](https://conda.io/miniconda.html) and refer to installation [instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).
You will need to accept the license agreement during installation and we recommend that you allow the Conda installer to prepend its path to your .bashrc file when asked.

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Check if the conda has successfully installed

```
conda -h
```

If conda has installed correctly, you should see the follow output.
If you do not see the below output, you may need to close and reopen your terminal.

```
$ conda
usage: conda [-h] [-V] command ...

conda is a tool for managing and deploying applications, environments and packages.

Options:

positional arguments:
  command
    clean        Remove unused packages and caches.
    config       Modify configuration values in .condarc. This is modeled
                 after the git config command. Writes to the user .condarc
                 file (/Users/prughani/.condarc) by default.
    create       Create a new conda environment from a list of specified
                 packages.
..............
```

---

#### Installation:

Clone *Pore-c Snakemake* git repository from https://git.oxfordnanolabs.local/genomic_apps_projects/pore-c-snakemake
A snakemake pipeline to drive the analysis of pore-C data. This workflow follows the format of those in the [Snakemake-Workflows](https://github.com/snakemake-workflows) project.

1. Create a new github repository using this workflow [as a template](https://help.github.com/en/articles/creating-a-repository-from-a-template).
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the newly created repository to your local system, into the place where you want to perform the data analysis.


`
git clone https://github.com/nanoporetech/Pore-C-Snakemake.git
cd pore-c-snakemake
## Creates environment and the dependencies will install automatically 
conda env create
conda activate pore_c_snakemake
`

--- 

#### Installing Pore-C tools

The [Pore-C tools](https://github.com/nanoporetech/pore-c) package required to run this pipeline. They are not yet available through conda so we have to create an environment manually.

`git clone https://github.com/nanoporetech/pore-c.git
cd pore-c
conda env create
pip install -e .
`
***********


# 3. Usage

#### Configure workflow:

Configure the workflow according to your needs via editing the file `config.yaml`. At the very least you need to provide a path to a reference genome.

#### Add basecall location:

Add the locations of your fastq files to `basecall.tsv`. The run IDs must *not* contain any spaces or special characters. 

#### Execute workflow:

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

### Testing:

Test cases are in the subfolder `.test`. To run the tests use

    snakemake --use-conda  -d .test all salsa_bed  juicebox


They are automtically executed via continuous integration with Travis CI and requires [git lfs](https://github.com/git-lfs/git-lfs/wiki/Installation) to be installed.

# 4. Output files
Once the pipeline has run successfully you should expect the following files in the output directoy:

*  **`align_table/`:**
    *  `*.at.catalog.yaml` - an [intake](https://intake.readthedocs.io/en/latest/) catalog containing metadata about the alignment table.
    *  `*.at.alignment.parquet` - records containing all alignment information in parquet format.
    *  `*.at.alignment_summary.csv` - alignment summary table.
    *  `*.at.overlap.parquet` - table of the overlaps between alignments and fragments.
    *  `*.at.read.parquet` - per-read alignment statistics (contacts per-read etc).
    *  `*.at.read_summary.csv` - summary stats of input reads, such as read N50.
*   **`basecall/`:** 
    *   `*.rd.catalog.yaml` - an [intake](https://intake.readthedocs.io/en/latest/) catalog containing read metadata.
    *   `*.rd.pass.fq.gz` - basecalls that have passed filtering.
    *   `*.rd.read_metadata.parquet` - a table of per-read statistics. 
    *   `*.rd.summary.csv` - a table of aggregate statistics for the reads.
*   **`mapping/`:**
    *   `*.read_sort.bam` - bam alignment file sorted by read name.
*   **`matrix/`**
    *   `*.matrix.catalog.yaml` - an [intake](https://intake.readthedocs.io/en/latest/) catalog containing metadata about the aggregate matrix.
    *   `*.matrix.coo.csv.gz` - aggregate read counts in the format 'bin1_id,bin2_id,count' - suitable for use with `cooler load` the bin width for this set by the `*base*` matrix resolution in the config file.
    *   `*.matrix.cool` - the aggregate contact counts in [cooler](https://mirnylab.github.io/cooler/) format
    *   `*.matrix.counts.mcool` - a multi-resolution cool file.
*   **`pairs/`:**
    *   `*.pairs.pairs.gz` - contains fragment position and fragment pairs in [pairs format](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md).
*   **`refgenome/`:**
    *   `*.rg.metadata.csv` - chromosome metadata in csv format.
    *   `*.rg.chromsizes`- reference genome chromosome lengths
    *   `*.rg.fa.gz` - reference genome compressed with bgzip
    *   `*.rg.fa.gz.fai` - samtools indexed reference genome
    *   `*..rg.fa.gz.bwt` - bwa index reference genome
*   **`virtual_digest/`:**
    *   `*.vd.fragments.parquet` - A table containing the intervals generated by the virtual digest.
    *   `*.vd.digest_stats.csv` - virtual digest aggregate statistics
*   **`assembly/`:**
    *   `*.salsa2.bed` - *optional* a bed file compatible with the [salsa2](https://github.com/marbl/SALSA) scaffolding tool.
*   **`juicebox/`:**
    *   `*.hicRef` - *optional* a [restriction site format](https://github.com/aidenlab/juicer/wiki/Pre#restriction-site-file-format) file.
    *   `*.hic.txt` - *optional* a [hic medium format](https://github.com/aidenlab/juicer/wiki/Pre#medium-format-most-common) file of pairwise contacts.

# 5. Help

#### Licence and Copyright:
© 2019 Oxford Nanopore Technologies Ltd.

Bioinformatics-Tutorials is distributed by Oxford Nanopore Technologies under the terms of the MPL-2.0 license.


