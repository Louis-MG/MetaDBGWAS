Metadbgwas
==========

* [Motivation](#motivation)
* [Requirements](#requirements)
* [Installation](#installation)
* [Usage](#usage)
* [Output](#output)
* [License](#license)

# Motivation

This tool expands the work of [DBGWAS](https://gitlab.com/leoisl/dbgwas) [(Jaillard et al, 2018)](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007758) and brings it to the metagenomic scale. It finds variants significantly associated with a given phenotype of interest, and output its findings in a web page. Thta web page shows the selected unitigs with their surounding graph component, and a user-provided annotation can be used.

# Overview

![Schematic of metadbgwas](/figures/figure_pipeline.png?raw=true "Pipeline Overview")

# Output

Here is an example of the output of Metadbgwas: significant components are shown in preview with annotation if provided, and a click on them will take you to the page with full information and interactive graph.

![Output example](/figures/suppl_comp_table.png?raw=true "Output example of Metadbgwas")

# Requirements

* [GCC](https://gcc.gnu.org/) >= 9.4
* [CMAKE](https://cmake.org/) > 3.10.0
* [zlib](https://en.wikipedia.org/wiki/Zlib)
* [pthreads](https://en.wikipedia.org/wiki/Pthreads)
* [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) The blast suite in your path
* [R](https://www.r-project.org/) >= 3.2.0 
* [Boost](https://www.boost.org/)

# Installation

1. Use the following command to download the repository :
    ```shell
    git clone --recursive https://github.com/Louis-MG/Metadbgwas.git
    ```
3. Complete the installation :
    ```shell
    cd Metadbgwas
    sh install.sh
    ```

# Usage

```
	* General
NOTE: path should be absolute.
--files <path> path to the directory containing the read files.
--output <path> path to the output folder. Default set to ./ .
--threads <int> number of threads to use. !! Default set to 4 !!
--verbose <int> level of verbosity. Default to 1, 0-1. 0 is equivalent to --quiet.
--clean removes intermediary files to save space if you are worried about your storage.
--skip1 skips the Lighter correction step. Corrected files are supposed to be in the output folder.
--skip2 skips the Lighter and Bcalm2 steps. Corrected files and unitigs folder are supposed to be in the output folder.
--skip3 skips the Lighter, Bcalm2 and REINDEER steps. Corrected files, unitigs and matrix folder are supposed to be in the output folder.

        * Lighter
NOTE: if your datset contains different bacterial genomes with very different size, it is better to choose --k option and provide the pick-rate (noted alpha).
--K <int> <int> kmer length and approximate genome size (in base). Recommended is 17 G.
        or
--k <int> <int> <float> kmer length and genome size (in base), alpha (probability of sampling a kmer). Recommended is 17 G alpha. alpha is best chosen at 70/coverage.

	* Bcalm2
--kmer <int> kmer length used for unitigs build. Default to 31.
--abundance-min <int> minimum number of occurence of a kmer to keep it in the union DBG. Default to 5, highly recommended to change to the 2.5% quantile of the Poisson law with lambda = coverage.

	* Reindeer
Reindeer uses kmer, threads, and output parameters. No others need to be specified.

        * DBGWAS
--strains A text file describing the strains containing 3 columns: 1) ID of the strain; 2) Phenotype (a real number or NA); 3) Path to a multi-fasta file containing the sequences of the strain. This fil>
--newick Optional path to a newick tree file. If (and only if) a newick tree file is provided, the lineage effect analysis is computed and PCs figures are generated.
--nc-db A list of Fasta files separated by comma containing annotations in a nucleotide alphabet format (e.g.: -nc-db path/to/file_1.fa,path/to/file_2.fa,etc). You can customize these files to work better with DBGWAS (see https://gitlab.com/leoisl/dbgwas/tree/master#customizing-annotation-databases).
--pt-db A list of Fasta files separated by comma containing annotations in a protein alphabet format (e.g.: -pt-db path/to/file_1.fa,path/to/file_2.fa,etc). You can customize these files to work better with DBGWAS (see https://gitlab.com/leoisl/dbgwas/tree/master#customizing-annotation-databases).
--threshold maximum value for which phenotype will be considered to be 0.

        * Miscellaneous
--license prints the license text in standard output.
--help displays help.
```

# Exemple

```bash
bash metadbgwas.sh --files /input --output /output --K 17 6000000
```

# Docker

An image is hosted on Docker hub. You can also build it localy using the dockerfile located in /docker. 

```bash
docker pull 007ptar007/metadbgwas:latest
docker run -v 'path/to/input/folder:/input' metadbgwas --files ./input --strains /input/strains --threads 40 --output ./output --K 17 G
```

# Singularity

You can also run the docker image with singularity:

```bash
singularity pull docker://007ptar007/metadbgwas
singularity run -H /path/to/input metadbgwas_latest.sif --files ./input --strains /input/strains --threads 40 --output ./output --K 17 G
```

# Output

User can find in the output folder :
- the corrected fasta files.
- unitigs folder with bcalm2 output, sample-wise and dataset-wise. 
- matrix folder with the matrix of absence/presence of kmers in the unitigs. Ouptut from Reindeer.
- step1, step2 and step3 that contains internal files of the modified DBGWAS
- visualisation contains visulatisation files.
- command_line.txt with the paremeters used for the execution

# How to reference :

Please cite this tool as :

Metadbgwas, Louis-Mael Gueguen, 2022.


# Issues :

You can post issues in the issue section of the github repository. You can also email me at `lm<dot>gueguen<at>orange<dot>fr` . I will do my best to resolve them.


# License

The work is available under the zlib license.
