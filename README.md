MetaDBGWAS
==========

* [Motivation](#motivation)
* [Requirements](#requirements)
* [Installation](#installation)
* [Usage](#usage)
* [Output](#output)
* [License](#license)

# Motivation

This tool expands the work of [DBGWAS](https://gitlab.com/leoisl/dbgwas) [(Jaillard et al, 2018)](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007758) and brings it to the metagenomic scale.

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
    git clone --recursive https://github.com/Louis-MG/MetaDBGWAS.git
    ```
2. Add the C++ libraries :
   ```shell
    sudo apt install libgatbcore-dev libhdf5-dev libboost-dev libpstreams-dev
    ```
3. Complete the installation :
    ```shell
    cd MetaDBGWAS
    sed -i "51i#include <limits>" ./REINDEER/blight/robin_hood.h #temporary fix for REINDEER compilation
    sh install.sh
    ```
4. Add the R libraries :
    ```R
    install.packages(c('ape', 'phangorn'))
    install.packages("https://raw.githubusercontent.com/sgearle/bugwas/master/build/bugwas_1.0.tar.gz", repos=NULL, type="source")
    ```

# Usage

```
	* General
NOTE: path should be absolute.
--files <path> path to the directory containing the read files.
--output <path> path to the output folder. Default set to ./ .
--threads <int> number of threads to use. !! Default set to 4 !!
--verbose <int> level of verbosity. Default to 1, 0-1. 0 is equivalent to --quiet.
--clean removes files from output directory if not empty.

        * Lighter
--K <int> <int> kmer length and genome size (in base). Recommended is 17 X.
        or
--k <int> <int> <float> kmer length and genome size (in base), alpha (probability of sampling a kmer). Recommended is 17 X X.

	* Bcalm2
--kmer <int> kmer length used for unitigs build. Default to 31.
--abundance-min <int> minimum number of occurence of a kmer to keep it in the union DBG. Default to 5, highly recommended to change it.

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
docker run -v 'path/to/input/folder:/input' metadbgwas
```

# Singularity

You can also run the docker image with singularity:

```bash
docker pull 007ptar007/metadbgwas
singularity run -H /path/to/input metadbgwas_latest.sif --files ./input --strains /input/strains --threads 40 --output ./output --K 17 X
```

# Output

User can find in the output folder :
- the corrected fasta files.
- unitigs folder with output from bcalm (intermediate .h5 file of gatb and the unitigs.fa)
- matrix folder with the matrix of absence/presence of kmers in the unitigs. Ouptut from Reindeer.
- fof.txt : file of files for input of bcalm
- fof_unitigs.txt : file of file for input of Reindeer.
- step1 folder which contains input files for gemma and bugwas.
- step2 folder folder which contains 
- step3 contains visulatisation files.

# How to reference :

Please cite this tool as :

Metadbgwas, L-M. Gueguen, 2022.


# Issues :

You can post issues in the issue section of the github repository. You can also email me at lm<dot>gueguen<at>orange<dot>fr . I will do my best to resolve them.


# License

The work is available under the zlib license.
