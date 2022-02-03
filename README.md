MetaDBGWAS
==========

* [Motivation](#motivation)
* [Requirements](#requirements)
* [Installation](#installation)
* [Usage](#usage)
* [Output](#output)
* [License](#license)

# Motivation

This tool expands thee work of [DBGWAS](https://gitlab.com/leoisl/dbgwas) [(Jaillard et al, 2018)](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007758) and bring to the metagenomic scale.

# Requirements

* [GCC](https://gcc.gnu.org/) >= 4.8
* [CMAKE](https://cmake.org/) > 3.10.0
* [zlib](https://en.wikipedia.org/wiki/Zlib)
* [pthreads](https://en.wikipedia.org/wiki/Pthreads)
* [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) The blast suite in your path
* [R](https://www.r-project.org/) >= 3.2.0 

# Installation

1. Use the following command to download the repository :
	```bash
	git clone --recursive https://github.com/Louis-MG/MetaDBGWAS.git
	```
2. Complete the installation :
	```bash
	cd MetaDBGWAS
	sh install.sh
	```

# Usage

```
 * General
--files <path> path to files
--output <path> path to the output folder, current directory by default
--threads <int> number of threads to use
--verbose <int> level of verbosity. Default to 1, 1-3
--clean removes files from output directory if not empty

        * Lighter
--K <kmer length (int)> <genome size (base, int)>
        or
--k <kmer length (int)> <genome size (in base, int)> <alpha (float)>

	* bcalm
--kmer <kmer length (int)> kmer length used for unitigs build.

        * DBGWAS
--strains A text file describing the strains containing 3 columns: 1) ID of the strain; 2) Phenotype (a real number or NA); 3) Path to a multi-fasta file containing the sequences of the strain. This fil>
--newick Optional path to a newick tree file. If (and only if) a newick tree file is provided, the lineage effect analysis is computed and PCs figures are generated.  [default '']

        * Miscellaneous
--license prints the license text in standard output
--help displays help\n"

```

# Output

# License

The work is available under the zlib license.
