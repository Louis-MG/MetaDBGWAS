#!/bin/bash

# Louis-Mael Gueguen lm.gueguen@orange.fr


#donner les valeurs par defaut aux parametres


#General
files=
output="./"
threads=4
verbose=0
clean=false

#Lighter
genome_size=
alpha=0
kmer_l=17

#bcalm2

kmer=31

#Reindeer

#dbgwas
strain=''
newick=''

#strains
#k
#newick



#miscealenous

Version()
{
	echo "\nMetadDBGWAS 0.1\n"
}

License()
{
	echo "Copyright (C) 2022 Louis-Maël Gueguen

This software is provided 'as-is', without any express or implied
warranty.  In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.

Louis-Maël Gueguen lm.gueguen@orange.fr\n"
}

Help()
{
   # Display Help
   echo "
	* General
--files <path> path to files.
--output <path> path to the output folder, current directory by default.
--threads <int> number of threads to use. Default to 4.
--verbose <int> level of verbosity. Default to 1, 1-3.

	* Lighter
--K <kmer length (int)> <genome size (base, int)>
	or
--k <kmer length (int)> <genome size (in base, int)> <alpha (float)>

	* DBGWAS
--strains A text file describing the strains containing 3 columns: 1) ID of the strain; 2) Phenotype (a real number or NA); 3) Path to a multi-fasta file containing the sequences of the strain. This file needs a header. Check the sample_example folder or https://gitlab.com/leoisl/dbgwas/raw/master/sample_example/strains for an example.
--newick Optional path to a newick tree file. If (and only if) a newick tree file is provided, the lineage effect analysis is computed and PCs figures are generated.

	* Miscellaneous
--license prints the license text in standard output.
--help displays help.\n"
}

#Parameters parsing

while [ $# > 1 ]
do
	case $1 in
	-f | --files) files="$2"
	shift 2;;
	-o | --output) output="$2"
	shift 2;;
	-c | --clean) clean=true
	shift;;
	-t | --threads) threads="$2"
	shift 2;;
	--K) kmer_l="$2" genome_size="$3"
	shift 3;;
	--k) kmer_l="$2" genome_size="$3" alpha="$4"
	shift 4;;
	-k | --kmer) kmer="#2"
	shift 2;;
	--version) Version; exit;;
	--license) License; exit;;
	-v | --verbose) verbose="$2"
	shift 2;;
	-h | --help) Help; exit;;
	-* | --*) echo "Unknown option"; exit;;
	*) break;
	esac
done


#creates output dir if it doesnt exist yet

if [ -d $output ]
then
	if [ "$(ls -A $output)" ] && [ $clean != true ]
	then
    		echo "$output is not Empty."
		exit 0
	else
		rm -r $output
	fi
else
	mkdir $output
fi


#Lighter
#else tells user that file is not found

if [ $verbose -ge 1 ]
then
	echo "Starting kmer corrections with Lighter ..."
fi

if [ -f $files ]
then
	echo $files > list_files
	if [ $alpha -gt 0 ]
	then
		for i in $files
		do
			./Lighter/lighter -r ${i} -od $output -t $threads -discard -k $kmer_l $genome_size $alpha
		done
	else
		for i in $files
		do
			./Lighter/lighter -r ${i} -od $output -t $threads -discard -K $kmer_l $genome_size
		done
	fi
elif [ -d $files  ]
then
	find $files -type f > list_files
	for i in $files/*
	do
		./Lighter/lighter -r ${i} -od $output -t $threads -discard -k $kmer_l $genome_size $alpha
	done
else
	echo "File not found, verify the path."
fi

#Bcalm 2

./bcalm/build/bcalm -in ./list_files -kmer-size $kmer -nb-cores $threads -out-dir $output
