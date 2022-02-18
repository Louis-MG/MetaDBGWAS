#!/bin/bash

# Louis-Mael Gueguen lm.gueguen@orange.fr


#donner les valeurs par defaut aux parametres


#General
files=
output="./"
threads=4
verbose=1
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
--files <path> path to one file or a directory containing the files.
--output <path> path to the output folder. Default set to ./ .
--threads <int> number of threads to use. Default set to 4.
--verbose <int> level of verbosity. Default to 1, 0-1. 0 is equivalent to --quiet.
--clean removes files from output directory if not empty.

        * Lighter
--K <int> <int> kmer length and genome size (in base). Recommended is 17 X.
        or
--k <int> <int> <float> kmer length and genome size (in base), alpha (probability of sampling a kmer). Recommended is 17 X X.

        * bcalm
--kmer <kmer length (int)> kmer length used for unitigs build. Default to 31.

        * Reindeer
Reindeer uses kmer, threads, and output parameters. No others need to be specified. 

        * DBGWAS
--strains A text file describing the strains containing 3 columns: 1) ID of the strain; 2) Phenotype (a real number or NA); 3) Path to a multi-fasta file containing the sequences of the strain. This fil>
--newick Optional path to a newick tree file. If (and only if) a newick tree file is provided, the lineage effect analysis is computed and PCs figures are generated.

        * Miscellaneous
--license prints the license text in standard output.
--help displays help.

        * Exemple
bash metadbgwas.sh --files /test/ --output ./output --threads 4 --verbose 1 --K 17 6000000\n
	"
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

# if verbose is set to 0 : silenceing of the commands (equivaluent to --quiet)
if [ $verbose -eq 0 ]
then
	exec 1>&1 &>/dev/null
fi


#Lighter
#else tells user that file is not found

if [ $verbose -ge 1 ]
then
	echo "Starting kmer corrections with Lighter ..."
fi

if [ -f $files ]
then
	if [ $alpha -gt 0 ]
	then
		./Lighter/lighter -r $files -od $output -t $threads -discard -k $kmer_l $genome_size $alpha
	else
		./Lighter/lighter -r $files -od $output -t $threads -discard -K $kmer_l $genome_size
	fi
elif [ -d $files  ]
then
	if [ $alpha -gt 0 ]
	then
		for i in $files/*.f*
		do
			./Lighter/lighter -r ${i} -od $output -t $threads -discard -k $kmer_l $genome_size $alpha
		done
	else
		for i in $files/*.f*
		do
			./Lighter/lighter -r ${i} -od $output -t $threads -discard -K $kmer_l $genome_size
		done
	fi
else
	echo "File/folder not found, verify the path."
fi

#Bcalm 2

find $output/*.fq* -type f > fof.txt
mv fof.txt $output
if [ $verbose -ge 1 ] #loop to silence the command if --verbose is at 0
then
	echo 'Starting bcalm2 ...'
	verbosity_level='-verbose $verbose'
else
	verbosity_level=''
fi

mkdir $output/unitigs
./bcalm/build/bcalm -in $output/fof.txt -kmer-size $kmer -nb-cores $threads -out-dir $output/unitigs $verbosity_level
mv ./fof.unitigs.fa $output/unitigs
echo "$output/unitigs/fof.unitigs.fa" > $output/fof_unitigs.txt #creates the file of file for reindeer with unitigs


# Reindeer
mkdir $output/matrix
./REINDEER/Reindeer -o $output/matrix -t $threads --nocount -k $kmer --index -f $output/fof_unitigs.txt

# MetaDBGWAS executable to get .edges and .nodes 

./src/MetaDBGWAS --files $output/fof.unitigs.fa --output $output --threads $threads --kmer $kmer
