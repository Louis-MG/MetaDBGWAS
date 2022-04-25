#!/bin/bash

# Louis-Mael Gueguen lm.gueguen@orange.fr


GREEN='\e[0;32m' # green color
NC='\e[0m' # No Color

#############################################
#
#       Parameters
#
#############################################


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
#kmer as bcalm

#dbgwas
strains='' #strais file
newick='' #phylo tree file
ncDB='' #nucleotide database file
ptDB='' #protein database file
keepNA=''
threshold=
#kmer as bcalm too

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
--strains A text file describing the strains containing 3 columns: 1) ID of the strain; 2) Phenotype (a real number or NA); 3) Path to a multi-fasta file containing the sequences of the strain. This file needs the header. Tab separated.
--newick Optional path to a newick tree file. If (and only if) a newick tree file is provided, the lineage effect analysis is computed and PCs figures are generated.
--nc-db Optional  A list of Fasta files separated by comma containing annotations in a nucleotide alphabet format (e.g.: -nc-db path/to/file_1.fa,path/to/file_2.fa,etc). You can customize these files to work better with DBGWAS (see https://gitlab.com/leoisl/dbgwas/tree/master#customizing-annotation-databases).
--pt-db Optionnal A list of Fasta files separated by comma containing annotations in a protein alphabet format (e.g.: -pt-db path/to/file_1.fa,path/to/file_2.fa,etc). You can customize these files to work better with DBGWAS (see https://gitlab.com/leoisl/dbgwas/tree/master#customizing-annotation-databases).
--keepNA Optionnal Keep strains with phenotype NA.
--threshold maximum value for which phenotype will be considered to be 0.

        * Miscellaneous
--license prints the license text in standard output.
--help displays help.

        * Exemple
bash metadbgwas.sh --files /test/ --output ./output --threads 4 --verbose 1 --K 17 6000000
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
	--strains) strains="$2"
	shift 2;;
  --newick) newick="$2"
  shift 2;;
	--nc-db) ncDB="-nc-db $2"
	shift 2;;
	--pt-db) ptDB="-pt-db $2"
	shift 2;;
	--keepNA) keepNA="-keepNA"
	shift;;
  --threshold) threshold="-threshold $2"
  shift;;
	--version) Version; exit;;
	--license) License; exit;;
	-v | --verbose) verbose="$2"
	shift 2;;
	-h | --help) Help; exit;;
	-* | --*) echo "Unknown option"; exit;;
	*) break;
	esac
done

#gets local directory of metadbgwas
metadbgwas_path=$(dirname $0)
metadbgwas_path=$(cd $metadbgwas_path && pwd)

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


#############################################
#
#       Lighter
#
#############################################

#else tells user that file is not found
if [ $verbose -ge 1 ]
then
        echo "${GREEN}Starting kmer corrections with Lighter ...${NC}"
fi
#checks folder existence
if [ -d $files ]
then
        if [ $alpha -gt 0 ]
        then
                for i in $files/*.f*
                do
                        $metadbgwas_path/Lighter/lighter -r ${i} -od $output -t $threads -discard -k $kmer_l $genome_size $alpha
                done
        else
                for i in $files/*.f*
                do
                        $metadbgwas_path/Lighter/lighter -r ${i} -od $output -t $threads -discard -K $kmer_l $genome_size
                done
	fi
else
        echo "Folder not found or is not a folder, verify the path."
	exit 0
fi

#############################################
#
#       Bcalm 2
#
#############################################

find $output/*.cor.f* -type f > $output/fof.txt
if [ $verbose -ge 1 ] #loop to silence the command if --verbose is at 0
then
	echo "${GREEN}Starting bcalm2 ...${NC}"
	verbosity_level='-verbose $verbose'
else
	verbosity_level=''
fi
mkdir $output/unitigs
#we create the de Bruijn Graph of the files we want to index
for i in $output/*.cor.f*
do
	echo $i
        $metadbgwas_path/bcalm/build/bcalm -in $i -kmer-size $kmer -nb-cores $threads  $verbosity_level -abundance-min 1 # TODO: add the -out option to give prefix and avoid moving files around
	mv *.unitigs.fa $output/unitigs
done
find $output/unitigs/*.unitigs.fa -type f > $output/unitigs/fof_unitigs_index.txt
#the option abundance min is used to keep all kmers: we already corrected them, and not keeping them all to build the unitigs would have 2 unfortunate consequences :
	# 1 some variation would be lost, and we want to analyse it !
	# 2 when mapping the kmers to the unitig graph, that causes a segfault (index > ULONG_MAX)
$metadbgwas_path/bcalm/build/bcalm -in $output/fof.txt -kmer-size $kmer -nb-cores $threads -out-dir $output/unitigs $verbosity_level -abundance-min 1
rm $output/fof.txt
mv ./fof.unitigs.fa ./unitigs.fa
mv ./unitigs.fa $output/unitigs


#############################################
#
#	Reindeer
#
#############################################

mkdir $output/step1
# first we index:
$metadbgwas_path/REINDEER/Reindeer --index -f $output/unitigs/fof_unitigs_index.txt -o $output/matrix -k $kmer -t 1
#then we query the unitigs on the index of kmers we built precendently:
$metadbgwas_path/REINDEER/Reindeer --query -l $output/matrix -q $output/unitigs/unitigs.fa -o $output/matrix -t 1

# MetaDBGWAS executable to get .edges and .nodes, gemman and bugwas input files, as well as the pheno files.
$metadbgwas_path/src/MetaDBGWAS --files $output/unitigs/unitigs.fa --output $output --threads $threads --kmer $kmer --strains $strains


#############################################
#
#       DBGWAS
#
#############################################


#creating the step 2 folder :
mv $output/graph.edges.dbg $output/graph.nodes $output/step1

#starting DBGWAS at step 2:

echo "${GREEN}Starting DBGWAS ...${NC}"
$metadbgwas_path/DBGWAS/bin/DBGWAS -k $kmer -strains $strains -keepNA -nb-cores $threads -output $output -skip1 $keepNA $ncDB $ptDB $threshold

