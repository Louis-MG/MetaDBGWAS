#!/bin/bash

# Louis-Mael Gueguen lm.gueguen@orange.fr


#donner les valeurs par defaut aux parametres


#General
files=
output="./"
threads=4
verbose=0

#Lighter
genome_size=
alpha=0
kmer_l=17

#Reindeer

#dbgwas
strain=''
newick=''
kmer=31

#strains
#k
#newick



#miscealenous
version="MetadDBGWAS 0.1"
license=false
license_text="Copyright (C) 2022 Louis-Maël Gueguen\n
\n
This software is provided 'as-is', without any express or implied\n
warranty.  In no event will the authors be held liable for any damages\n
arising from the use of this software.\n
\n
Permission is granted to anyone to use this software for any purpose,\n
including commercial applications, and to alter it and redistribute it\n
freely, subject to the following restrictions:\n
\n
1. The origin of this software must not be misrepresented; you must not\n
   claim that you wrote the original software. If you use this software\n
   in a product, an acknowledgment in the product documentation would be\n
   appreciated but is not required.\n
2. Altered source versions must be plainly marked as such, and must not be\n
   misrepresented as being the original software.\n
3. This notice may not be removed or altered from any source distribution.\n
\n
Louis-Maël Gueguen lm.gueguen@orange.fr\n"

help=false
help_text="To_be_written"

#boucle while pour iterer sur les parametres donnes
#case pour assignation

while [ $# > 1 ]
do
	case $1 in
	-f | --files) files="$2"
	shift 2;;
	-o | --output) output="$2"
	shift 2;;
	-t | --threads) threads="$2"
	shift 2;;
	--K) kmer_l="$2" genome_size="$3"
	shift 3;;			#should add a true/false variable to select the type of run for Lighter
	--k) kmer_l="$2" genome_size="$3" alpha="$4"
	shift 4;;
	-k | --kmer) kmer="#2"
	shift 2;;
	--version) echo $version; exit 0;;
	--license) echo $license_text; exit 0 ;;
	-v | --verbose) verbose="$2"
	shift 2;;
	-h | --help) echo $help_text; exit 0
	shift 2;;
	-* | --*) echo "Unknown option"; exit 1
	shift 2;;
	*) break;
	esac
done


echo ${kmer_l} ${threads} ${genome_size} ${output} $files

#if the file exists, then the runs for Lighter kmer correction are differentiated by alpha value

if [ $verbose -ge 1 ]
then
	echo "Starting kmer corrections with Lighter ..."
fi

if [ -f $files ]
then
	if [ $alpha -gt 0 ]
	then
		for i in $files
		do
			echo "lighter petit k"
			../Lighter/lighter -r ${i} -od $output -t $threads -discard -k $kmer_l $genome_size $alpha
		done
	else
		for i in $files
		do
			echo 'lighter K'
			../Lighter/lighter -r ${i} -od $output -t $threads -discard -K $kmer_l $genome_size
		done
	fi
fi
