#!/bin/bash

#donner les valeurs par defaut aux parametres


#General
files=false
output="./"
threads=4

#Lighter
genome_size=0
alpha=0.1
kmer=17

#Reindeer

#dbgwas



#miscealenous
version=0.1
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

#boucle while pour iterer sur les parametres donnes
#case pour assignation

while [ $# > 1 ]
do
	case $1 in
	--files) file="$2";;
	--output) output="$2";;
	--threads) threads="$2";;
	--K) kmer="$2" genome_size="$3";;
	--k) kmer="$2" genome_size="$3" alpha="$4";;
	--version) echo $version;;
	--license) echo -e $license_text;;
	*) break;
	esac
	shift
done

if [ -e files ]
then
	for i in $files
	do
		./Lighter/lighter -r ${i} -od $output -t $threads -discard -K $kmer $genome_size
	done

	for i in $files
	do
		./Lighter/lighter -r ${i} -od $output -t $threads -discard -k $kmer $genome_size $alpha
	done
fi
