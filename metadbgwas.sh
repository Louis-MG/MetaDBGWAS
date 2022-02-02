#!/bin/bash

#donner les valeurs par defaut aux parametres


#General
files=NONE
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


#boucle while pour iterer sur les parametres donnes
#case pour assignation

while [ $# > 1 ]
do
	case $1 in
	--file) file="$2";;
	--output) output="$2";;
	--threads) threads="$2";;
	--K) kmer="$2" genome_size="$3";;
	--k) kmer="$2" genome_size="$3" alpha="$4";;
	--version) echo $version;;
	*) break;
	esac
	shift
done

#if [ files==NONE ]
#do
#	echo "Please submit a file with the --file flag."
#	break
#done

for i in $files
do
	./lighter -r ${i} -od $output -t $threads -discard -K $kmer $genome_size
done

for i in $files
do
        ./lighter -r ${i} -od $output -t $threads -discard -k $kmer $genome_size $alpha
done


echo $kmer
echo $genome_size
echo $alpha
