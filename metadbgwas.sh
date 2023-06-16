#!/bin/bash

# Louis-Mael Gueguen lm.gueguen@orange.fr


GREEN=$(tput setaf 2) # green color
NC=$(tput sgr0) # No Color

#############################################
#
#       Parameters
#
#############################################


#General
files=
output="./"
declare -i threads=4
declare -i verbose=1
clean=false
skip1=false
skip2=false
skip3=false

#Lighter
declare -i genome_size=
alpha=0
declare -i kmer_l=17

#bcalm2
declare -i kmer=31
decalre -i abundance_min=5

#Bifrost
#kmer as bcalm

#dbgwas
strains='' #strains file
newick='' #phylo tree file
ncDB='' #nucleotide database file
ptDB='' #protein database file
keepNA=''
threshold=
#kmer as bcalm too

#miscealenous

Version()
{
	echo "\nMetadDBGWAS 1.0\n"
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
--clean removes intermediary files to save space if you are worried about your storage.
--skip1 skips the Lighter correction step. Corrected files are supposed to be in the output folder.
--skip2 skips the Lighter and Bcalm2 steps. Corrected files and unitigs folder are supposed to be in the output folder.
--skip3 skips the Lighter, Bcalm2 and REINDEER steps. Corrected files, unitigs and matrix folder are supposed to be in the output folder.

        * Lighter
NOTE: if your datset contains different bacterial genomes with very different size, it is better to choose --k option and provide the pick-rate (noted alpha).
--K <int> <int> kmer length and genome size (in base). Recommended is 17 G.
        or
--k <int> <int> <float> kmer length and genome size (in base), alpha (probability of sampling a kmer). Recommended is 17 G alpha. alpha is best chosen at 7/coverage.

        * bcalm
--kmer <int> kmer length used for unitigs build. Default to 31.
--abundance-min <int> Minimum number of occurence of a kmer to keep it in the Union DBG. Default to 5, highly recommended to change to the 2.5% quantile of the Poisson law with lambda = coverage.

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
bash metadbgwas.sh --files ./test/ --output ./output --strains ./strains --threads 4 --verbose 1 --K 17 6000000
	"
}

#Parameters parsing

if [[ $# -eq 0 ]]
then
	Help
	exit
fi

while [[ $# -gt 0 ]]
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
	--skip1) skip1=true
	shift;;
	--skip2) skip2=true skip1=true
	shift;;
	--skip3) skip3=true skip2=true skip1=true
	shift;;
	--K) kmer_l="$2" genome_size="$3"
	shift 3;;
	--k) kmer_l="$2" genome_size="$3" alpha="$4"
	shift 4;;
	-k | --kmer) kmer="$2"
	shift 2;;
	--abundance-min) abundance_min="$2"
	shift 2;;
	--strains) strains="$2"
	shift 2;;
	--newick) newick="-newick $2"
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
	-* | --*) unknown="$1" ;echo "Unknown option: ${unknown}"; Help; exit;;
	*) ;;
	esac
done

#gets local directory of metadbgwas
metadbgwas_path=$(dirname $0)
metadbgwas_path=$(cd $metadbgwas_path && pwd)

#creates output dir if it doesnt exist yet
#if folder exists and no step must be skipped:
if [ -d $output ]
then
	if [ $skip1 = false ]
	then
		if [ "$(ls $output)" ]
		then
	    		echo "Warning: $output is not Empty." #folder should be emptied or removed for a new run
		fi
	fi
#else the folder does not exists and it is created
else
	mkdir $output
fi

#saves command line:
echo -e "--files $files\n--strains $strains\n--output $output\n--kmer $kmer\n--abundance-min $abundance_min\n--k/K $kmer_l $genome_size $alpha\n--newick $newick\n--nc-DB $ncDB\n--pt-db $ptDB\n--keepNA $keepNA\n--threshold $threshold" > $output/command_line.txt


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

if [ $verbose -ge 1 ] && [ $skip1 = false ]
then
        echo -e "${GREEN}Starting kmer corrections with Lighter ...${NC}"
fi

if [ $skip1 = false ]
then
	if [ -d $files ] #if the argument is a path to the files
	then
		counter_total=$(ls -1 $files/*.f* | wc -l)
		counter=1
	        if [ $alpha -gt 0 ]
	        then
			for i in $files/*.f*
			do
				echo -ne "${counter}/${counter_total} Processing files ...\n"
				$metadbgwas_path/Lighter/lighter -r ${i} -od $output -t $threads -discard -k $kmer_l $genome_size $alpha
				tput cuu 14
				tput ed
				counter=$(($counter+1))
			done
		else
			for i in $files/*.f*
			do
				echo -ne "${counter}/${counter_total} Processing files ...\n"
				$metadbgwas_path/Lighter/lighter -r ${i} -od $output -t $threads -discard -K $kmer_l $genome_size
				tput cuu 14
				tput ed
				counter=$(($counter+1))
			done
		fi
	else #if the argument if a fof
		counter_total=$(cat $files | wc -l)
                counter=1
                if [ $alpha -gt 0 ]
                then
                        while read line :
                        do
                                echo -ne "${counter}/${counter_total} Processing files ...\n"
                                $metadbgwas_path/Lighter/lighter -r ${line} -od $output -t $threads -discard -k $kmer_l $genome_size $alpha
                                tput cuu 14
                                tput ed
                                counter=$(($counter+1))
                        done < $files
                else
                        while read line:
                        do
                                echo -ne "${counter}/${counter_total} Processing files ...\n"
                                $metadbgwas_path/Lighter/lighter -r ${line} -od $output -t $threads -discard -K $kmer_l $genome_size
                                tput cuu 14
                                tput ed
                                counter=$(($counter+1))
                        done < $files
                fi
	fi
else
	echo -e "${GREEN}Skipping Lighter step ...${NC}"
fi


#############################################
#
#       Bcalm 2
#
#############################################

if [ $skip2 = false ]
then
	find $output/*.cor.f* -type f > $output/fof.txt
	if [ $verbose -ge 1 ] #loop to silence the command if --verbose is at 0
	then
		echo -e "${GREEN}Starting Bcalm2 ...${NC}"
		verbosity_level='-verbose $verbose'
	else
		verbosity_level=''
	fi
	mkdir $output/unitigs
	counter_total=$(ls -1 $output/*.cor.f* | wc -l)
	declare -i counter; counter=1
	#we create the de Bruijn Graph of the files we want to index
	for i in $output/*.cor.f*
	do
		echo "${counter}/${counter_total} Processing files ..."
	        $metadbgwas_path/bcalm/build/bcalm -in $i -kmer-size $kmer -nb-cores $threads  $verbosity_level -abundance-min 1 # TODO: add the -out option to give prefix and avoid moving files around
		mv *.unitigs.fa $output/unitigs
		counter+=1 
		tput cuu 4
		tput ed
	done
	find $output/unitigs/*.unitigs.fa -type f > $output/unitigs/fof_unitigs_index.txt
	#the option abundance min is used to keep all kmers: we already corrected them, and not keeping them all to build the unitigs would have 2 unfortunate consequences :
		# 1 some variation would be lost, and we want to analyse it !
		# 2 when mapping the kmers to the unitig graph, that causes a segfault (index > ULONG_MAX)
	$metadbgwas_path/bcalm/build/bcalm -in $output/fof.txt -kmer-size $kmer -nb-cores $threads -out-dir $output/unitigs $verbosity_level -abundance-min $abundance_min
	rm $output/fof.txt
	mv ./fof.unitigs.fa ./unitigs.fa
	mv ./unitigs.fa $output/unitigs
	if [ $verbose -ge 2 ] #loop to silence the command if --verbose is at 0
	then
		echo "Cleaning temporary files ..."
	fi
	rm $output/unitigs/fof.h5
else
	echo -e "${GREEN}Skipping Bcalm2 step ...${NC}"
fi

#############################################
#
#	Bifrost
#
#############################################


#TODO: make the fof, use the correct folder structure
if [ $skip3 = false ]
then
	if [ ! -d $output/step1 ]
	then
		mkdir $output/step1
	fi
	echo -e "${GREEN}Starting Bifrost${NC}"
	# first we build the colored dbg:
	$metadbgwas_path/bifrost/build/src/Bifrost build -t $threads --colors --input-ref-file $output/unitigs/fof_unitigs_index.txt -o $output/step1/bifrost_colored_dbg
	#then we query the unitigs on the bdg of kmers we built precendently:
	$metadbgwas_path/bifrost/build/src/Bifrost query -t $threads -e 1 --input-graph-file $output/step1/bifrost_colored_dbg.gfa.gz --input-query-file $output/unitigs/unitigs.fa --input-color-file $output/step1/bifrost_colored_dbg.color.bfg -o $output/step1/result_genomes_bifrost_query
	header=$(cut -f1 $strains | sed -z 's/\n/\t/g')
	header=${header/ID/query_name}
	sed -i "1s/.*/$header/g" $output/step1/result_genomes_bifrost_query.tsv
else
	echo -e "${GREEN}Skipping Bifrost step ...${NC}"
fi

#############################################
#
#	MetaDBGWAS
#
#############################################

#conversion of the information of unitigs.fa to a GFA (Graphical fragment assembly http://gfa-spec.github.io/GFA-spec/GFA1.html)
python3 $metadbgwas_path/bcalm/scripts/convertToGFA.py $output/unitigs/unitigs.fa $output/step1/graph.gfa $kmer
#see if I can change that script to cpp

#fix for pahntomjs
export OPENSSL_CONF=/dev/null
# MetaDBGWAS executable that generates input ifles for bugwas and gemma, runs statistical tests, then finally generates output
$metadbgwas_path/tools/src/MetaDBGWAS --files $output/unitigs/unitigs.fa --output $output --threads $threads --kmer $kmer --strains $strains $keepNA --threads $threads --output $output $ncDB $ptDB $newick $threshold --no-preview

if [ $clean = true ]
then
	if [ $verbose -ge 1 ]
	then
		echo "Compressing step1 files ..."
	fi
	gzip $output/step1/unitigs2PhenoCounter
	gzip $output/step1/graph.gfa
	gzip $output/step1/bugwas_input.all_rows.binary
	gzip $output/step1/bugwas_input.unique_rows.binary
fi
