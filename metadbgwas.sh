#!/bin/bash
while getopts o1:o2:o3:o4:o5: flag 
do
	case '$flag' in 
		o1) option1=${OPTARG}
		o2) option1=${OPTARG}
		o3) option1=${OPTARG}
		o4) option1=${OPTARG}
		o5) option1=${OPTARG}
	esac
done
