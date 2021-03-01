#!/bin/bash


while getopts s:o:k:d:p:f: flag
do
	case "${flag}" in
		s) trnascan_SSfile=${OPTARG};;
		o) output_directiry=${OPTARG};;
		k) kingdom=${OPTARG};;
		d) database=${OPTARG};;
		p) proteome_file=${OPTARG};;
		f) hmm_file=${OPTARG};;
	esac
done


python tRNA-MAP.py -s $trnascan_SSfile -o $output_directiry -k $kingdom -d $database &
python HMMprogram.py -k $kingdom -p $proteome_file -f $hmm_file
