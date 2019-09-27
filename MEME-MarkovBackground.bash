#!/bin/bash
#PBS -k o
#PBS -M tlicknac@asu.edu
#PBS -m abe
#PBS -l nodes=1:ppn=4,walltime=08:00:00,vmem=5gb

input_dir="/N/dc2/scratch/tlicknac/Paramecium_-200to0/Newest_All-Aurelias_Upstream"	
	
for filename in $input_dir/*; do 
	if [[ $filename == *[0-9][0-9][0-9][0-9][0-9]* ]]; then
		output_file=$(basename "$filename"); 
		output_dir="/N/dc2/scratch/tlicknac/Paramecium_-200to0/Newest_All-Aurelias_MEME_Results-Markov/"; 
		output="$output_dir""$output_file";  
		/N/u/tlicknac/Carbonate/meme/bin/meme $filename -dna -oc $output -objfun de -test mrs -bfile /N/dc2/scratch/tlicknac/Paramecium_-200to0/concatenated-upstream-markov.txt -nostatus -time 18000 -maxsize 60000  -nmotifs 5 -minw 6 -maxw 15
	fi
done
