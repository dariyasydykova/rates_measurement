#!/bin/bash
num_sim=10

if [ ! -d "./hyphy/rates" ]; then
	mkdir ./hyphy/rates
fi

if [ ! -d "./hyphy/rates/raw_rates" ]; then
	mkdir ./hyphy/rates/raw_rates
fi

for br_len in ${br_len_arr[*]} 
do	
	for n in $(seq 1 $num_sim) 
	do
		aln_tree_file=n2_bl${br_len}_${n}_aln_tree.txt
		rates_file=n2_bl${br_len}_${n}_rates.txt
		echo "INFILE = aln_tree_input_files/${aln_tree_file}" #> hyphy/setup.txt
		echo "OUTFILE = rates/raw_rates/${rates_file}" >> hyphy/setup.txt
		echo "MODELFILE = JC_aa.mdl" >> hyphy/setup.txt	
		#HYPHYMP hyphy/fitrates.bf
	done
done
