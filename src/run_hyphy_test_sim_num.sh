#!/bin/bash
num_sim=30
site_dupl_lst=(1000000) ##1000 100000 1000000)
br_len=0.4

if [ ! -d "./site_num_test/rates/" ]; then
	mkdir ./site_num_test/rates/
fi

if [ ! -d "./site_num_test/rates/raw_rates" ]; then
	mkdir ./site_num_test/rates/raw_rates
fi


for site_dupl in ${site_dupl_lst[*]}
do
	for n in $(seq 1 $num_sim) 
	do
		aln_tree_file=n2_bl${br_len}_n${site_dupl}_${n}_aln_tree.txt
		rates_file=n2_bl${br_len}_n${site_dupl}_${n}_rates.txt
		echo "INFILE = ../site_num_test/aln_tree_input_files/${aln_tree_file}" > hyphy/setup.txt
		echo "OUTFILE = ../site_num_test/rates/raw_rates/${rates_file}" >> hyphy/setup.txt
		HYPHYMP hyphy/fitrates.bf
	done
done
