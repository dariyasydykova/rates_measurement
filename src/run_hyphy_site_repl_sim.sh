#!/bin/bash
site_dupl_arr=(10 100 1000 10000 100000)
num_sim=50

tree_file=n2_bl0.24.tre

for site_dupl in ${site_dupl_arr[*]}
do	
	for i in $(seq 1 $num_sim) 
	do
		aln_tree_file=n2_bl0.24_site_dupl_${site_dupl}_${i}_aln_tree.txt
		aln_tree_dir=aln_tree_files/site_dupl
		rates_file=n2_bl0.24_site_dupl_${site_dupl}_${i}_rates.txt
		echo "INFILE = ${aln_tree_dir}/${aln_tree_file}" > hyphy/setup.txt
		echo "OUTFILE = ../inferred_rates/raw_rates/site_dupl/${rates_file}" >> hyphy/setup.txt
		echo "SITE_DUPL = ${site_dupl}" >> hyphy/setup.txt
		HYPHYMP hyphy/fitrates_JC_equalf.bf
	done
done
