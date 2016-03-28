#!/bin/bash
num_sim=30
site_dupl_lst=(1000 10000 100000 1000000)
br_len=0.4

if [ ! -d "./site_num_test/aln_tree_input_files/" ]; then
	mkdir ./site_num_test/aln_tree_input_files/
fi

for site_dupl in ${site_dupl_lst[*]}
do
	for n in $(seq 1 $num_sim) 
	do
		tree_file=n2_bl${br_len}.tre
		aln_file=n2_bl${br_len}_n${site_dupl}_${n}.fa
		cat site_num_test/aln/${aln_file} trees/${tree_file} > site_num_test/aln_tree_input_files/n2_bl${br_len}_n${site_dupl}_${n}_aln_tree.txt
	done
done

