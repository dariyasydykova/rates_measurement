#!/bin/bash
br_len_arr=(0 0.25 0.5 0.75 1 1.25 1.5 1.75 2)
num_sim=10

if [ ! -d "./hyphy/aln_tree_input_files/" ]; then
	mkdir ./hyphy/aln_tree_input_files/
fi

for br_len in ${br_len_arr[*]} 
do	
	for n in $(seq 1 $num_sim) 
	do
		tree_file=n2_bl${br_len}.tre
		aln_file=n2_bl${br_len}_${n}.fa
		cat aln/${aln_file} trees/${tree_file} > hyphy/aln_tree_input_files/n2_bl${br_len}_${n}_aln_tree.txt
	done
done

