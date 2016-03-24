#!/bin/bash
num_sim=30

if [ ! -d "./hyphy/aln_tree_input_files/" ]; then
	mkdir ./hyphy/aln_tree_input_files/
fi

for br_len in `seq 0.02 0.02 0.52` 
do	
	for n in $(seq 1 $num_sim) 
	do
		tree_file=n2_bl${br_len}.tre
		aln_file=n2_bl${br_len}_${n}.fa
		cat aln/${aln_file} trees/${tree_file} > hyphy/aln_tree_input_files/n2_bl${br_len}_${n}_aln_tree.txt
	done
done

