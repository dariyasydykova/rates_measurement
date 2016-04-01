#!/bin/bash
num_sim=30
site_dupl_lst=100000

if [ ! -d "./hyphy/aln_tree_files/" ]; then
	mkdir ./hyphy/aln_tree_files/
fi

for br_len in `seq -f %.2f 0.02 0.02 0.52` 
do	
	for n in $(seq 1 $num_sim) 
	do
		tree_file=n2_bl${br_len}.tre
		aln_file=n2_bl${br_len}_${n}.fa
		cat aln/${aln_file} trees/${tree_file} > hyphy/aln_tree_files/n2_bl${br_len}_${n}_aln_tree.txt
	done
done

