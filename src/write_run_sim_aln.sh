#!/bin/bash
br_len_arr=(0 0.25 0.5 0.75 1 1.25 1.5 1.75 2)

if [ ! -d "./aln" ]; then
	mkdir aln
fi

if [ -f "./src/run_sim_aln.sh" ]; then
	rm ./src/run_sim_aln.sh 
fi

for br_len in ${br_len_arr[*]} 
do	
	tree_file=n2_bl${br_len}.tre
	echo "python src/simulate_aln.py ../therm_constraints_rate_variation/ddG_calculations/foldX/foldx_ddG/132L_A_foldx_ddG.txt trees/${tree_file}" >> ./src/run_sim_aln.sh	
done

chmod +x ./src/run_sim_aln.sh