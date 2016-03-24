#!/bin/bash
br_len_arr=(0 0.25 0.5 0.75 1 1.25 1.5 1.75 2)
num_sim=10

if [ ! -d "./aln" ]; then
	mkdir aln
fi

if [ -f "./src/run_sim_aln.sh" ]; then
	rm ./src/run_sim_aln.sh 
fi

for br_len in ${br_len_arr[*]} 
do	
	for n in $(seq 1 $num_sim) 
	do
		tree_file=n2_bl${br_len}.tre
		aln_file=n2_bl${br_len}_${n}.tre
		echo "python src/simulate_aln.py ../therm_constraints_rate_variation/ddG_calculations/foldX/foldx_ddG/132L_A_foldx_ddG.txt trees/${tree_file} aln/${aln_file}" >> ./src/run_sim_aln.sh	
	done
done

chmod +x ./src/run_sim_aln.sh