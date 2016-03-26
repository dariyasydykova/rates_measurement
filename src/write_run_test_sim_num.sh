#!/bin/bash
num_sim=30
site_dupl_lst=(1000 10000 100000 1000000)
br_len=0.4
num_model=1
##we use site1 for the test

if [ ! -d "./site_num_test/aln" ]; then
	mkdir .site_num_test/aln
fi

if [ -f "./src/run_test_sim_num.sh" ]; then
	rm ./src/run_test_sim_num.sh
fi

for site_dupl in ${site_dupl_lst[*]}
do
	for n in $(seq 1 $num_sim) 
	do
		tree_file=n2_bl${br_len}.tre
		aln_file=n2_bl${br_len}_n${site_dupl}_${n}.fa
		echo "python $HOME/substitution_matrices_in_pheno_models/src/simulate_aln.py $HOME/132L_A_foldx_ddG.txt $HOME/substitution_matrices_in_pheno_models/trees/${tree_file} $HOME/substitution_matrices_in_pheno_models/site_num_test/aln/${aln_file} ${site_dupl} ${num_model}" >> ./src/run_test_sim_num.sh	
	done
done

chmod +x ./src/run_test_sim_num.sh