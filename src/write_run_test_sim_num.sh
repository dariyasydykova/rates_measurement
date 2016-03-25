#!/bin/bash
num_sim=30
site_dupl_lst=(1000 10000 100000 1000000 10000000)
br_len=0.4
##we use site1 for the test

if [ ! -d "./aln" ]; then
	mkdir aln
fi

if [ ! -d "./aln/test_site_dupl" ]; then
	mkdir aln/test_site_dupl
fi

if [ -f "./src/run_test_sim_num.sh" ]; then
	rm ./src/run_test_sim_num.sh
fi

for site_dupl in ${site_dupl_lst[*]}
do
	for n in $(seq 1 $num_sim) 
	do
		tree_file=n2_bl${br_len}.tre
		aln_file=n2_bl${br_len}_${n}.fa
		echo "python src/simulate_aln.py ../therm_constraints_rate_variation/ddG_calculations/foldX/foldx_ddG/132L_A_foldx_ddG.txt trees/${tree_file} aln/test_site_dupl/${aln_file} ${site_dupl}" >> ./src/run_test_sim_num.sh	
	done
done

chmod +x ./src/run_test_sim_num.sh