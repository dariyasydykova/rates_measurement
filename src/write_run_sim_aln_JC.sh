#!/bin/bash
num_sim=30
site_dupl=100000

if [ -f "run_sim_aln_JC.sh" ]; then
	rm run_sim_aln_JC.sh 
fi

aln_dir=../../substitution_matrices_in_pheno_models_data/aa_aln/JC
for br_len in `seq -f %.2f 0.02 0.02 0.52` 
do	
	for n in $(seq 1 $num_sim) 
	do
		tree_file=n2_bl${br_len}.tre
		aln_file=n2_bl${br_len}_${n}.fa
		
		echo python simulate_aln_true_JC.py ../trees/${tree_file} ${aln_dir}/${aln_file} ${site_dupl} >> run_sim_aln_JC.sh	 
	done
done

chmod +x run_sim_aln_JC.sh