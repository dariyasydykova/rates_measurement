#!/bin/bash
num_sim=30
site_dupl=100000
total_sites=11
sim_space='codon'
#sim_space='aa'

if [ ! -d "./aln" ]; then
	mkdir aln
fi

if [ -f "./src/run_sim_aln.sh" ]; then
	rm ./src/run_sim_aln.sh 
fi

for br_len in `seq -f %.2f 0.02 0.02 0.52` 
do	
	for n in $(seq 1 $num_sim) 
	do
		tree_file=n2_bl${br_len}.tre
		aln_file=n2_bl${br_len}_${n}.fa
		
		if [ $sim_space = 'aa' ]; then
			echo "python src/simulate_aa_aln.py 132L_A_foldx_ddG.txt trees/${tree_file} ../substitution_matrices_in_pheno_models_data/aa_aln/${aln_file} ${site_dupl} ${total_sites}" >> ./src/run_sim_aln.sh	
		fi
		
		if [ $sim_space = 'codon' ]; then 
			echo "python src/simulate_codon_aln.py 132L_A_foldx_ddG.txt trees/${tree_file} ../substitution_matrices_in_pheno_models_data/codon_aln/${aln_file} ${site_dupl} ${total_sites}" >> ./src/run_sim_aln.sh	
		fi 
	done
done

chmod +x ./src/run_sim_aln.sh