#!/bin/bash
num_sim=30

if [ -f ./src/run_translate.sh ]; then
	rm ./src/run_translate.sh  
fi

for br_len in `seq -f %.2f 0.02 0.02 0.52` 
do	 
	for i in $(seq 1 $num_sim) 
	do
		aln=n2_bl${br_len}_${i}.fa
		echo "python ./src/translate_aln_codon_to_aa.py ../substitution_matrices_in_pheno_models_data/codon_aln/${aln} ../substitution_matrices_in_pheno_models_data/translated_aln/${aln}" >> ./src/run_translate.sh 
	done
done

chmod +x ./src/run_translate.sh