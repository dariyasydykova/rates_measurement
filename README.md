# Theory of measurement for site-specific evolutionary rates in amino-acid sequences

This repository contains all the scripts and data to reproduce the results of:

D. K. Sydykova, C. O. Wilke (2018). Theory of measurement for site-specific evolutionary rates in amino-acid sequences

## Contents

* `analytical_rates` contains rates that were derived analytically and were calculated  
* `inferred_rates`
* `q_matrix` contains site-wise substitution matrices (*Q*) used for simulating alignments. Files that start with `132L_A` indicate substitution matrices that were calculated using data from Echave et al (2015) for egg white lysozyme (PDB ID: 132L). For example, file `132L_A_site79_q_matrix.npy` corresponds to the substitution matrix calculated for site 79. The site position corresponds to the site position given by the file `132L_A_foldx_ddG.txt`, which was directly copied from the git repository for Echave et al (2015) [https://github.com/wilkelab/therm_constraints_rate_variation](https://github.com/wilkelab/therm_constraints_rate_variation/tree/master/ddG_calculations/foldX/foldx_ddG)
    * `amino_acid` contains amino acid substitution matrices
    * `codon` contains codon substitution matrice
* `trees`
* `hyphy`
* `plots`
* `src`

## Analysis
