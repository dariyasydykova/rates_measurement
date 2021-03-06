
Theory of measurement for site-specific evolutionary rates in amino-acid sequences
==================================================================================

This repository contains all the scripts and data to reproduce the results of:

D. K. Sydykova, C. O. Wilke (2018). Theory of measurement for site-specific evolutionary rates in amino-acid sequences

Contents
--------

-   `analytical_rates` contains rates that were calculated using analytical derivations. The following list is the files contained in this direcotry and their descriptions.
    -   `all_sites_aa.csv` contians site-wise rates for every site in egg white lysozyme (PDB ID: 132L) calculated for different times (column `site` in `all_sites_aa.csv` directly corresponds to column `SITE` in `132L_A_foldx_ddG.txt`). These rates were calculated assuming that the true model is a mutation-selection model, and the inference model is Jukes-Cantor (equations 3-5). This file was generated with the command

        `python analytical_rate_aa.py -m 125 -q q_matrix/amino_acid/ -o all_sites_aa.csv`.

    -   `ten_sites_aa.csv` contains site-wise rates for the first ten sites in egg white lysozyme (PDB ID: 132L) calculated for different times (column `site` in `ten_sites_aa.csv` directly corresponds to column `SITE` in `132L_A_foldx_ddG.txt`). These rates were calculated assuming that the true model is a mutation-selection model, and the inference model is Jukes-Cantor (equations 3-5). This file was generated with the command

        `python analytical_rate_aa.py -m 10 -q q_matrix/amino_acid/ -o ten_sites_aa.csv`.

    -   `ten_sites_aa_true_JC.csv` contains site-wise rates for ten sites that were calculated under the assumption that both the true model and the inference model are Jukes-Cantor. This file was generated with the command

        `python analytical_rate_aa_true_JC.py -m 10 -o ten_sites_aa_true_JC.csv`.

    -   `ten_sites_aa_QM.csv` contains site-wise rates for when rate is measured with an arbitrary *Q*<sub>M</sub> matrix and for when rate is measured with a Jukes-Cantor matrix (equation 1).
    -   `ten_sites_codon.csv` contians site-wise rates for every site in egg white lysozyme (PDB ID: 132L) calculated for different times (column `site` in `all_sites_aa.csv` directly corresponds to column `SITE` in `132L_A_foldx_ddG.txt`). These rates were calculated assuming that the true model is a codon mutation-selection model, and the inference model is an amino acid Jukes-Cantor (equation 6 and equations 22S and 24S). This file was generated with the command

        `python analytical_rate_codon.py -m 10 -q q_matrix/codon/ -o ten_sites_codon.csv`.

-   `inferred_rates` contains files with site-wise rates inferred with HyPhy. There are two directories in `inferred rates`, `raw_rates` and `processed_rates`. `raw_rates` contains individual files for a simulated alignment (one file per alignment), and `processed_rates` contains concatenated files from `raw_rates`.The following list describes the directories contained in `raw_rates`:
    -   `JC`, inferred rates when the true and the inference models are both Jukes-Cantor-like.
    -   `all_sites`, inferred rates when the true model is MutSel, and the inference model is either Jukes-Cantor-like (JC), WAG, JTT, or LG.
    -   `site_dupl`, rates inferred with JC for alignments with different number of site duplicates.
    -   `ten_sites`, inferred rates when the true model is MutSel, and the inference model is JC.
    -   `translated`, inferred rates when the true model is a codon MutSel model, and the inference model is amino acid JC.
-   `q_matrix` contains site-wise substitution matrices *Q* used for simulating alignments and for calculating site-wise rates. There are two directories in `q_matrix`, `amino_acid` and `codon`. `amino acid` contains amino acid substitution matrices, and `codon` contains codon substitution matrices.

    Files that start with `132L_A` indicate substitution matrices that were calculated using data from Echave et al (2015) for egg white lysozyme (PDB ID: 132L). For example, file `132L_A_site79_q_matrix.npy` corresponds to the substitution matrix calculated for site 79. The site positions here correspond to the site positions given by the file `132L_A_foldx_ddG.txt`, which was directly copied from the git repository for Echave et al (2015) [https://github.com/wilkelab/therm\_constraints\_rate\_variation](https://github.com/wilkelab/therm_constraints_rate_variation/tree/master/ddG_calculations/foldX/foldx_ddG). These matrices were calculated according to the mutation-selection (MutSel) theory by Halpern and Bruno (1998). The script to calculate amino acid MutSel matrices is `src/calculate_aa_mutsel_Q.py`, and the script to calculate codon MutSel matrices is `src/calculate_codon_mutsel_Q.py`.

    Files that start with `site0_JC` indicate substitution matrices *Q* defined as *Q* = *r*<sup>(*k*)</sup>*Q*<sub>JC</sub>. Here, *Q*<sub>JC</sub> is the Jukes-Cantor-like matrix, and *r*<sup>(*k*)</sup> is the true rate at site *k*. True rates were generated by the script `src/analytical_rate_aa_true_JC.py` and were stored in `analytical_rates/ten_sites_aa_true_JC.csv`. The script to calculate amino acid substitution matrices *Q* defined as *Q* = *r*<sup>(*k*)</sup>*Q*<sub>JC</sub> is `src/calculate_aa_true_JC_Q.py`.

-   `trees` contains tree files that were used for simulating alignments. Each file name stores the number of branches (`n`) in the tree and the branch lengths (`bl`). For example, file `n2_bl0.005.tre` describes a tree with 2 branches of lengths 0.005 each.
-   `hyphy` contains all scripts and files to run HyPhy.
-   `plots` contains plots generated for the manuscript.
-   `src` contains code to run the analysis.

Analysis
--------
