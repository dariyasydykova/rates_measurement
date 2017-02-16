library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

setwd("substitution_matrices_in_pheno_models/")
t <- read.table("analytically_derived_rates/codon_aa_equil_freq.txt",header=T)

t %>% group_by(site,amino_acid) %>% 
  summarize(eq_freq_aa=eq_freq_aa[1],eq_freq_codon_sum=sum(eq_freq_codon)) %>%
  mutate(eq_freq_diff=eq_freq_aa-eq_freq_codon_sum) -> eq_freq_diff

eq_freq_diff %>% summarize(diff_sum=sum(abs(eq_freq_diff))) -> sum_diff

p <- ggplot(eq_freq_diff,aes(x=eq_freq_diff))+
  geom_histogram(binwidth=0.01)+
  xlab("Differences in equilibrium frequencies")+
  facet_wrap(~site)
