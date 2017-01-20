library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(readr)

setwd("substitution_matrices_in_pheno_models/")
r_inf <- read_csv( "inferred_rates/processed_rates/rates_site_dupl.csv")

p <- ggplot(r,aes(site_dupl,inf_rate)) +
  #background_grid("xy")+
  geom_point(size=0.9,alpha=0.8) + 
  geom_hline(aes(yintercept=1),color="red")+
  ylab("Rate") +
  xlab("Site Dupl") +
  coord_cartesian(ylim=c(0.01,100))+
  scale_y_log10(breaks=c(0.01,0.1,1,10,100),label=c("0.01","0.1","1","10","100")) +
  scale_x_log10(breaks=c(100,1000,10000,100000,1000000),label=c("100","1000","10000","100000","1000000"))+
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 10))

save_plot("plots/inf_rate_accuracy_v_site_dupl.png",plot=p)
