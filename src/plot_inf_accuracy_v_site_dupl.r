library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(readr)

setwd("substitution_matrices_in_pheno_models/")
r_inf <- read.csv( "inferred_rates/processed_rates/rates_site_dupl.csv")
r_an <- read_tsv("analytically_derived_rates/rates_all_sites.txt")

true_r <- r_an %>%
  filter(site==r_an$site[1],time==0.480002) %>%
  mutate(site=1)

r <- r_inf %>%
  group_by(site_dupl, rep) %>% 
  mutate(inf_rate_norm=inf_rate/mean(inf_rate)) %>%
  filter(site==1) %>%
  left_join(true_r,by="site")
  
p <- ggplot(r,aes(site_dupl,inf_rate_norm)) +
  #background_grid("xy")+
  geom_point(size=0.9,alpha=0.8) + 
  geom_hline(aes(yintercept=r_tilde_ms),color="red")+
  ylab("Rate") +
  xlab("Site Duplicates") +
  coord_cartesian(ylim=c(0.001,1000))+
  scale_y_log10(breaks=c(0.001,0.01,0.1,1,10,100,1000),label=c("0.001","0.01","0.1","1","10","100","1,000")) +
  scale_x_log10(breaks=c(10,100,1000,10000,100000),label=c("10","100","1,000","10,000","100,000"))+
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 10))

save_plot("plots/inf_rate_accuracy_v_site_dupl.png",plot=p)
