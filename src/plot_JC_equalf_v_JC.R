library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

setwd("substitution_matrices_in_pheno_models/")
r_an <- read.table("analytically_derived_rates/rates_all_sites_aa.txt",header=T)
r_inf <- read.csv("inferred_rates/processed_rates/rates_all_sites.csv")

#get the mean inferred rate at each time point
r_an_temp <- r_an %>% 
  group_by(site) %>% 
  summarise(r_tilde_small_t=r_tilde_small_t[1]) 
r_an_temp$site <- 1:length(r_an_temp$site)

d_label <- data.frame(time=c(0.0009,0.009,0.09,0.9),
                      time_label=c('0.0009','0.009','0.09','0.9'))
r <- r_inf %>% 
  group_by(time,model,num_taxa) %>% 
  mutate(inf_rate_mean=mean(rate),inf_rate_norm=rate/inf_rate_mean)  %>%
  filter(num_taxa==512,model=="JC" | model=="JC_equalf") %>%
  group_by(site,model,time) %>%
  summarise(inf_rate_norm_mean=mean(inf_rate_norm)) %>%
  left_join(r_an_temp,by='site') %>%
  left_join(d_label) 

r_JC <- r %>% ungroup() %>%
  filter(model=="JC") %>% 
  mutate(rate_JC=inf_rate_norm_mean) %>%
  select(rate_JC,site,time)

r_JC_equalf <-  r %>% ungroup() %>%
  filter(model=="JC_equalf") %>% 
  mutate(rate_JC_equalf=inf_rate_norm_mean) %>%
  select(rate_JC_equalf,site,time,time_label)

r_combined <- r_JC %>% left_join(r_JC_equalf,by=c("site","time"))

p <- ggplot(r_combined,aes(rate_JC_equalf,rate_JC)) +
  #background_grid("xy")+
  geom_point(size=0.9,alpha=0.8) + 
  geom_abline(color="red")+
  xlab("Mean inferred rate (equal frequencies)") +
  ylab("Mean inferred rate (observed frequencies)") +
  facet_wrap( ~ time_label) +
  coord_cartesian(ylim=c(0,2.5), xlim=c(0,2.5))+
  scale_y_continuous(breaks=seq(0,2.5,1)) +
  scale_x_continuous(breaks=seq(0,2.5,1)) +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))+
  panel_border()

ggsave(paste0("plots/rates_inf_diff_JC.png"))


