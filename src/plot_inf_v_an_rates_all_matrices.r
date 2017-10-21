library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

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
  filter(num_taxa==512) %>%
  group_by(site,model,time) %>%
  summarise(inf_rate_norm_mean=mean(inf_rate_norm)) %>%
  left_join(r_an_temp,by='site') %>%
  left_join(d_label)

p <- ggplot(r,aes(r_tilde_small_t,inf_rate_norm_mean)) +
    #background_grid("xy")+
    geom_point(size=0.9,alpha=0.8) + 
    geom_abline(color="red")+
    ylab("Mean inferred rate") +
    xlab("Analytically derived rate") +
    #xlab(expression(paste("Analytically Derived Rates (", hat(r)^(k), "for small t)"))) +
    facet_grid(model ~ time_label) +
    coord_cartesian(ylim=c(0,2.5), xlim=c(0,2.5))+
    scale_y_continuous(breaks=seq(0,2.5,1)) +
    scale_x_continuous(breaks=seq(0,2.5,1)) +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12))+
    panel_border()
  
ggsave(paste0("plots/inf_v_an_rates_all_matrices.png"))


