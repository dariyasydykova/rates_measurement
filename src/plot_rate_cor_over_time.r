library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

t1 <- read.table("numerically_derived_rates/all_rates.txt",header=T)
t2 <- read.csv("inferred_rates/processed_rates/all_rates_all_sites.csv")

r_inf <- group_by(t2,time,model) %>% mutate(rate_mean=mean(rate)) %>% mutate(rate_norm = rate / rate_mean) 
r_an <- group_by(t1,site) %>% mutate(r_tilde_ms_norm_small_t=r_tilde_ms_norm[1])

an_r_lst <- unique(r_an$r_tilde_ms_norm_small_t)
r_inf <- mutate(r_inf, an_r=rep(an_r_lst,30))

r_inf <- group_by(r_inf,site,time,model) %>% mutate(inf_mean_rate_per_an_r=mean(rate_norm))
  
num_taxa_lst <- unique(r_inf$num_taxa)
for (nt in num_taxa_lst){
  r_inf_filtered <- r_inf %>% filter(num_taxa==nt,model!="JC")
  
  ##formating time to avoid "e" in facet_grid labeling
  if (nt==64) r_inf_filtered$time_str <- c(rep("0.0006",11160),rep("0.006",11160),rep("0.06",11160),rep("0.6",11160))
  if (nt==128) r_inf_filtered$time_str <- c(rep("0.0007",11160),rep("0.007",11160),rep("0.07",11160),rep("0.7",11160))
  if (nt==256) r_inf_filtered$time_str <- c(rep("0.0008",11160),rep("0.008",11160),rep("0.08",11160),rep("0.8",11160))
  if (nt==512) r_inf_filtered$time_str <- c(rep("0.0009",11160),rep("0.009",11160),rep("0.09",11160),rep("0.9",11160))
  
  p_rates <- ggplot(r_inf_filtered,aes(an_r,inf_mean_rate_per_an_r)) +
    background_grid("xy")+
    geom_point(alpha=0.03) + 
    #geom_smooth(method = "lm", se = FALSE)+
    geom_abline(color="red")+
    ylab("Mean inferred Rate") +
    xlab(expression(paste("Numerically Derived Rates (", tilde(r)^(k), "for small t)"))) +
    facet_grid(model ~ time_str) +
    coord_cartesian(ylim=c(0,2.5), xlim=c(0,2.5))+
    scale_y_continuous(breaks=seq(0,2.5,1)) +
    scale_x_continuous(breaks=seq(0,2.5,1))
  #     theme(axis.title = element_text(size = 16),
  #           axis.text = element_text(size = 16))
  
  ggsave(paste0("plots/n",nt,"_an_v_num_all_sites_facet_grid.png"))
}

