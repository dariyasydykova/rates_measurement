library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

setwd("substitution_matrices_in_pheno_models/")
t1 <- read.table("analytically_derived_rates/rates_ten_sites.txt",header=T)
t2<- read.csv("inferred_rates/processed_rates/rates_ten_sites.csv")

t_hyphy <- t2 %>% group_by(time) %>% 
  mutate(rate_mean=mean(rate), rate_norm = rate / rate_mean)

sites_to_plot <- c(1,2,4,5,7,9)
plot_lst <- list()

for (i in sites_to_plot){
  r_an <- filter(t1,site==i+1)
  r_inf <- filter(t_hyphy,site==i) 
  
  p_rates <- ggplot(r_an,aes(x=time)) +
    background_grid("xy")+
    geom_line(aes(y=r_tilde_ms),color="black",size=0.8) + 
    geom_line(aes(y=r_tilde_ms_small_t), color="blue",size=0.8) + 
    stat_summary(data=r_inf,
                 inherit.aes=FALSE,
                 aes(x=time,y=rate_norm),
                 color="red",
                 fun.y = mean,
                 fun.ymin = function(x) mean(x) - sd(x)/sqrt(length(x)), 
                 fun.ymax = function(x) mean(x) + sd(x)/sqrt(length(x)), 
                 geom = "pointrange",
                 size=0.25)+
    xlab("Time") +
    ylab("Relative rate") +
    coord_cartesian(ylim=c(0,2.5),xlim=c(0,1))+
    scale_y_continuous(breaks=seq(0,2.5,0.5),label=c("0","0.5","1.0","1.5","2.0","2.5")) +
    scale_x_continuous(breaks=seq(0,1,0.2),expand = c(0.01, 0),label=c("0","0.2","0.4","0.6","0.8","1.0")) +
    geom_hline(yintercept=1)+
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 15),
          legend.position="none")
  
  if (i==2 | i==4 | i==7 | i==9) p_rates <- p_rates+theme(axis.title.y = element_blank())
  plot_lst[[length(plot_lst)+1]] <- p_rates
}

prow <- plot_grid(plotlist=plot_lst,
                  labels="AUTO",
                  align = 'vh',
                  hjust = -1,
                  ncol=3,
                  nrow=2)

save_plot("plots/rates_true_MutSel_inf_JC.png", prow,
          ncol = 3, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3)

r <- all_r %>% 
  mutate(cor=cor(rate_norm,r_tilde_ms_norm,method="spearman",use="pairwise.complete.obs"),
         rmsd=sqrt(mean((rate_norm - r_tilde_ms_norm)^2)))
