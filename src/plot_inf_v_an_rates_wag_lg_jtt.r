library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

setwd("substitution_matrices_in_pheno_models/")
t1 <- read.table("numerically_derived_rates/all_rates.txt",header=T)
t2 <- read.csv("inferred_rates/processed_rates/all_rates_all_sites.csv")

r_inf <- group_by(t2,time,model,num_taxa) %>% mutate(rate_mean=mean(rate)) %>% mutate(rate_norm = rate / rate_mean) 
r_an <- group_by(t1,site) %>% mutate(r_tilde_ms_norm_small_t=r_tilde_ms_norm[1])

an_r_lst <- unique(r_an$r_tilde_ms_norm_small_t)
r_inf <- mutate(r_inf, an_rate_norm_small_t=rep(an_r_lst,30))

r <- r_inf %>% mutate(cor=cor(rate_norm,an_rate_norm_small_t,method="spearman",use="pairwise.complete.obs"),
                   rmsd=sqrt(mean((rate_norm - an_rate_norm_small_t)^2)),
                   bias=mean(rate_norm - an_rate_norm_small_t)) 

model_lst <- unique(r_inf$model)
for (m in model_lst){
  r_filtered <- r %>% filter(model==m)
  
  p_cor <- ggplot(r_filtered,aes(br_len,cor,colour=factor(num_taxa)))+
    stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x)/sqrt(length(x)), 
                 fun.ymax = function(x) mean(x) + sd(x)/sqrt(length(x)), 
                 geom = "pointrange",
                 size=0.4)+
    scale_x_log10(breaks=c(0.00005,0.0005,0.005,0.05),labels=c("0.00005","0.0005","0.005","0.05")) +
    guides(col = guide_legend(title="Number of Taxa",reverse = TRUE)) +
    stat_summary(fun.y = mean,geom = "line",aes(color=factor(num_taxa)),size=0.6)+
    xlab("Branch Length") +
    ylab("Correlation (spearman)") +
    coord_cartesian(ylim=c(0, 1))+
    scale_y_continuous(breaks=seq(0,1,0.2))
    #+theme(axis.title = element_text(size = 14),
     #     axis.text = element_text(size = 12),
      #    legend.text = element_text(size = 11),
       #   legend.title = element_text(size = 12))  
  
  p_rmsd <- ggplot(r_filtered,aes(br_len,rmsd,colour=factor(num_taxa))) + 
    stat_summary(fun.y = mean,
                 fun.ymin = function(x) mean(x) - sd(x)/sqrt(length(x)), 
                 fun.ymax = function(x) mean(x) + sd(x)/sqrt(length(x)), 
                 geom = "pointrange",
                 size=0.4)+
    scale_x_log10(breaks=c(0.00005,0.0005,0.005,0.05),labels=c("0.00005","0.0005","0.005","0.05")) +
    guides(col = guide_legend(title="Number of Taxa",reverse = TRUE)) +
    stat_summary(fun.y = mean,geom = "line",aes(color=factor(num_taxa)),size=0.6)+
    xlab("Branch Length") +
    ylab("RMSD") +
    coord_cartesian(ylim=c(0,30))
    #scale_y_continuous(breaks=seq(0,15,2.5))
    #theme(axis.title = element_text(size = 14),
    #      axis.text = element_text(size = 12),
    #      legend.text = element_text(size = 11),
    #      legend.title = element_text(size = 12))
  
  p_bias <- ggplot(r_filtered,aes(br_len,bias,colour=factor(num_taxa))) + 
    stat_summary(fun.y = mean,
                 fun.ymin = function(x) mean(x) - sd(x)/sqrt(length(x)), 
                 fun.ymax = function(x) mean(x) + sd(x)/sqrt(length(x)), 
                 geom = "pointrange",
                 size=0.4)+
    scale_x_log10(breaks=c(0.00005,0.0005,0.005,0.05),labels=c("0.00005","0.0005","0.005","0.05")) +
    guides(col = guide_legend(title="Number of Taxa",reverse = TRUE)) +
    stat_summary(fun.y = mean,geom = "line",aes(color=factor(num_taxa)),size=0.6)+
    xlab("Branch Length") +
    ylab("Bias") +
    geom_hline(yintercept=0,color="grey")
    #coord_cartesian(ylim=c(-2e-16,5e-16),xlim=c(0.0023,0.66))+
    #scale_y_continuous(breaks=seq(-2e-16,5e-16,1e-16),labels=fancy_scientific)+ 
    #theme(axis.title = element_text(size = 14),
     #     axis.text = element_text(size = 12),
      #    legend.text = element_text(size = 11),
       #   legend.title = element_text(size = 12))
  
  grobs <- ggplotGrob(p_bias)$grobs
  legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
  
  p_row <- plot_grid(p_cor+theme(legend.position="none"),
                           p_rmsd+theme(legend.position="none"),
                           p_bias+theme(legend.position="none"),
                           labels=c("A","B","C"),
                           align = 'vh',
                           hjust = -1,
                           ncol=3,
                           nrow=1)
  final_p <- plot_grid( p_row, legend, rel_widths = c(2, .3))
  save_plot(paste0("plots/",m,"_cor_rmsd_bias.png"), final_p,
            ncol = 3, # we're saving a grid plot of 2 columns
            nrow = 1, # and 2 rows
            # each individual subplot should have an aspect ratio of 1.3
            base_aspect_ratio = 1.3)
  
}