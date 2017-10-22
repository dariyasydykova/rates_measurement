library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(readr)

####Figure: Analytically derived rate and inferred rate when true model is MutSel and inference model is JC#####
r_an <- read_csv("../analytical_rates/ten_sites_aa.csv")
r_inf<- read_csv("../inferred_rates/processed_rates/rates_ten_sites_aa.csv")

#normalize inferred rates
r_inf %>% group_by(time,rep) %>% mutate(rate_norm = rate / mean(rate)) -> r_norm_inf

#sites to plot
sites_to_plot <- c(1,2,4,5,7,9)
plot_lst <- list()

for (i in sites_to_plot){
  r_an_filtered <- filter(r_an,site==i+1)
  r_inf_filtered <- filter(r_norm_inf,site==i) 
  
  p_rates <- ggplot(r_an_filtered,aes(x=time)) +
    background_grid("xy")+
    geom_line(aes(y=r_tilde),color="black",size=1.2) + 
    geom_line(aes(y=r_tilde_small_t), color="royalblue1",size=1.2) + 
    geom_line(aes(y=r_tilde_large_t), color="green3",size=1.2) + 
    stat_summary(data= r_inf_filtered,
                 inherit.aes=FALSE,
                 aes(x=time,y=rate_norm),
                 color="orangered",
                 fun.y = mean,
                 fun.ymin = function(x) mean(x) - sd(x)/sqrt(length(x)), 
                 fun.ymax = function(x) mean(x) + sd(x)/sqrt(length(x)), 
                 geom = "pointrange",
                 size=0.5)+
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

save_plot("../plots/rates_true_MutSel_inf_JC.png", prow,
          ncol = 3, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3)

####Figure: Accuracy of the inferred rate#####
r_inf <- read.csv( "../inferred_rates/processed_rates/rates_site_dupl.csv")
r_an <- read_csv("../analytical_rates/all_sites_aa.csv")

true_r <- r_an %>%
  filter(site==4,time==0.480002)

r <- r_inf %>%
  group_by(site_dupl, rep) %>% 
  mutate(rate_norm=rate/mean(rate)) %>%
  filter(site==3) 

p <- ggplot(r,aes(site_dupl,rate_norm)) +
  geom_point(size=0.9,alpha=0.8) + 
  geom_hline(aes(yintercept=true_r$r_tilde),color="orangered",size=1.2)+
  ylab("Relative rate") +
  xlab("Site duplicates") +
  coord_cartesian(ylim=c(0.001,1000),xlim=c(10,100000))+
  scale_y_log10(breaks=c(0.001,0.01,0.1,1,10,100,1000),label=c("0.001","0.01","0.1","1","10","100","1,000")) +
  scale_x_log10(breaks=c(10,100,1000,10000,100000),label=c("10","100","1,000","10,000","100,000"))+
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 10))

save_plot("../plots/inf_rate_accuracy_v_site_dupl.png",plot=p)

####Figure: Analytically derived rate and the true rate when true model and inference model is JC#####
r_an <- read_csv("../analytical_rates/ten_sites_aa_true_JC.csv")

sites_to_plot <- c(1,4,5,6,7,9)
plot_lst <- list()

for (i in sites_to_plot){
  r_an_filtered <- filter(r_an,site==i)
  
  p_rates <- ggplot(r_an_filtered) +
    background_grid("xy")+
    geom_line(aes(time,true_r),color="black",size=1.2) + 
    geom_line(aes(time,r_tilde+0.03),color="dodgerblue",size=1.2) + 
    xlab("Time") +
    ylab("Relative rate") +
    coord_cartesian(ylim=c(0,2.5),xlim=c(0,1))+
    scale_y_continuous(breaks=seq(0,2.5,0.5),label=c("0","0.5","1.0","1.5","2.0","2.5")) +
    scale_x_continuous(breaks=seq(0,1,0.2),expand = c(0.01, 0),label=c("0","0.2","0.4","0.6","0.8","1.0")) +
    geom_hline(yintercept=1)+
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 15),
          legend.position="none")
  
  if (i==4 | i==5 | i==7 | i==9) p_rates <- p_rates+theme(axis.title.y = element_blank())
  plot_lst[[length(plot_lst)+1]] <- p_rates
}

prow <- plot_grid(plotlist=plot_lst,
                  labels="AUTO",
                  align = 'vh',
                  hjust = -1,
                  ncol=3,
                  nrow=2)

save_plot("../plots/rates_true_JC_inf_JC.png", prow,
          ncol = 3, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3)

####Figure: a comparison of rates inferred with Jukes-Cantor-like matrix with observed frequencies and equilibrium frequencies #####
r_an <- read_csv("../analytical_rates/all_sites_aa.csv")
r_inf <- read_csv("../inferred_rates/processed_rates/rates_all_sites.csv")

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
  geom_point(size=1.2,alpha=0.8) + 
  geom_abline(color="red")+
  xlab("Inferred rate (equal frequencies)") +
  ylab("Inferred rate (observed frequencies)") +
  facet_wrap( ~ time_label) +
  coord_cartesian(ylim=c(0,2.1), xlim=c(0,2.1))+
  scale_y_continuous(breaks=seq(0,2.2,1)) +
  scale_x_continuous(breaks=seq(0,2.2,1)) +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))+
  panel_border()


save_plot("../plots/rates_inf_diff_JC.png", p,
          ncol = 1, # we're saving a grid plot of 2 columns
          nrow = 1, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3)

####Figure: a comparison of rates inferred with JC, JTT, WAG, LG#####

r_an <- read_csv("../analytical_rates/all_sites_aa.csv")
r_inf <- read_csv("../inferred_rates/processed_rates/rates_all_sites.csv")

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
  filter(num_taxa==512,model!="JC_equalf") %>%
  group_by(site,model,time) %>%
  summarise(inf_rate_norm_mean=mean(inf_rate_norm)) %>%
  left_join(r_an_temp,by='site') %>%
  left_join(d_label)

p <- ggplot(r,aes(r_tilde_small_t,inf_rate_norm_mean)) +
  #background_grid("xy")+
  geom_point(size=1.2,alpha=0.8) + 
  geom_abline(color="red")+
  ylab("Mean inferred rate") +
  xlab("Analytically derived rate") +
  #xlab(expression(paste("Analytically Derived Rates (", hat(r)^(k), "for small t)"))) +
  facet_grid(model ~ time_label) +
  coord_cartesian(ylim=c(0,2.5), xlim=c(0,2.5))+
  scale_y_continuous(breaks=seq(0,2,1)) +
  scale_x_continuous(breaks=seq(0,2,1)) +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))+
  panel_border()

save_plot("../plots/inf_v_an_rates_all_matrices.png", p,
          ncol = 1, # we're saving a grid plot of 2 columns
          nrow = 1, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_height=8,
          base_width=8)

####Figure: Analytically derived rate and inferred rate when true model is MutSel and inference model is JC#####
r_an <- read_csv("../analytical_rates/ten_sites_codon.csv")
r_inf<- read_csv("../inferred_rates/processed_rates/rates_translated.csv")

#normalize inferred rates
r_inf %>% group_by(time,rep) %>% mutate(rate_norm = rate / mean(rate)) -> r_norm_inf

#sites to plot
sites_to_plot <- c(1,2,4,5,7,9)
plot_lst <- list()

for (i in sites_to_plot){
  r_an_filtered <- filter(r_an,site==i+1)
  r_inf_filtered <- filter(r_norm_inf,site==i) 
  
  p_rates <- ggplot(r_an_filtered,aes(x=time)) +
    background_grid("xy")+
    geom_line(aes(y=r_tilde),color="black",size=1.2) + 
    geom_line(aes(y=r_tilde_small_t), color="royalblue1",size=1.2) + 
    geom_line(aes(y=r_tilde_large_t), color="green3",size=1.2) + 
    stat_summary(data= r_inf_filtered,
                 inherit.aes=FALSE,
                 aes(x=time,y=rate_norm),
                 color="orangered",
                 fun.y = mean,
                 fun.ymin = function(x) mean(x) - sd(x)/sqrt(length(x)), 
                 fun.ymax = function(x) mean(x) + sd(x)/sqrt(length(x)), 
                 geom = "pointrange",
                 size=0.5)+
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

save_plot("../plots/rates_true_codon_MutSel_inf_JC.png", prow,
          ncol = 3, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3)


