library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

setwd("substitution_matrices_in_pheno_models/")
r_file_lst=list.files(path = "analytically_derived_rates/",pattern="rates_with_mu",full.names=T)
t1 <- read.table("analytically_derived_rates/rates_ten_sites.txt",header=T)

d <- data.frame()
for (f in r_file_lst){
  r <- read.table(f,header=T,colClasses = c("numeric","numeric","numeric","character","numeric"))
  if (nrow(d)==0) {
    d <- r
  } else d <- rbind(d,r)
}

sites_to_plot <- c(1,7,4,6)
#sites_to_plot <- unique(t1$site)
plot_lst <- list()

for (i in sites_to_plot){
  # r_mu_A1 <- filter(d,site==i+1,mu_nuc=='A',mu_rate==1)
  # r_mu_A3 <- filter(d,site==i+1,mu_nuc=='A',mu_rate==3)
  # r_mu_A5 <- filter(d,site==i+1,mu_nuc=='A',mu_rate==5)
  # r_mu_T3 <- filter(d,site==i+1,mu_nuc=='T',mu_rate==3)
  # r_mu_T5 <- filter(d,site==i+1,mu_nuc=='T',mu_rate==5)
  # r_mu_C3 <- filter(d,site==i+1,mu_nuc=='C',mu_rate==3)
  # r_mu_C5 <- filter(d,site==i+1,mu_nuc=='C',mu_rate==5)
  # r_mu_G3 <- filter(d,site==i+1,mu_nuc=='G',mu_rate==3)
  # r_mu_G5 <- filter(d,site==i+1,mu_nuc=='G',mu_rate==5)
  r_mu1 <- filter(d,site==i+1,mu_rate==1)
  #r_mu3 <- filter(d,site==i+1,mu_rate==3)
  #r_mu5 <- filter(d,site==i+1,mu_rate==5)
  r_an <- filter(t1,site==i+1)
  
  p_rates <- ggplot() +
    background_grid("xy")+
    # geom_line(data=r_mu_A3,aes(x=time,y=r_tilde),color="blue",size=0.8) +
    # geom_line(data=r_mu_A5,aes(x=time,y=r_tilde),color="blue",size=0.8) +
    # geom_line(data=r_mu_T3,aes(x=time,y=r_tilde),color="green",size=0.8) +
    # geom_line(data=r_mu_T5,aes(x=time,y=r_tilde),color="green",size=0.8) +
    # geom_line(data=r_mu_C3,aes(x=time,y=r_tilde),color="red",size=0.8) +
    # geom_line(data=r_mu_C5,aes(x=time,y=r_tilde),color="red",size=0.8) +
    # geom_line(data=r_mu_G3,aes(x=time,y=r_tilde),color="brown",size=0.8) +
    # geom_line(data=r_mu_G5,aes(x=time,y=r_tilde),color="brown",size=0.8) +
    # geom_line(data=r_mu_A1,aes(x=time,y=r_tilde),color="yellow",size=0.8) +
    geom_line(data=r_mu1,aes(x=time,y=r_tilde),color="red",size=0.8) +
    #geom_line(data=r_mu3,aes(x=time,y=r_tilde,group=mu_nuc),color="green",size=0.8) +
    #geom_line(data=r_mu5,aes(x=time,y=r_tilde,group=mu_nuc),color="blue",size=0.8) +
    geom_line(data=r_an,aes(x=time,y=r_tilde_ms),color="black",size=0.8) +
    xlab("Time") +
    ylab("Relative rate") +
    coord_cartesian(ylim=c(0,3.5),xlim=c(0,2))+
    scale_y_continuous(breaks=seq(0,3.5,0.5)) + #,label=c("0","0.5","1.0","1.5","2.0","2.5")) +
    scale_x_continuous(breaks=seq(0,2,0.5),expand = c(0.01, 0)) + #,label=c("0","0.2","0.4","0.6","0.8","1.0")) +
    #scale_y_continuous(breaks=seq(0.0,1.6,0.1))+ #,label=c("0","0.5","1.0","1.5","2.0","2.5")) +
    #scale_x_continuous(breaks=seq(0,0.04,0.01),expand = c(0.01, 0)) + #,label=c("0","0.2","0.4","0.6","0.8","1.0")) +
    geom_hline(yintercept=1)+
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 15),
          legend.position="none")
  
  if (i==2) p_rates <- p_rates+theme(axis.title.y = element_blank())
  plot_lst[[length(plot_lst)+1]] <- p_rates
}

prow <- plot_grid(plotlist=plot_lst,
                  labels="AUTO",
                  align = 'vh',
                  hjust = -1,
                  ncol=2,
                  nrow=2)

save_plot("plots/rates_true_MutSel_inf_JC_with_mu.png", prow,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3)

