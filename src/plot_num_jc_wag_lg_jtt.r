library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

t <- read.table("numerically_derived_rates/all_rates.txt",header=T)

for (i in c(1:10)){
  r <- filter(t,site==i)

  #true_r <- r$true_r[1]/r$mean_true_r[1]
  p_rates <- ggplot() +
    background_grid("xy")+
    geom_line(data=r,aes(time,r_tilde_jc_norm),color="black") + 
    xlab("Time") +
    ylab("Rate") +
    coord_cartesian(ylim=c(0, 2.5),xlim=c(0,1))+
    scale_y_continuous(breaks=seq(0,2.5,0.5)) +
    scale_x_continuous(breaks=seq(0,1,0.2),expand = c(0.01, 0)) +
    geom_hline(yintercept=1)+
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 16))
  print(p_rates)
  #ggsave(paste0("plots/site",i,"_true_model_JC.png"))
}