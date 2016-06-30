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
    geom_line(data=r,aes(time,r_tilde_jc_norm),color="red") + 
    geom_line(data=r,aes(time,r_tilde_wag_norm),color="blue") + 
    xlab("Time") +
    ylab("Rate") +
    scale_y_continuous(breaks=seq(0,4.0,0.5), limits = c(0,4.0)) +
    scale_x_continuous(breaks=seq(0,1,0.2),limits = c(0,1.01),expand = c(0.01, 0)) +
    geom_hline(yintercept=1)+
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 16))
  print(p_rates)
  ggsave(paste0("plots/site",i,"_empirical_q_matrices.png"))
}