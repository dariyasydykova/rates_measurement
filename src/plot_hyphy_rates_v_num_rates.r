library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

hyphy_rates_file <- read.csv("hyphy/rates/processed_rates/all_rates.csv")
num_rates_file<-read.table("numerically_derived_rates/processed_rates/all_numerical_rates.txt",header=T)

for (i in c(1:5)) {
  t1 <- filter(hyphy_rates_file,site==i)
  t2 <- filter(num_rates_file,site==i)
  
  mean_lst <- c()
  for (b in unique(t1$branch_len)){
    a <- filter(t1,branch_len==b)
    mean <- mean(a$rate)
    mean_lst <- c(mean_lst,rep(mean,length(a$site)))
  }
  t1$mean <- mean_lst

  p_hyphy_rates <- ggplot() +
    background_grid("xy")+
    geom_line(data=t2,aes(time,r),color="red") + 
    geom_line(data=t2,aes(time,r_large_t),color="grey")+
    geom_line(data=t2,aes(time,r_small_t),color="grey")+
    geom_boxplot(data=t1,aes(x=time,group=time,y=rate),color="blue", width=.01) +
    xlab("Time") +
    ylab("Rate") +
    scale_y_continuous(breaks=seq(-1,4.0,0.5), limits = c(-1,4.0)) +
    scale_x_continuous(breaks=seq(0,1,0.1),limits = c(0,1.01),expand = c(0.01, 0)) +
    geom_hline(yintercept=1)+
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 12))
  print(p_hyphy_rates)
  ggsave(paste0("plots/site",i,"_num_rate_v_hyphy_rates.png"))
}