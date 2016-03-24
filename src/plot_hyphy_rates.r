library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

hyphy_rates_file <- read.csv("hyphy/rates/processed_rates/all_rates.csv")
rate_file_names <- list.files("numerically_derived_rates",full.names=T)

for (i in c(1:5)) {
  t1 <- filter(hyphy_rates_file,site==i)
  t2 <- read.table(rate_file_names[i],header=T)
  
  mean_lst <- c()
  for (b in unique(t1$branch_len)){
    a <- filter(t1,branch_len==b)
    mean <- mean(a$rate)
    mean_lst <- c(mean_lst,rep(mean,length(a$site)))
  }
  t1$mean <- mean_lst

  p_hyphy_rates <- ggplot() +
    geom_line(data=t2,aes(time,r),color="red") + 
    geom_point(data=t1,aes(branch_len,rate),color="blue") +
    geom_line(data=t1,aes(branch_len,mean),color="blue") +
    xlab("Time") +
    ylab("Rate") +
    scale_y_continuous(breaks=seq(-1,6,0.5), limits = c(-1,6)) +
    scale_x_continuous(breaks=seq(0,2,0.5),limits = c(0,2),expand = c(0.01, 0)) +
    geom_hline(yintercept=1)+
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 12))
  ggsave(paste0("plots/site",i,"num_rate_v_hyphy_rates.png"))
}