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
    geom_boxplot(data=t1,aes(x=branch_len,group=branch_len,y=rate),color="blue", width=.01) +
    xlab("Time") +
    ylab("Rate") +
    scale_y_continuous(breaks=seq(-1,4.5,0.5), limits = c(-1,4.5)) +
    scale_x_continuous(breaks=seq(0,0.5,0.1),limits = c(0,0.5),expand = c(0.01, 0)) +
    geom_hline(yintercept=1)+
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 12))
  ggsave(paste0("plots/site",i,"_num_rate_v_hyphy_rates.png"))
}