library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

t <- read.csv("site_num_test/rates/processed_rates/all_rates.csv")
temp1 <- read.table("numerically_derived_rates/raw_rates/site1_rates_132L_A.txt",header=T)
temp2 <- filter(temp1,time=="0.800001")
num_rate <- temp2$r
t$diff <- num_rate-t$rate

p_hyphy_rates <- ggplot() +
    background_grid("xy")+
    geom_line(data=t,aes(time,r),color="red") + 
    xlab("Time") +
    ylab("Rate") +
    scale_y_continuous(breaks=seq(-1,4.5,0.5), limits = c(-1,4.5)) +
    scale_x_continuous(breaks=seq(0,1,0.1),limits = c(0,1.01),expand = c(0.01, 0)) +
    geom_hline(yintercept=1)+
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 12))
  #ggsave(paste0("plots/site",i,"_num_rate_v_hyphy_rates.png"))