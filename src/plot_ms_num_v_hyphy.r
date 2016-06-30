library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

t1 <- read.table("test_numpy.txt",header=T)
t2<- read.csv("hyphy/rates/processed_rates/all_rates.csv")
t_hyphy <- group_by(t2,time) %>% mutate(rate_mean=mean(rate))

for (i in c(1:10)){
  r <- filter(t1,site==i)
  r_hyphy <- filter(t_hyphy,site==i)
  
  #true_r <- r$true_r[1]/r$mean_true_r[1]
  p_rates <- ggplot() +
    background_grid("xy")+
    geom_line(data=r,aes(time,r_tilde_ms_norm),color="red") + 
    geom_boxplot(data=r_hyphy,aes(x=time,group=time,y=rate/rate_mean),color="blue", width=.01) +
    xlab("Time") +
    ylab("Rate") +
    scale_y_continuous(breaks=seq(0,4.0,0.5), limits = c(0,4.0)) +
    scale_x_continuous(breaks=seq(0,1,0.2),limits = c(0,1.01),expand = c(0.01, 0)) +
    geom_hline(yintercept=1)+
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 16))
  print(p_rates)
  ggsave(paste0("plots/site",i,"_ms_num_v_hyphy.png"))
}