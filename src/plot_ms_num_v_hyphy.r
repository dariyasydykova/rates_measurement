library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

t1 <- read.table("numerically_derived_rates/all_rates.txt",header=T)
t2<- read.csv("hyphy/rates/processed_rates/all_rates.csv")
t_hyphy <- group_by(t2,time,model) %>% mutate(rate_mean=mean(rate))

for (i in c(1:10)){
  r_an <- filter(t1,site==i) 
  r_jc <- filter(t_hyphy,site==i,model=="JC")
  r_wag <- filter(t_hyphy,site==i,model=="WAG")
  r_jtt <- filter(t_hyphy,site==i,model=="JTT")
  r_lg <- filter(t_hyphy,site==i,model=="LG")

  p_rates <- ggplot() +
    background_grid("xy")+
    geom_line(data=r_an,aes(time,r_tilde_ms_norm),color="black") + 
    geom_boxplot(data=r_jc,aes(x=time,y=rate/rate_mean,group=time),color="red", width=.01) +
    geom_boxplot(data=r_wag,aes(x=time,y=rate/rate_mean,group=time),color="blue", width=.01) +
    geom_boxplot(data=r_jtt,aes(x=time,y=rate/rate_mean,group=time),color="green", width=.01) +
    geom_boxplot(data=r_lg,aes(x=time,y=rate/rate_mean,group=time),color="yellow", width=.01) +
    xlab("Time") +
    ylab("Rate") +
    scale_y_continuous(breaks=seq(0,4.0,0.5), limits = c(0,4.0)) +
    scale_x_continuous(breaks=seq(0,1,0.2),limits = c(0,1.01),expand = c(0.01, 0)) +
    geom_hline(yintercept=1)+
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 16))
  print(p_rates)
  ggsave(paste0("plots/site",i,"_ms_an_v_num.png"))
}