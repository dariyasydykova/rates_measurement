library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

t1 <- read.table("numerically_derived_rates/all_rates.txt",header=T)
t2<- read.csv("hyphy/rates/processed_rates/all_rates.csv")
t_hyphy <- group_by(t2,time,model) %>% mutate(rate_mean=mean(rate))

colfunc <- colorRampPalette("set1")
for (i in c(1:10)){
  r_an <- filter(t1,site==i)
  r_inf <- filter(t_hyphy,site==i) 
  r_inf$rate_norm <- r_inf$rate/r_inf$rate_mean
  # offset <- c(-0.005, -0.0025, 0, 0.0025, 0.005)
#    r_jc <- filter(t_hyphy,site==i,model=="JC") 
#    r_jc_ef <- filter(t_hyphy,site==i,model=="JC_equalf")
#    r_wag <- filter(t_hyphy,site==i,model=="WAG")
#    r_jtt <- filter(t_hyphy,site==i,model=="JTT")
#    r_lg <- filter(t_hyphy,site==i,model=="LG")
  
  p_rates <- ggplot() +
    background_grid("xy")+
    geom_line(data=r_an,aes(time,r_tilde_ms_norm),color="black") + 
#    stat_summary(data=r_inf,fun.y = mean,geom = "line",aes(x=time,y=rate/rate_mean,color=factor(model)),size=0.6)+
#      geom_boxplot(data=r_jc_ef,aes(x=time,y=rate/rate_mean,group=time),color="grey", width=.01,alpha = 0.5) +
#      geom_boxplot(data=r_jc,aes(x=time,y=rate/rate_mean,group=time),color="grey", width=.01,alpha = 0.5) +
#      geom_boxplot(data=r_wag,aes(x=time,y=rate/rate_mean,group=time),color="grey", width=.01,alpha = 0.5) +
#      geom_boxplot(data=r_jtt,aes(x=time,y=rate/rate_mean,group=time),color="grey", width=.01,alpha = 0.5) +
#      geom_boxplot(data=r_lg,aes(x=time,y=rate/rate_mean,group=time),color="grey", width=.01,alpha = 0.5) +
    stat_summary(data=r_inf,
                 aes(x=time,y=rate_norm,color=factor(model)),
                 fun.y = mean,
                 fun.ymin = function(x) mean(x) - sd(x)/sqrt(length(x)), 
                 fun.ymax = function(x) mean(x) + sd(x)/sqrt(length(x)), 
                 geom = "pointrange",
                 size=0.3)+
    xlab("Time") +
    ylab("Rate") +
    coord_cartesian(ylim=c(0, 2.5),xlim=c(0,1))+
    scale_y_continuous(breaks=seq(0,2.5,0.5)) +
    scale_x_continuous(breaks=seq(0,1,0.2),expand = c(0.01, 0)) +
    geom_hline(yintercept=1)+
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 16))+
    scale_color_discrete(guide = guide_legend(title = "Model"),labels = c("JC","JC equal freq","JTT","LG","WAG"))
  print(p_rates)
  ggsave(paste0("plots/site",i,"_ms_an_v_num.png"))
}