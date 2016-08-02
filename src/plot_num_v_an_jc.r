library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

t1 <- read.table("numerically_derived_rates/all_rates.txt",header=T)
t2<- read.csv("hyphy/rates/processed_rates/all_rates_ten_sites.csv")

t_an <- group_by(t1[1:1000,],time) %>% mutate(rate_mean=mean(r_tilde_ms))
t_hyphy <- group_by(t2,time) %>% mutate(rate_mean=mean(rate))

site_lst_t1 <- unique(t_an$site)
site_lst_t2 <- unique(t2$site)

for (i in site_lst_t2){
  r_an <- filter(t_an,site==site_lst_t1[i])
  r_inf <- filter(t_hyphy,site==i) 
  r_an <- mutate(r_an, r_tilde_ms_norm_ten_sites = r_tilde_ms / rate_mean)
  r_inf <- mutate(r_inf, rate_norm = rate / rate_mean)
  
  p_rates <- ggplot() +
    background_grid("xy")+
    geom_line(data=r_an,aes(time,r_tilde_ms_norm_ten_sites),color="black") + 
    stat_summary(data=r_inf,
                 aes(x=time,y=rate_norm),
                 color="red",
                 fun.y = mean,
                 fun.ymin = function(x) mean(x) - sd(x)/sqrt(length(x)), 
                 fun.ymax = function(x) mean(x) + sd(x)/sqrt(length(x)), 
                 geom = "pointrange",
                 size=0.3)+
    xlab("Time") +
    ylab("Rate") +
    coord_cartesian(ylim=c(0,2.5),xlim=c(0,1))+
    scale_y_continuous(breaks=seq(0,2.5,0.5)) +
    scale_x_continuous(breaks=seq(0,1,0.2),expand = c(0.01, 0)) +
    geom_hline(yintercept=1)+
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 16))
  print(p_rates)
  ggsave(paste0("plots/site",i,"_ms_an_v_num_JC_equalf.png"))
}