library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

t1 <- read.table("numerically_derived_rates/all_rates.txt",header=T)
t2<- read.csv("hyphy/rates/processed_rates/all_rates.csv")
t_hyphy <- group_by(t2,time,model) %>% mutate(rate_mean=mean(rate))
t_hyphy$rate_norm <- t_hyphy$rate/t_hyphy$rate_mean

time_lst <- unique(t_hyphy$time)
t1$temp_time <- t1$time-0.000002
d_temp <- data.frame()
for (t in time_lst) {
  temp <- filter(t1,temp_time==t)
  d_temp <- rbind(d_temp, temp)
}

ind <- match(t_hyphy$time,d_temp$temp_time)
t_hyphy$true_rate <- d_temp$r_tilde_ms_norm[ind]

##test <-  t_hyphy %>% mutate(cor=cor(rate_norm,true_rate,method="spearman"))

r <- filter(t_hyphy,time==0.04 | time==0.2 | time==0.6 | time==0.8 | time==1)
p <- ggplot(r,aes(true_rate,rate_norm)) +
  geom_point()+
  facet_grid(model ~ time) +
  xlab("True rate") +
  ylab("Inferred rate")+ 
  coord_cartesian(ylim=c(0,5)) ##,xlim=c(-1,3))
#theme(axis.title = element_text(size = 16),
#     axis.text = element_text(size = 16))
print(p)
#ggsave(paste0("plots/site",i,"_ms_an_v_num.png"))
