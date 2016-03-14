library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

rate_file_names <- list.files("numerically_derived_rates",full.names=T)

for (rate_table_name in rate_file_names) {
  t <- read.csv(rate_table_name,sep="\t")
  
  p_rates <- ggplot(t,aes(time,r))+
    geom_line(color="#FF0000")+
    geom_line(aes(y=r_large_t),color="#3399FF")+
    geom_line(aes(y=r_small_t),color="#33CC00")+
    scale_y_continuous(breaks=seq(-1,4,0.5), limits = c(-1,4))+
    scale_x_continuous(breaks=seq(0,2,0.5), limits = c(0,2),expand = c(0, 0))+
    geom_hline(yintercept=1)+
    xlab("Time") +
    ylab("Site-specific rate") 
  
  plot_file_name <- gsub(".txt",".png",substr(rate_table_name,27,90))
  ggsave(paste0("plots/",plot_file_name))
} 
#+
 # geom_line(aes(y=r_small_t,color="red"))+
  #geom_line(aes(y=r_large_t,color="blue"))

            
