library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

file_lst <- list.files("hyphy/rates/raw_rates",full.names=T)

d <- data.frame()
for (file_name in file_lst) {
  t <- read.table(file_name,header=T)
  
  str <- regexpr("bl\\d+",file_name)[1]
  end <- regexpr("_\\d+_rates.txt",file_name)[1]
  bl <- as.numeric(substr(file_name,str+2,end-1))
  t$branch_len <- rep(bl,length(t$site))
  
  t$site <- c(1:5)
  d <- rbind(d, t)
}

ordered_d <- arrange(d,site)
write.csv(ordered_d,file="hyphy/rates/processed_rates/all_rates.csv",quote=F,row.names=F)


  
