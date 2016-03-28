library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

input_dir <- "site_num_test/rates/raw_rates"
output_dir <- "site_num_test/rates/processed_rates/"

file_lst <- list.files(input_dir,full.names=T)

d <- data.frame()
for (file_name in file_lst) {
  t <- read.table(file_name,header=T)
  
  str <- regexpr("bl\\d+",file_name)[1]
  end <- regexpr("_n\\d+_\\d+_rates.txt",file_name)[1]
  time <- as.numeric(substr(file_name,str+2,end-1))*2
  t$time <- rep(time,length(t$site))
  
  str <- regexpr("_n\\d+",file_name)[1]
  end <- regexpr("_\\d+_rates.txt",file_name)[1]
  site_dupl <- as.numeric(substr(file_name,str+2,end-1))
  t$site_dupl <- rep(site_dupl,length(t$site))
  
  d <- rbind(d, t)
}

ordered_d <- arrange(d,site)
write.csv(ordered_d,file=paste0(output_dir,"all_rates.csv"),quote=F,row.names=F)


  
