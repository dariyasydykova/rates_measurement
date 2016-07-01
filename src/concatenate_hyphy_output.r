library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

input_dir <- "hyphy/rates/raw_rates"
output_dir <- "hyphy/rates/processed_rates/"

file_lst <- list.files(input_dir,full.names=T)

d <- data.frame()
for (file_name in file_lst) {
  t <- read.table(file_name,header=T)
  
  str <- regexpr("bl\\d+",file_name)[1]
  end <- regexpr("_\\d+_\\w+_rates.txt",file_name)[1]
  time <- as.numeric(substr(file_name,str+2,end-1))*2
  t$time <- rep(time,length(t$site))
  
  str <- regexpr("[[:upper:]]+",file_name)[1]
  end <- regexpr("_rates.txt",file_name)[1]
  model <- substr(file_name,str,end-1)
  t$model <- rep(model,length(t$site))
  
  t$site <- c(1:10)
  d <- rbind(d, t)
}

write.csv(d,file=paste0(output_dir,"all_rates.csv"),quote=F,row.names=F)


  
