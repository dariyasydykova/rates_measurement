library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

setwd("substitution_matrices_in_pheno_models/")
input_dir <- "inferred_rates/raw_rates/all_sites"
output_file <- "inferred_rates/processed_rates/rates_all_sites.csv"

file_lst <- list.files(input_dir,full.names=T)

d <- data.frame()
for (file_name in file_lst) {
  t <- read.table(file_name,header=T)
  
  str <- regexpr("n\\d+",file_name)[1]
  end <- regexpr("_bl\\d+",file_name)[1]
  num <- as.numeric(substr(file_name,str+1,end-1))
  t$num_taxa <- rep(num,length(t$site)) 
  
  str <- regexpr("bl\\d+",file_name)[1]
  end <- regexpr("_\\d+_[[:upper:]]+_rates.txt",file_name)[1]
  time <- as.numeric(substr(file_name,str+2,end-1))*logb(num,2)*2
  t$time <- rep(time,length(t$site))
  br_len <- as.numeric(substr(file_name,str+2,end-1))
  t$br_len <- rep(br_len,length(t$site))
  
  if (grepl("site_dupl",file_name)){
    str <- regexpr("site_dupl_\\d+",file_name)[1]
    end <- regexpr("_\\d+_rates.txt",file_name)[1]
    s <- as.numeric(substr(file_name,str+10,end-1))
    t$site_dupl <- rep(s,length(t$site))
    
    str <- regexpr("_\\d+_rates.txt",file_name)[1]
    end <- regexpr("_rates.txt",file_name)[1]
    sim_num <- as.numeric(substr(file_name,str+1,end-1))
    t$rep <- rep(sim_num,length(t$site))
    
   } else {
    str <- regexpr("[[:upper:]]+",file_name)[1]
    end <- regexpr("_rates.txt",file_name)[1]
    model <- substr(file_name,str,end-1)
    t$model <- rep(model,length(t$site))
    
    str <- regexpr("_\\d+_\\w+_rates.txt",file_name)[1]
    end <- regexpr("[[:upper:]]+",file_name)[1]
    sim_num <- as.numeric(substr(file_name,str+1,end-2))
    t$rep <- rep(sim_num,length(t$site))
  }
  t$site <- c(1:length(t$site))
  d <- rbind(d, t)
}

write.csv(d,file=output_file,quote=F,row.names=F)



