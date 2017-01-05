library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

input_dir <- "analytically_derived_rates/raw_rates"
output_dir <- "analytically_derived_rates/processed_rates/"

file_lst <- list.files(input_dir,full.names=T)

d <- data.frame()
for (rate_file in file_lst){
  r <- read.table(rate_file,header=T)
  d <- rbind(d,r)
}

write.table(d,paste0(output_dir,"all_numerical_rates.txt"),sep="\t",row.names=F,quote=F)