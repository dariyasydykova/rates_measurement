# R code for generating a balanced tree
library(ape)

setwd("substitution_matrices_in_pheno_models/")
##Generating trees for simulations with 10 sites for JC with equal equilib frequencies. 
##These simulations pair with analytical calculations 
tree <- stree(2, type = "balanced") # generated binary tree w/ num of taxa.
br_len <- seq(0.02,0.50,by=0.02)
for (bl in br_len) {
  tree$edge.length <- rep(bl,nrow(tree$edge)) # sets all branch lengths to bl
  
  bl_str <- sprintf("%.2f",bl)
  f = paste0("trees/n2_bl",bl_str,".tre")
  write.tree(tree, file=f)
}

##Generating trees for simulations with all sites for all empirical matrices 
num_taxa <- c(128,256,512,1024,2048)
br_len <- c(0.0001,0.001,0.01,0.1)
for (num in num_taxa) {
  tree <- stree(num, type = "balanced") # generated binary tree w/ num of taxa.
  for (bl in br_len) {
    tree$edge.length <- rep(bl,nrow(tree$edge)) # sets all branch lengths to bl
  
    if (bl==0.0001) {
      bl_str <- sprintf("%.4f",bl)
    } else bl_str <- bl
    
    f = paste0("trees/n",num,"_bl",bl_str,".tre")
    write.tree(tree, file=f)
  }
}  