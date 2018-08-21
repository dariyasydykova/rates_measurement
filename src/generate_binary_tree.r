# R code for generating a balanced tree
library(ape)

##Generating trees for simulations with 10 sites for JC with equal equilib frequencies. 
##These simulations pair with analytical calculations 
tree <- stree(2, type = "balanced") # generated binary tree w/ num of taxa.
br_len <- seq(0.02,0.50,by=0.02)
for (bl in br_len) {
  tree$edge.length <- rep(bl,nrow(tree$edge)) # sets all branch lengths to bl
  
  bl_str <- sprintf("%.2f",bl)
  f = paste0("../trees/n2_bl",bl_str,".tre")
  write.tree(tree, file=f)
}

##Generating trees for simulations with all sites for all empirical matrices 
num_taxa <- c(64,128,256,512)
br_len <- c(0.00005,0.0005,0.005,0.05)
for (num in num_taxa) {
  tree <- stree(num, type = "balanced") # generated binary tree w/ num of taxa.
  for (bl in br_len) {
    tree$edge.length <- rep(bl,nrow(tree$edge)) # sets all branch lengths to bl
  
    if (bl==0.00005) {
     ≈
    } else if(bl==0.0005) {
      bl_str <- sprintf("%.4f",bl)
    } else bl_str <- bl
    
    f = paste0("../trees/n",num,"_bl",bl_str,".tre")
    write.tree(tree, file=f)
  }
}  

##Generating trees for simulations with small time 't' to test arbitrary QM
num <- 2
br_len <- c(0.005,0.01)

tree <- stree(num, type = "balanced") # generated binary tree w/ num of taxa.
for (bl in br_len) {
  tree$edge.length <- rep(bl,nrow(tree$edge)) # sets all branch lengths to bl
  
  if (bl==0.005) {
    bl_str <- sprintf("%.3f",bl)
  } else bl_str <- sprintf("%.2f",bl)
  f = paste0("../trees/n2_bl",bl_str,".tre")
  write.tree(tree, file=f)
  }