# R code for generating a balanced tree
library(ape)

tree <- stree(2, type = "balanced") # generated binary tree w/ num of taxa.
br_len <- seq(0.02,0.52,by=0.02)
for (bl in br_len) {
  tree$edge.length <- rep(bl,nrow(tree$edge)) # sets all branch lengths to bl
  
  bl_str <- sprintf("%.2f",bl)
  f = paste0("trees/n2_bl",bl_str,".tre")
  write.tree(tree, file=f)
}
