# R code for generating a balanced tree
library(ape)

tree <- stree(2, type = "balanced") # generated binary tree w/ num of taxa.
br_len <- seq(0.02,0.5,by=0.02)
for (bl in br_len) {
  tree$edge.length <- rep(bl,nrow(tree$edge)) # sets all branch lengths to bl
  
  f = paste0("trees/n2_bl",bl,".tre")
  write.tree(tree, file=f)
}
##tree$edge.length<-runif(n=nrow(tree$edge),min=0,max=0.1) # draws branch lengths from uniform distribution with min0, max2
