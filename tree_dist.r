#!/usr/bin/Rscript

library(phangorn)
library(plyr)
path = "/Users/roblanfear/Dropbox/Projects_Current/systematic_bias/processed_data/IQtree/"


random.path.dist = function(N){
  
  a = rtree(N, rooted = FALSE, tip.label = NULL, br = NULL)
  b = rtree(N, rooted = FALSE, tip.label = NULL, br = NULL)

  return(path.dist(a, b))
}

normalised.path.dist = function(t1, t2){
  
    observed = path.dist(t1, t2)
    
    N = length(t1$tip.label)
    
    mean.pd = mean(replicate(n = 1000, expr = random.path.dist(100)))
  
    norm.pd = observed / mean.pd

    return(norm.pd)    
}


get_dists = function(treefile){
    print(treefile)
    print(readLines(treefile))
    trees = read.tree(treefile)
     
    if(length(trees)==3){

        # check for missing tips and drop them if necessary
        m1 = setdiff(trees[[1]]$tip.label, trees[[2]]$tip.label)
	      print(m1)
	      m2 = setdiff(trees[[1]]$tip.label, trees[[3]]$tip.label)
	      print(m2)
	      drop = union(m1, m2)
	      print(drop)        
        
	      if(length(drop)>0){
	        all_bad = all_not = bad_not = NA
	      }else{
        
            all_bad = normalised.path.dist(trees[[1]], trees[[2]])
            all_not = normalised.path.dist(trees[[1]], trees[[3]])
            bad_not = normalised.path.dist(trees[[2]], trees[[3]])    
	      }
    }else{
        all_bad = all_not = bad_not = NA    
    }    
    a = data.frame(t1 = 'all', t2 = 'bad', dist = all_bad)    
    b = data.frame(t1 = 'all', t2 = 'not', dist = all_not)    
    c = data.frame(t1 = 'bad', t2 = 'not', dist = bad_not)    
    return(rbind(a, b, c))    
}

# get the tree files
basename.matches = list.files(path, pattern='trees.nex', recursive=TRUE, full.names=TRUE)
treefiles = grep(pattern='S/trees.nex', basename.matches, value = TRUE)

print(treefiles)

# get the distances
r = lapply(treefiles, get_dists)
rd = ldply(r, data.frame)

# get the names and tests
t = strsplit(treefiles, '/')
u = lapply(t, function (x) tail(x, n=3))
dataset = unlist(lapply(u, function (x) x[1]))
test = unlist(lapply(u, function (x) x[2]))

# we rep by 3 each since each output gives three consequtive lines in the dataframe
rd$dataset = rep(dataset, each = 3)
rd$test = rep(test, each = 3)

write.csv(rd, file = "~/Dropbox/Projects_Current/systematic_bias/tables/tree_distances.csv")
