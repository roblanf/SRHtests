#!/usr/bin/Rscript

library(phangorn)
library(plyr)
path = "~/Desktop/IQtree/"

get_dists = function(treefile){
    trees = read.tree(treefile)
    if(length(trees)==3){

        # check for missing tips and drop them if necessary
        m1 = setdiff(trees[[1]]$tip.label, trees[[2]]$tip.label)
        m2 = setdiff(trees[[1]]$tip.label, trees[[3]]$tip.label)
        drop = union(m1, m2)
        
        if(length(drop)>0){
            trees[[1]] = drop.tip(trees[[1]], drop)   
            trees[[2]] = drop.tip(trees[[2]], drop)   
            trees[[3]] = drop.tip(trees[[3]], drop)   
        }
        
        all_bad = path.dist(trees[[1]], trees[[2]])
        all_not = path.dist(trees[[1]], trees[[3]])
        bad_not = path.dist(trees[[2]], trees[[3]])    
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

# get the distances
r = lapply(treefiles, get_dists)
rd = ldply(r, data.frame)

# get the names and tests
t = strsplit(treefiles, '/')
u = lapply(t, function (x) tail(x, n=3))
dataset = unlist(lapply(u, function (x) x[1]))
test = unlist(lapply(u, function (x) x[2]))

rd$dataset = dataset
rd$test = test

write.csv(rd, file = "~/Desktop/IQtree/treedistances.csv")