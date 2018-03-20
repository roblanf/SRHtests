library(ggplot2)
library(reshape2)
library(cowplot)

d = read.csv("~/Dropbox/Projects_Current/systematic_bias/tables/summary_trees.csv")
t = read.csv("~/Dropbox/Projects_Current/systematic_bias/tables/tree_distances.csv")

# add p values to t
t$d1.t2.pWSH = NA
t$t1.d2.pWSH = NA


for(i in 1:nrow(t)){
    
    print(t[i,])
    
    row = t[i,]
    if(!is.na(row$norm.dist)){
        t1 = row$t1
        t2 = row$t2
        if(t1=='all'){ t1 = "All"}
        if(t2=='all'){ t2 = "All"}
        if(t1=='bad'){ t1 = "Bad"}
        if(t2=='bad'){ t2 = "Bad"}
        if(t1=='not'){ t1 = "Not_Bad"}
        if(t2=='not'){ t2 = "Not_Bad"}
        
        d1t2 = subset(d, partition == t1)
        d1t2 = subset(d1t2, dataset == as.character(row$dataset))
        d1t2 = subset(d1t2, test == as.character(row$test))
        d1t2 = subset(d1t2, tree == t2)
        
        t$d1.t2.pWSH[i] = d1t2$pWSH
    
        t1d2 = subset(d, partition == t2)
        t1d2 = subset(t1d2, dataset == as.character(row$dataset))
        t1d2 = subset(t1d2, test == as.character(row$test))
        t1d2 = subset(t1d2, tree == t1)
        
        t$t1.d2.pWSH[i] = t1d2$pWSH
    }
    
}


# tree distances
t$comparison = paste(t$t1, "vs.", t$t2)
t$comparison = as.character(lapply(t$comparison, function(x) { gsub("bad", "fail", x) }))
t$comparison = as.character(lapply(t$comparison, function(x) { gsub("not", "pass", x) }))

# dt is left dataset rejects right tree. td is right dataset rejects left tree 
t$dt = 0
t$dt[which(t$d1.t2.pWSH<0.05)] = 1

t$td = 0
t$td[which(t$t1.d2.pWSH<0.05)] = 1

ggplot(t, aes(x = dataset, y = norm.dist)) + 
    geom_point(aes(colour = comparison), position = position_dodge(width = 0.5), alpha = 0.75, size = 2, na.rm=TRUE) + 
    geom_point(aes(group = comparison), position = position_dodge(width = 0.5), alpha = 0.75, size = (t$dt-0.5), na.rm=TRUE) + 
    geom_point(aes(group = comparison), position = position_dodge(width = 0.5), alpha = 0.75, size = 4*(t$td-0.5), na.rm=TRUE, shape = 21) + 
    facet_wrap(~test, ncol=1) + 
    theme(axis.text.x=element_text(angle=90,hjust=1)) +
    ylim(0, 2)

t = t[c(2:8)]
t1 = melt(t, id.vars = c("t1", "t2", "norm.dist", "dataset", "test"))
head(t1)

t1$t1 = as.character(lapply(t1$t1, function(x) { gsub("bad", "fail", x) }))
t1$t1 = as.character(lapply(t1$t1, function(x) { gsub("not", "pass", x) }))

t1$t2 = as.character(lapply(t1$t2, function(x) { gsub("bad", "fail", x) }))
t1$t2 = as.character(lapply(t1$t2, function(x) { gsub("not", "pass", x) }))



t1$data = as.character(t1$t1)
t1$tree = as.character(t1$t1)

t1$data[which(t1$variable=="t1.d2.pWSH")] = as.character(t1$t2[which(t1$variable=="t1.d2.pWSH")])
t1$tree[which(t1$variable=="d1.t2.pWSH")] = as.character(t1$t2[which(t1$variable=="d1.t2.pWSH")])

t1$data = paste("alignment:", t1$data)
t1$tree = paste("tree:", t1$tree)


# Figure 4: datasets are rows, trees are columns. Shows that (as expected) the ALL dataset tends to reject the other trees. The 'not' datasets rejects the 'bad' tree too. 
ggplot(subset(t1, test=="MPTIS"), aes(x = value)) + geom_histogram(bins = 11) + 
    facet_grid(data~tree) + 
    geom_vline(xintercept = 0.05, color = 'red', linetype = 'dashed') + 
    xlim(NA, 1)

ggplot(subset(t1, test=="MPTMS"), aes(x = value)) + geom_histogram(bins = 11) + 
    facet_grid(data~tree) + 
    geom_vline(xintercept = 0.05, color = 'red', linetype = 'dashed') + 
    xlim(NA, 1)

ggplot(subset(t1, test=="MPTS"), aes(x = value)) + geom_histogram(bins = 11) + 
    facet_grid(data~tree) + 
    geom_vline(xintercept = 0.05, color = 'red', linetype = 'dashed') + 
    xlim(NA, 1)

