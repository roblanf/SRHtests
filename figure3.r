library(ggplot2)

d = read.csv("~/Dropbox/Projects_Current/systematic_bias/tables/table3.csv")
t = read.csv("~/Dropbox/Projects_Current/systematic_bias/tables/tree_distances.csv")

# add p values to t
t$d1.t2.pSH = NA
t$t1.d2.pSH = NA


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
        
        d1t2 = subset(d, Partition == t1)
        d1t2 = subset(d1t2, Dataset == as.character(row$dataset))
        d1t2 = subset(d1t2, Test == as.character(row$test))
        d1t2 = subset(d1t2, Tree == t2)
        
        t$d1.t2.pSH[i] = d1t2$pSH
    
        
        t1d2 = subset(d, Partition == t2)
        t1d2 = subset(t1d2, Dataset == as.character(row$dataset))
        t1d2 = subset(t1d2, Test == as.character(row$test))
        t1d2 = subset(t1d2, Tree == t1)
        
        t$t1.d2.pSH[i] = t1d2$pSH
    }
    
}


# tree distances
t$comparison = paste(t$t1, "vs.", t$t2)
ggplot(t, aes(x=dataset, y = norm.dist)) + 
    geom_point(aes(colour = comparison), position = position_dodge(width = 0.5), alpha = 0.75, size = 2, na.rm=TRUE) + 
    facet_wrap(~test, ncol=1) + 
    theme(axis.text.x=element_text(angle=90,hjust=1)) + scale_y_log10()


# TO DO: combine d and t in such a way that we can colour each bar of the plot by a factor that indicates which of four cases are true
# t1 does not reject t2, and vice versa
# t1 rejects t2, but not vice versa
# t2 rejects t1, but not vice versa
# both trees reject each other
# all according to pSH from d
# make this the point shapes. Different shapes for different combinations.