library(ggplot2)

d = read.csv("~/Dropbox/Projects_Current/systematic_bias/tables/table3.csv")

d1 = subset(d, Tree=='Not_Bad')
d1 = subset(d1, Partition=='All')

# overall p values from binomial tests
ggplot(d1, aes(Dataset)) + geom_bar(aes(weight = pSH)) + facet_wrap(~Test, ncol=1) + theme(axis.text.x=element_text(angle=90,hjust=1))
