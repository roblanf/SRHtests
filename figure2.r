library(ggplot2)

d = read.csv("~/Dropbox/Projects_Current/systematic_bias/tables/table2.csv")

# overall p values from binomial tests
ggplot(d, aes(Dataset)) + geom_bar(aes(weight = Bad_Charsets..)) + facet_wrap(~Test, ncol=1) + theme(axis.text.x=element_text(angle=90,hjust=1))
