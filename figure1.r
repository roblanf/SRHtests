library(ggplot2)
library(plyr)

d = read.csv("~/Dropbox/Projects_Current/systematic_bias/tables/table1.csv")
d = d[,-c(1,2, 10, 11, 12)]
names(d) = c("charset", "test", "fail", "pass", "NA", "p_binomial", "dataset", "sites")

# overall p values from binomial tests
ggplot(d, aes(x=p_binomial)) + geom_histogram(bins=40) + facet_wrap(~test, ncol=1) + geom_vline(xintercept = 0.05, colour = 'red')
