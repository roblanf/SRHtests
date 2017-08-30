library(ggplot2)

d = read.csv("~/Dropbox/Projects_Current/systematic_bias/tables/table3.csv")
t = read.csv("~/Dropbox/Projects_Current/systematic_bias/tables/tree_distances.csv")


d1 = subset(d, Tree=='Not_Bad')
d1 = subset(d1, Partition=='All')

# overall p values from binomial tests
# not a great plot, eh...
ggplot(d1, aes(Dataset)) + geom_bar(aes(weight = 1-pSH)) + facet_wrap(~Test, ncol=1) + theme(axis.text.x=element_text(angle=90,hjust=1))


# tree distances
t$comparison = paste(t$t1, "vs.", t$t2)


ggplot(t, aes(x=dataset)) + geom_bar(aes(weight = dist)) + facet_grid(test~comparison) + theme(axis.text.x=element_text(angle=90,hjust=1))

#rescale the big bars to 1, and then we can look at the detail a bit better
t$dist[t$dist>1] = 1
ggplot(t, aes(x=dataset)) + geom_bar(aes(weight = dist)) + facet_grid(test~comparison) + theme(axis.text.x=element_text(angle=90,hjust=1)) + ylim(0, 1)


# TO DO: combine d and t in such a way that we can colour each bar of the plot by a factor that indicates which of four cases are true
# t1 does not reject t2, and vice versa
# t1 rejects t2, but not vice versa
# t2 rejects t1, but not vice versa
# both trees reject each other
# all according to pSH from d