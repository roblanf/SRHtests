library(ggplot2)

d = read.csv("~/Dropbox/Projects_Current/systematic_bias/tables/summary_charsets.csv")



# is the number of bad subsets linked to the 
# we should really do this with tree lengths...
ggplot(d, aes(x = Avg.Bad.Charset.length., y = Avg.Not.bad.Charset.length.)) + geom_point()

ggplot(d, aes(x = Delta.length..Bad.NB.)) + geom_histogram() + xlim(-100, 100)

ggplot(d, aes(x = Delta.length..Bad.NB.)) + geom_density() + xlim(-100, 100)


# sign test
diffs = d$Delta.length..Bad.NB.
pos = sum(diffs>0, na.rm = T)
neg = sum(diffs<0, na.rm = T)
