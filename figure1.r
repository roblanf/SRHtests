library(ggplot2)
library(plyr)
library(scales)

d = read.csv("~/Dropbox/Projects_Current/systematic_bias/tables/is_charset_bad.csv")
d = d[,-c(1)]
names(d) = c("charset", "test", "fail", "pass", "undetermined", "p_binomial", "dataset", "bad")

# we only count charsets where we could measure at least 20 pairs of sequences for a test
d$p = d$p_binomial
d$N = d$fail+d$pass
d$p[which(d$N<20)] = NA

# now we can categorise pass/fail/NA for the p-values
d$result = "undetermined"
d$result[which(d$p<0.05)] = 'fail'
d$result[which(d$p>=0.05)] = 'pass'

d$result = as.factor(d$result)
d$result = addNA(d$result)
d$result = factor(d$result, levels = c("pass", "undetermined", "fail"))



# Figure 1A overall p values from binomial tests
ggplot(d, aes(x=p)) + 
    geom_histogram(bins=40) + 
    facet_wrap(~test, ncol=1) + 
    geom_vline(xintercept = 0.05, colour = 'red', linetype = 'dashed') +
    facet_wrap(~test, ncol=1)

# dataset-by-test how many are pass/fail/NA
gg_color_hue <- function(n) {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
}
mycols <- gg_color_hue(length(unique(d$result)))
cols = mycols
cols[2] = 'grey'
cols[1] = mycols[3]
cols[3] = mycols[1]

# Figure 1B - proportion of each dataset that fails...
ggplot(d, aes(x = dataset)) + 
    geom_bar(aes(fill = result), position = "fill", na.rm = T) + 
    facet_wrap(~test, ncol=1) + 
    scale_fill_manual(values=cols) + 
    ylab("Proportion of character sets") +
    theme(axis.text.x=element_text(angle=90,hjust=1))
