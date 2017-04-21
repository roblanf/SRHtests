library(ggplot2)
library(plyr)


proportion_p <- function(pvals){
  
  # ignore NAs
  pvals = pvals[complete.cases(pvals)]
  
  if(length(pvals)>9){
    sig = sum(pvals<0.05)
    return(sig / length(pvals))
  }else{
    return(NA)
  }
}

d = read.csv("~/Desktop/test.csv")

# this shows the proportion of p values <0.2 for each charset (where there were at least 10 p values)
d_prop = ddply(d, .(Dataset, Charset, Test), summarise, proportion = proportion_p(p.value))
p = ggplot(d_prop, aes(x = proportion, group = Test))
p + geom_histogram() + geom_vline(xintercept = 0.05, colour='red', alpha = 0.5) + facet_wrap(~Dataset, scales = 'free_y')

# we can summarise this further as a simple number of pass / fail charsets for each alignment
pass <- function(proportions){ proportions = proportions[complete.cases(proportions)]; return(sum(proportions<=0.05))}
fail <- function(proportions){ proportions = proportions[complete.cases(proportions)]; return(sum(proportions>0.05))}
na <- function(proportions){ return(sum(is.na(proportions)))}
d_passfail = ddply(d_prop, .(Dataset), summarise, fail = fail(proportion), pass = pass(proportion), no_data = na(proportion))

# we can look at an individual dataset / charset like this
dnum=6
name = levels(d$Dataset)[dnum]
propsub = subset(d_prop, Dataset==name)

# worst charset in the dataset
cs = propsub$Charset[which(propsub$proportion==max(propsub$proportion, na.rm=T))][1]
c = ggplot(subset(d, Dataset==name & Charset==as.character(cs)), aes(x = p.value))
c + geom_histogram() + ggtitle(paste(cs, 'from', name))


# whole dataset
p = ggplot(subset(d, Dataset==levels(d$Dataset)[dnum]), aes(x=p.value, colour = Charset))

p + geom_histogram(aes(fill=Charset)) + facet_wrap(~Charset)
p + geom_histogram(aes(fill=Charset)) + facet_grid(Charset~Test, scales='free_y')

