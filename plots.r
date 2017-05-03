library(ggplot2)
library(plyr)


sig_p <- function(pvals){
	pvals = pvals[complete.cases(pvals)]
	return(sum(pvals<0.05))
}

nonsig_p <- function(pvals){
	pvals = pvals[complete.cases(pvals)]
	return(sum(pvals>=0.05))
}

na_p <- function(pvals){
	return(sum(is.na(pvals)))
}

binom_p <- function(sig, nonsig){
	total = sig + nonsig
	binom.test(sig, total, p = 0.05, alternative = 'greater')
}

#read in the data
d = read.csv("Documents/github/SRHtests/ouput.csv")
d$p.value = as.numeric(as.character(d$p.value)) # Oh R

# make a dataframe with one row per charset, p values that are <0.05, >=0.05, NA, binomial
charsets = ddply(d, .(Dataset, Charset, Test), summarise, sig = sig_p(p.value), non = nonsig_p(p.value), na = na_p(p.value))
charsets$binomial = apply(charsets, 1, function(x)  binom.test(c(as.integer(x[4]), as.integer(x[5])), p = 0.05, alternative = 'greater')$p.value)




# we can summarise this further as a simple number of pass / fail charsets for each alignment
pass <- function(pvals){ return(sum(pvals>=0.05))}
fail <- function(pvals){ return(sum(pvals<0.05))}
passfail = ddply(charsets, .(Dataset), summarise, fail = fail(binomial), pass = pass(binomial))



# we can look at an individual dataset / charset like this
dnum=6
name = levels(d$Dataset)[dnum]
propsub = subset(charsets, Dataset==name)

# worst charset in the dataset
cs = propsub$Charset[which(propsub$binomial==min(propsub$binomial, na.rm=T))][1]
c = ggplot(subset(d, Dataset==name & Charset==as.character(cs)), aes(x = p.value))
c + geom_histogram() + ggtitle(paste(cs, 'from', name))


# plot a whole dataset (gets tricky with lots of lines...)
p = ggplot(subset(d, Dataset==levels(d$Dataset)[dnum]), aes(x=p.value))

p + geom_histogram() + facet_wrap(~Charset)
p + geom_histogram() + facet_grid(Charset~Test, scales='free_y')
p + geom_line(aes(color=Charset), stat="density", size=1, alpha=0.4)

# for lots of charsets


# here's a bad dataset
p = ggplot(subset(d, Dataset=="Borowiec_2015"), aes(x=p.value))
p + geom_line(stat="density", alpha=0.25, aes(group=Charset, y = ..scaled..))
