library(ggplot2)
library(lme4)
library(MuMIn)

# load it and get what we want
d = read.csv("~/Documents/github/SRHtests/summary_substitutions.csv")
d = subset(d, partition %in% c('pass', 'fail'))

# number of substitutions
d$n_subs = d$Speed * d$tree_length * d$sites

# do each test independently
mptis = subset(d, test == "MPTIS")
mpts = subset(d, test == "MPTS")
mptms = subset(d, test == "MPTMS")

# visualise it
ggplot(mptis, aes(y = log(n_subs), x = partition)) + geom_violin() + geom_jitter(size = 0.1, alpha = 0.4) + facet_wrap(~dataset)
ggplot(mptms, aes(y = log(n_subs), x = partition)) + geom_violin() + geom_jitter(size = 0.1, alpha = 0.4) + facet_wrap(~dataset)
ggplot(mpts, aes(y = log(n_subs), x = partition)) + geom_violin() + geom_jitter(size = 0.1, alpha = 0.4) + facet_wrap(~dataset)

# Fit models with n_subs as a fixed effect, and dataset as a random effect
fit.mptis = glmer(partition ~ log(n_subs) + (1 | dataset), family = binomial, data = mptis)
summary(fit.mptis)
r.squaredGLMM(fit.mptis)

fit.mpts = glmer(partition ~ log(n_subs) + (1 | dataset), family = binomial, data = mpts)
summary(fit.mpts)
r.squaredGLMM(fit.mpts)

fit.mptms = glmer(partition ~ log(n_subs) + (1 | dataset), family = binomial, data = mptms)
summary(fit.mptms)
r.squaredGLMM(fit.mptms)