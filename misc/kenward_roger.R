## lots of miscellaneous stuff for exploring properties of K-R ddf, especially
## for GLMMs
library(glmmTMB)
devtools::load_all("glmmTMB")
library(lmerTest)
library(parameters)
library(MASS)  
library(dplyr)  

source("misc/kenward_roger_funs.R")
data("sleepstudy")

comp_df <- function(form, data) {
    t1 <- system.time(m1 <- lmer(form, data))
    t2 <- system.time(m2 <- glmmTMB(form, data, REML = TRUE))
    aa <- function(val, time) { attr(val, "time") <- time; return(val) }
    list(lmer =    aa(coef(summary(m1))[,"df"], t1),
         glmmTMB = aa(c(glmmTMB::dof_kenward(m2)), t2))
}

comp_df(Reaction ~ Days + (Days|Subject), sleepstudy)
comp_df(Reaction ~ Days + (1|Subject), sleepstudy)

## ?? glmmTMB model convergence problem
comp_df(Reaction ~ Days + (Days||Subject), sleepstudy)

################ Testing K-R dfs ##################
set.seed(444) 

###### Testing K-R dfs for normally distributed data #################

nsim <- 1000
z_stat_normal <- rep(NA, nsim)
pb <- txtProgressBar(max=nsim, style = 3)

sim_data <- function(n_genes = 500, n_samples = 4, n_treatments = 2,
                     beta = c(5, 0),
                     theta = c(0, 0),
                     betadisp = 0,
                     family = gaussian) {
    n_obs <- n_samples * n_treatments  # Total observations per gene
    
    ## Simulating data for genes and samples
    gene <- rep(1:n_genes, each = n_obs)  # Replicate genes for each treatment (same genes measured in both treatments)
    sample <- rep(1:n_samples, times = n_genes * n_treatments)  # Replicate samples for each gene-treatment pair
    treatment <- rep(c("Control", "Treated"), times = n_genes * n_samples)  # Two treatments per sampl
    dd <- data.frame(Gene = factor(gene), Sample = factor(sample),
                     Treatment = factor(treatment))
    dd$expression <- simulate_new(~ Treatment + (1|Gene) + (1|Sample:Treatment),
                                  newdata = dd,
                                  newparams = list(beta = beta, theta = theta,
                                                   betadisp = betadisp),
                                  family = family)[[1]]
    return(dd)
}

set.seed(101)
for (i in seq.int(nsim)) {
    setTxtProgressBar(pb, i)
    
    ## Create the data frame
    simulated_data <- sim_data()
    model_normal <- glmmTMB(expression ~ Treatment + (1|Gene) + (1|Sample:Treatment),
                            data = simulated_data, REML=TRUE)
    
    z_stat_normal[i] <- coef(summary(model_normal))$cond["TreatmentTreated","z value"]
}

close(pb)

z_stat_normal <- z_stat_normal[is.finite(z_stat_normal)]

mean(abs(z_stat_normal)>1.96)

set.seed(101)
simulated_data <- sim_data()
comp_df(expression ~ Treatment + (1|Gene) + (1|Sample:Treatment),
        simulated_data)


system.time(
    dof_kenward(model_normal) # K-R degrees of freedom
)

m <- MASS::fitdistr(z_stat_normal[is.finite(z_stat_normal)], "t") # Fit t-distribution to z-statistics

print(m)

## Plot fitted t-distribution over histogram of z-statistics
param.t <- m$estimate
dh <- hist(z_stat_normal, prob=TRUE, breaks=100)
##to get the scaling correct, then do
par(new = TRUE)
ss=seq(range(dh$mids)[1],range(dh$mids)[2],length.out = 1000)
x=((ss-param.t["m"])/param.t["s"])
plot(x,100*dt(x=x,df=param.t["df"]), type="l",
     col="red", lwd=3,xlab="",axes=F,ylab="")


theoretical_quantiles <- with(m,
                              qt(ppoints(length(z_stat_normal)),
                                 df = estimate["df"],
                                 ncp = (estimate["m"] / estimate["s"]))*estimate["s"])

## Create the Q-Q plot
qqplot(theoretical_quantiles, z_stat_normal,
       main = "Q-Q plot: Fitted t-distribution",
       xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
       pch = 19, col = "black")
abline(0, 1, col = "black", lty = 1)

####### Now for negative binomial data ###########



pb <- txtProgressBar(max=nsim, style = 3)
z_stat_nbinom <- rep(NA, nsim)
set.seed(101)
for (i in seq.int(nsim)) {
    setTxtProgressBar(pb, i)

    ## Create the data frame
    simulated_data <- sim_data(n_genes=100, n_samples = 4, n_treatments = 2,
                               betadisp = log(1000),
                               beta = c(6, 0),
                               theta = c(2,2),
                               family = nbinom2)
    
    model_nbinom <- glmmTMB(expression ~ Treatment + (1|Gene) + (1|Sample:Treatment),
                         data = simulated_data,
                         family = nbinom2,
                         REML=TRUE)
    
    z_stat_nbinom[i] <- coef(summary(model_nbinom))$cond["TreatmentTreated","z value"]
}
close(pb)

ggplot(simulated_data_nbinom, aes(x=Treatment, y=expression)) +
    geom_point()+
    facet_wrap(~Gene, scales="free")

glmmTMB::dof_kenward(model_nbinom) # K-R degrees of freedom

length(z_stat_nbinom[abs(z_stat_nbinom)>1.96])/length(z_stat_nbinom)
length(z_stat_nbinom)
summary(model_nbinom)


glmmTMB(y ~ 1, family = t_family,
        data = data.frame(y = z_stat_nbinom[is.finite(z_stat_nbinom)]))
        
m <- MASS::fitdistr(z_stat_nbinom[is.finite(z_stat_nbinom)],
                    "t", start = list(m=0, s=1,df = 100)) # fit t-distribution to nbinom z-stat
print(m)

## Plot fitted t-distribution over histogram of z-statistics
param.t <- m$estimate
dh <- hist(z_stat_nbinom, prob=TRUE, breaks=200)
##to get the scaling correct, then do
par(new = TRUE)
ss <- seq(range(dh$mids)[1],range(dh$mids)[2],length.out = 1000)
x <- ((ss-param.t["m"])/param.t["s"])
plot(x,100*dt(x=x,df=param.t["df"]), type="l",
     col="red", lwd=3,xlab="",axes=F,ylab="")

theoretical_quantiles <- qt(ppoints(length(z_stat_nbinom)), df  = m$estimate[["df"]],
                            ncp = (m$estimate["m"] / m$estimate["s"]))*m$estimate["s"]

## Create the Q-Q plot
qqplot(theoretical_quantiles, z_stat_nbinom,
       main = "Q-Q plot: Fitted t-distribution",
       xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
       pch = 19, col = "black")
abline(0, 1, col = "black", lty = 1)



sim_data <- function(n_genes = 500, n_samples = 4, n_treatments = 2,
                     beta = c(5, 0),
                     theta = c(0, 0),
                     betadisp = 0,
                     family = gaussian,
                     seed = NULL) {
    n_obs <- n_samples * n_treatments  # Total observations per gene
    
    ## Simulating data for genes and samples
    gene <- rep(1:n_genes, each = n_obs)  # Replicate genes for each treatment (same genes measured in both treatments)
    sample <- rep(1:n_samples, times = n_genes * n_treatments)  # Replicate samples for each gene-treatment pair
    treatment <- rep(c("Control", "Treated"), times = n_genes * n_samples)  # Two treatments per sampl
    dd <- data.frame(Gene = factor(gene), Sample = factor(sample),
                     Treatment = factor(treatment))
    dd$expression <- simulate_new(~ Treatment + (1|Gene) + (1|Sample:Treatment),
                                  newdata = dd,
                                  newparams = list(beta = beta, theta = theta,
                                                   betadisp = betadisp),
                                  family = family,
                                  seed = seed)[[1]]
    return(dd)
}

s1 <- sim_data(n_genes=100, beta = c(3, 2),
               betadisp = log(100),
               family = nbinom2,
               seed = 101)
m1 <- glmmTMB(expression ~ Treatment + (1|Gene) + (1|Sample:Treatment),
              family = nbinom2, data = s1)
glmmTMB::dof_kenward(m1)


nsim <- 400
set.seed(101)
rr <- replicate(nsim, {
    cat(".")
    s1 <- sim_data(n_genes=100, beta = c(3, 2),
                   betadisp = log(100),
                   family = nbinom2)
    m1 <- glmmTMB(expression ~ Treatment + (1|Gene) + (1|Sample:Treatment),
                  family = nbinom2, data = s1)
    ddf <- try(glmmTMB::dof_kenward(m1), silent = TRUE)
    if (inherits(ddf, "try-error")) ddf <- NA
    c(coef(summary(m1))$cond["TreatmentTreated", "z value"],
          ddf[2])
})
rr <- t(rr)

df_fit <- glmmTMB(y~1, family = t_family,
                  data = data.frame(y = rr[,1]),
                  start = list(psi = -3))
cc <- coef(summary(df_fit))
qdist <- \(x) { qt(x, df = family_params(df_fit), ncp = (cc$cond[3]))*cc$cond[2] }
ddist <- \(x) { dt(x/cc$cond[2], df = family_params(df_fit), ncp = (cc$cond[3]))/cc$cond[2] }

qqplot(qdist(ppoints(nsim)), rr[,1], conf.level = 0.95)
qqline(rr[,1], distribution = qdist, probs = c(0.1, 0.6), col = 2)

hist(rr[,1] , breaks = 100, prob = TRUE)
curve(ddist(x), add = TRUE, col = 2, lwd = 2)

krdf <- rr[,2]
krdf <- krdf[is.finite(krdf) & krdf < 100]
length(krdf)
hist(log(krdf), breaks = 100)
abline(v = log(family_params(df_fit)), col = 2)
