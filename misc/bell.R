library(glmmTMB)
library(GLMMadaptive)
library(bellreg)
set.seed(1234)

W <- gsl::lambert_W0

##user defined function for GLMMadaptive package##
my_bell<- function (link = "log") {
    stats <- make.link(link)
    log_dens_mu <- function (y, eta, mu_fun, phis, eta_zi) {
                                        # the log density function
                                        # you link log(mu) to covariates
        mu <- as.matrix(mu_fun(eta))
        comp1 <- y*log(W(mu))
        comp2 <- -exp(W(mu))
        out <- comp1 + comp2
        attr(out, "mu_y") <- mu
        out
    }

    log_dens <- function (y, eta, mu_fun, phis, eta_zi){
                                        # the derivative of the log density w.r.t. mu
        mu <- as.matrix(mu_fun(eta))
        theta <- W(mu)
        
        comp1 <- -exp(theta)
        comp3 <- log(theta)*(y)
        
        out <- (comp1 + comp3)
        attr(out, "mu_y") <- mu
        out
    }
    score_eta_fun <- function (y, mu, phis, eta_zi) {
                                        # the derivative of the log density w.r.t. mu
        mu <- as.matrix(mu)
        theta <- W(mu)
                                        # the derivative of the log density w.r.t. theta
        comp1 <-  -exp(theta)
        comp3 <- (1/theta)*(y)
                                        # the derivative of theta w.r.t mu
        tht_mu <- W(mu)/(mu*W(mu)+mu)
                                        # the derivative of mu w.r.t. eta log link
        mu_eta <- mu
        (comp1 + comp3) * tht_mu * mu_eta

    }

    score_phis_fun <- function (y, mu, phis, eta_zi) {
                                        # the derivative of the log density w.r.t. phis

    }
    structure(list(family = "user bell", link = stats$name, linkfun = stats$linkfun,
                   linkinv = stats$linkinv, 
                   log_dens = log_dens, 
                   score_eta_fun = score_eta_fun),
              class = "family")
}
######################
## simulated homogeneous case

set.seed(101)
dd <- data.frame(y = bellreg::rbell(1000, theta = 1))
## summary(dd$y)
table(dd$y)
m0 <- bellreg(y ~ 1, data = dd)
m1 <- glmmTMB(y ~ 1, data = dd, family = bell, start = list(beta=0.5))

## checking gradients

parvec <- seq(0, 2, length = 101)
likvec <- vapply(parvec, m1$obj$fn, numeric(1))
likvec2 <- vapply(parvec, \(t) -sum(dbell(dd$y, theta = W(exp(t)), log = TRUE)), numeric(1))

par(mfrow=c(2,1))
plot(parvec, likvec, type = "l", col = "red", lwd = 2, log="y")
points(parvec, likvec2, col = "blue", lwd = 2)

stopifnot(max(abs(likvec-likvec2)) < 1e-5)

grvec <- vapply(parvec, m1$obj$gr, numeric(1))
plot(parvec, grvec, type = "l", col = "red", lwd = 2)
dx <- diff(parvec[1:2])
lines(parvec[-1]+dx/2, diff(likvec)/dx)
abline(h=0, col="grey")

m1$fit$objective + m0$fit$value
m1$fit$objective
m1$obj$fn(m0$fit$par)
-sum(dbell(dd$y, theta = W(exp(m0$fit$par)), log = TRUE))
## something slightly weird with bellreg?

################
## Leprosy dataset

leprosy <- read.table("leprosy.txt", header = FALSE)
names(leprosy) <- c("drug", "y1", "y2")
## Placebo as control
leprosy$drug <- relevel(factor(leprosy$drug), "C")

## leprosy <- read.dta("leprosy.dta")

## To long format
leprosy.long <- reshape(data= leprosy,
                        varying= paste0("y",1:2),
                        ##v.names   = "",
                        timevar= "time",
                        idvar= "id",
                        direction= "long",
                        sep= ""
                        )
## Change time 1, 2 to 0, 1
leprosy.long$time <- leprosy.long$time - 1
## Order by id and then time
leprosy.long <- leprosy.long[with(leprosy.long, order(id, time)),]

fm1 <- mixed_model(fixed = y ~ time*drug, random = ~ 1| id,data=leprosy.long,
                   n_phis = 0,
                   family =my_bell(),max_coef_value = 50,initial_values = list("betas" = poisson()))# Laplace approximation
summary(fm1)
#####glmmTMB code
glmm.model <- glmmTMB(y ~  time*drug+ (1|id),
                      data = leprosy.long,
                      na.action = na.omit,
                      family = glmmTMB::bell())

all.equal(fixef(glmm.model)$cond, fixef(fm1))

fabric <- read.table("fabric.txt", header = TRUE)
plot(faults ~ roll_length, data = fabric)
glmmTMB(faults ~ log(roll_length), data = fabric, family = bell)

glmmadapt <- function() mixed_model(fixed = y ~ time*drug,
                                    random = ~ 1| id,data=leprosy.long,
                                    n_phis = 0,
                                    family =my_bell(),max_coef_value = 50,
                                    initial_values = list("betas" = poisson()))

tmbfun <- function() glmmTMB(y ~  time*drug+ (1|id),
                      data = leprosy.long,
                      na.action = na.omit,
                      family = glmmTMB::bell())

library(rbenchmark)
benchmark(glmmadapt(), tmbfun())

bench::mark(glmmadapt(), tmbfun(), check = FALSE)

fabric <- read.table("fabric.txt", header = TRUE)
plot(faults ~ roll_length, data = fabric)

pois_fit <- glm(faults ~ roll_length, data = fabric, family = poisson)
bell_fit <- glmmTMB(faults ~ roll_length, data = fabric, family = bell)
## agrees with Castellares et al. (2018
printCoefmat(coef(summary(bell_fit))$cond, digits = 3)
