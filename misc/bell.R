library(glmmTMB)
library(GLMMadaptive)
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
fixef(fm1)
summary(fm1)
#####glmmTMB code
glmm.model <- glmmTMB(y ~  time*drug+ (1|id),
                      data = leprosy.long,
                      na.action = na.omit,
                      family = glmmTMB::bell())

fixef(glmm.model)$cond
glmm.model$obj$fn()
tmb_pars <- with(glmm.model$obj$env, last.par.best[-random])
fm1_pars <- c(fixef(fm1), logsd=log(attr(fm1$D, "L")))
glmm.model$obj$fn(c(fixef(fm1), log(attr(fm1$D, "L"))))
logLik(fm1)
summary(glmm.model)

## OK, why doesn't this work? Test on something easier ...
## set.seed(102)  ## neg log-likelihood problem
dd <- data.frame(f = factor(rep(1:10, each = 20)))
dd$y <- simulate_new( ~ 1 + (1|f),
                     family = bell,
                     newdata = dd,
                     newparams = list(theta = 0, beta = 0))[[1]]

boxplot(y~f, data = dd)
glmmTMB(y ~ 1 + (1|f),
        data = dd,
        family = bell) |> summary()
