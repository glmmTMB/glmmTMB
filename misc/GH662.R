library(tidyverse)
library(glmmTMB)
library(broom)
library(broom.mixed)
library(tidyverse)
library(future)
library(furrr)

## set number of cores appropriately, and/or use plan("sequential")
##  to do non-parallel evaluation
options(mc.cores=10)
plan("multisession")
f_pkgs <- c("broom","broom.mixed","glmmTMB")

### DATA

#' simulation function
#' @param N number of obs (should be even)
#' @param seed random-number seed
sim_fun <- function(N=25000,seed=NULL, do_gamma=FALSE) {
    ## sample from Gamma *or* uniform dist
    rval <- if (!do_gamma) sample(0:100000, N/2, replace=TRUE) else rgamma(N/2,shape=2,scale=1)
    if (!is.null(seed)) set.seed(seed)
    ## sensible ranges, but null model (no relationship between predictors & response
    df <- data.frame(
        id=1:N
      , disease=sample(0:1, N, replace=TRUE)
      , age=sample(18:88, N, replace=TRUE)
      , gender=sample(0:1, N, replace=TRUE)
      , cost=c(rval, rep(0, N/2))
      , time=sample(30:3287, N, replace=TRUE))
    df$cost_binary <- as.numeric(df$cost>0)
    return(df)
}

#' summary function: fit GLM via glm() and glmmTMB()
#' @param df data frame from sim_fun() 
sum_fun <- function(df,
                    form = cost ~ disease + gender + age + offset(log(time)),
                    do_subset=TRUE,
                    family=Gamma(link="log")) {
    if (do_subset) df <- subset(df, cost>0)
    glm_fit <- glm(form,data=df,family=family)
    glmmTMB_fit <- glmmTMB(form,data=df,family=family)
    ## extract estimates/SEs and combine
    res <- (purrr::map_dfr(list(glm=glm_fit,glmmTMB=glmmTMB_fit),
                           tidy, .id="pkg")
        %>% dplyr::select(pkg,term,estimate,std.error)
    )
    return(res)
}

## compute mean SE and sd(estimate) across the ensemble, by pkg & term
collapse_fun <- function(all_sims) {
    a0 <- (all_sims
        %>% group_by(pkg,term)
        %>% summarise(mean_se=mean(std.error),
                      mean_est=mean(estimate),
                      calc_se=sd(estimate),.groups="drop")
        %>% ungroup()
    )
    a1 <- a0 %>% tidyr::pivot_wider(id_cols=c(pkg,term),names_from=pkg,values_from=mean_se:calc_se)
    return(a1)
}

## pretty-print (leave out mean estimates)
pp <- function(x) x %>% collapse_fun() %>% select(-starts_with("mean_est_"))

### 
set.seed(101)
all_sims <- furrr::future_map_dfr(1:1000,
                                  ~sum_fun(sim_fun()),
                                  .progress=TRUE,
                                  .options=furrr_options(packages=f_pkgs, seed=TRUE))


## equivalent, but for no offset
set.seed(101)
all_sims_nooff <- furrr::future_map_dfr(1:1000,
                                 ~sum_fun(sim_fun(),
                                          form = cost ~ disease + gender + age
                                          ),
                                 .progress=TRUE,
                                 .options=furrr_options(packages=f_pkgs, seed=TRUE))


## response actually Gamma-distributed ...
set.seed(101)
all_sims_gamma <- furrr::future_map_dfr(1:1000,
                                  ~sum_fun(sim_fun(do_gamma=TRUE)),
                                  .progress=TRUE,
                                  .options=furrr_options(packages=f_pkgs, seed=TRUE))


## response actually Gamma-distributed, no offset ...
set.seed(101)
all_sims_gamma_nooff <- furrr::future_map_dfr(1:1000,
              ~sum_fun(sim_fun(do_gamma=TRUE)),
              form = cost ~ disease + gender + age,
                                  .progress=TRUE,
                                  .options=furrr_options(packages=f_pkgs, seed=TRUE))

pp(all_sims_gamma)

all_sims_binom <- furrr::future_map_dfr(1:1000,
                                 ~sum_fun(sim_fun(),
                                          form = cost_binary ~ disease + gender + age,
                                          do_subset=FALSE,
                                          family=binomial),
                                 .progress=TRUE,
                                 .options=furrr_options(packages=f_pkgs, seed=TRUE))

## RESULTS
pp(all_sims)
   term        mean_se_glm mean_se_glmmTMB calc_se_glm calc_se_glmmTMB
   <chr>             <dbl>           <dbl>       <dbl>           <dbl>
 1 (Intercept)     0.0655         0.0340       0.0669          0.0669 
 2 age             0.00103        0.000533     0.00104         0.00104
 3 disease         0.0421         0.0219       0.0427          0.0427 
 4 gender          0.0421         0.0219       0.0421          0.0421 

pp(all_sims_binom)
  term        mean_se_glm mean_se_glmmTMB calc_se_glm calc_se_glmmTMB
  <chr>             <dbl>           <dbl>       <dbl>           <dbl>
 1 (Intercept)    0.0394          0.0394      0.0384          0.0384  
 2 age            0.000617        0.000618    0.000630        0.000630
 3 disease        0.0253          0.0253      0.0263          0.0263  
 4 gender         0.0253          0.0253      0.0253          0.0253  

pp(all_sims_nooff)
  term        mean_se_glm mean_se_glmmTMB calc_se_glm calc_se_glmmTMB
  <chr>             <dbl>           <dbl>       <dbl>           <dbl>
1 (Intercept)    0.0161          0.0209      0.0166          0.0166  
2 age            0.000252        0.000327    0.000266        0.000266
3 disease        0.0103          0.0134      0.0104          0.0104  
4 gender         0.0103          0.0134      0.0104          0.0104  

pp(all_sims_gamma_nooff)
## standard error too SMALL??

## simpler/standalone
## this seems OK (based on one replicate?)

N <- 5e4
set.seed(101)
library(glmmTMB)
library(bbmle)
dd <- data.frame(y=rgamma(N,shape=2,scale=1))
m1 <- glm(y~1,data=dd,family=Gamma(link="log"))
m2 <- mle2(y~dgamma(shape=exp(log_shape),
                    scale=exp(log_mu-log_shape)),
                    data=dd,
                    start=list(log_mu=0,log_shape=0))
m3 <- glmmTMB(y~1,data=dd,family=Gamma(link="log"))

cc <- c("Estimate", "Std. Error")
coef(summary(m1))[,cc]
coef(summary(m2))["log_mu",cc]
coef(summary(m3))$cond[,cc]


rgamma2 <- function(n, shape, mu) {
    rgamma(n, shape=shape, scale=mu/shape)
}
dgamma2 <- function(x, shape, mu, log=FALSE) {
    dgamma(x, shape=shape, scale=mu/shape, log=log)
}

## also fine with a single covariate
set.seed(101)
dd <- data.frame(x=rnorm(N))
dd$y <- with(dd,rgamma2(N,shape=2,mu=exp(1+2*x)))
m1 <- glm(y~x,data=dd,family=Gamma(link="log"))
m2 <- mle2(y~dgamma2(shape=exp(log_shape), mu=exp(log_mu)),
                    data=dd,
                    parameters=list(log_mu~x),
                    start=list(log_mu=0,log_shape=0))
m3 <- glmmTMB(y~x,data=dd,family=Gamma(link="log"))

cc <- c("Estimate", "Std. Error")
coef(summary(m1))[,cc]
coef(summary(m2))[1:2,cc]
coef(summary(m3))$cond[,cc]
