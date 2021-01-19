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
sim_fun <- function(N=25000,seed=NULL) {
    if (!is.null(seed)) set.seed(seed)
    df <- data.frame(
        id=1:N
      , disease=sample(0:1, N, replace=TRUE)
      , age=sample(18:88, N, replace=TRUE)
      , gender=sample(0:1, N, replace=TRUE)
      , cost=c(sample(0:100000, (N/2), replace=TRUE), rep(0, (N/2)))
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
    res <- (purrr::map_dfr(list(glm=glm_fit,glmmTMB=glmmTMB_fit),
                           tidy, .id="pkg")
        %>% dplyr::select(pkg,term,estimate,std.error)
    )
    return(res)
}

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

set.seed(101)
all_sims <- furrr::future_map_dfr(1:1000,
                                  ~sum_fun(sim_fun()),
                                  .progress=TRUE,
                                  .options=furrr_options(packages=f_pkgs, seed=TRUE))

pp <- function(x) x %>% collapse_fun() %>% select(-starts_with("mean_est_"))

## Gamma results
pp(all_sims)

all_sims_binom <- furrr::future_map_dfr(1:1000,
                                 ~sum_fun(sim_fun(),
                                          form = cost_binary ~ disease + gender + age,
                                          do_subset=FALSE,
                                          family=binomial),
                                 .progress=TRUE,
                                 .options=furrr_options(packages=f_pkgs, seed=TRUE))

pp(all_sims_binom)
##  term        mean_se_glm mean_se_glmmTMB calc_se_glm calc_se_glmmTMB
##  <chr>             <dbl>           <dbl>       <dbl>           <dbl>
## 1 (Intercept)    0.0394          0.0394      0.0384          0.0384  
## 2 age            0.000617        0.000618    0.000630        0.000630
## 3 disease        0.0253          0.0253      0.0263          0.0263  
## 4 gender         0.0253          0.0253      0.0253          0.0253  

