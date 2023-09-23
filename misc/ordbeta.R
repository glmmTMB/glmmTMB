library(glmmTMB)
devtools::load_all("~/R/pkgs/glmmTMB/glmmTMB")
library(ordbetareg)
data(pew)
library(tidyverse)
library(cmdstanr)
library(broom.mixed)

do_slow <- TRUE

model_data <- select(pew,therm,age="F_AGECAT_FINAL",
                        sex="F_SEX_FINAL",
                        income="F_INCOME_FINAL",
                        ideology="F_IDEO_FINAL",
                        race="F_RACETHN_RECRUITMENT",
                        education="F_EDUCCAT2_FINAL",
                     region="F_CREGION_FINAL",
                        approval="POL1DT_W28",
                       born_again="F_BORN_FINAL",
                       relig="F_RELIG_FINAL",
                        news="NEWS_PLATFORMA_W28") %>% 
    mutate(across(c("race","ideology","income","approval","sex","education","born_again","relig"), function(c) {
      factor(c, exclude=levels(c)[length(levels(c))])
    })) %>% 
    # need to make these ordered factors for BRMS
    mutate(education=ordered(education),
           income=ordered(income))

## note glmmTMB doesn't automatically normalize
TMB_fit <- glmmTMB(formula=therm/100 ~ education + income +
                       (1|region), 
                   data=model_data,
                   family = ordbeta)


y <- model.frame(TMB_fit)$therm/100

model0 <- update(TMB_fit, . ~ 1)
ss <- simulate(model0)
hist(ss[[1]])
hist(y)

ss2 <- simulate(TMB_fit)
par(mfrow = c(1,2))
hist(ss2[[1]], breaks = 20)
hist(y, breaks = 20)

vcov(TMB_fit, full = TRUE)

simfun <- function(form, data, pars, show_pars = FALSE, ...) {
    ## for now assume response variable is in data
    r1 <- glmmTMB(form,
              data = data,
              ...,
              doFit = FALSE)
## construct TMB object, but don't fit it
    r2 <- fitTMB(r1, doOptim = FALSE)
    if (show_pars) return(r2$env$last.par)
    pars <- do.call("make_pars",
                    c(list(r2$env$last.par), pars))
    r2$simulate(par = pars)$yobs
}

make_pars <- function(pars, ...) {
    ## FIXME: check for name matches, length matches etc.
    L <- list(...)
    for (nm in names(L)) {
        pars[names(pars) == nm] <- L[[nm]]
    }
    return(pars)
}


simfun(y ~ 1, data=data.frame(fake = 1:1000),
       family = ordbeta, show_pars = TRUE)

ss <- simfun(y ~ 1, data=data.frame(y = rep(0.5, 1000), fake = 1:1000),
       family = ordbeta,
       pars = make_pars(list(beta=0, betad = 3,  psi = c(-2, 2))))
hist(ss, breaks = 100, freq = FALSE)
midprob <- plogis(-2, lower.tail = FALSE) * plogis(2)
shape1 <- plogis(0)*exp(3); shape2 <- (1-plogis(0))*exp(3)
curve(midprob*dbeta(x, shape1, shape2), add = TRUE, col = 2, add = TRUE)

ord_fit_mean <- ordbetareg(formula=therm ~ education + income +
                               (1|region), 
                           data=model_data,
                           cores=12,
                           chains=2,
                           iter=1000,
                           refresh=0,
                           threads = threading(5),
                           backend  = "cmdstanr")

res <- (list(glmmTMB = TMB_fit, ordbetareg = ord_fit_mean)
    |> purrr::map_dfr(tidy, effects = "fixed", .id = "pkg", conf.int = TRUE)
    |> select(pkg, term, estimate, lwr = conf.low, upr = conf.high)
    |> mutate(across(term, ~ gsub("E([0-9])","^\\1", .)))
)

if (do_slow) {
    res2 <- (tidy(TMB_fit, effects = "fixed", conf.int = TRUE, conf.method = "profile")
        |> select(term, estimate, lwr = conf.low, upr = conf.high)
        |> mutate(pkg = "glmmTMB", method = "profile", .before = 1)
    )
    res <- (res
        |> mutate(method = ifelse(pkg == "glmmTMB", "wald", "marginal"), .after = 1)
        |> bind_rows(res2)
    )
}
ggplot(res, aes(estimate, term, colour = pkg, shape = method)) +
    geom_pointrange(aes(xmin = lwr, xmax = upr), position = position_dodge(width = 0.2))

res |> arrange(term, method) |> View()

### explore profiles: confirm
pp <- profile(TMB_fit, trace = 10)
## compute and plot signed square root
pp2 <- (pp
    |> na.omit()
    |> full_join(tt, by = ".par")
    |> group_by(.par)
    |> mutate(.mid = .focal[value == 0][1],
              .zeta = sqrt(value)*sign(.focal-.mid))
)
library(ggplot2); theme_set(theme_bw())
ggplot(pp2, aes(.focal,value)) + geom_point() + geom_line() + facet_wrap(~.par, scale = "free_x")
ggplot(pp2, aes(.focal,.zeta)) + geom_point() + geom_line() + facet_wrap(~.par, scale = "free_x") +
    geom_smooth(method="lm", colour = adjustcolor("red", alpha.f = 0.5), se  = FALSE) +
    geom_hline(yintercept = 0, colour = "blue", lty = 2)



###
library(dplyr)
library(glmmTMB)

nobs <- 100
ngroup <- 5
dat <- tibble(
  x = sample(1:3, size = nobs, replace = TRUE, prob = c(0.15, 0.5, 0.35)),
  y = runif(nobs),
  z = case_when(
    x == 1 ~ 0,
    x == 2 ~ y,
    x == 3 ~ 1
  ),
  g = rep(seq_len(ngroup), length.out = nobs)
)
m1 <- glmmTMB(z ~ (1|g), data = dat, family = ordbeta)
confint(m1)
confint(m1, method = "profile")

