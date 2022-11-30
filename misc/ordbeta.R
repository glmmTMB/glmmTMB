library(glmmTMB)
## devtools::load_all("glmmTMB")
library(ordbetareg)
data(pew)
library(tidyverse)
library(cmdstanr)
library(broom.mixed)

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


TMB_fit <- glmmTMB(formula=therm/100 ~ education + income +
                       (1|region), 
                   data=model_data,
                   family = ordbeta,
                   start = list(psi = c(-1, 1)))

vcov(TMB_fit, full = TRUE)
View(tidy(TMB_fit, conf.int = TRUE))

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
)

ggplot(res, aes(estimate, term, colour = pkg)) +
    geom_pointrange(aes(xmin = lwr, xmax = upr))

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
glmmTMB(z ~ (1|g), data = dat, family = ordbeta,
        start = list(psi = c(-1, 1)))

