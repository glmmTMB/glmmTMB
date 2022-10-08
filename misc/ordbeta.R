## cutpoints k1, k2
## alpha = 1- g(eta -k1)
## delta = g(eta-k1) -



  ## if(y==0) {
  ##     return log1m_inv_logit(mu - thresh[1]);
  ##   } else if(y==1) {
  ##     return log_inv_logit(mu  - thresh[2]);
  ##   } else {
  ##     return log_diff_exp(log_inv_logit(mu   - thresh[1]), log_inv_logit(mu - thresh[2])) +
  ##               beta_lpdf(y|exp(log_inv_logit(mu) + log(phi)),exp(log1m_inv_logit(mu) + log(phi)));
  ##   }

library(glmmTMB)
library(ordbetareg)
data(pew)
library(dplyr)

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
                   start = list(thetaf = c(-1, 1)))

ord_fit_mean <- ordbetareg(formula=therm ~ mo(education)*mo(income) +
                               (1|region), 
                           data=model_data,
                           cores=2,chains=2,iter=1000,
                           refresh=0)
## NOTE: to do parallel processing within chains
## add the options below
##threads=threading(5),
##backend="cmdstanr"
##where threads is the number of cores per chain
## you must have cmdstanr set up to do so
## see https://mc-stan.org/cmdstanr/

## TO DO
## check vs ordbetareg (w/o mo())
## initialize() check in family (0<= y <= 1)
