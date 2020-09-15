## experimenting with a mildly overfitted model; testing out confint/profile robustness

## add block info to existing R data set
InsectSprays$block <- rep(rep(1:6,each=2),6)
library(glmmTMB)
library(broom.mixed)
library(tidyverse)
library(ggstance)
## devtools::load_all("~/R/pkgs/glmmTMB/glmmTMB")
## debug(confint.profile.glmmTMB)
m1 <- glmmTMB(count~spray + (1|block/spray), family=poisson, data=InsectSprays)
## confint(m1,method="profile",parm="theta_")
m <- c("Wald","profile","uniroot")
names(m) <- m ## for map_dfr id
tt <- purrr::map_dfr(m, ~tidy(x=m1, effects="ran_pars", conf.int=TRUE, conf.method=.),
              .id="method")
tt2 <- (tt
    ## back-transform uniroot/profile results by hand
    %>% mutate_at(c("conf.low","conf.high"), ~ifelse(.<0,exp(.),.))
    ## unbounded CIs below
    %>% mutate_at("conf.low",replace_na,0)
)

ggplot(tt2,aes(x=estimate,xmin=conf.low,xmax=conf.high,y=group, colour=method)) +
    geom_pointrange(position=position_dodgev(height=0.5))+
    geom_vline(xintercept=0,lty=2)
