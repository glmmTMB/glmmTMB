library(HMMpa)
library(LaplacesDemon)
library(testthat)
library(glmmTMB)

xs <- c(0, 1, 5, 10, 20)
expect_equal(
  glmmTMB::dgenpois(xs, lambda1 = 3, lambda2 = 0.4),
  HMMpa::dgenpois(xs, lambda1 = 3, lambda2 = 0.4),
  tolerance = 1e-7
)

## Ntzoufras et al 2005 parameterization: lambda = lambda1/(1-lambda2)
expect_equal(
  glmmTMB::dgenpois(xs, lambda1 = 3, lambda2 = 0.4),
  LaplacesDemon::dgpois(xs, lambda = 3/(1-0.4), omega = 0.4),
  tolerance = 1e-7
)

gendat <- data.frame(y = c(11,10,9,10,9,8,11,7,9,9,9,8,11,10,11,9,10,7,13,9))
gen1 <- glmmTMB(y ~ 1, family = genpois(), data = gendat)
simfun <- function(n=100) simulate_new(~1, family = genpois, newdata = data.frame(x = rep(0, n)),
                                  newparams = list(beta = 2.25, betadisp = log(0.25)))[[1]]
simfun2 <- function(n=100) rgenpois(n, lambda1 = exp(2.25)*(1/sqrt(0.25)), lambda2 = 1-1/sqrt(0.25))
all.equal(var(simfun2(10000)), var(simfun(10000)), tolerance = 2e-2)
residuals(gen1, type = "dunn-smyth")

phi <- predict(gen1, type = "disp")
x <- sort(unique(gendat$y))
pgenpois_mu(x-1, 9.5, phi = rep(phi[1], length(x)))
pgenpois_mu(x-1, 9.5, phi = rep(phi[1], length(x)), clamp = FALSE)




set.seed(123)
dat_gp <- data.frame( y = rpois(300, lambda = 3), x =  rnorm(300), id = factor(rep(1:30, each = 10)) )
m_gp <- glmmTMB( y ~ x +  (1 | id), family = genpois, data = dat_gp )
summary(m_gp)
 
r_gp <- residuals(m_gp, type = "dunn-smyth")
plot(fitted(m_gp), r_gp)

## 1- 1/sqrt(phi) > -1 -> 1/sqrt(phi) < 2 -> sqrt(phi) < 1/2 -> phi 
