devtools::load_all("~/R/pkgs/glmmTMB/glmmTMB")
data("sleepstudy", package = "lme4")
pp <- list(beta = c(280, 3),
           betad = 1,
           theta = c(3, 0.5, 2))


## basic sim (checking for sensible parameters)
s1 <- simulate_new( ~ Days + (Days|Subject),
                   seed = 101,
                   newdata = sleepstudy,
                   family = gaussian,
                   newparams = pp)[[1]]
plot(sleepstudy$Days, s1, ylim = c(180,450))
with(sleepstudy, points(Days, Reaction, col = 2))

## debug(simulate_new)
s2 <- simulate_new( ~ Days + (Days|Subject),
             newdata = sleepstudy,
             family = gaussian,
             newparams = pp,
             return_val = "object")
fixef(s2)
VarCorr(s2) ## SDs are exp(3), exp(0.5), cor is 2/(sqrt(2^2 +1)
getME(s2, "b")

s3 <- simulate_new( ~ Days + (Days|Subject),
             newdata = sleepstudy,
             family = gaussian,
             newparams = pp,
             return_val = "pars")

### now work on fixing particular elements ...
pp2 <- c(pp, list(b = rep(c(-1, 0, 1), length.out = 36)))

s4 <- simulate_new( ~ Days + (Days|Subject),
                   seed = 101,
                   newdata = sleepstudy,
                   family = gaussian,
                   newparams = pp2)[[1]]
plot(sleepstudy$Days, s4, ylim = c(180,450))
with(sleepstudy, points(Days, Reaction, col = 2))
points(sleepstudy$Days, s1, col = 4)

s5 <- simulate_new( ~ Days + (Days|Subject),
                   seed = 101,
                   newdata = sleepstudy,
                   family = gaussian,
                   newparams = pp2,
                   return_val = "pars")

print(s5)

s6 <- simulate_new( ~ Days + (Days|Subject),
                   seed = 101,
                   newdata = sleepstudy,
                   family = gaussian,
                   newparams = pp2,
                   return_val = "object")

getME(s6, "b")  ## correct, although names are wrong?

pp3 <- c(pp, list(b = list("Days|Subject" = pp2$b)))

s7 <- simulate_new( ~ Days + (Days|Subject),
                   seed = 101,
                   newdata = sleepstudy,
                   family = gaussian,
                   newparams = pp3,
                   return_val = "pars")
## same as before

## try some examples with multiple b terms ...

pp4 <- list(beta = c(280),
           betad = 1,
           theta = c(-1, 1, 0))


## basic sim (checking for sensible parameters)
s8 <- simulate_new( ~ 1 + (1|Subject) + ar1(0+factor(Days)|Subject),
                   seed = 101,
                   newdata = sleepstudy,
                   family = gaussian,
                   newparams = pp4)[[1]]

s9 <- simulate_new( ~ 1 + (1|Subject) + ar1(0+factor(Days)|Subject),
                   seed = 101,
                   newdata = sleepstudy,
                   family = gaussian,
                   newparams = pp4,
                   return_val = "pars")

nb <- sum(names(s9) == "b")
pp5 <- c(pp4, list(b = rep(c(-1, 0, 1), length.out = nb)))
## basic sim (checking for sensible parameters)
s10 <- simulate_new( ~ 1 + (1|Subject) + ar1(0+factor(Days)|Subject),
                   seed = 101,
                   newdata = sleepstudy,
                   family = gaussian,
                   newparams = pp5)[[1]]


ns <- length(unique(sleepstudy$Subject))
pp6 <- c(pp4, list(b = list("0+factor(Days)|Subject" =
                                rep(c(-1, 0, 1), length.out = nb - ns))))

s11 <- simulate_new( ~ 1 + (1|Subject) + ar1(0+factor(Days)|Subject),
                   seed = 101,
                   newdata = sleepstudy,
                   family = gaussian,
                   newparams = pp6)[[1]]

s12 <- simulate_new( ~ 1 + (1|Subject) + ar1(0+factor(Days)|Subject),
                   seed = 101,
                   newdata = sleepstudy,
                   family = gaussian,
                   newparams = pp6,
                   return_val = "pars")


