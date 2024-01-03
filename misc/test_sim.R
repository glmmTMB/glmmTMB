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
s1 <- simulate_new( ~ Days + (Days|Subject),
             newdata = sleepstudy,
             family = gaussian,
             newparams = pp,
             return_val = "object")
fixef(s1)
VarCorr(s1) ## SDs are exp(3), exp(0.5), cor is 2/(sqrt(2^2 +1)
getME(s1, "b")

debug(simulate_new)

s1 <- simulate_new( ~ Days + (Days|Subject),
             newdata = sleepstudy,
             family = gaussian,
             newparams = pp,
             return_val = "pars")

### now work on fixing particular elements ...
pp2 <- c(pp, list(b = rep(c(-1, 0, 1), length.out = 36)))

s1 <- simulate_new( ~ Days + (Days|Subject),
                   seed = 101,
                   newdata = sleepstudy,
                   family = gaussian,
                   newparams = pp2)[[1]]
plot(sleepstudy$Days, s1, ylim = c(180,450))
with(sleepstudy, points(Days, Reaction, col = 2))

s1 <- simulate_new( ~ Days + (Days|Subject),
                   seed = 101,
                   newdata = sleepstudy,
                   family = gaussian,
                   newparams = pp2,
                   return_val = "pars")
