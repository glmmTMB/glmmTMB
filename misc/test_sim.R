devtools::load_all("~/R/pkgs/glmmTMB/glmmTMB")
debug(simulate_new)

pp <- list(beta = c(280, 3),
           betad = 1,
           theta = c(3, 0.5, 2))
data("sleepstudy", package = "lme4")
s1 <- simulate_new( ~ Days + (Days|Subject),
                   seed = 101,
                   newdata = sleepstudy,
                   family = gaussian,
                   newparams = pp)[[1]]
plot(sleepstudy$Days, s1, ylim = c(180,450))
with(sleepstudy, points(Days, Reaction, col = 2))

s1 <- simulate_new( ~ Days + (1|Subject),
             newdata = sleepstudy,
             family = gaussian,
             newparams = list(beta = c(200, 1),
                              betad = 1,
                              theta = rep(1,3)),
             return_val = "object")

s1$obj$env$last.par
