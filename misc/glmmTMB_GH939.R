library(glmmTMB)
sim_count <- simulate_new(~ mined + (1|site),
                  newdata = Salamanders,
                  family = betabinomial,
                  weights = rep(10, nrow(Salamanders)),
                  newparams = list(beta = c(2, 1),
                                   betad = log(2), ## log(NB dispersion)
                                   theta = log(1)) ## log(among-site SD))
)
