## example of setting simcodes
library(glmmTMB)
library(tidyverse)
library(patchwork)
data("sleepstudy", package = "lme4")
m1 <- glmmTMB(Reaction ~ Days + (Days|Subject), data = sleepstudy)

## modify sim code and simulate
sim_fun <- function(model, simcode = "zero") {
    set_simcodes(model$obj, val = simcode, terms = "ALL")

    nvec <- seq(nrow(sleepstudy))
    ss0 <- sleepstudy |> mutate(n = nvec) |> select(-Reaction)
    ss <- (simulate(m1, nsim = 100, seed = 101)
        |> bind_rows()
        |> mutate(n = nvec)
        |> pivot_longer(-n, values_to = "Reaction")
        |> full_join(ss0, by = "n")
    )
    return(ss)
}

sim_zero <- sim_fun(m1, "zero")
sim_fix <- sim_fun(m1, "fix")
sim_random <- sim_fun(m1, "random")

gg0 <- (ggplot(sim_zero,
               aes(Days, Reaction, colour = Subject))
    + geom_point()
    + lims(y = c(100, 550))
    + facet_wrap(~Subject)
)

gg0 + (gg0 %+% sim_fix) + (gg0 %+% sim_random)
