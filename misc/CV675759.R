library(glmmTMB)
dd <- expand.grid(g = 1:50, time = 1:10) |>
    transform(ftime = factor(time)) |> 
    subset(! time %in% 7:9)

dd$y <- simulate_new( ~ 1 + ar1(0 + ftime | g),
               newdata = dd,
               newparams = list(beta = 0, theta = c(1,1), betadisp = -2),
               family = gaussian,
               seed = 101)[[1]]
devtools::load_all("glmmTMB")
debug(glmmTMB)
fit <- glmmTMB( y ~ 1 + ar1(0 + ftime | g), family = gaussian,
               data = dd)
str(model.frame(fit)$ftime)
fitC <- update(fit, control = glmmTMBControl(drop_unused_levels = FALSE))
str(model.frame(fitC)$ftime)
