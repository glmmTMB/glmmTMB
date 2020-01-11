remotes::install_local(".")
library(glmmTMB)
set.seed(1)
nt <- min(parallel::detectCores(),5)
N <- 5e5
xdata <- rnorm(N, 1, 2)
groups <- 200
data_use <- data.frame(obs = 1:N)
data_use <- within(data_use,
{
  
  group_var <- rep(seq(groups), times = nrow(data_use) / groups)
  group_intercept <- rnorm(groups, 0, 0.1)[group_var]
  xdata <- xdata
  ydata <- 0.3 + group_intercept + 0.5*xdata + rnorm(N, 0, 0.25)
})

(t_serial <- system.time(
     model3 <- glmmTMB(formula = ydata ~ 1 + xdata + (1 | group_var),
                       data = data_use,
                       control = glmmTMBControl(parallel = 1))
 )
)
(t_parallel <- system.time(
     update(model3,  control = glmmTMBControl(parallel = nt))
 )
)

print(t_serial[["elapsed"]]/t_parallel[["elapsed"]])
