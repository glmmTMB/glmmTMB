b0 <- 5
b1 <- 1
sigma_resid <- 2
sigma_group <- 1

# simulate data
set.seed(54321)
x <- rnorm(30)
group <- rep(letters[1:6], each = 5)
resid <- rnorm(30, 0, sigma_resid)
group_effect <- rnorm(6, 0, sigma_group)
names(group_effect) <- letters[1:6]
y <- b0 + b1 * x + group_effect[group] + resid
dat <- data.frame(x = x, group = group, y = y)

# fit model with & without REML
fit_no_REML <- glmmTMB::glmmTMB(y ~ x + (1|group), data = dat, REML = FALSE)
fit_REML <- glmmTMB::glmmTMB(y ~ x + (1|group), data = dat, REML = TRUE)

# extract fixed effect coefs
coef_no_REML <- glmmTMB::fixef(fit_no_REML)$cond
coef_REML <- glmmTMB::fixef(fit_REML)$cond

# build a design matrix (DM) and a prediction data set with x = 0.5
pred_DM <- matrix(c(1, 0.5), nrow = 1)
pred_dat <- data.frame(x = 0.5, group = NA)

pred_DM %*% coef_no_REML                       # gives fixed effect prediction
predict(fit_no_REML, pred_dat, re.form = NA)   # returns fixed effect prediction
predict(fit_no_REML, pred_dat, re.form = ~0)   # returns fixed effect prediction
predict(fit_no_REML, pred_dat, re.form = NULL)

debug(predict.glmmTMB)
predict(fit_REML, se.fit = TRUE)
pred_DM %*% coef_REML                       # gives fixed effect prediction
predict(fit_REML, pred_dat, re.form = NA)   # doesn't return fixed effect prediction
predict(fit_REML, pred_dat, re.form = ~0)   # doesn't return fixed effect prediction
predict(fit_REML, pred_dat, re.form = NULL) # does return fixed effect prediction

mean(getME(fit_REML, "b"))
