library(ape)
library(glmmTMB)
library(reformulas)
data(carni70, package = "ade4")

carnidat <- carni70$tab

# load tree, make tip labeling match
tree <- read.tree(text = carni70$tre)
tree$tip.label <- gsub('\\.', '_', tree$tip.label)

# make ordering match
carnidat <- carnidat[tree$tip.label, ]

# phylogenetic variance-covariance matrix
## row/column names of phylo_varcov must match factor levels in data
phylo_varcov <- vcv(tree)

# convert species to factor (in order)
carnidat$species <- factor(rownames(carnidat), levels = rownames(carnidat))

# add dummy variable to specify a single cluster
carnidat$dummy <- factor(1)

## The propto() term fixes Sigma = lambda * K for a *known* matrix K
## (here, Brownian-motion covariance from a fixed tree), with lambda
## the only free (fitted) parameter -- so the "covariance function"
## for kriging is known analytically for any pair of tips, not fit
## from an empirical variogram. This is simple kriging with a known
## covariance function and an estimated linear trend, i.e. the
## approach of Garland & Ives (2000, Am Nat 155:346-364), and the
## prediction SE below follows their two-term decomposition (kriging
## variance of the phylogenetic effect + GLS uncertainty in the fixed
## effects), plus residual variance for a predicted observation
## rather than a predicted species mean.

##' Phylogenetic imputation for a glmmTMB propto() model
##'
##' @param formula a glmmTMB formula containing exactly one
##'   propto(<species> | <group>, <matrix>) term; the matrix argument
##'   given in the formula is not used (it is replaced internally by
##'   the training-species submatrix of \code{phylo_vcov}) but must be
##'   syntactically present
##' @param phylo_vcov full phylogenetic variance-covariance matrix
##'   (all species, training + missing), with row/column names
##'   matching the levels of the species variable in \code{data}
##' @param data data frame of tip-level data; rows with \code{NA} in
##'   the response column are the species to be imputed
##' @return \code{list(fit, se)}: named vectors (by row name of the
##'   missing rows in \code{data}) of predicted response values and
##'   prediction SEs. The SE combines all three sources of prediction
##'   uncertainty: the phylogenetic (kriging) variance of the missing
##'   species' random effect given the fitted ones, GLS uncertainty
##'   in the fixed-effect estimates, and residual variance.
phylo_impute <- function(formula, phylo_vcov, data) {

    resp_name <- deparse(formula[[2]])
    y <- data[[resp_name]]
    test_idx  <- is.na(y)
    train_idx <- !test_idx
    if (!any(test_idx)) stop("no missing (NA) values found in the response column")

    ## locate the (single) propto() term
    bars <- findbars_x(formula, specials = "propto", default.special = NULL)
    is_propto <- vapply(bars, function(x) is.call(x) && identical(x[[1]], quote(propto)), logical(1))
    if (sum(is_propto) != 1) stop("formula must contain exactly one propto() term")
    propto_call <- bars[[which(is_propto)]]
    bar_expr    <- propto_call[[2]]                 # e.g. 0 + species | dummy
    species_var <- all.vars(bar_expr[[2]])
    if (length(species_var) != 1)
        stop("expected a single species grouping variable in the propto() term")
    group_var <- deparse(bar_expr[[3]])

    species <- as.character(data[[species_var]])
    test_sp <- unique(species[test_idx])

    train_data <- data[train_idx, , drop = FALSE]
    train_data[[species_var]] <- factor(train_data[[species_var]])
    train_sp <- levels(train_data[[species_var]])

    if (!all(c(train_sp, test_sp) %in% rownames(phylo_vcov)))
        stop("not all species in 'data' are found in 'phylo_vcov'")

    K11 <- phylo_vcov[train_sp, train_sp, drop = FALSE]
    K21 <- phylo_vcov[test_sp,  train_sp, drop = FALSE]
    K22 <- phylo_vcov[test_sp,  test_sp,  drop = FALSE]

    ## rebuild the formula with the propto() matrix replaced by the training
    ## submatrix; glmmTMB() always resets environment(formula) <- parent.frame()
    ## (i.e. this function's own frame), so a plain local binding suffices
    fixed_formula <- nobars(formula)
    new_re_term   <- call("propto", bar_expr, quote(.phylo_vcov_train))
    fit_formula   <- addForm(fixed_formula, new_re_term)
    .phylo_vcov_train <- K11

    fit <- glmmTMB(fit_formula, data = train_data)

    vc         <- VarCorr(fit)$cond[[group_var]]
    lambda_hat <- unname(vc[1, 1] / K11[1, 1])

    b1 <- unlist(ranef(fit)$cond[[group_var]][1, ])
    names(b1) <- sub(paste0("^", species_var), "", names(b1))
    b1 <- b1[train_sp]

    beta_hat    <- fixef(fit)$cond
    vcov_beta   <- vcov(fit)$cond
    sigma_resid <- sigma(fit)

    ## kriging (BLUP) equations for the missing species' random effects
    K11inv  <- solve(K11)
    b2_hat  <- as.vector(K21 %*% K11inv %*% b1)
    names(b2_hat) <- test_sp
    condvar   <- lambda_hat * (K22 - K21 %*% K11inv %*% t(K21))
    krige_var <- diag(condvar)
    names(krige_var) <- test_sp

    ## fixed-effect part (+ its uncertainty) for the missing rows
    fe_terms_obj <- delete.response(terms(fixed_formula))
    Xtest <- model.matrix(fe_terms_obj, data = data[test_idx, , drop = FALSE])

    row_species  <- species[test_idx]
    fixed_part   <- as.vector(Xtest %*% beta_hat)
    beta_var_row <- diag(Xtest %*% vcov_beta %*% t(Xtest))

    fitted_vals <- fixed_part + b2_hat[row_species]
    se_vals     <- sqrt(krige_var[row_species] + beta_var_row + sigma_resid^2)

    row_names <- rownames(data)[test_idx]
    if (is.null(row_names) || anyDuplicated(row_names)) row_names <- which(test_idx)
    names(fitted_vals) <- names(se_vals) <- row_names

    list(fit = fitted_vals, se = se_vals)
}

## --- usage on the original (real) data / GH1309 example ---

carnidat_na <- carnidat
carnidat_na$log_range <- log(carnidat_na$range)
set.seed(101)
test_rows <- sample(nrow(carnidat_na), 20)
carnidat_na$log_range[test_rows] <- NA

imputed <- phylo_impute(log_range ~ log(size) + propto(0 + species | dummy, phylo_varcov),
                         phylo_vcov = phylo_varcov,
                         data = carnidat_na)
data.frame(fit = imputed$fit, se = imputed$se,
           actual = log(carnidat[names(imputed$fit), "range"]))

## --- validation with simulated data (known true parameters) ---
##
## Simulate from the same model structure via simulate_new(), using
## "pretty" rounded parameters loosely based on the real-data fit
## (beta ~ (1.3, 0.3), phylo SD ~ 0.7) but with residual SD reduced
## well below the real-data value, so the phylogenetic signal
## dominates and lambda (and hence the kriging predictions) should be
## well identified from a 50-species training set.

beta_target     <- c(1.3, 0.3)  # (Intercept), log(size)
phylo_sd_target <- 0.7
resid_sd_target <- 0.05
diagval <- mean(diag(phylo_varcov))  # tree is ultrametric: constant root-to-tip variance
lambda_target <- phylo_sd_target^2 / diagval

## theta has a single free entry (log(lambda)) since the rest of the
## propto covariance structure is fixed by the supplied matrix;
## no "b" component -> species-level random effects are resampled
newparams <- list(beta = beta_target, betadisp = log(resid_sd_target),
                   theta = log(lambda_target))

sim_y <- simulate_new(~ log(size) + propto(0 + species | dummy, phylo_varcov),
                       newdata = carnidat, family = gaussian,
                       newparams = newparams, seed = 1, return_val = "sim")[[1]]

carnidat_sim <- transform(carnidat, sim_logrange = sim_y)

set.seed(107)
test_rows_sim <- sample(nrow(carnidat_sim), 20)
carnidat_sim_na <- carnidat_sim
carnidat_sim_na$sim_logrange[test_rows_sim] <- NA

imputed_sim <- phylo_impute(sim_logrange ~ log(size) + propto(0 + species | dummy, phylo_varcov),
                             phylo_vcov = phylo_varcov,
                             data = carnidat_sim_na)
res_sim <- data.frame(fit = imputed_sim$fit, se = imputed_sim$se,
                       actual = carnidat_sim[names(imputed_sim$fit), "sim_logrange"])
print(res_sim)
cat("cor(fit, actual):", cor(res_sim$fit, res_sim$actual), "\n")
cat("RMSE:", sqrt(mean((res_sim$fit - res_sim$actual)^2)), "\n")
z <- (res_sim$fit - res_sim$actual) / res_sim$se
cat("coverage of 95% interval:", mean(abs(z) < 1.96), "\n")
