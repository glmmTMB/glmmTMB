## ------------------------------------------------------------------------
## Type-I error comparison on heteroscedastic data:
##   lmerTest::lmer               (Satterthwaite ddf, but NO heteroscedasticity
##                                  accommodation -- included as a baseline)
##   nlme::lme + varIdent          (heteroscedasticity via variance function,
##                                  ordinary t/F tests, no ddf correction)
##   glmmTMB, dispformula = ~trt, ddf = "asymptotic"     (LRT, no ddf correction)
##   glmmTMB, dispformula = ~trt, car::Anova Type II Wald chisq (no ddf correction;
##                                  the original blog post's uncorrected glmmTMB test)
##   glmmTMB, dispformula = ~trt, ddf = "satterthwaite"
##   glmmTMB, dispformula = ~trt, ddf = "kenward-roger"
##
## Design and data-generating process adapted from
## https://lgraz.com/posts/lmm-heteroskedastic/#type-i-error-rate-on-heteroscedastic-data
##
## Requires the glmmTMB version from https://github.com/glmmTMB/glmmTMB/pull/1302
## (branch "anova_ddf"), which adds a `ddf` argument to anova.glmmTMB for
## Kenward-Roger/Satterthwaite F-ratio tests comparing two nested models:
##   remotes::install_github("glmmTMB/glmmTMB", ref = "anova_ddf", subdir = "glmmTMB")
## ------------------------------------------------------------------------

library(glmmTMB)
library(lmerTest)
library(nlme)
library(car)
library(future)
library(future.apply)

if (!"ddf" %in% names(formals(getS3method("anova", "glmmTMB")))) {
    stop("this script requires the 'anova_ddf' branch of glmmTMB (PR #1302); install via\n",
         "remotes::install_github('glmmTMB/glmmTMB', ref = 'anova_ddf', subdir = 'glmmTMB')")
}

nsim   <- 5000
ncores <- 24
## multicore forks the current R process (fine on Linux/macOS in a
## non-interactive session); use plan(multisession, workers = ncores) instead
## on Windows or inside RStudio, where forking isn't supported.
plan(multicore, workers = ncores)

## ---- fixed design matrix (unbalanced, as in the original post) --------
## 40 ids x 10 obs, 4 trt levels (10 ids each); then randomly drop 130 of
## the 200 rows belonging to trt B/C so those groups are under-represented.
## This matrix is generated once and reused, unchanged, across all simulations.
set.seed(123)
X <- data.frame(
    id  = factor(rep(1:40, each = 10)),
    trt = factor(rep(c("A", "B", "C", "D"), each = 100))
)
X <- X[-sample(101:300, 130, replace = FALSE), ]

## ---- heteroscedastic null response -------------------------------------
## fixef = 0 for every trt level (i.e. we're simulating under H0: no trt
## effect); random intercept sd = 3 (Var = 9), residual sd = 1, matching
## the makeLmer()/simulate() setup in the original post but done directly
## with rnorm() so the script has no dependency on the 'simr' package.
## Extra N(0, 4) noise is added on top for trt groups C and D, so the
## conditional variance differs by treatment (dispformula = ~trt is needed
## to correctly model this; lmerTest and plain glmmTMB cannot).
get_hetero_data <- function(X) {
    D <- X
    b_id <- rnorm(nlevels(D$id), 0, 3)
    D$y  <- b_id[as.integer(D$id)] + rnorm(nrow(D), 0, 1)
    grpCD <- D$trt %in% c("C", "D")
    D$y[grpCD] <- D$y[grpCD] + rnorm(sum(grpCD), 0, 4)
    D
}

## last-column p-value, robust to the differing column layouts of
## Wald-z / Satterthwaite-t / Kenward-Roger-t / chisq / F tables
last_col <- function(tab, row) tab[row, ncol(tab)]

## ---- one fitting function per "batch"; each is run over nsim replicates ----
## with the SAME sequence of simulated datasets (set.seed(i) below), so all
## methods are compared on matched data.  Each returns one row per method
## (glmmTMB satterthwaite/kenward-roger share a single REML fit, since the
## ddf calculation is the only thing that differs between them).

fit_lmerTest <- function(D) {
    fit <- lmerTest::lmer(y ~ trt + (1 | id), data = D)
    A <- anova(fit)
    C <- coef(summary(fit))
    data.frame(method = "lmerTest", anova = last_col(A, "trt"), contrast = last_col(C, "trtC"))
}

fit_lme_varIdent <- function(D) {
    fit <- nlme::lme(y ~ trt, random = ~1 | id, data = D,
                      weights = varIdent(form = ~1 | trt))
    A <- anova(fit)
    C <- coef(summary(fit))
    data.frame(method = "lme_varIdent", anova = last_col(A, "trt"), contrast = last_col(C, "trtC"))
}

## anova.glmmTMB(ddf = ) compares two nested models rather than testing terms
## within a single model, so we need an explicit null (no trt) vs full (trt) pair.
## ddf = "asymptotic" is a likelihood-ratio test and needs ML (REML = FALSE)
## fits, since REML log-likelihoods aren't comparable across differing
## fixed-effect structures.
##
## Also reports "glmmTMB_wald": car::Anova(full, type = "II"), a single-model
## Wald chi-square test (no model comparison, no ddf correction at all) --
## this is the same uncorrected test used for glmmTMB in the original blog
## post (via `anova.glmmTMB <- glmmTMB:::Anova.glmmTMB`), included here so the
## LRT-based "no correction" and the post's original Wald-based "no
## correction" can both be compared against the ddf-corrected versions below.
## Both share the same full/null ML fits, so this costs no extra model fits.
fit_glmmTMB_ml <- function(D) {
    full <- glmmTMB(y ~ trt + (1 | id), dispformula = ~trt, data = D, REML = FALSE)
    null <- glmmTMB(y ~ 1   + (1 | id), dispformula = ~trt, data = D, REML = FALSE)
    A_lrt  <- suppressMessages(anova(null, full, ddf = "asymptotic"))
    A_wald <- car::Anova(full, type = "II")
    C <- coef(summary(full))$cond
    rbind(
        data.frame(method = "glmmTMB_none", anova = last_col(A_lrt, 2),       contrast = last_col(C, "trtC")),
        data.frame(method = "glmmTMB_wald", anova = last_col(A_wald, "trt"),  contrast = last_col(C, "trtC"))
    )
}

## Satterthwaite/Kenward-Roger F-tests are computed from a single shared REML
## fit of the null/full pair (kenward-roger requires REML; satterthwaite is
## conventionally computed on REML fits too, matching lmerTest's default).
fit_glmmTMB_reml <- function(D) {
    full <- glmmTMB(y ~ trt + (1 | id), dispformula = ~trt, data = D, REML = TRUE)
    null <- glmmTMB(y ~ 1   + (1 | id), dispformula = ~trt, data = D, REML = TRUE)
    out <- lapply(c(satterthwaite = "satterthwaite", `kenward-roger` = "kenward-roger"), function(ddf) {
        A <- suppressMessages(anova(null, full, ddf = ddf))
        C <- suppressMessages(coef(summary(full, ddf = ddf))$cond)
        data.frame(method = paste0("glmmTMB_", ddf),
                   anova = last_col(A, 2), contrast = last_col(C, "trtC"))
    })
    do.call(rbind, out)
}

method_specs <- list(
    list(name = "lmerTest",     fun = fit_lmerTest,       methods = "lmerTest"),
    list(name = "lme_varIdent", fun = fit_lme_varIdent,   methods = "lme_varIdent"),
    list(name = "glmmTMB_ml",   fun = fit_glmmTMB_ml,     methods = c("glmmTMB_none", "glmmTMB_wald")),
    list(name = "glmmTMB_reml", fun = fit_glmmTMB_reml,   methods = c("glmmTMB_satterthwaite", "glmmTMB_kenward-roger"))
)

## each replicate is timed individually (elapsed seconds for data generation +
## fitting), and the whole batch's wall-clock time is timed separately, so we
## can report per-package fit cost and check parallel scaling (batch wall time
## vs. naive serial time = sum of per-replicate times).
run_batch <- function(spec) {
    batch_t0 <- Sys.time()
    res <- future_lapply(seq_len(nsim), function(i) {
        rep_t0 <- Sys.time()
        set.seed(i)
        D <- get_hetero_data(X)
        out <- tryCatch(spec$fun(D), error = function(e) {
            data.frame(method = spec$methods, anova = NA_real_, contrast = NA_real_)
        })
        out$sim <- i
        out$elapsed <- as.numeric(difftime(Sys.time(), rep_t0, units = "secs"))
        out
    }, future.seed = TRUE)
    wall <- as.numeric(difftime(Sys.time(), batch_t0, units = "secs"))
    res <- do.call(rbind, res)
    cat(sprintf("[%s] nsim = %d, wall = %.1fs, mean/rep = %.2fs, naive serial = %.1fs\n",
                spec$name, nsim, wall, mean(res$elapsed[!duplicated(res$sim)]),
                sum(res$elapsed[!duplicated(res$sim)])))
    attr(res, "wall") <- wall
    res
}

script_t0 <- Sys.time()
batch_results <- lapply(method_specs, run_batch)
results_wide <- do.call(rbind, batch_results)
rownames(results_wide) <- NULL

## timing summary: per-batch wall time and mean per-replicate fit time
timing <- data.frame(
    batch          = vapply(method_specs, `[[`, "", "name"),
    nsim           = nsim,
    ncores         = ncores,
    wall_sec       = vapply(batch_results, attr, numeric(1), "wall")
)
per_rep <- do.call(rbind, lapply(seq_along(batch_results), function(k) {
    r <- batch_results[[k]]
    r <- r[!duplicated(r$sim), ]
    data.frame(batch = method_specs[[k]]$name, mean_sec_per_rep = mean(r$elapsed))
}))
timing <- merge(timing, per_rep, by = "batch")
cat("\n---- timing ----\n")
print(timing)
cat(sprintf("\ntotal script wall time: %.1fs\n",
            as.numeric(difftime(Sys.time(), script_t0, units = "secs"))))

## long form: one row per (method, sim, test type, p-value)
results_long <- reshape(
    results_wide,
    varying   = c("anova", "contrast"),
    v.names   = "p_value",
    timevar   = "test",
    times     = c("anova", "contrast"),
    direction = "long"
)
results_long$id <- NULL
rownames(results_long) <- NULL

## Type-I error summary (fraction of p < .05, with binomial test against 0.05)
tab <- aggregate(p_value ~ method + test, data = results_long,
                  FUN = function(p) mean(p < 0.05, na.rm = TRUE))
names(tab)[3] <- "type1_error"
tab$n_na <- aggregate(p_value ~ method + test, data = results_long,
                       FUN = function(p) sum(is.na(p)))$p_value
tab$binom_p <- mapply(function(m, t) {
    p <- results_long$p_value[results_long$method == m & results_long$test == t]
    p <- p[!is.na(p)]
    binom.test(sum(p < 0.05), length(p), 0.05)$p.value
}, tab$method, tab$test)

print(tab)

saveRDS(results_long, "hetero_ddf_sim_results.rds")
