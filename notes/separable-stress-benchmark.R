## Stress benchmark for separable member x AR(1) covariance structures.
##
## Run from the repository root with:
##   Rscript notes/separable-stress-benchmark.R
##
## The script supports both APIs:
## - current generic syntax:
##     separable(homcs(0 + member) %x% ar1(0 + time) | group)
##     separable(us(0 + member) %x% ar1(0 + time) | group)
## - archived specialized syntax:
##     homcsxar1(mt + 0 | group)
##     unxar1(mt + 0 | group)
##
## Optional environment variables:
##   BENCH_CASES: comma-separated case names
##   BENCH_GROUP_SCALE: multiplier for the number of simulated groups
##   BENCH_OUT: optional CSV output path

load_local_glmmTMB <- function() {
    pkg_dir <- file.path(getwd(), "glmmTMB")
    if (dir.exists(pkg_dir) && requireNamespace("devtools", quietly = TRUE)) {
        devtools::load_all(pkg_dir, quiet = TRUE)
    } else {
        library(glmmTMB)
    }
}

load_local_glmmTMB()

rho_to_theta <- function(rho) rho / sqrt(1 - rho^2)

case_specs <- list(
    homcs_2x12 = list(struc = "homcs", n_member = 2L, n_time = 12L,
                      n_group = 400L, sd = 1.4, phi = 0.55, rho = 0.35,
                      seed = 1001L),
    homcs_4x10 = list(struc = "homcs", n_member = 4L, n_time = 10L,
                      n_group = 300L, sd = 1.2, phi = 0.45, rho = 0.20,
                      seed = 1002L),
    us_2x12 = list(struc = "us", n_member = 2L, n_time = 12L,
                   n_group = 400L, sd = c(0.9, 1.5), phi = 0.50,
                   R_member = matrix(c(1, -0.30, -0.30, 1), 2, 2),
                   seed = 1003L),
    us_4x10 = list(struc = "us", n_member = 4L, n_time = 10L,
                   n_group = 300L, sd = c(0.8, 1.0, 1.25, 1.55),
                   phi = 0.40,
                   R_member = outer(seq_len(4), seq_len(4),
                                    function(i, j) 0.45^abs(i - j)),
                   seed = 1004L)
)

case_filter <- Sys.getenv("BENCH_CASES", unset = "")
if (nzchar(case_filter)) {
    keep <- strsplit(case_filter, ",", fixed = TRUE)[[1]]
    case_specs <- case_specs[keep]
    if (any(vapply(case_specs, is.null, logical(1)))) {
        stop("Unknown BENCH_CASES entry.")
    }
}

group_scale <- as.numeric(Sys.getenv("BENCH_GROUP_SCALE", unset = "1"))
if (!is.finite(group_scale) || group_scale <= 0) {
    stop("BENCH_GROUP_SCALE must be a positive number.")
}

member_corr <- function(spec) {
    if (spec$struc == "homcs") {
        R <- matrix(spec$rho, spec$n_member, spec$n_member)
        diag(R) <- 1
        R
    } else {
        spec$R_member
    }
}

make_dat <- function(spec) {
    set.seed(spec$seed)
    M <- spec$n_member
    T <- spec$n_time
    G <- as.integer(round(spec$n_group * group_scale))
    members <- factor(paste0("m", seq_len(M)), levels = paste0("m", seq_len(M)))
    times <- factor(seq_len(T), levels = seq_len(T))
    base <- expand.grid(member = members, time = times)

    R_member <- member_corr(spec)
    S_member <- diag(spec$sd, M, M) %*% R_member %*% diag(spec$sd, M, M)
    R_time <- outer(seq_len(T), seq_len(T),
                    function(i, j) spec$phi^abs(i - j))
    Sigma <- kronecker(R_time, S_member)
    y <- matrix(rnorm(G * M * T), nrow = G) %*% chol(Sigma)

    dd <- base[rep(seq_len(nrow(base)), G), , drop = FALSE]
    dd$group <- factor(rep(seq_len(G), each = nrow(base)))
    dd$member_i <- as.integer(dd$member)
    dd$time_i <- as.integer(dd$time)
    dd$mt <- numFactor(dd$member_i, dd$time_i)
    dd$y <- as.vector(t(y))
    dd
}

api_type <- function() {
    vc <- names(glmmTMB:::.valid_covstruct)
    if ("separable" %in% vc) return("generic")
    if (all(c("homcsxar1", "unxar1") %in% vc)) return("specialized")
    stop("Loaded glmmTMB does not expose either separable() or xar1 structures.")
}

case_formula <- function(spec, api) {
    if (api == "generic") {
        if (spec$struc == "homcs") {
            y ~ 1 + separable(homcs(0 + member) %x% ar1(0 + time) | group)
        } else {
            y ~ 1 + separable(us(0 + member) %x% ar1(0 + time) | group)
        }
    } else {
        if (spec$struc == "homcs") {
            y ~ 1 + homcsxar1(mt + 0 | group)
        } else {
            y ~ 1 + unxar1(mt + 0 | group)
        }
    }
}

extract_recovery <- function(fit, spec) {
    vc <- VarCorr(fit)$cond[[1]]
    sd_full <- as.numeric(attr(vc, "stddev"))
    R_full <- unname(attr(vc, "correlation"))
    M <- spec$n_member
    T <- spec$n_time
    idx <- matrix(seq_len(M * T), nrow = M, ncol = T)

    sd_hat <- rowMeans(matrix(sd_full, nrow = M, ncol = T))

    phi_vals <- c()
    for (m in seq_len(M)) {
        for (tt in seq_len(T - 1L)) {
            phi_vals <- c(phi_vals, R_full[idx[m, tt], idx[m, tt + 1L]])
        }
    }
    phi_hat <- mean(phi_vals)

    R_member_hat <- diag(1, M)
    for (i in seq_len(M - 1L)) {
        for (j in (i + 1L):M) {
            vals <- vapply(seq_len(T), function(tt) {
                R_full[idx[i, tt], idx[j, tt]]
            }, numeric(1))
            R_member_hat[i, j] <- R_member_hat[j, i] <- mean(vals)
        }
    }

    true_R <- member_corr(spec)
    list(
        sd_max_abs_err = max(abs(sd_hat - spec$sd)),
        phi_abs_err = abs(phi_hat - spec$phi),
        member_cor_max_abs_err = max(abs(R_member_hat - true_R)),
        sd_hat = paste(signif(sd_hat, 4), collapse = ";"),
        phi_hat = phi_hat,
        member_cor_hat = paste(signif(R_member_hat[lower.tri(R_member_hat)], 4),
                               collapse = ";")
    )
}

run_case <- function(name, spec, api) {
    dd <- make_dat(spec)
    form <- case_formula(spec, api)
    timing <- system.time({
        fit <- suppressWarnings(
            glmmTMB(form, data = dd, dispformula = ~0)
        )
    })
    rec <- extract_recovery(fit, spec)
    data.frame(
        api = api,
        case = name,
        struc = spec$struc,
        n_member = spec$n_member,
        n_time = spec$n_time,
        n_group = as.integer(round(spec$n_group * group_scale)),
        n_obs = nrow(dd),
        elapsed_sec = unname(timing[["elapsed"]]),
        convergence = fit$fit$convergence,
        objective = fit$fit$objective,
        sd_max_abs_err = rec$sd_max_abs_err,
        phi_abs_err = rec$phi_abs_err,
        member_cor_max_abs_err = rec$member_cor_max_abs_err,
        sd_hat = rec$sd_hat,
        phi_hat = rec$phi_hat,
        member_cor_hat = rec$member_cor_hat,
        row.names = NULL
    )
}

api <- api_type()
cat("API:", api, "\n")
cat("Cases:", paste(names(case_specs), collapse = ", "), "\n")
cat("Group scale:", group_scale, "\n")

res <- do.call(rbind, Map(run_case, names(case_specs), case_specs,
                          MoreArgs = list(api = api)))
print(res[, c("api", "case", "n_obs", "elapsed_sec", "convergence",
              "sd_max_abs_err", "phi_abs_err", "member_cor_max_abs_err",
              "sd_hat", "phi_hat", "member_cor_hat")],
      row.names = FALSE)

out <- Sys.getenv("BENCH_OUT", unset = "")
if (nzchar(out)) {
    write.csv(res, out, row.names = FALSE)
    cat("Wrote:", out, "\n")
}
