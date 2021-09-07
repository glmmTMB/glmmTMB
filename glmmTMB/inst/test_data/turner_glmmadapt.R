##

if (require("GLMMadaptive")) {

    beta.binomial <- function (link = "logit") {
        .link <- link
        env <- new.env(parent = .GlobalEnv)
        assign(".link", link, envir = env)
        stats <- make.link(link)
        dbbinom <- function (x, size, prob, phi, log = FALSE) {
            A <- phi * prob
            B <- phi * (1 - prob)
            log_numerator <- lbeta(x + A, size - x + B)
            log_denominator <- lbeta(A, B)
            fact <- lchoose(size, x)
            if (log) {
                fact + log_numerator - log_denominator
            } else {
                exp(fact + log_numerator - log_denominator)
            }
        }
        log_dens <- function (y, eta, mu_fun, phis, eta_zi) {
            phi <- exp(phis)
            eta <- as.matrix(eta)
            mu_y <- mu_fun(eta)
            out <- if (NCOL(y) == 2L) {
                       dbbinom(y[, 1L], y[, 1L] + y[, 2L], mu_y, phi, TRUE)
                   } else {
                       dbbinom(y, rep(1L, length(y)), mu_y, phi, TRUE)
                   }
            attr(out, "mu_y") <- mu_y
            out
        }
        score_eta_fun <- function (y, mu, phis, eta_zi) {
            phi <- exp(phis)
            mu <- as.matrix(mu)
            if (NCOL(y) == 2L) {
                size <- y[, 1L] + y[, 2L]
                y <- y[, 1L]
            } else {
                size <- rep(1L, length(y))
            }
            phi_mu <- phi * mu
            phi_1mu <- phi * (1 - mu)
            comp1 <- (digamma(y + phi_mu) - digamma(size - y + phi_1mu)) * phi
            comp2 <- (digamma(phi_mu) - digamma(phi_1mu)) * phi
            mu.eta <- switch(.link,
                             "logit" = mu - mu * mu,
                             "cloglog" = - (1 - mu) * log(1 - mu))
            out <- (comp1 - comp2) * mu.eta
            out
        }
        score_phis_fun <- function (y, mu, phis, eta_zi) {
            phi <- exp(phis)
            mu <- as.matrix(mu)
            if (NCOL(y) == 2L) {
                size <- y[, 1L] + y[, 2L]
                y <- y[, 1L]
            } else {
                size <- rep(1L, length(y))
            }
            mu1 <- 1 - mu
            phi_mu <- phi * mu
            phi_1mu <- phi * mu1
            comp1 <- digamma(y + phi_mu) * mu + digamma(size - y + phi_1mu) * mu1 - 
                digamma(size + phi)
            comp2 <- digamma(phi_mu) * mu + digamma(phi_1mu) * mu1 - digamma(phi)
            out <- (comp1 - comp2) * phi
            out
        }
        structure(list(family = "beta binomial", link = stats$name, 
                       linkfun = stats$linkfun, linkinv = stats$linkinv, log_dens = log_dens,
                       score_eta_fun = score_eta_fun,
                       score_phis_fun = score_phis_fun),
                  class = "family")
    } ## beta.binomial family

    X <- readRDS(system.file("test_data","turner_bb.rds",
                             package="glmmTMB"))
    
    nvec <- 1:12
    pnames <- c(colnames(model.matrix(~ (Trt + 0) / Dose, data=X)),
                "nll")
    res <- matrix(NA,nrow=length(pnames),ncol=length(nvec),
                  dimnames=list(param=pnames,
                                nAGQ=nvec))

    for (i in seq_along(nvec)) {
        cat(nvec[i],"\n")
        fit <- mixed_model(cbind(Dead, Alive) ~ (Trt + 0) / Dose, data = X,
                                random = ~ Dose | Rep, 
                                family = beta.binomial(link = "cloglog"),
                                n_phis = 1,
                           nAGQ = nvec[i])
        res[,i] <- c(fixef(fit),c(-logLik(fit)))
    }
    par(las=1,bty="l")
    parmat <- res[-nrow(res),]
    matplot(nvec,t(parmat),type="l")
    saveRDS(res,file="turner_bb_GLMMadaptive.rds", version=2)
}
