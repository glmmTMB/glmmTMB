## Test cases for the RTMB ordered beta family
## Fit each model with RTMB and legacy TMB, then compare likelihoods,
## parameters, predictions, and covariance estimates where applicable.

context("RTMB ordered beta backend")

skip_if_not_installed("RTMB")

tol_logLik <- 1e-5
tol_fixef <- 1e-5
tol_varcorr <- 1e-4

set.seed(2001)
ordbeta_dat <- expand.grid(
  obs = seq_len(5),
  g = factor(seq_len(30))
)
ordbeta_dat$x <- rnorm(nrow(ordbeta_dat))
ordbeta_dat$z <- rnorm(nrow(ordbeta_dat))
ordbeta_dat$off <- runif(nrow(ordbeta_dat), -0.15, 0.15)
ordbeta_dat$w <- runif(nrow(ordbeta_dat), 0.5, 1.5)
ordbeta_effect <- rnorm(nlevels(ordbeta_dat$g), sd = 0.25)
ordbeta_eta <- -0.2 + 0.35 * ordbeta_dat$x + ordbeta_effect[ordbeta_dat$g]
ordbeta_mu <- plogis(ordbeta_eta)
ordbeta_phi <- 8
ordbeta_cutpoints <- c(-1, 1)
ordbeta_u <- runif(nrow(ordbeta_dat))
ordbeta_p0 <- plogis(ordbeta_cutpoints[1L] - ordbeta_eta)
ordbeta_p1 <- plogis(ordbeta_eta - ordbeta_cutpoints[2L])
ordbeta_dat$y <- ifelse(
  ordbeta_u < ordbeta_p0,
  0,
  ifelse(
    ordbeta_u > 1 - ordbeta_p1,
    1,
    rbeta(
      nrow(ordbeta_dat),
      shape1 = ordbeta_mu * ordbeta_phi,
      shape2 = (1 - ordbeta_mu) * ordbeta_phi
    )
  )
)

ordbeta_start <- list(psi = ordbeta_cutpoints)
ordbeta_map <- list(psi = factor(c(NA, NA)))

expect_ordbeta_matches <- function(m_rtmb, m_tmb, check_varcorr = FALSE) {
  expect_equal(
    as.numeric(logLik(m_rtmb)),
    as.numeric(logLik(m_tmb)),
    tolerance = tol_logLik
  )
  expect_equal(
    fixef(m_rtmb)$cond,
    fixef(m_tmb)$cond,
    tolerance = tol_fixef
  )
  if (check_varcorr) {
    expect_equal(
      as.numeric(VarCorr(m_rtmb)$cond$g),
      as.numeric(VarCorr(m_tmb)$cond$g),
      tolerance = tol_varcorr
    )
  }
}

test_that("ordbeta: RTMB density matches C++ parameterization", {
  x <- c(0, 0.25, 0.75, 1)
  eta <- c(-0.5, -0.2, 0.2, 0.5)
  mu <- plogis(eta)
  phi <- rep(6, length(x))
  cutpoints <- c(-1, 1)

  expected <- c(
    log(plogis(cutpoints[1L] - eta[1L])),
    log(plogis(eta[2L] - cutpoints[1L]) -
          plogis(eta[2L] - cutpoints[2L])) +
      dbeta(x[2L], mu[2L] * phi[2L], (1 - mu[2L]) * phi[2L], log = TRUE),
    log(plogis(eta[3L] - cutpoints[1L]) -
          plogis(eta[3L] - cutpoints[2L])) +
      dbeta(x[3L], mu[3L] * phi[3L], (1 - mu[3L]) * phi[3L], log = TRUE),
    log(plogis(eta[4L] - cutpoints[2L]))
  )

  expect_equal(
    glmmTMB:::dordbeta_rtmb(
      x, eta = eta, mean = mu, phi = phi, cutpoints = cutpoints,
      log = TRUE
    ),
    expected,
    tolerance = 1e-12
  )
})

test_that("ordbeta: fixed conditional effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    family = ordbeta(),
    data = ordbeta_dat,
    start = ordbeta_start,
    map = ordbeta_map,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    family = ordbeta(),
    data = ordbeta_dat,
    start = ordbeta_start,
    map = ordbeta_map,
    se = FALSE
  )

  expect_ordbeta_matches(m_rtmb, m_tmb)
})

test_that("ordbeta: offsets, weights, and dispersion effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x + offset(off),
    dispformula = ~ z,
    weights = w,
    family = ordbeta(),
    data = ordbeta_dat,
    start = ordbeta_start,
    map = ordbeta_map,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x + offset(off),
    dispformula = ~ z,
    weights = w,
    family = ordbeta(),
    data = ordbeta_dat,
    start = ordbeta_start,
    map = ordbeta_map,
    se = FALSE
  )

  expect_ordbeta_matches(m_rtmb, m_tmb)
  expect_equal(fixef(m_rtmb)$disp, fixef(m_tmb)$disp, tolerance = tol_fixef)
})

test_that("ordbeta: conditional random intercept", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x + (1 | g),
    family = ordbeta(),
    data = ordbeta_dat,
    start = ordbeta_start,
    map = ordbeta_map,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x + (1 | g),
    family = ordbeta(),
    data = ordbeta_dat,
    start = ordbeta_start,
    map = ordbeta_map,
    se = FALSE
  )

  expect_ordbeta_matches(m_rtmb, m_tmb, check_varcorr = TRUE)
})

test_that("ordbeta: estimated cutpoints", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    family = ordbeta(),
    data = ordbeta_dat,
    start = ordbeta_start,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    family = ordbeta(),
    data = ordbeta_dat,
    start = ordbeta_start,
    se = FALSE
  )

  expect_ordbeta_matches(m_rtmb, m_tmb)
  rtmb_psi <- m_rtmb$obj$env$last.par.best[
    grepl("^psi", names(m_rtmb$obj$env$last.par.best))
  ]
  tmb_psi <- m_tmb$obj$env$last.par.best[
    grepl("^psi", names(m_tmb$obj$env$last.par.best))
  ]
  expect_equal(
    unname(rtmb_psi),
    unname(tmb_psi),
    tolerance = tol_fixef
  )
  expect_equal(
    family_params(m_rtmb),
    family_params(m_tmb),
    tolerance = tol_fixef
  )
})

test_that("ordbeta: zero inflation", {
  set.seed(2002)
  zi_dat <- ordbeta_dat
  zi_dat$y[rbinom(nrow(zi_dat), 1, plogis(-2 + 0.5 * zi_dat$z)) == 1] <- 0

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    ziformula = ~ z,
    family = ordbeta(),
    data = zi_dat,
    start = ordbeta_start,
    map = ordbeta_map,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    ziformula = ~ z,
    family = ordbeta(),
    data = zi_dat,
    start = ordbeta_start,
    map = ordbeta_map,
    se = FALSE
  )

  expect_ordbeta_matches(m_rtmb, m_tmb)
  expect_equal(fixef(m_rtmb)$zi, fixef(m_tmb)$zi, tolerance = tol_fixef)
})

test_that("ordbeta: prediction and simulation", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    family = ordbeta(),
    data = ordbeta_dat,
    start = ordbeta_start,
    map = ordbeta_map,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    family = ordbeta(),
    data = ordbeta_dat,
    start = ordbeta_start,
    map = ordbeta_map,
    se = FALSE
  )

  expect_equal(
    predict(m_rtmb, type = "response"),
    predict(m_tmb, type = "response"),
    tolerance = tol_fixef
  )
  expect_equal(
    predict(m_rtmb, type = "disp"),
    predict(m_tmb, type = "disp"),
    tolerance = tol_fixef
  )

  set.seed(2003)
  sim_y <- simulate(m_rtmb, nsim = 1)[[1]]
  expect_equal(length(sim_y), nrow(ordbeta_dat))
  expect_true(all(sim_y >= 0 & sim_y <= 1))
})

## Ordered-beta regressions for GH #1359, item 1.
ordbeta_link <- function(name) {
    if (name != "lambertW") return(name)
    ## W(z)*exp(W(z)) = z.
    ## With W(mu)=exp(eta), invert and differentiate this identity for the test link.
    structure(list(name="lambertW",
        linkfun=function(mu) log(.Call("lambertW_R",as.double(mu),PACKAGE="glmmTMB")),
        linkinv=function(eta) exp(eta)*exp(exp(eta)),
        mu.eta=function(eta) {
            exp(eta)*exp(exp(eta))*(1+exp(eta))
        },
        valideta=function(eta) all(is.finite(eta))),class="link-glm")
}
## Check E[Y^2]-E[Y]^2 independently of the centered implementation.
test_that("ordered-beta variance includes the two endpoint masses and precision", {
    mu <- c(0.2, 0.4, 0.7)
    phi <- c(2, 8, 12)
    cut <- c(-1.3, 0.9)
    links <- c("logit", "probit", "cloglog", "identity", "log", "inverse", "sqrt")
    for (link in links) {
        eta <- make.link(link)$linkfun(mu)
        p0 <- plogis(cut[1] - eta)
        p1 <- plogis(eta - cut[2])
        pc <- 1 - p0 - p1
        M <- p1 + pc * mu
        expected <- p1 + pc * (mu^2 + mu * (1 - mu)/(phi + 1)) - M^2
        v <- ordbeta(link)$variance
        expect_equal(v(mu, phi, cut), expected, tolerance = 1e-12)
        expect_equal(v(mu, phi, cut, eta), expected, tolerance = 1e-12)
    }
    link <- ordbeta_link("lambertW")
    eta <- link$linkfun(mu)
    p0 <- plogis(cut[1] - eta)
    p1 <- plogis(eta - cut[2])
    pc <- 1 - p0 - p1
    M <- p1 + pc * mu
    expected <- p1 + pc * (mu^2 + mu * (1 - mu)/(phi + 1)) - M^2
    expect_equal(ordbeta(link)$variance(mu, phi, cut), expected, tolerance=1e-12)
    expect_error(ordbeta()$variance(0.4), "phi|cutpoints")
})
test_that("ordered-beta variance remains finite when endpoint probabilities saturate", {
    eta <- c(-1000, 1000)
    expect_equal(ordbeta()$variance(plogis(eta), 8, c(-1, 1), eta), c(0, 0))
})
ordbeta_data <- function() {
    data.frame(y = c(rep(0, 200), rep(1, 160), qbeta(ppoints(440), 5, 5)))
}
ordbeta_fit <- function(link="logit", map=NULL) {
    link <- ordbeta_link(link)
    eta <- if (is.character(link)) make.link(link)$linkfun(0.5) else link$linkfun(0.5)
    glmmTMB(y ~ 1, data=ordbeta_data(), family=ordbeta(link),
        start=list(beta=eta, betadisp=log(10),
                   psi=eta+c(qlogis(0.25), -qlogis(0.2))), map=map)
}
ordbeta_moments <- function(mu, phi, eta, cut) {
    p0 <- plogis(cut[1]-eta)
    p1 <- plogis(eta-cut[2])
    pc <- 1-p0-p1
    M <- p1+pc*mu
    ## Raw second moment: the one adds p1; the beta adds pc*(mu^2+Var[B]).
    V <- p1+pc*(mu^2+mu*(1-mu)/(phi+1))-M^2
    list(p0=p0,p1=p1,pc=pc,mean=M,variance=V)
}
for (rtmb in c(FALSE,TRUE)) local({
    backend <- rtmb
    label <- if (backend) "RTMB" else "TMB"
    test_that(paste("ordered-beta simulation matches all mixture moments with",label), {
        skip_on_cran()
        local_useRTMB(backend)
        m <- ordbeta_fit()
        expect_equal(m$fit$convergence,0)
        expect_true(m$sdr$pdHess)
        eta <- predict(m,type="link")[[1L]]
        mu <- predict(m,type="conditional")[[1L]]
        cut <- glmmTMB:::getParList(m)$psi
        z <- ordbeta_moments(mu,sigma(m),eta,cut)
        s <- unlist(simulate(m,nsim=40,seed=101),use.names=FALSE)
        ## Indicator means have SE sqrt(p*(1-p)/N); response means have SE sqrt(V/N).
        ## The multipliers and slack below are simulation-test tolerances.
        N <- length(s)
        expect_equal(mean(s==0),z$p0,tolerance=6*sqrt(z$p0*(1-z$p0)/N)+2/N,scale=1)
        expect_equal(mean(s==1),z$p1,tolerance=6*sqrt(z$p1*(1-z$p1)/N)+2/N,scale=1)
        expect_equal(mean(s),z$mean,tolerance=6*sqrt(z$variance/N),scale=1)
        expect_equal(var(s),z$variance,tolerance=6/sqrt(N),scale=1)
    })
    test_that(paste("ordered-beta prediction, SE and residuals use full moments with",label), {
        skip_on_cran()
        local_useRTMB(backend)
        for (link in c("logit","probit","cloglog","identity","log","inverse","sqrt","lambertW")) {
            m <- ordbeta_fit(link)
            expect_equal(m$fit$convergence,0)
            expect_true(m$sdr$pdHess)
            eta <- predict(m,type="link")[[1L]]
            mu <- family(m)$linkinv(eta)
            cut <- glmmTMB:::getParList(m)$psi
            z <- ordbeta_moments(mu,sigma(m),eta,cut)
            expect_equal(predict(m,type="conditional"),rep(mu,nobs(m)),tolerance=1e-9)
            expect_equal(predict(m,type="response"),rep(z$mean,nobs(m)),tolerance=1e-9)
            expect_equal(predict(m,type="disp"),rep(sigma(m),nobs(m)),tolerance=1e-9)
            d <- ordbeta_data()
            logdensity <- ifelse(d$y==0,log(z$p0),
                ifelse(d$y==1,log(z$p1),log(z$pc)+dbeta(d$y,mu*sigma(m),(1-mu)*sigma(m),log=TRUE)))
            expect_equal(as.numeric(logLik(m)),sum(logdensity),tolerance=1e-8)
            full_mean <- function(pp) {
                pl <- m$obj$env$parList(pp,m$fit$parfull)
                e <- pl$beta[1]
                ordbeta_moments(family(m)$linkinv(e),exp(pl$betadisp[1]),e,pl$psi)$mean
            }
            ## Check gradient^T V gradient by differentiating the full mean numerically.
            gradient <- numDeriv::grad(full_mean,m$fit$par)
            expected_se <- sqrt(as.numeric(t(gradient)%*%m$sdr$cov.fixed%*%gradient))
            pred <- predict(m,newdata=data.frame(i=1:3),type="response",se.fit=TRUE)
            expect_equal(unname(pred$se.fit),rep(expected_se,3),tolerance=1e-5)
            expect_equal(unname(residuals(m,type="pearson")),
                (d$y-z$mean)/sqrt(z$variance),tolerance=1e-8)
            derivative <- numDeriv::grad(function(e)
                ordbeta_moments(family(m)$linkinv(e),sigma(m),e,cut)$mean,eta)
            expect_equal(unname(residuals(m,type="working")),
                (d$y-z$mean)/derivative,tolerance=1e-6)
        }
    })
    test_that(paste("ordered-beta cutoff estimates and SEs respect fixed and reordered maps with",label), {
        skip_on_cran()
        local_useRTMB(backend)
        maps <- list(NULL,
            list(psi=factor(c(NA,1))),
            list(psi=factor(c(1,NA))),
            list(psi=factor(c(NA,NA))),
            list(psi=factor(c(2,1),levels=1:2)))
        for (map in maps) {
            m <- ordbeta_fit(map=map)
            expect_equal(m$fit$convergence,0)
            expect_true(m$sdr$pdHess)
            expected <- plogis(glmmTMB:::getParList(m)$psi)
            expect_equal(unname(family_params(m)),expected)
            table <- summary(m)$cutoffs
            expect_identical(rownames(table),c("lower cutoff","upper cutoff"))
            ## Propagate free-parameter covariance through the mapped cutoffs: J V J^T.
            jac <- numDeriv::jacobian(function(pp)
                plogis(m$obj$env$parList(pp,m$fit$parfull)$psi),m$fit$par)
            se <- sqrt(diag(jac%*%m$sdr$cov.fixed%*%t(jac)))
            expect_equal(unname(table[,"Estimate"]),expected)
            expect_equal(unname(table[,"Std. Error"]),se,tolerance=1e-6)
            expect_equal(unname(family_params(m)),c(0.25,0.8),tolerance=1e-5)
            grad_mean <- numDeriv::grad(function(pp) {
                pl <- m$obj$env$parList(pp,m$fit$parfull)
                ordbeta_moments(plogis(pl$beta[1]),exp(pl$betadisp[1]),
                                pl$beta[1],pl$psi)$mean
            },m$fit$par)
            variance_mean <- as.numeric(t(grad_mean)%*%m$sdr$cov.fixed%*%grad_mean)
            pred <- predict(m,newdata=data.frame(i=1:3),type="response",cov.fit=TRUE)
            expect_equal(unname(pred$se.fit),rep(sqrt(variance_mean),3),tolerance=1e-5)
            expect_equal(pred$cov.fit,matrix(variance_mean,3,3),tolerance=1e-5)
            expect_output(print(m),"lower cutoff")
            expect_output(print(m),"upper cutoff")
            expect_output(print(summary(m)),"Ordered beta cutoffs")
            if (is.null(map)) {
                m0 <- update(m, data = ordbeta_data(), family = ordbeta(),
                             start = NULL, map = NULL, se = FALSE)
                expect_equal(suppressWarnings(summary(m0)$cutoffs), table,
                             tolerance = 1e-5)

            }
        }
    })
    test_that(paste("ordered-beta shared cutoffs retain their common uncertainty with",label), {
        skip_on_cran()
        local_useRTMB(backend)
        set.seed(2718)
        d <- data.frame(x=rep(seq(-1,1,length.out=5),each=80))
        d$y <- rbinom(nrow(d),1,plogis(0.6*d$x+0.2))
        m <- glmmTMB(y~0+x,data=d,family=ordbeta(),
            start=list(beta=0.6,psi=c(-0.2,-0.2),betadisp=log(8)),
            map=list(psi=factor(c(1,1)),betadisp=factor(NA)))
        expect_equal(m$fit$convergence,0)
        expect_true(m$sdr$pdHess)
        table <- summary(m)$cutoffs
        expect_equal(table[1,],table[2,])
        expect_true(all(table[,"Std. Error"]>0))
        expect_equal(unname(family_params(m)),rep(table[1,"Estimate"],2))
    })
    test_that(paste("ordered-beta mixture residuals include zero inflation with",label), {
        skip_on_cran()
        local_useRTMB(backend)
        d <- data.frame(y=c(rep(0,371),rep(1,248),qbeta(ppoints(381),5,5)))
        m <- glmmTMB(y~1,data=d,family=ordbeta(),ziformula=~1,
            start=list(beta=0,betadisp=log(10),psi=c(-1.3,0.8),betazi=qlogis(0.2)),
            map=list(psi=factor(c(NA,NA))))
        expect_equal(m$fit$convergence,0)
        expect_true(m$sdr$pdHess)
        e <- predict(m,type="link")[[1L]]
        mu <- predict(m,type="conditional")[[1L]]
        pi <- predict(m,type="zprob")[[1L]]
        expect_true(pi>0.1 && pi<0.3)
        cut <- glmmTMB:::getParList(m)$psi
        z <- ordbeta_moments(mu,sigma(m),e,cut)
        M <- (1-pi)*z$mean
        ## Condition first on the extra zero: within-group variance plus variance of group means.
        V <- (1-pi)*z$variance+pi*(1-pi)*z$mean^2
        expect_equal(predict(m,type="response"),rep(M,nrow(d)),tolerance=1e-9)
        expect_equal(unname(residuals(m,type="pearson")),(d$y-M)/sqrt(V),tolerance=1e-8)
        derivative <- (1-pi)*numDeriv::grad(function(eta)
            ordbeta_moments(plogis(eta),sigma(m),eta,cut)$mean,e)
        expect_equal(unname(residuals(m,type="working")),(d$y-M)/derivative,tolerance=1e-6)
    })



    test_that(paste("ordered-beta cutoff SEs handle an empty dispersion coefficient block with", label), {
        skip_on_cran()
        local_useRTMB(backend)
        d <- ordbeta_data()
        d$x <- rep(c(-1, 1), length.out = nrow(d))
        d$logphi <- log(10)
        m <- glmmTMB(y ~ x, data = d, family = ordbeta(),
            dispformula = ~0 + offset(logphi),
            start = list(beta = c(0, 0), psi = qlogis(c(0.25, 0.8))),
            map = list(beta = factor(c(1, 1))))
        expect_equal(m$fit$convergence, 0)
        expect_true(m$sdr$pdHess)
        table <- summary(m)$cutoffs
        expect_true(all(is.finite(table)))
        jac <- numDeriv::jacobian(function(pp)
            plogis(m$obj$env$parList(pp, m$fit$parfull)$psi), m$fit$par)
        expected_se <- sqrt(diag(jac %*% m$sdr$cov.fixed %*% t(jac)))
        expect_equal(unname(table[, "Std. Error"]), expected_se,
                     tolerance = 1e-6)
        expect_equal(table[, "Estimate"], family_params(m))
    })

})
