## RTMB objective, link transformations, priors, and reporting.

cmb <- function(f, d) function(p) f(p, d)

osa_keep <- function(x) {
  if (inherits(x, "osa")) {
    as.vector(x@keep[, 1L])
  } else {
    rep(1, length(x))
  }
}

osa_value <- function(x) {
  if (inherits(x, "osa")) x@x else x
}

## Translated from logit_inverse_linkfun(), glmmTMB.cpp:213-232.
## The binomial likelihood uses logit(probability), not necessarily eta.
logit_inverse_linkfun_rtmb <- function(eta, link) {
  link_name <- link_name_rtmb(link)
  switch(
    link_name,
    logit = eta,
    probit = RTMB::pnorm(eta, log.p = TRUE) -
      RTMB::pnorm(eta, lower.tail = FALSE, log.p = TRUE),
    cloglog = RTMB::logspace_sub(exp(eta), 0),
    {
      mu <- switch(
        link_name,
        log = exp(eta),
        identity = eta,
        sqrt = eta * eta,
        inverse = 1 / eta,
        lambertW = exp(eta) * exp(exp(eta)),
        stop("link not yet implemented for binomial: ", link_name)
      )
      log(mu) - log1p(-mu)
    }
  )
}

## Translated from log_inverse_linkfun(), glmmTMB.cpp:234-249.
## Negative-binomial likelihoods use log(mu) for robust density evaluation.
log_inverse_linkfun_rtmb <- function(eta, link) {
  link_name <- link_name_rtmb(link)
  switch(
    link_name,
    log = eta,
    logit = -RTMB::logspace_add(0, -eta),
    {
      mu <- switch(
        link_name,
        probit = RTMB::pnorm(eta),
        cloglog = -expm1(-exp(eta)),
        identity = eta,
        sqrt = eta * eta,
        inverse = 1 / eta,
        lambertW = exp(eta) * exp(exp(eta)),
        stop("link not yet implemented for log inverse-link: ", link_name)
      )
      log(mu)
    }
  )
}

logit_mu_rtmb <- function(eta, link) {
  logit_inverse_linkfun_rtmb(eta, link)
}

log_mu_rtmb <- function(eta, link) {
  log_inverse_linkfun_rtmb(eta, link)
}

apply_zi_prediction <- function(mu, eta, etazi, ziPredictCode) {
  zi_predict_name <- zipredict_name_rtmb(ziPredictCode)
  if (zi_predict_name == "corrected") {
    pz <- 1 / (1 + exp(-etazi))
    mu <- mu * (1 - pz)
  } else if (zi_predict_name == "uncorrected") {
    ## leave mu and eta unchanged
  } else if (zi_predict_name == "prob") {
    mu <- 1 / (1 + exp(-etazi))
    eta <- etazi
  } else if (zi_predict_name == "disp") {
    ## handled separately by caller
  } else {
    stop("Invalid ziPredictCode: ", ziPredictCode)
  }

  list(mu = mu, eta = eta)
}

linkfun_rtmb <- function(mu, link) {
  link_name <- link_name_rtmb(link)
  switch(
    link_name,
    log = log(mu),
    identity = mu,
    sqrt = sqrt(mu),
    logit = log(mu) - log1p(-mu),
    probit = stats::qnorm(mu),
    cloglog = log(-log1p(-mu)),
    inverse = 1 / mu,
    lambertW = stop("linkfun for lambertW not yet implemented"),
    stop("link not yet implemented for prediction aggregation: ", link_name)
  )
}

family_name_rtmb <- function(family) {
  if (is.character(family) && length(family) == 1L) {
    return(family)
  }
  if (is.list(family) && !is.null(family$family)) {
    return(family$family)
  }

  family_name <- names(family)
  if (length(family_name) == 0L) {
    family_name <- names(.valid_family)[match(family, .valid_family)]
  }
  family_name
}

link_name_rtmb <- function(link) {
  if (is.character(link) && length(link) == 1L) {
    return(link)
  }

  link_name <- names(link)
  if (length(link_name) == 0L) {
    link_name <- names(.valid_link)[match(link, .valid_link)]
  }
  link_name
}

zipredict_name_rtmb <- function(ziPredictCode) {
  if (is.character(ziPredictCode) && length(ziPredictCode) == 1L) {
    return(ziPredictCode)
  }

  zi_predict_name <- names(ziPredictCode)
  if (length(zi_predict_name) == 0L) {
    zi_predict_name <- names(.valid_zipredictcode)[
      match(ziPredictCode, .valid_zipredictcode)
    ]
  }
  zi_predict_name
}

prior_nll <- function(beta, betazi, betadisp, theta, thetazi, psi,
                      prior_distrib, prior_whichpar, prior_distrib_name,
                      prior_whichpar_name, prior_elstart, prior_elend,
                      prior_npar, prior_params) {
  nll <- 0
  par_ind <- 1L

  if (length(prior_distrib_name) == 0L && length(prior_distrib) > 0L) {
    prior_distrib_name <- names(.valid_prior)[
      match(prior_distrib, .valid_prior)
    ]
  } else {
    prior_distrib_name <- as.character(prior_distrib_name)
  }
  if (length(prior_whichpar_name) == 0L && length(prior_whichpar) > 0L) {
    prior_whichpar_name <- names(.valid_vprior)[
      match(prior_whichpar, .valid_vprior)
    ]
  } else {
    prior_whichpar_name <- as.character(prior_whichpar_name)
  }

  for (i in seq_along(prior_distrib)) {
    parvec <- switch(
      prior_whichpar_name[i],
      beta = beta,
      betazi = betazi,
      betadisp = betadisp,
      theta = theta,
      thetazi = thetazi,
      psi = psi,
      stop("Unknown prior parameter vector name: ", prior_whichpar_name[i])
    )

    par_start <- prior_elstart[i] + 1L
    par_end <- prior_elend[i] + 1L
    if (par_start > par_end) {
      par_ind <- par_ind + prior_npar[i]
      next
    }
    if (par_start < 1L || par_end > length(parvec)) {
      stop(
        "Bad prior index for prior ", i, ": requested elements ",
        prior_elstart[i], ":", prior_elend[i],
        " in a parameter vector of length ", length(parvec)
      )
    }

    if (prior_distrib_name[i] == "lkj") {
      corpars <- parvec[par_start:par_end]
      nll <- nll - dlkj_rtmb(
        corpars,
        prior_params[par_ind],
        log = TRUE
      )
    } else {
      for (j in par_start:par_end) {
        parval <- parvec[j]
        logpriorval <- switch(
          prior_distrib_name[i],
          normal = RTMB::dnorm(
            parval,
            mean = prior_params[par_ind],
            sd = prior_params[par_ind + 1L],
            log = TRUE
          ),
          t = {
            location <- prior_params[par_ind]
            scale <- prior_params[par_ind + 1L]
            df <- prior_params[par_ind + 2L]
            RTMB::dt((parval - location) / scale, df = df, log = TRUE) -
              log(scale)
          },
          cauchy = dcauchy_rtmb(
            parval,
            location = prior_params[par_ind],
            scale = prior_params[par_ind + 1L],
            log = TRUE
          ),
          gamma = {
            shape <- prior_params[par_ind + 1L]
            scale <- prior_params[par_ind] / prior_params[par_ind + 1L]
            RTMB::dgamma(exp(parval), shape = shape, scale = scale,
                         log = TRUE)
          },
          stop("Prior distribution not implemented: ", prior_distrib_name[i])
        )
        nll <- nll - logpriorval
      }
    }
    par_ind <- par_ind + prior_npar[i]
  }
  nll
}

## Variables injected into rtmb_tpl() by RTMB::getAll()
utils::globalVariables(c(
  "X", "XS", "Z", "offset", "terms", "family", "link", "weights", "size",
  "beta", "b", "theta",
  "Xzi", "XziS", "Zzi", "zioffset", "termszi",
  "betazi", "bzi", "thetazi",
  "Xdisp", "XdispS", "Zdisp", "dispoffset", "termsdisp",
  "betadisp", "bdisp", "thetadisp",
  "psi", "combinom_disp_link", "ziPredictCode", "doPredict",
  "whichPredict", "aggregate",
  "prior_distrib", "prior_whichpar", "prior_elstart", "prior_elend",
  "prior_npar", "prior_params"
))

rtmb_tpl <- function(parameters, data) {
  RTMB::getAll(data, parameters)
  family_name <- data$family_name %||% family_name_rtmb(family)
  link_name <- data$link_name %||% link_name_rtmb(link)
  zi_predict_name <- zipredict_name_rtmb(ziPredictCode)
  ## Keep the original response for NA and structural-zero checks; OBS() may
  ## replace yobs with a simulation or OSA reference. During OSA calculations
  ## yobs is moved from data into parameters, so data$yobs may be absent.
  yobs_obs <- if (!is.null(data$yobs)) data$yobs else osa_value(yobs)
  yobs <- RTMB::OBS(yobs)

  nll <- 0

  ## Random-effects contribution; translated from glmmTMB.cpp:900-903
  cond_re <- allterms_nll(b, theta, terms)
  zi_re <- allterms_nll(bzi, thetazi, termszi)
  disp_re <- allterms_nll(bdisp, thetadisp, termsdisp)
  nll <- nll + cond_re$nll + zi_re$nll + disp_re$nll
  b <- cond_re$u
  bzi <- zi_re$u
  bdisp <- disp_re$u

  ## Conditional linear predictor and inverse link; adapted from
  ## glmmTMB.cpp:833, 911-918, and 934-937
  sparseX <- nrow(X) == 0 && ncol(X) == 0
  Xc <- if (sparseX) XS else X
  eta <- Xc %*% beta + Z %*% b + offset
  eta <- as.vector(eta)

  mu <- switch(
    link_name,
    log = exp(eta),
    identity = eta,
    sqrt = eta * eta,
    logit = 1 / (1 + exp(-eta)),
    probit = RTMB::pnorm(eta),
    cloglog = -expm1(-exp(eta)),
    inverse = 1 / eta,
    lambertW = exp(eta) * exp(exp(eta)),
    stop("link not yet implemented: ", link_name)
  )

  ## Zero-inflation linear predictor; adapted from
  ## glmmTMB.cpp:836, 880, and 919-925
  has_zi <- length(betazi) > 0 || length(bzi) > 0
  if (has_zi) {
    sparseXzi <- nrow(Xzi) == 0 && ncol(Xzi) == 0
    Xzic <- if (sparseXzi) XziS else Xzi
    etazi <- Xzic %*% betazi + Zzi %*% bzi + zioffset
    etazi <- as.vector(etazi)
  }

  ## Dispersion linear predictor; adapted from
  ## glmmTMB.cpp:839, 926-932, and 939
  sparseXdisp <- nrow(Xdisp) == 0 && ncol(Xdisp) == 0
  Xdispc <- if (sparseXdisp) XdispS else Xdisp
  etadisp <- Xdispc %*% betadisp + Zdisp %*% bdisp + dispoffset
  etadisp <- as.vector(etadisp)
  phi <- exp(etadisp)
  if (family_name == "combinomial" && combinom_disp_link == 1L) {
    phi <- etadisp
  }

  ## Observation likelihoods; adapted from glmmTMB.cpp:961-978,
  ## 1095-1101, and 1180-1199
  i <- !is.na(yobs_obs) | inherits(yobs, "simref")
  yobs_i <- yobs[i]
  keep <- osa_keep(yobs_i)
  eta_zi <- if (has_zi) etazi[i] else NULL
  logit_mu <- function() logit_mu_rtmb(eta, link_name)
  log_mu <- function() log_mu_rtmb(eta, link_name)
  log_var_minus_mu <- function() {
    log_var_minus_mu_rtmb(family_name, log_mu(), etadisp, psi)
  }

  tmp_loglik <- switch(
    family_name,
    poisson = dZI(RTMB::dpois)(yobs_i, lambda = mu[i], eta_zi = eta_zi, log = TRUE,
                               is_zero = yobs_obs[i] == 0),
    truncated_poisson = dZI(dtruncated_poisson_rtmb)(
      yobs_i, lambda = mu[i], eta_zi = eta_zi, log = TRUE,
      is_zero = yobs_obs[i] == 0
    ),
    gaussian = dZI(RTMB::dnorm)(
      yobs_i, mean = mu[i], sd = phi[i], eta_zi = eta_zi,
      log = TRUE, is_zero = yobs_obs[i] == 0
    ),
    ## Translated from the Gamma_family case in glmmTMB.cpp:991-996.
    Gamma = dZI(dgamma_rtmb)(
      yobs_i, mean = mu[i], shape = phi[i], eta_zi = eta_zi,
      log = TRUE, is_zero = yobs_obs[i] == 0
    ),
    ## Translated from the beta_family case in glmmTMB.cpp:997-1002.
    beta = dZI(dbeta_rtmb)(
      yobs_i, mean = mu[i], phi = phi[i], eta_zi = eta_zi,
      log = TRUE, is_zero = yobs_obs[i] == 0
    ),
    ## Translated from the ordbeta_family case in glmmTMB.cpp:1004-1031.
    ordbeta = dZI(dordbeta_rtmb)(
      yobs_i, eta = eta[i], mean = mu[i], phi = phi[i], cutpoints = psi,
      eta_zi = eta_zi, log = TRUE, is_zero = yobs_obs[i] == 0
    ),
    ## Translated from the lognormal_family case in glmmTMB.cpp:1164-1179.
    lognormal = dZI(dlognormal_rtmb)(
      yobs_i, mean = mu[i], sd = phi[i], eta_zi = eta_zi,
      log = TRUE, is_zero = yobs_obs[i] == 0
    ),
    ## Translated from the t_family case in glmmTMB.cpp:1182-1190.
    t = dZI(dt_rtmb)(
      yobs_i, mean = mu[i], scale = phi[i], df = exp(psi[1L]),
      eta_zi = eta_zi, log = TRUE, is_zero = yobs_obs[i] == 0
    ),
    ## Translated from the skewnormal_family case in glmmTMB.cpp:975-980.
    skewnormal = dZI(dskewnormal_rtmb)(
      yobs_i, mean = mu[i], sd = phi[i], alpha = psi[1L],
      eta_zi = eta_zi, log = TRUE, is_zero = yobs_obs[i] == 0
    ),
    ## Translated from the tweedie_family case in glmmTMB.cpp:1155-1163.
    tweedie = dZI(RTMB::dtweedie)(
      yobs_i, mu = mu[i], phi = phi[i],
      p = 1 / (1 + exp(-psi[1L])) + 1,
      eta_zi = eta_zi, log = TRUE, is_zero = yobs_obs[i] == 0
    ),
    ## Translated from the bell_family case in glmmTMB.cpp:1191-1200.
    bell = dZI(dbell_rtmb)(
      yobs_i, mean = mu[i], eta_zi = eta_zi,
      log = TRUE, is_zero = yobs_obs[i] == 0
    ),
    ## Translated from the binomial_family case in glmmTMB.cpp:979-983.
    binomial = dZI(dbinom_robust_rtmb)(
      yobs_i, size = size[i], logit_p = logit_mu()[i], eta_zi = eta_zi,
      log = TRUE, is_zero = yobs_obs[i] == 0),
    ## Translated from betabinomial_family, glmmTMB.cpp:1037-1047.
    betabinomial = {
      logit_p <- logit_mu()[i]
      dZI(dbetabinom_robust_rtmb)(
        yobs_i,
        log_shape1 = -RTMB::logspace_add(0, -logit_p) + etadisp[i],
        log_shape2 = -RTMB::logspace_add(0, logit_p) + etadisp[i],
        size = size[i], eta_zi = eta_zi, log = TRUE,
        is_zero = yobs_obs[i] == 0
      )
    },
    ## Translated from combinomial_family, glmmTMB.cpp:1049-1067.
    combinomial = dZI(dcombinom2_rtmb)(
      yobs_i, size = size[i], mean = mu[i] * size[i], nu = phi[i],
      eta_zi = eta_zi, log = TRUE, is_zero = yobs_obs[i] == 0
    ),
    ## Translated from the nbinom1_family case in glmmTMB.cpp:1042-1056.
    nbinom1 = dZI(dnbinom_robust_rtmb)(
      yobs_i, log_mu = log_mu()[i], log_var_minus_mu = log_var_minus_mu()[i],
      eta_zi = eta_zi, log = TRUE, is_zero = yobs_obs[i] == 0),
    ## Translated from truncated_nbinom1_family, glmmTMB.cpp:1042-1064.
    truncated_nbinom1 = dZI(dtruncated_nbinom1_rtmb)(
      yobs_i, log_mu = log_mu()[i], log_var_minus_mu = log_var_minus_mu()[i],
      log_phi = etadisp[i], eta_zi = eta_zi, log = TRUE,
      is_zero = yobs_obs[i] == 0),
    ## Translated from the nbinom2_family case in glmmTMB.cpp:1066-1075.
    nbinom2 = dZI(dnbinom_robust_rtmb)(
      yobs_i, log_mu = log_mu()[i], log_var_minus_mu = log_var_minus_mu()[i],
      eta_zi = eta_zi, log = TRUE, is_zero = yobs_obs[i] == 0),
    ## Translated from the nbinom12_family case in glmmTMB.cpp:1084-1094.
    nbinom12 = dZI(dnbinom_robust_rtmb)(
      yobs_i, log_mu = log_mu()[i], log_var_minus_mu = log_var_minus_mu()[i],
      eta_zi = eta_zi, log = TRUE, is_zero = yobs_obs[i] == 0),
    ## Translated from the genpois_family case in glmmTMB.cpp:1128-1133.
    genpois = dZI(dgenpois_rtmb)(
      yobs_i, theta = mu[i] / sqrt(phi[i]),
      lambda = 1 - 1 / sqrt(phi[i]), eta_zi = eta_zi,
      log = TRUE, is_zero = yobs_obs[i] == 0),
    ## Translated from truncated_genpois_family, glmmTMB.cpp:1134-1139.
    truncated_genpois = dZI(dtruncated_genpois_rtmb)(
      yobs_i, theta = mu[i] / sqrt(phi[i]),
      lambda = 1 - 1 / sqrt(phi[i]), eta_zi = eta_zi,
      log = TRUE, is_zero = yobs_obs[i] == 0),
    ## Translated from the compois_family case in glmmTMB.cpp:1115-1119.
    compois = dZI(dcompois2_rtmb)(
      yobs_i, mean = mu[i], nu = 1 / phi[i], eta_zi = eta_zi,
      log = TRUE, is_zero = yobs_obs[i] == 0),
    ## Translated from truncated_compois_family, glmmTMB.cpp:1121-1127.
    truncated_compois = dZI(dtruncated_compois2_rtmb)(
      yobs_i, mean = mu[i], nu = 1 / phi[i], eta_zi = eta_zi,
      log = TRUE, is_zero = yobs_obs[i] == 0),
    ## Translated from truncated_nbinom2_family, glmmTMB.cpp:1066-1081.
    truncated_nbinom2 = dZI(dtruncated_nbinom2_rtmb)(
      yobs_i, log_mu = log_mu()[i], log_var_minus_mu = log_var_minus_mu()[i],
      log_size = etadisp[i], eta_zi = eta_zi, log = TRUE,
      is_zero = yobs_obs[i] == 0),
    stop(
      "distribution not implemented yet for use with RTMB backend: ",
      family_name
    )
  )

  nll <- nll - sum(keep * weights[i] * tmp_loglik)

  ## Prior contribution; translated from glmmTMB.cpp:1203-1267
  nll <- nll + prior_nll(
    beta = beta,
    betazi = betazi,
    betadisp = betadisp,
    theta = theta,
    thetazi = thetazi,
    psi = psi,
    prior_distrib = prior_distrib,
    prior_whichpar = prior_whichpar,
    prior_distrib_name = data$rtmb_prior_distrib_name %||% character(0),
    prior_whichpar_name = data$rtmb_prior_whichpar_name %||% character(0),
    prior_elstart = prior_elstart,
    prior_elend = prior_elend,
    prior_npar = prior_npar,
    prior_params = prior_params
  )

  ## Prediction output; translated from glmmTMB.cpp:1353-1379
  mu_pred_all <- mu
  eta_pred_all <- eta

  ## Convert untruncated mean to the conditional mean of truncated distribution
  ## translated from glmmTMB.cpp:1331-1334
  if (family_name == "truncated_poisson") {
    mu_vector <- mu[seq_along(mu)]
    log_nzprob_pred <- log_nzprob_truncated_poisson_rtmb(mu_vector)
    mu_pred_all <- mu_pred_all / exp(log_nzprob_pred)
  } else if (family_name == "truncated_nbinom1") {
    mu_vector <- mu[seq_along(mu)]
    etadisp_vector <- etadisp[seq_along(etadisp)]
    log_nzprob_pred <- log_nzprob_truncated_nbinom1_rtmb(
      mu_vector,
      etadisp_vector
    )
    mu_pred_all <- mu_pred_all / exp(log_nzprob_pred)
  } else if (family_name == "truncated_nbinom2") {
    log_mu_value <- log_mu()
    log_mu_vector <- log_mu_value[seq_along(log_mu_value)]
    etadisp_vector <- etadisp[seq_along(etadisp)]
    log_nzprob_pred <- log_nzprob_truncated_nbinom2_rtmb(
      log_mu_vector,
      etadisp_vector
    )
    mu_pred_all <- mu_pred_all / exp(log_nzprob_pred)
  } else if (family_name == "truncated_genpois") {
    mu_vector <- mu[seq_along(mu)]
    phi_vector <- phi[seq_along(phi)]
    log_nzprob_pred <- log_nzprob_truncated_genpois_rtmb(
      mu_vector / sqrt(phi_vector)
    )
    mu_pred_all <- mu_pred_all / exp(log_nzprob_pred)
  } else if (family_name == "truncated_compois") {
    mu_vector <- mu[seq_along(mu)]
    phi_vector <- phi[seq_along(phi)]
    log_nzprob_pred <- log_nzprob_truncated_compois_rtmb(
      mu_vector,
      1 / phi_vector
    )
    mu_pred_all <- mu_pred_all / exp(log_nzprob_pred)
  }

  if (has_zi || zi_predict_name == "prob") {
    zi_pred <- apply_zi_prediction(
      mu = mu_pred_all,
      eta = eta_pred_all,
      etazi = etazi,
      ziPredictCode = zi_predict_name
    )
    mu_pred_all <- zi_pred$mu
    eta_pred_all <- zi_pred$eta
  }

  if (zi_predict_name == "disp") {
    mu_pred_all <- if (family_name == "Gamma") 1 / sqrt(phi) else phi
    eta_pred_all <- etadisp
  }

  mu_predict <- mu_pred_all[whichPredict]
  eta_predict <- eta_pred_all[whichPredict]

  if (length(aggregate) > 0) {
    if (length(aggregate) != length(mu_predict)) {
      stop(
        "'aggregate' wrong size; got length ", length(aggregate),
        " but prediction length is ", length(mu_predict)
      )
    }

    "[<-" <- RTMB::ADoverload("[<-")
    n_aggregate <- max(as.integer(aggregate))
    tmp <- rep(mu_predict[1L] * 0, n_aggregate)
    for (j in seq_along(mu_predict)) {
      tmp[as.integer(aggregate[j])] <- tmp[as.integer(aggregate[j])] +
        mu_predict[j]
    }

    mu_predict <- tmp
    eta_predict <- linkfun_rtmb(mu_predict, link_name)
  }

  corr <- cond_re$corr
  sd <- cond_re$sd
  corrzi <- zi_re$corr
  sdzi <- zi_re$sd
  corrdisp <- disp_re$corr
  sddisp <- disp_re$sd
  fact_load <- cond_re$fact_load

  REPORT(corr)
  REPORT(sd)
  REPORT(corrzi)
  REPORT(sdzi)
  REPORT(corrdisp)
  REPORT(sddisp)
  REPORT(fact_load)
  REPORT(b)
  REPORT(bzi)
  REPORT(bdisp)
  REPORT(mu_predict)
  REPORT(eta_predict)

  if (doPredict == 1) {
    ADREPORT(mu_predict)
  } else if (doPredict == 2) {
    ADREPORT(eta_predict)
  } else if (doPredict == 3) {
    ADREPORT(b)
    ADREPORT(bzi)
    ADREPORT(bdisp)
  }

  nll
}
