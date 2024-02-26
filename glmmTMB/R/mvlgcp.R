#' Main wrapper function for fitting uni- or multivariate, 2D spatial point process models via \code{glmmTMB}
#'
#' @description This function acts as a wrapper for \code{glmmTMB()}, permitting the modelling point patterns (in two dimensions) as a log-Gaussian Cox process (or Poisson process) using a similar approach to Dovers et al. (2024). The appropriate likelihood(s) are approximated using the approach of Berman and Turner (1992). Latent field(s) included in the Cox process are approximated via spatial basis functions while multivariate versions involve an additional approximation using a factor-analytic approach.
#'
#' @param formula a formula for the fixed & random effects, \emph{excluding} latent fields. The response of the formula must be a binary point event/quadrature identifier.
#' @param data a \code{data.frame} comprising of the response and predictors supplied within \code{formula}, as well as \code{coord.names}. The 'response' must be a binary that indicates whether a datum is a point event (1) or quadrature point (0) - see Dovers et al. 2024 for details.
#' @param weights a vector of quadrature weights with length equal to \code{nrow(data)} OR the name of the column within \code{data} where these can be found. These must be supplied for a point process model approximated in this way.
#' @param basis.functions a basis functions object. See \code{make_basis()}.
#' @param coord.names a vector of two character strings describing the column names of the coordinates (horizontal and vertical resp.) within \code{data}.
#' @param response.id a vector of length \code{nrow(data)} that categorises data into the various responses in the multivariate case. 
#' @param bf.matrix.type one of \code{"sparse"} or \code{"dense"} indicating whether to use a sparse/dense matrix to represent the basis function approximation to the latent field(s). The default of \code{"sparse"} is usually preferred, although in some instances (e.g. with single-response TPRS basis) using a dense matrix may be more efficient.
#' @param n_factors an integer describing the number of latent factors to be used to characterise the multivariate latent fields within the multivariate log-Gaussian Cox process.
#' @param ... other terms to be passed to \code{glmmTMB()}.
#'
#' @details
#' An approximate log-Gaussian Cox process will be fitted to the data when \code{basis.functions} are supplied. Otherwise, a Poisson process framework is used.
#' Additionally if \code{response.id} is supplied a multivariate LGCP is fitted with \code{n_factor} independent latent fields that provide correlated, response-specific fields via factor loadings (using the existing \code{rr()} framework in \code{glmmTMB}).
#' 
#' Note that for a Poisson process \code{response.id} is not used and a multivariate form should be specified via the formula, e.g. \code{response ~ (1|response.id) + (0 + covariate|response.id)} will fit random intercepts and slopes (for the "\code{covariate}") according to \code{response.id}.
#'
#'\code{basis.functions} should be an object created via \code{make_basis()} which currently supports those from \code{mgcv}, \code{FRK}, and \code{scampr} packages. There is scope to include any kind of basis functions, this would require writing an appropriate method within \code{bf_matrix()}.
#' 
#' @return a fitted \code{glmmTMB} model object with some additional list items, specific to the point process model.
#' @export
#' 
#' @references
#' Berman, M. and Turner, T. R. (1992), Approximating Point Process Likelihoods with GLIM, \emph{Journal of the Royal Statistical Society. Series C} (Applied Statistics), \bold{41}, 31–38.
#'
#' Dovers, E and Stoklosa, J and Warton, D. I. (2024) Fitting log-Gaussian Cox processes using generalized additive model software, \emph{The American Statistician}, DOI: 10.1080/00031305.2024.2316725 
#'
#' @examples
#' # using a subset of the Lansing Wood data:
#' dat1 <- lansing[lansing$tree == "blackoak", ]
#' dat2 <- lansing[lansing$tree %in% c("blackoak", "hickory", "maple"), ]
#' dat2$tree <- factor(as.character(dat2$tree)) # remove unused tree species
#' 
#' # set up the basis functions with dimension k = 100
#' bfs <- make_basis(k = 100, dat1)
#' 
#' # fit univariate LGCP to the blackoak point pattern data
#' m1 <- mvlgcp(pt ~ 1, data = dat1, weights = dat1$wt, basis.functions = bfs)
#' # fit a multivariate LGCP to the point patterns of "blackoak", "hickory" and "maple" trees
#' m2 <- mvlgcp(pt ~ (1 | tree), data = dat2, weights = dat2$wt, basis.functions = bfs, response.id = dat2$tree)
#' 
#' # compare the latent fields for blackoaks with/without accounting for between species correlations
#' 
#' # first we need to specify the domain (in this example the quadrature points of the data are a regular grid of the domain of interest)
#' domain <- dat1[dat1$pt == 0, ]
#' 
#' # extract the fields
#' fld1 <- get_field(m1, newdata = domain)
#' fld1_mv <- get_field(m2, newdata = domain, which.response = 1)
#' fld2_mv <- get_field(m2, newdata = domain, which.response = 2)
#' fld3_mv <- get_field(m2, newdata = domain, which.response = 3)
#' 
#' # plot the fields
#' par(mfrow = c(2, 2), mar = c(0,0,1,1.5))
#' plot(vec2im(fld1, domain$x, domain$y), main = "Univariate LGCP\nBlack Oak")
#' points(dat1[dat1$pt == 1, c("x", "y")])
#' plot(vec2im(fld1_mv, domain$x, domain$y), main = "Multivariate LGCP\nBlack Oak")
#' points(dat2[dat2$pt == 1 & dat2$tree == "blackoak", c("x", "y")])
#' plot(vec2im(fld2_mv, domain$x, domain$y), main = "Multivariate LGCP\nHickory")
#' points(dat2[dat2$pt == 1 & dat2$tree == "hickory", c("x", "y")], pch = 2)
#' plot(vec2im(fld3_mv, domain$x, domain$y), main = "Multivariate LGCP\nMaple")
#' points(dat2[dat2$pt == 1 & dat2$tree == "maple", c("x", "y")], pch = 3)
#' par(mfrow = c(1, 1), mar = c(5.1,4.1,4.1,2.1))
mvlgcp <- function(formula, data, weights = NULL, basis.functions, coord.names = c("x", "y"), response.id, bf.matrix.type = c("sparse", "dense"), n_factors = 2, ...) {
  
  # check the basis function matrix type
  bf.matrix.type <- match.arg(bf.matrix.type)
  
  mc <- match.call() # gets the arguments (must be updated for a point process model as below)
  call.list <- as.list(mc)
  
  # check the form of the weights
  object.supplied <- tryCatch(!is.null(weights), error = function(e) FALSE)
  if (!object.supplied) {
    weight.name <- deparse(substitute(weights))
    if (weight.name == "NULL") {
      stop("weights must be supplied for fitting a LGCP.\n\nThese must be quadrature weights in rows where the formula response = 0,\nThese must be some small number (e.g. 1e-6) in rows where the formula response = 1.")
    } else {
      wt.vec <- as.vector(data[,weight.name])
    }
  } else {
    wt.vec <- weights
  }
  
  # checks
  if ((!all(coord.names %in% colnames(data)))) {
    stop(paste0("At least one of 'coord.names', ", paste(coord.names, collapse = " or "), ", not found 'data'"))
  }
  if (!missing(basis.functions)) {
    if ((!all(coord.names == attr(basis.functions, "coord.names")))) {
      stop(paste0("At least one of the supplied 'coord.names', ", paste(coord.names, collapse = " or "), ", do not match the coord.names associated with the 'basis.functions' provided."))
    }
  }
  # if (!missing(family)) { # TODO: missing doesn't work for ... arguments
  #   warning("'family' argument supplied will be forced to 'poisson()' as needed for these LGCP or Poisson process")
  # }
  
  # get the response variable out of the formula
  resp <- all.vars(formula[[2]])
  data$new.response <- data[ , resp] / wt.vec
  
  # alter the call according to requirements for a point process model
  call.list$family <- poisson()
  call.list$data <- data
  call.list$weights <- wt.vec
  call.list$formula <- update.formula(formula, "new.response ~ .")
  # remove the function of the call
  call.list[[1]] <- NULL
  
  # check for basis functions (this indicates whether fitting a Poisson or Cox process)
  if (missing(basis.functions)) {
    
    # remove the unecessary arguments for glmmTMB
    call.list$basis.functions <- NULL
    call.list$coord.names <- NULL
    call.list$response.id <- NULL
    call.list$bf.matrix.type <- NULL
    call.list$n_factors <- NULL
    
    # NOTE: there is no need to differentiate single or mutliple responses for a Poisson process (this should be specified via random slopes/intercepts via the formula)
    tmp.time <- system.time(assign("mod", do.call("glmmTMB", call.list)))
    # add in some extra info
    mod$call <- mc
    mod$timing <- tmp.time
    mod$n_responses <- 1
    mod$n_basis <- 0
    mod$n_factors <- 0
    mod$basis.functions <- NULL
    mod$coord.names <- coord.names
    mod$weights <- wt.vec
    mod$response.id <- rep(1, nrow(data))

  } else {
    
    # check for response id (this indicates whether a multivariate LGCP is required)
    if (missing(response.id)) {
      
      ## Univariate LGCP ##

      # determine the number of basis functions, k, depending on type
      if (is(basis.functions, "mgcv.bf")) {
        k <- ncol(basis.functions$re$rand$Xr)
      } else if (is(basis.functions, "scampr.bf")) {
        k <- nrow(basis.functions)
      } else if (is(basis.functions, "Basis")) {
        k <- basis.functions@n
      } else {
        stop("supplied 'basis.functions' are not a recognised class")
      }

      # set the basis function matrix
      newZ <- bf_matrix(basis.functions, point.locations = data[ , coord.names], bf.matrix.type = bf.matrix.type)
      if (is(basis.functions, "mgcv.bf")) {
        call.list$data$basis.fixed <- as.matrix(newZ[, 1:sum(basis.functions$re$pen.ind == 0)])
        newZ <- newZ[, (sum(basis.functions$re$pen.ind == 0) + 1):ncol(newZ)]
      }
      
      # set the dummy factor for the basis functions and smoothing penalty
      call.list$data$basis.functions <- factor(rep(1:k, length.out = nrow(data)))
      call.list$data$b <- rnorm(nrow(data)) # this needs to be any continuous variable

      # update the formula to include basis functions (i.e. approx. latent field)
      if (is(basis.functions, "mgcv.bf")) {
        call.list$formula <- update(call.list$formula, ~ . + basis.fixed + (0 + b|basis.functions))
      } else {
        call.list$formula <- update(call.list$formula, ~ . + (0 + b|basis.functions))
      }
      
      # create an indicator for whether the model should be fitted based on the original user-setting
      not_fitting <- FALSE
      if (!is.null(call.list$doFit)) {
        if (!call.list$doFit) {
          not_fitting <- TRUE
        }
      }
      # ensure model is not fitted before we replace the basis function matrix
      call.list$doFit = FALSE
      # remove the unecessary arguments for glmmTMB
      call.list$basis.functions <- NULL
      call.list$coord.names <- NULL
      call.list$response.id <- NULL
      call.list$bf.matrix.type <- NULL
      call.list$n_factors <- NULL
      
      # set the reduced rank model structure
      m.init <- do.call("glmmTMB", call.list)
      
      # check that the dimensions of the other random effects match:
      if (any(names(m.init$condReStruc) != "0 + b | basis.functions")) {
        if (!sum(sapply(m.init$condReStruc[names(m.init$condReStruc) != "0 + b | basis.functions"], function(x){x[[1]] * x[[2]]})) == ncol(m.init$data.tmb$Z) - ncol(newZ)){
          stop("dimension of additional random effects are incorrect, please check the model spec.")
        }
      } else {
        if (!sapply(m.init$condReStruc, function(x){x[[1]] * x[[2]]}) == ncol(newZ)){
          stop("internal error in the dimension of random effects approximating the latent field, please check the model spec.")
        }
      }
      
      # replace the random effect matrix
      if (any(names(m.init$condReStruc) != "0 + b | basis.functions")) {
        m.init$data.tmb$Z <- cbind(m.init$data.tmb$Z[,1:(ncol(m.init$data.tmb$Z) - ncol(newZ))], newZ)
      } else {
        m.init$data.tmb$Z <- newZ
      }
      
      # fit the model
      if (not_fitting) {
        mod <- m.init
        tmp.time <- NA
      } else {
        tmp.time <- system.time(assign("mod", glmmTMB::fitTMB(m.init)))
      }
      # add in some extra info
      mod$call <- mc
      mod$timing <- tmp.time
      mod$n_responses <- 1
      mod$n_basis <- 0
      mod$n_factors <- 0
      mod$basis.functions <- basis.functions
      mod$coord.names <- coord.names
      mod$weights <- wt.vec
      mod$response.id <- rep(1, nrow(data))

    } else {
      
      # check that there are more than 1 response.id provided
      if (length(unique(response.id)) == 1) {
        warning("only one response provided in 'response.id'. If a univariate LGCP is desired, omit 'response.id' argument.")
      }
      
      ## Multivariate LGCP ##
      
      # re-order the data according to the response.id
      call.list$data <- call.list$data[order(response.id), ]
      call.list$weights <- call.list$weights[order(response.id)]
      call.list$data$response.id <- factor(sort(response.id))
      responses <- levels(call.list$data$response.id)

      # get the number of basis functions and responses
      if (is(basis.functions, "mgcv.bf")) {
        k <- ncol(basis.functions$re$rand$Xr)
        # if (attr(basis.functions, "rm.fixed")) {
        #   k <- basis.functions$sm[[1]]$bs.dim - ncol(basis.functions$re$Xf)
        # } else {
        #   if (attr(basis.functions, "rm.intercept")) {
        #     k <- basis.functions$sm[[1]]$bs.dim - 1
        #   } else {
        #     k <- basis.functions$sm[[1]]$bs.dim
        #   }
        # }
      } else if (is(basis.functions, "scampr.bf")) {
        k <- nrow(basis.functions)
      } else if (is(basis.functions, "Basis")) {
        k <- basis.functions@n
      } else {
        stop("supplied 'basis.functions' are not a recognised class")
      }
      m <- length(responses)
      
      # set up the basis function matrix
      newZ <- as(matrix(0, nrow = nrow(data), ncol = m * k), "sparseMatrix") # working in sparse for efficiency (or will this be more inefficient if the matrix ends up dense?)
      for (i in 1:m) {
        # calculate the response specific basis function matrix component
        newZ[call.list$data$response.id == responses[i], (1:k) + (k * (i-1))] <- bf_matrix(basis.functions, point.locations = call.list$data[call.list$data$response.id == responses[i], coord.names], bf.matrix.type = bf.matrix.type)
      }
      gc() # heavy work above so garbage collect may help
      
      # set up column indices that reflect the order of coefficients by responses (for using reduced rank structure within glmmTMB)
      col.idx <- rep(1:k, each = m)
      for (i in 2:m) {
        col.idx[seq(i, m*k, by = m)] <- (1:k) + (k*(i-1))
      }
      # adjust the column order to match that required by glmmTMB
      newZ <- newZ[,col.idx]
      
      # set the dummy factor to initialise the correct dimension for Z
      call.list$data$basis.functions <- factor(rep(1:k, length.out = nrow(call.list$data)))
      
      # update the formula to include basis functions (for some reason update() seems to add unwanted brackets # MAEVE MAY WANT TO KNOW THIS)
      # form <- update(formula, ~ . + rr(0 + response.id|dummy_k, d))
      call.list$formula <- as.formula(paste0(paste(deparse(call.list$formula), collapse = ""), " + rr(0 + response.id|basis.functions, n_factors)"))
      
      # create an indicator for whether the model should be fitted based on the original user-setting
      not_fitting <- FALSE
      if (!is.null(call.list$doFit)) {
        if (!call.list$doFit) {
          not_fitting <- TRUE
        }
      }
      # ensure model is not fitted before we replace the basis function matrix
      call.list$doFit = FALSE
      # remove the unnecessary arguments for glmmTMB
      call.list$basis.functions <- NULL
      call.list$coord.names <- NULL
      call.list$response.id <- NULL
      call.list$bf.matrix.type <- NULL
      call.list$n_factors <- NULL
      
      # set the reduced rank model structure
      rrm.init <- do.call("glmmTMB", call.list)
      
      # check that the dimensions of the other random effects match:
      if (any(names(rrm.init$condReStruc) != "0 + response.id | basis.functions")) {
        if (!sum(sapply(rrm.init$condReStruc[names(rrm.init$condReStruc) != "0 + response.id | basis.functions"], function(x){x[[1]] * x[[2]]})) == ncol(rrm.init$data.tmb$Z) - ncol(newZ)){
          stop("dimension of additional random effects are incorrect, please check the model spec.")
        }
      } else {
        if (!sapply(rrm.init$condReStruc, function(x){x[[1]] * x[[2]]}) == ncol(newZ)){
          stop("internal error in the dimension of random effects approximating the latent field, please check the model spec.")
        }
      }
      
      # replace the random effect matrix
      if (any(names(rrm.init$condReStruc) != "0 + response.id | basis.functions")) {
        rrm.init$data.tmb$Z <- cbind(rrm.init$data.tmb$Z[,1:(ncol(rrm.init$data.tmb$Z) - ncol(newZ))], newZ)
      } else {
        rrm.init$data.tmb$Z <- newZ
      }
      
      # set the starting pars to be calculated (gets initial factor loadings from DS residuals within rr)
      rrm.init$control$start_method$method = "res"
      
      # free up some extra space
      rm(newZ)
      gc()
      
      # fit the model
      if (not_fitting) {
        mod <- rrm.init
        tmp.time <- NA
      } else {
        tmp.time <- system.time(assign("mod", glmmTMB::fitTMB(rrm.init)))
      }
      # add in some extra info
      mod$call <- mc
      mod$timing <- tmp.time
      mod$n_responses <- m
      mod$n_basis <- k
      mod$n_factors <- n_factors
      mod$basis.functions <- basis.functions
      mod$coord.names <- coord.names
      mod$weights <- weights
      mod$response.id <- response.id
      mod$col.idx <- col.idx

    }
  }
  
  return(mod)
}

#' Create 2D basis functions for use in \code{glmmTMB::mvlgcp()}
#'
#' @description Uses data provided to place points across the range of the provided coordinates. Options include basis functions as in the packages \code{"scampr"} (see Dovers et al. 2023), \code{"mgcv"} (see Wood 2017), or \code{"FRK"} (see Zammit-Mangion and Cressie 2021). Note that the final basis dimension may not match \code{k}, particularly on unequal coordinate ranges - trial-and-error may be needed to get your desired basis dimension.
#'
#' @param k integer describing the desired basis dimension.
#' @param data a data frame containing the two coordinates described by \code{coord.names}.
#' @param from.package a character string that describing the package from which to create basis functions.
#' @param coord.names vector of character strings describing the names of the coordinates in 'data'. Ordered horizontal then vertical axes
#' @param longlat specific to \code{scampr} basis functions, a logical indicating whether the coordinates are in Longitude and Latitude so that geodesic distances are used (the radius is also extended to the maximum geodesic distance between nodes). Defaults to \code{FALSE}.
#' @param raw specific to \code{mgcv} basis functions, a logical of whether to return raw or transformed values (mgcv does the latter by default, assists fitting algorithms to converge). Defaults to \code{FALSE}.
#' @param ... additional arguments to be added for the underlying basis construction. See \code{mgcv::smoothCon()} and \code{FRK::auto_basis()} for example.
#' 
#' @return an object of class specific to the basis function type required. To be passed to bf_matrix()
#' @export
#'
#' @importFrom sp spDists
#' @importFrom mgcv smoothCon smooth2random
#' @importFrom FRK auto_basis
#'
#' @references
#' Wood, S. N. (2017), Generalized additive models: an introduction with R, \emph{CRC Press}.
#' 
#' Zammit-Mangion, A. and Cressie, N. (2021), FRK: An R package for spatial and spatio-temporal prediction with large datasets, \emph{Journal of Statistical Software}, \bold{98}, 1–48.
#' 
#' Dovers, E., Brooks, W., Popovic, G. C., and Warton, D. I. (2023), Fast, Approximate Maximum Likelihood Estimation of Log-Gaussian Cox Processes, \emph{Journal of Computational and Graphical Statistics}, \bold{32}, 1660–1670.
#'
#' @examples
#' # Create basis function nodes on the locations of presence records and quadrature
#' bfs <- make_basis(k = 100, data = gorillas)
make_basis <- function(k, data, from.package = c("scampr", "mgcv", "FRK"), coord.names = c("x", "y"), longlat = FALSE, raw = FALSE, ...) {
  if (!all(coord.names %in% colnames(data))) {
    stop("at least one of 'coord.names' not found in the data provided")
  }
  from.package <- match.arg(from.package)
  
  if (from.package == "scampr") {
    
    nodes.on.long.edge <- round(sqrt(k))
    
    # Get the data ranges
    xrange <- range(data[ , coord.names[1]], na.rm = T)
    yrange <- range(data[ , coord.names[2]], na.rm = T)
    # Get the axis sizes
    dx <- diff(xrange)
    dy <- diff(yrange)
    # Get the largest axis identifier
    big.axis.id <- which.max(c(dx, dy))
    small.axis.id <- c(1, 2)[!c(1, 2) %in% big.axis.id]
    # Jitter the first node placement by 1/4 the distance between nodes (so basis fn doesn't land on an edge point)
    big.axis.jitter <- c(dx, dy)[big.axis.id] / (nodes.on.long.edge * 4)
    big.axis.nodes <- seq(c(xrange[1], yrange[1])[big.axis.id] - big.axis.jitter, c(xrange[2], yrange[2])[big.axis.id] + big.axis.jitter, length.out = nodes.on.long.edge)
    # Create small axis nodes using the same distance between nodes
    small.axis.nodes <- seq(c(xrange[1], yrange[1])[small.axis.id] - big.axis.jitter, c(xrange[2], yrange[2])[small.axis.id] + big.axis.jitter, by = diff(big.axis.nodes)[1])
    node.list <- list()
    node.list[[big.axis.id]] <- big.axis.nodes
    node.list[[small.axis.id]] <- small.axis.nodes
    basis.functions <- as.data.frame(expand.grid(node.list[[1]], node.list[[2]]))
    colnames(basis.functions) <- coord.names
    basis.functions$scale <- diff(big.axis.nodes)[1] * 1.5
    
    # if dealing with longitudes/latitudes, change the radius to reflect km distances on a curved surface
    if (longlat) {
      max_diag_dists <- NULL
      for (i in 1:nrow(basis.functions)) {
        if (i %% nodes.on.long.edge != 0) {
          tmp_bloc <- basis.functions[c(i:(i+1), i:(i+1) + nodes.on.long.edge), 1:2]
          tmp_dists <- sp::spDists(as.matrix(na.omit(tmp_bloc)), as.matrix(na.omit(tmp_bloc)), longlat = longlat)
          max_diag_dists[i] <- max(tmp_dists)
        } else {
          max_diag_dists[i] <- 0
        }
      }
      basis.functions$scale <- max(max_diag_dists) * 1.5
      attr(basis.functions, "longlat") <- TRUE
    } else {
      attr(basis.functions, "longlat") <- FALSE
    }
    basis.functions$res <- 1
    class(basis.functions) <- c(class(basis.functions), "scampr.bf")
    
  } else if (from.package == "mgcv") {
    
    # set up the bivariate smooth
    eval(parse(text = paste0("sm <- mgcv::smoothCon(s(", coord.names[1], ", ", coord.names[2], ", k = k, ...), data = data)")))
    re <- mgcv::smooth2random(sm[[1]], "", type = "2")
    basis.functions <- list(sm = sm, re = re, coord.names = coord.names)
    class(basis.functions) <- "mgcv.bf"
    attr(basis.functions, "raw") <- raw
    
  } else if (from.package == "FRK") {
    
    # convert the data to a spatial points dataframe
    sp::coordinates(data) <- coord.names
    basis.functions <- FRK::auto_basis(data = data, max_basis = k, ...)
    # class(basis.functions) <- c(class(basis.functions), "mgcv.bf") # unfortunately this messes up the FRK object interfacing
    
  }
  
  # retain coord.names of the basis functions (attribute maintains the functionality of FRK bfs)
  attr(basis.functions, "coord.names") <- coord.names
  
  return(basis.functions)
}

#' Function to evaluate basis functions (created via \code{make_basis()}) at point locations provided and return them in a matrix
#' 
#' @description Calculates a matrix of basis function values at the point locations supplied (rows are locations and columns are basis functions). Primarily used internally for model fitting, also used to predict the basis at new locations.
#'
#' @param basis.functions an object created by \code{make_basis()}
#' @param point.locations either a matrix or data.frame describing the point locations at which the basis functions are to be evaluated.
#' @param bf.matrix.type a character string, one of 'sparse' or 'dense' indicating whether to create sparse or dense matrix.
#'
#' @return a matrix (sparse or dense depending on \code{bf.matrix.type})
#' @export
#'
#' @importFrom FRK eval_basis
#' @importFrom methods as
#' @importFrom sp spDists
#' @importFrom mgcv PredictMat
#'
#' @examples
#' data("gorillas")
#' # Create basis function nodes on the locations of presence records and quadrature
#' bfs <- make_basis(k = 100, data = gorillas)
#' # Calculate the basis function matrix at the point events in gorillas
#' X <- bf_matrix(bfs, gorillas[gorillas$pres == 1, c("x", "y")])
bf_matrix <- function(basis.functions, point.locations, bf.matrix.type = c("sparse", "dense")) {
  
  ## WANT TO ADD IN FUNCTIONALITY TO MAKE BF MATRIX ORTHOGONAL TO FIXED EFFECT MATRIX AS IN MGCV
  # Uses: PX = X(X′X)−1X′ and QX = I − PX
  
  bf.matrix.type <- match.arg(bf.matrix.type)
  
  # alter approach based on the type of basis functions
  if (is(basis.functions, "Basis")) { # NOTE THESE ARE FRK BASIS FUNCTIONS
    
    # evaluate the basis function matrix using FRK::eval_basis
    bf.mat <- FRK::eval_basis(basis = basis.functions, as.matrix(point.locations))
    # since these default to sparse we need to convert to dense if needed
    if (bf.matrix.type == "dense") {
      bf.mat <- methods::as(bf.mat, "matrix")
    }

  } else if (is(basis.functions, "scampr.bf")) {
    
    # initialise the basis function matrix
    bf.mat <- NULL
    # loop through each resolution of basis functions
    for (res in unique(basis.functions$res)) {
      # obtain the radius
      radius <- basis.functions$scale[basis.functions$res == res][1]
      # calculate the distances from points to basis function nodes
      if (nrow(point.locations) == 0) { # allows the edge case where there are no points (could arise in cross-validation for example)
        dist.mat <- matrix(NA, nrow = 0, ncol = nrow(as.matrix(basis.functions[,1:2][basis.functions$res == res, ])))
      } else {
        dist.mat <- sp::spDists(as.matrix(point.locations), as.matrix(basis.functions[,1:2][basis.functions$res == res, ]), longlat = attr(basis.functions, "longlat"))
      }
      # We want to get these matrices into sparse format ASAP to save memory and computation times
      
      # record instances of any point exactly on the nodes of the basis functions
      pts.on.nodes <- methods::as(dist.mat == 0, "sparseMatrix")
      # set distances beyond the radius to zero
      dist.mat[dist.mat > radius] <- 0
      # if sparse we can save computation time here
      if (bf.matrix.type == "sparse") {
        dist.mat <- methods::as(dist.mat, "sparseMatrix")
      }
      # calculate the bi-square basis function matrix
      Z <- ((dist.mat != 0 | pts.on.nodes) - (dist.mat / radius)^2)^2
      # add resolution to matrix via columns
      bf.mat <- cbind(bf.mat, Z)
      rm(Z)
    }
    # ensure the matrix is sparse if required
    if (bf.matrix.type == "sparse") { # TODO: check the current form of the bf.mat
      bf.mat <- methods::as(bf.mat, "sparseMatrix")
    }
  } else  if (is(basis.functions, "mgcv.bf")) {
    
    if (length(basis.functions$re$rand) > 1) {
      warning("Multiple random effects from 'mgcv::smoothCon()' detected.\nNote 'make_basis(..., type = 'mgcv')' builds a basis function matrix for a single random effect - only the first from the list provided will be used.")
    }

    # set the smooth information
    sm <- basis.functions$sm[[1]]
    # set the random effect information
    re <- basis.functions$re
    
    # Return either the raw (back-transformed basis functions) or otherwise
    if (attr(basis.functions, "raw")) {
      ##############################################################################
      # THE FOLLOWING IS TAKEN DIRECTLY FROM mgcv::smoothCon() HELP FILES:
      X <- mgcv::PredictMat(sm, point.locations) ## get prediction matrix for new data
      ## transform to r.e. parameterization
      if (!is.null(re$trans.U)) X <- X%*%re$trans.U
      X <- t(t(X)*re$trans.D)
      ## re-order columns according to random effect re-ordering...
      X[,re$rind] <- X[,re$pen.ind!=0] 
      ## re-order penalization index in same way  
      pen.ind <- re$pen.ind; pen.ind[re$rind] <- pen.ind[pen.ind>0]
      ## start return object...
      r <- list(rand=list(),Xf=X[,which(re$pen.ind==0),drop=FALSE])
      for (i in 1:length(re$rand)) { ## loop over random effect matrices
        r$rand[[i]] <- X[,which(pen.ind==i),drop=FALSE]
        attr(r$rand[[i]],"s.label") <- attr(re$rand[[i]],"s.label")
      }
      names(r$rand) <- names(re$rand)
      ##############################################################################
      
      # set the fixed and random components of the basis
      X_fixed <- r$Xf
      X_random <- r$rand[[1]]
      
    } else {
      
      X <- mgcv::PredictMat(sm, point.locations) ## get prediction matrix for new data
      
      # set the fixed and random components of the basis
      X_fixed <- X[ , re$pen.ind == 0]
      X_random <- X[ , re$pen.ind != 0]
      
    }
    
    # construct the full basis matrix
    bf.mat <- unname(cbind(X_fixed, X_random))
    # ensure the matrix is sparse if required
    if (bf.matrix.type == "sparse") { # TODO: check the current form of the bf.mat
      bf.mat <- methods::as(bf.mat, "sparseMatrix")
    }
    attr(bf.mat, "pen.ind") <- c(rep(0, ncol(X_fixed)), rep(1, ncol(X_random)))
    
  } else {
    
    stop("Basis functions provided are not of the correct type. See documentation for details")
    
  }
  
  return(bf.mat)
}

#' Prune basis functions (created via \code{make_basis()}) that do not intersect the data supplied
#' 
#' @description Computes the basis function matrix and checks for columns that do not have at least \code{n_non_zero} non-zero entires. If \code{prune=TRUE} basis functions not meeting the requirement will be removed.
#'
#' @param basis.functions an object created by \code{make_basis()}
#' @param data a data frame containing the coordinates associated with the supplied \code{basis.functions}
#' @param n_non_zero an integer describing the number of non-zero points required for pruning. Default is 1.
#' @param prune a logcial indicating whether to prune the basis functions, or simply identify problem ones. Default is \code{TRUE}.
#'
#' @return a basis functions object as in \code{make_basis} (or a vector of indices describing problem basis functions if \code{prune=FALSE})
#' @export
#'
#' @importFrom FRK eval_basis
#' @importFrom methods as
#' @importFrom sp spDists
#' @importFrom mgcv PredictMat
#'
#' @examples
#' data("gorillas")
#' # Create basis function nodes on the locations of presence records and quadrature
#' bfs <- make_basis(k = 100, data = gorillas)
#' # prune the basis functions according to a subset of the gorillas data
#' new.bfs <- prune(bfs, gorillas[1:1000, ])
prune_basis <- function(basis.functions, data, n_non_zero = 1, prune = TRUE) {
  
  bf.mat <- bf_matrix(basis.functions, data[,attr(basis.functions, "coord.names")])
  
  # if we only require a single non-zero entry we can use colSums which is computationally more efficient
  if (n_non_zero == 1) {
    check <- data.frame(basis.fn = 1:ncol(bf.mat), prune = colSums(as.matrix(bf.mat)) == 0)
  } else {
    check <- data.frame(basis.fn = 1:ncol(bf.mat), prune = apply(bf.mat, 2, function(x){sum(x != 0)}) < n_non_zero)
  }
  
  if (prune) {
    if (!is(basis.functions, "scampr.bf")) {
      stop("pruning is currently only available for basis.functions built with from.package='scampr'.\nEither change these or set prune=FALSE to identify problem basis functions.")
    } else {
      check <- basis.functions[!check$prune, ]
    }
  }

  return(check)
}

#' Compute the influence matrix of a fitted \code{glmmTMB} model
#'
#' @description This function computes the influence matrix from a fitted \code{glmmTMB} model along the lines of the package \code{mgcv} (Wood 2017).
#'
#' @param model a fitted \code{glmmTMB} model object.
#' @param ... additionally arguments for the S3 generic
#' @param by.obs a logical indicating whether the influence should be calculated by observations (\code{=TRUE}) or by terms (\code{=FALSE}). Since the function extracts the diagonal of \eqn{A}, with the goal to compute its trace (for example to calculate effective degrees of freedom), when \code{by.obs=FALSE} the formula for \eqn{A} (below) is re-arranged to avoid some additional matrix multiplication - this effectively calculates the influence by terms rather than observations.
#' 
#' @details
#' The Laplace approximation to a marginalised likelihood function is essentially a penalised likelihood so that there is an equivalence between GAMs (as fitted via \code{mgcv}) and GLMMs as fitted within \code{glmmTMB}.
#' This means we can use the approach in \code{mgcv} to obtain the influence (hat) matrix (\eqn{\boldsymbol{A}}) of the fitted model. \deqn{\boldsymbol{A} = \boldsymbol{X}(\boldsymbol{X}^{T}\boldsymbol{W}\boldsymbol{X} + \boldsymbol{S})^{-1}\boldsymbol{X}^{T}\boldsymbol{W}}
#' where \eqn{\boldsymbol{X}} is the \eqn{n\times(p + k)} fixed (\eqn{p}) and random (\eqn{k}) effect design matrices and \eqn{\boldsymbol{W}} is the "weights" matrix which is diagonal with \deqn{[\boldsymbol{W}]_{i,i}=\frac{\text{weights}}{V(\mu_{i})g'(\mu_{i})^{2}}}
#' where \code{weights} are the observation weights of the fitted model.
#' 
#' \eqn{\boldsymbol{S}} is the \eqn{(p + k)\times(p + k)} penalty matrix comprising blocks for each of the fixed and random components where the fixed effect block is zero everywhere (as these are unpenalised terms). The random effects component has a block diagonal structure where blocks represent independent sets of random effects. For the \eqn{j^{th}} set of random terms, the block is equal to \eqn{\lambda_{j}\boldsymbol{C}^{-1}_{b_{j}}}
#' where \eqn{\boldsymbol{C}_{b_{j}}} is the correlation matrix of the random effects (within set \eqn{j}) while the smoothing penalty is the dispersion parameter (\eqn{\phi}) squared, divided by the \eqn{j^{th}} random effect variance (\eqn{\sigma^{2}_{b_{j}}}), i.e. \deqn{\lambda_{j} = \frac{\phi^{2}}{\sigma^{2}_{b_{j}}}}
#' 
#' 
#' @return a vector of influence of length observations in the fitted model (or terms if \code{by.obs=FALSE})
#' @export
#'
#' @importFrom Matrix bdiag Diagonal
#'
#' @references
#' Wood, S. N. (2017), Generalized additive models: an introduction with R, \emph{CRC Press}.
#'
#' @examples
#' # fit a GLMM
#' m <- glmmTMB(count ~ spp + (1 | site), data = Salamanders, family = poisson)
#' # fit the equivalent GAM
#' library(mgcv)
#' m_gam <- gam(count ~ spp + s(site, bs = "re"), data = Salamanders, family = poisson, method = "ML")
#' # compare influence
#' cbind(glmmTMB = influence(m), mgcv = m_gam$hat)
#' # compare effective degrees of freedom
#' c(glmmTMB = sum(influence(m)), mgcv = sum(m_gam$hat))
influence.glmmTMB <- function(model, ..., by.obs = TRUE) {
  # get the call
  call.list <- as.list(model$call)
  # indicate not to fit the model
  call.list$doFit <- FALSE
  # get the model structure (by re-running the model without fitting)
  mod_str <- eval(as.call(call.list)) # NOTE: this gets the corrected Z for smoothers since it evaluates mvlgcp()
  
  # FOR MVLGCP: if the model is from mvlgcp() and has a rr() structure we need to re-order the Z matrix
  if (call.list[[1]] == "mvlgcp" & !is.null(model$col.idx)) {
    # separate out the basis functions from the other random effects
    other_re <- mod_str$data.tmb$Z[ , 1:(ncol(mod_str$data.tmb$Z) - length(mod_str$col.idx))]
    newZ <- mod_str$data.tmb$Z[ , ((ncol(mod_str$data.tmb$Z) - length(mod_str$col.idx)) + 1):ncol(mod_str$data.tmb$Z)]
    # re-order the basis function component according to the column indexing
    mod_str$data.tmb$Z <- cbind(other_re, newZ[ , order(model$col.idx)])
  }
  
  # check for any dropped fixed effects (due to, e.g., duplicate intercepts from mgcv basis)
  n_fixed <- length(fixef.glmmTMB(model)$cond)
  if (n_fixed != sum(names(model$fit$par) == "beta")) {
    term.dropped <- TRUE
    dropped.term.idx <- unname(is.na(fixef.glmmTMB(model)$cond))
    n_fixed <- sum(!dropped.term.idx)
    mod_str$data.tmb$X <- mod_str$data.tmb$X[ , !dropped.term.idx]
  }
  
  # combine the design matrices (fixed and random) # TODO: could include zi and dispersion components down the track?
  X <- as.matrix(cbind(mod_str$data.tmb$X, mod_str$data.tmb$Z))

  # set up the penalty matrix S #
  
  # check for the presence of random effects
  if (ncol(mod_str$data.tmb$Z) == 0) {
    S <- matrix(data = 0, nrow = n_fixed, ncol = n_fixed) # set S as an unpenalised matrix for fixed effects only
  } else {
    # get the random effect covariance matrices
    Sigmas_theta <- VarCorr.glmmTMB(model)
    # get the dispersion parameter
    disp.par <- sigma.glmmTMB(model)^2 # TODO: Squaring here will work for Gaussian and (trivially) Poisson/Binomial. This may be incorrect for other families
    
    # initialise a list to store covariance matrices
    Sigmas <- list()
    
    # loop through the random effect components
    for (i_re in 1:length(model$modelInfo$reStruc$condReStruc)) {

      # check for a reduced rank covariance structure
      if (names(model$modelInfo$reStruc$condReStruc[[i_re]]$blockCode) != "rr") {
        
        # perform a Kronecker on the between effects covariance structure (Sigmas_theta) and the within effect covariance structure
        Sigmas[[i_re]] <- kronecker(Sigmas_theta[[1]][[i_re]], diag(model$modelInfo$reStruc$condReStruc[[i_re]]$blockReps))
        # TODO: pretty sure this is the direction of the Kronecker... but could try:
        # Sigmas[[i_re]] <- kronecker(diag(model$modelInfo$reStruc$condReStruc[[i_re]]$blockReps), Sigmas_theta[[1]][[i_re]])
      
      } else {
        
        # extract the factor loadings
        Lambda <- model$obj$env$report(model$fit$parfull)$fact_load[[i_re]]
        
        # latent factors are independent, standard normal so use identity penalty matrix
        
        # also need to adjust the random effect design matrix by multiplying by the factor loadings

        # expand the factor loadings by blockReps to align with Z matrix
        IxLambda <- kronecker(diag(model$modelInfo$reStruc$condReStruc[[i_re]]$blockReps), Lambda)
        
        # get the current random effect indices within Z
        n_re <- unlist(lapply(model$modelInfo$reStruc$condReStruc, function(x){ x[["blockReps"]] * x[["blockSize"]] }))
        if (i_re == 1) {
          # for the first random effects start at 1
          re_idx <- 1:n_re[[i_re]]
        } else {
          # for other random effects start at the previous index + 1 and up to the cumulative sum of random effect indices
          re_idx <- (n_re[[i_re - 1]] + 1):cumsum(n_re)[[i_re]]
        }
        # compute the new Z component
        newZ_comp <- mod_str$data.tmb$Z[ , re_idx] %*% IxLambda
        # put the current Z matrix back together again
        X <- as.matrix(cbind(mod_str$data.tmb$X,
                             mod_str$data.tmb$Z[ , 1:ncol(mod_str$data.tmb$Z) < min(re_idx)],
                             newZ_comp,
                             mod_str$data.tmb$Z[ , 1:ncol(mod_str$data.tmb$Z) > max(re_idx)])
        )
        
        # need to adjust the model$modelInfo$reStruc$condReStruc blockReps/Sizes etc. so that subsequent rr() will be adjusted correctly

        
      }
      
      # perform a Kronecker on the between effects covariance structure (Sigmas_theta) and the within effect covariance structure
      Sigmas[[i_re]] <- kronecker(Sigmas_theta[[1]][[i_re]], diag(model$modelInfo$reStruc$condReStruc[[i_re]]$blockReps))
      # TODO: pretty sure this is the direction of the Kronecker... but could try:
      # Sigmas[[i_re]] <- kronecker(diag(model$modelInfo$reStruc$condReStruc[[i_re]]$blockReps), Sigmas_theta[[1]][[i_re]])
    }
    
    # invert the covariance matrices, multiply by dispersion parameter, and create block diagonal form 
    S_random <- Matrix::bdiag(lapply(Sigmas, function(x){disp.par * solve(x)}))
    
    # combine with unpenalised fixed effect components
    S <- Matrix::bdiag(matrix(data = 0, nrow = n_fixed, ncol = n_fixed), S_random)
  }
  
  # calculate the weights matrix for the hat matrix which are, according to mgcv: 1 / ( V(mu) * g'(mu)^2 )
  eta <- predict.glmmTMB(model)
  mu <- model$modelInfo$family$linkinv(eta)
  g_dash_squared <- 1 / model$modelInfo$family$mu.eta(eta)^2 # family$mu.eta() gets the derivative of mu w.r.t. eta which we invert for the required term
  V <- model$modelInfo$family$variance(mu)
  
  # multiply by the observation weights
  if (!is.null(model$frame$`(weights)`)) {
    w <- model$frame$`(weights)` / (V * g_dash_squared)
  } else {
    w <- 1 / (V * g_dash_squared)
  }

  # convert the weights into a diagonal matrix
  W <- Matrix::Diagonal(x = w)

  # perform matrix calculations according to formula in Wood S. N. 2017 Section 4.3
  XtW <- t(X) %*% W
  XtWX <- XtW %*% X
  XtWXplusS <- XtWX + S
  XtWXplusS_inv <- solve(XtWXplusS)
  if (by.obs) {
    A <- (X %*% XtWXplusS_inv) %*% XtW
  } else {
    A <- XtWXplusS_inv %*% XtWX # can re-arrange to this simpler multiplication if we are just interested in the trace
  }
  # TODO: I'm sure there is a more efficient way to do most of this (e.g. using a QR decomp of X?)
  
  return(A[cbind(1:min(dim(A)), 1:min(dim(A)))])
}

#' Predict the latent field components from a model fitted via \code{mvlgcp()}
#'
#' @description This function sets up the appropriate basis function matrix for the \code{newdata} provided, then multiplies by the random coefficient terms (at the mode of the Laplace approximation).
#'
#' @param object a fitted \code{glmmTMB} model object.
#' @param newdata a data frame containing the point locations (i.e. coordinates used in the original model fit) at which to calculate the field(s). If \code{newdata} is a two-column matrix, these are assumed to be the horizontal and vertical coordinates respectively.
#' @param which.response (optional) a single integer (or vector thereof) indicating the which of the response-specific fields to return. Applies only to multivaraite LGCP. Response numbers correspond to the response in order of \code{unique(object$response.id)}.
#' 
#' 
#' @details
#' The function first calculates basis functions (associated with the fitted \code{object}) at the coordinates within \code{newdata}. The resulting basis function matrix is then multiplied by the random coefficients at their mode (as used in the Laplace approximation of the marginal likelihood). If a specific response field is required, the \code{object$n_factor} latent fields are then multiplied by the response-specific factor loadings. 
#' 
#' @return a matrix with \code{nrow(newdata)} rows and k columns - where k is the number of latent fields. k = 1 in the univariate case, k = \code{length(which.response)} when supplied otherwise, k = \code{object$n_factors}.
#' @export
#'
#' @examples
#' # using a subset of the Lansing Wood data:
#' dat <- lansing[lansing$tree == "blackoak", ]
#' 
#' # set up the basis functions with dimension k = 100
#' bfs <- make_basis(k = 100, dat)
#' 
#' # fit univariate LGCP to the blackoak point pattern data
#' m <- mvlgcp(pt ~ 1, data = dat, weights = dat$wt, basis.functions = bfs)
#' 
#' # extract the field
#' fld <- get_field(m)
#' 
#' # get the domain
#' domain <- dat[dat$pt == 0, ]
#' 
#' plot(vec2im(fld, domain$x, domain$y), main = "Black Oak Latent Field")
get_field <- function(object, newdata, which.response = NULL) {
  
  # check that the fitted model has bsais functions attached to it
  if (is.null(object$basis.functions)) {
    stop("Supplied fitted model does not include basis functions")
  }
  # if newdata is missing using the model data
  if (missing(newdata)) {
    warning("'newdata' is not provided. Fitted values will be used...")
    newdata <- eval(as.list(object$call)$data)
  }
  # if more than two columns are supplied, check that these contain the correct coordinates
  if (ncol(newdata) > 2) {
    if (!all(attr(object$basis.functions, "coord.names") %in% colnames(newdata))) {
      stop("'coord.names' associated with the fitted model's basis functions are not found within the newdata supplied")
    } else{
      newdata <- newdata[,attr(object$basis.functions, "coord.names")]
    }
  }
  
  # get the basis function matrix on the new data
  X <- bf_matrix(object$basis.functions, newdata)
  
  # get the basis function coefficients ########################################
  
  # get all of the random effect coefficients from the model
  all_b <- getME.glmmTMB(object, "b")
  
  # get information on the other random effect terms (assumes "basis.functions" is a term reserved for these approximate field models)
  other_re <- object$modelInfo$reStruc$condReStruc[-which(object$modelInfo$grpVar == "basis.functions")]
  
  # determine the number of random effects to  skip (assumes that the field(s) are last in the formula, as is the setup for mvlgcp())
  skip_n <- do.call(sum, lapply(other_re, function(x){as.numeric(x[1]) * as.numeric(x[2])}))

  # get the full set of random coefficients for the latent fields (will include zeros in the reduced rank case)
  if (skip_n > 0) {
    tmp.b_flds <- all_b[((skip_n)+1):length(all_b)]
  } else {
    tmp.b_flds <- all_b
  }
  
  # check if the model is multivariate and adjust the coefficient vector accordingly
  b_flds <- list()
  if (object$n_factors > 0) {
    for (i in 1:object$n_factors) {
      b_flds[[i]] <- all_b[seq(i,length(tmp.b_flds), by = object$n_responses)]
    }
    # while here, extract the factor loadings for the rank reduced fields
    Lambda <- object$obj$env$report(object$fit$parfull)$fact_load[[which(object$modelInfo$grpVar == "basis.functions")]]
  } else {
    b_flds <- tmp.b_flds
    Lambda <- NULL
  }
  b <- sapply(b_flds, cbind)
  
  # compute the field
  fld <- X %*% b
  
  # adjust according to "which.response"
  if (!is.null(which.response)) {
    if (object$n_responses == 1) {
      
      warning("'which.response' is not applicable to a univariate model... returning the single latent field")
      
    } else if (object$n_responses < length(which.response)) {
      
      stop("more fields requested through 'which.response' than responses in the fitted model. Please check and retry")
      
    } else {
      
      # extract the factor loadings for the rank reduced fields
      Lambda <- object$obj$env$report(object$fit$parfull)$fact_load[[which(object$modelInfo$grpVar == "basis.functions")]]
      
      # calculate the response-specific field
      if (length(which.response) > 1) {
        fld <- fld %*% t(Lambda[which.response, ])
      } else {
        fld <- fld %*% Lambda[which.response, ]
      }
      
    }
  }
  
  attr(fld, "coefficients") <- b
  attr(fld, "loadings") <- Lambda
  attr(fld, "response.names") <- object$response.id
  
  return(fld)
}

#' Create a biplot from a multivarite LGCP
#'
#' @description This function creates a biplot of basis function coefficients (scores) forming the \code{n_factors} shared latent fields approximating the multivarate latent fields in a LGCP. Arrows/lines represent response-specific factor loadings.
#'
#' @param x a fitted model object from \code{mvlgcp}.
#' @param ... additional parameters as in \code{stats::biplot}
#' @param alpha a numeric in [0,1] indicating the scaling of arrows/lines representing factor loading lengths. Defaults to 0.5
#' @param load.names (optional) a character vector of length \code{object$n_reponses} of labels to be used on response-specific loadings. Default is \code{NULL} which means the fitted model's \code{response.id} will be used to determine labels.
#' @param score.col (optional) a vector of length \code{object$n_basis} of colours to be used for scores.
#' @param load.col (optional) a vector of length \code{object$n_reponses} of colours to be used on response-specific loadings.
#' @param show.top.n (optional) a numeric in [0,\code{object$n_reponses}] indicating a subset of the response-specific loadings (from largest to smallest).
#' @param vmax.rotation a logical indicating whether to perform a varimax rotation of the loadings matrix prior to plotting.
#' @param rm.bias a logical indicating whether to remove bias from both the scores and factor loadings by scaling them #NEEDS WORK
#' 
#' @return See \code{stats::biplot}
#' @exportS3Method stats::biplot glmmTMB
#'
#' @examples
#' dat <- lansing[lansing$tree %in% c("blackoak", "hickory", "maple"), ]
#' dat$tree <- factor(as.character(dat$tree)) # remove unused tree species
#'
#' # set up the basis functions with dimension k = 100
#' bfs <- make_basis(k = 100, dat)
#'
#' # fit a multivariate LGCP to the point patterns of "blackoak", "hickory" and "maple" trees
#' m <- mvlgcp(pt ~ (1 | tree), data = dat, weights = dat$wt, basis.functions = bfs, response.id = dat$tree)
#' 
#' biplot(m)
biplot.glmmTMB <- function(x, ..., alpha = 0.5, load.names, score.col, load.col, show.top.n, vmax.rotation = FALSE, rm.bias = T) {

  # use get_field() to extract loadings and scores
  suppressWarnings(assign("object", get_field(x)))
  
  fact_loads <- attr(object, "loadings")
  bf.coeffs <- attr(object, "coefficients")
  
  if (vmax.rotation) {
    rot <- varimax(fact_loads)
    Lambda <- rot$loadings
    b <- bf.coeffs %*% rot$rotmat
  } else {
    if (rm.bias) {
      Lambda <- scale(fact_loads)
      b <- scale(bf.coeffs)
    } else {
      Lambda <- fact_loads
      b <- bf.coeffs
    }
  }
  if (ncol(Lambda) > 2) {
    warning("More than two latent factors found. Only the first two will be used")
  }
  
  # order the loadings by arrow lengths
  arrow_lengths <- apply(Lambda, 1, function(x){sqrt(x[1]^2 + x[2]^2)})
  floads <- Lambda[order(arrow_lengths, decreasing = T), ]
  # get arrow names
  if (missing(load.names)) {
    resp.names <- attr(object, "response.names")
    load.names <- if (is(resp.names, "factor")) {levels(resp.names)} else {as.character(unique(resp.names))}
  }
  # subset the factor loading if required
  if (missing(show.top.n)) {
    show.top.n <- nrow(fact_loads)
  }
  sub.floads <- floads[1:show.top.n, ]
  sub.names <- load.names[order(arrow_lengths, decreasing = T)][1:show.top.n]
  # determine colors
  if (missing(score.col)) {
    score.col <- "black"
  } 
  if (missing(load.col)) {
    load.col <- "black"
  } else {
    load.col <- load.col[order(arrow_lengths, decreasing = T)][1:show.top.n]
  }
  
  ### Create the biplot ########################################################
  {
    # par(mar = c(3.1,2.1,0.5,0),mgp=c(1.75,0.75,0))
    plot(b, xlim = range(b[,1]), ylim = range(b[,2]),
         pch = 16, col = score.col, ...)
    text(x = (sub.floads[,1] * alpha) + 0.1 * sign(sub.floads[,1]), y = (sub.floads[,2] * alpha) + 0.1 * sign(sub.floads[,2]), labels = sub.names, cex = 0.75, offset = 0.8, col=load.col)
    arrows(x0 = rep(0, nrow(sub.floads)), y0 = rep(0, nrow(sub.floads)), x1 = sub.floads[,1] * alpha, y1 = sub.floads[,2] * alpha, length = 0.0, angle = 30, lwd = 0.5, col=load.col)
  }
}

#' Function to convert a vector (with corresponding coordinates) to spatstat image.
#'
#' @description This function converts the supplied vector into a matrix of vertical (rows) and horizontal (columns) points according to the supplied x, y locations to create a pixel image. Default is to create a \code{spatstat} image object.
#'
#' @param vec a vector of field values of equal length to locations x.loc and y.loc
#' @param x.loc a vector of horizontal locations in 2D
#' @param y.loc a vector of vertical locations in 2D
#' @param return.matrix a logical indicating whether to return the vector in a matrix of y.loc (rows) and x.loc (columns) instead of the default \code{im} object. Useful if other image plotting is preferred.
#'
#' @return a spatstat image (im) object (or matrix if \code{return.matrix = TRUE})
#' @export
#'
#' @importFrom spatstat.geom im
#'
#' @examples
#' # Get the domain from the gorillas data
#' domain <- gorillas[gorillas$pres == 0, ]
#' elevation <- vec2im(vec = domain$elevation, x.loc = domain$x, y.loc = domain$y)
#' plot(elevation)
vec2im <- function(vec, x.loc, y.loc, return.matrix = FALSE) {
  ux <- sort(unique(x.loc))
  uy <- sort(unique(y.loc))
  nx <- length(ux)
  ny <- length(uy)
  row.ref <- match(x.loc, ux)
  col.ref <- match(y.loc, uy)
  Vec <- rep(NA, max(row.ref)*max(col.ref))
  vec.ref <- (col.ref - 1)*max(row.ref) + row.ref
  Vec[vec.ref] <- vec
  var.mat <- matrix(Vec, max(row.ref), max(col.ref), dimnames = list(ux, uy))
  var.im <- spatstat.geom::im(t(var.mat), xcol = sort(unique(x.loc)), yrow = sort(unique(y.loc)))
  if (return.matrix) {
    ret.obj <- var.mat
  } else {
    ret.obj <- var.im
  }
  return(ret.obj)
}

#' Gorilla Nesting Sites
#'
#' This is the `gorillas` dataset from the package `spatstat` (also featured in `inlabru`), formatted for use with `scampr` models. Describes the nesting locations of gorillas, and associated predictors, in Kagwene Gorilla Sanctuary in Cameroon.
#'
#' @format A data frame with 26020 rows and 11 variables:
#' \describe{
#'   \item{x}{Horizontal coordinate in units Easting.}
#'   \item{y}{Vertical coordinate in units Northing.}
#'   \item{quad.size}{Size of the quadrat in Northing/Easting units squared. Set to 1e-6 at the nesting sites.}
#'   \item{aspect}{Compass direction of the terrain slope. Categorical, with levels N, NE, E, SE, S, SW, W and NW, which are coded as integers 1 to 8.}
#'   \item{elevation}{Digital elevation of terrain, in metres.}
#'   \item{heat}{Heat Load Index at each point on the surface (Beer's aspect), discretised. Categorical with values Warmest (Beer's aspect between 0 and 0.999), Moderate (Beer's aspect between 1 and 1.999), Coolest (Beer's aspect equals 2). These are coded as integers 1, 2 and 3, in that order.}
#'   \item{slopeangle}{Terrain slope, in degrees.}
#'   \item{slopetype}{Type of slope. Categorical, with values Valley, Toe (toe slope), Flat, Midslope, Upper and Ridge. These are coded as integers 1 to 6.}
#'   \item{vegetation}{Vegetation type: a categorical variable with 6 levels coded as integers 1 to 6 (in order of increasing expected habitat suitability).}
#'   \item{waterdist}{Euclidean distance from nearest water body, in metres.}
#'   \item{pres}{Binary indicating whether the row is a presence record (1) or a quadrat (0).}
#' }
#' @source \url{Library `inlabru`}
"gorillas"

#' Lansing Woods Point Pattern (as a data.frame)
#'
#' Locations and botanical classification of trees in Lansing Woods. The data come from an investigation of a 924 ft x 924 ft (19.6 acre) plot in Lansing Woods, Clinton County, Michigan USA by D.J. Gerrard. The data give the locations of 2251 trees and their botanical classification (into hickories, maples, red oaks, white oaks, black oaks and miscellaneous trees). The original plot size (924 x 924 feet) has been rescaled to the unit square. This dataset comes from the \code{spatstat} package, convert to a \code{data.frame} here for use in \code{glmmTMB::mvlgcp()}. Note that a regular grid of quadrature points has been added (for each species), as needed for the multivariate model. See \code{spatstat.data::lansing}. 
#'
#' @usage data(lansing)
#'
#' @format An object of class \code{"data.frame"} with 2251 rows and 3 variables:
#' \describe{
#'   \item{x}{Horizontal coordinate.}
#'   \item{y}{Vertical coordinate.}
#'   \item{tree}{Botanical classification of trees: hickories; maples; red, white and black oaks; and miscellaneous.}
#'   \item{pt}{Binary indicating whether the row is a point event (1) or a quadrature point (0).}
#'   \item{wt}{The weighting for each quadrature point. Set to 1e-10 at the point events.}
#'  }
#'
#' @references
#' Gerrard, D.J. (1969) Competition quotient: a new measure of the competition affecting individual forest trees. Research Bulletin 20, Agricultural Experiment Station, Michigan State University.
#'
#' Turner, R. and Baddeley, A. (2005), SPATSTAT: an R package for analyzing spatial point patterns, \emph{Journal of Statistical Software}, \bold{12}.
#'
"lansing"