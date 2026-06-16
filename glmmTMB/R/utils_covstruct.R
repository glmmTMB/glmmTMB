## Workaround to associate numeric values with factor levels in a way
## that survives through the lme4 machinery.

##' Create a factor with numeric interpretable factor levels.
##'
##' Some \code{glmmTMB} covariance structures require extra
##' information, such as temporal or spatial
##' coordinates. \code{numFactor} allows to associate such extra
##' information as part of a factor via the factor levels. The
##' original numeric coordinates are recoverable without loss of
##' precision using the function \code{parseNumLevels}.  Factor levels
##' are sorted coordinate wise from left to right: first coordinate is
##' fastest running.
##'
##' \code{sepgrid} is similar to \code{numFactor}, but creates levels for
##' the complete Cartesian product of the supplied coordinate levels.  Factor
##' inputs preserve unused levels; non-factor inputs use sorted observed
##' values.  Use factors with explicit levels when globally unobserved cells
##' are part of the intended separable grid, for example an unobserved day in
##' an AR(1) time series.
##' @title Factor with numeric interpretable levels.
##' @param x Vector, matrix or data.frame that constitute the
##'     coordinates.
##' @param ... Additional vectors, matrices or data.frames that
##'     constitute the coordinates.
##' @return Factor with specialized coding of levels.
##' @examples
##' ## 1D example
##' numFactor(sample(1:5,20,TRUE))
##' ## 2D example
##' coords <- cbind( sample(1:5,20,TRUE), sample(1:5,20,TRUE) )
##' (f <- numFactor(coords))
##' parseNumLevels(levels(f)) ## Sorted
##' ## Used as part of a model.matrix
##' model.matrix( ~f )
##' ## parseNumLevels( colnames(model.matrix( ~f )) )
##' ## Error: 'Failed to parse numeric levels: (Intercept)'
##' parseNumLevels( colnames(model.matrix( ~ f-1 )) )
##' @export
numFactor <- function(x, ...) {
    y <- data.frame(x, ...)
    if( !all( sapply(y, is.numeric) | sapply(y, is.factor)) )
        stop("All arguments to 'numFactor' must be numeric or factor.")
    asChar <- function(y) {
        y <- lapply(y, as.character)
        ans <- do.call("paste", c(y, list(sep=",")))
        paste0("(", ans, ")")
    }
    fac <- asChar(y)
    ndup <- !duplicated(fac)
    y0 <- y[ndup, , drop=FALSE]
    for (col in seq_along(y0) ) {
        y0 <- y0[ order( y0[[col]] ), , drop=FALSE]
    }
    facLevels <- asChar(y0)
    factor( fac, levels = facLevels )
}

##' @rdname numFactor
##' @export
sepgrid <- function(x, ...) {
    ## `sepgrid()` is deliberately close to `numFactor()`: it returns an
    ## ordinary factor whose levels encode coordinate pairs.  The important
    ## difference is that `sepgrid()` always creates the *complete Cartesian
    ## product* of the coordinate levels, even if some cells are not observed in
    ## the data.  That is essential for AR(1) margins: if day 2 is globally
    ## missing but still part of the study design, days 1 and 3 must be two time
    ## steps apart, not adjacent.
    ##
    ## The current prototype uses only two coordinates:
    ##   sepgrid(member, time)
    ## but the helper is written for a general number of coordinates because the
    ## representation is naturally extensible.  The separable covariance parser
    ## currently rejects anything other than two dimensions.
    y <- data.frame(x, ...)

    ## We accept character/integer variables as a convenience and convert them to
    ## stable integer coordinates below.  Factors preserve their full level set,
    ## including levels that are absent from the observed data.  Non-factors use
    ## sorted observed values, so users should use factors when globally absent
    ## levels are meaningful.
    ok <- vapply(y, function(z) is.numeric(z) || is.factor(z) ||
                   is.character(z) || is.integer(z), logical(1))
    if (!all(ok))
        stop("All arguments to 'sepgrid' must be numeric, factor, integer, or character.")

    ## `parseNumLevels()` expects numeric coordinates in the level labels.  We
    ## therefore store the coordinate *indices* rather than the raw labels:
    ##   member A/B -> 1/2
    ##   time levels 1:30 -> 1:30
    ## This makes downstream C++ code independent of whether member was coded as
    ## "A"/"B", "male"/"female", 0/1, or 1/2.  Human-readable labels can be
    ## added later for polished VarCorr output.
    levs <- lapply(y, function(z) {
        if (is.factor(z)) levels(z) else sort(unique(z[!is.na(z)]))
    })
    vals <- Map(function(z, lev) match(if (is.factor(z)) as.character(z) else z, lev),
                y, levs)
    vals <- as.data.frame(vals)

    asChar <- function(y) {
        is_na <- !stats::complete.cases(y)
        y <- lapply(y, as.character)
        ans <- do.call("paste", c(y, list(sep=",")))
        ans <- paste0("(", ans, ")")
        ans[is_na] <- NA_character_
        ans
    }

    ## `expand.grid()` varies its first argument fastest.  This is exactly the
    ## array convention we use in C++:
    ##   linear index = member_index + n_member * time_index
    ## so the random-effect block is ordered as:
    ##   member1_time1, member2_time1, member1_time2, member2_time2, ...
    grid <- do.call(expand.grid, c(lapply(levs, seq_along),
                                   list(KEEP.OUT.ATTRS = FALSE)))
    factor(asChar(vals), levels = asChar(grid))
}

##' @rdname numFactor
##' @param levels Character vector to parse into numeric values.
##' @importFrom stats complete.cases
##' @export
parseNumLevels <- function(levels) {
    ## Strip initial (irrelevant) characters:
    tmp <- sub("^.*(\\(.+\\))$", "\\1", levels)
    ## Now tmp must have the form ([0-9]*,[0-9]*,...)
    ## Otherwise it's an error
    tmp <- sub("^\\(", "", tmp)
    tmp <- sub("\\)$", "", tmp)
    ## Split string and convert to numeric
    ans <- lapply( strsplit(tmp, ","), as.numeric )
    ans <- t( do.call("cbind", ans) )
    ## if(any(is.na(ans))) stop("Failed to parse numeric levels.")
    if(any(is.na(ans))) {
        stop("Failed to parse numeric levels: ",
             levels[!complete.cases(ans)])
    }
    ans
}

## The helpers below are intentionally small and explicit.  They are a local
## parser layer for the public product syntax
##
##   separable(us(0 + member) %x% ar1(0 + time) | group,
##             scale = us(0 + member))
##
## `reformulas::splitForm()` can understand `separable(grid + 0 | group)` as a
## special random-effect term, but it does not currently return a structured
## object for the product metadata.  To avoid a larger formula-parser change,
## we rewrite the public call above into
##
##   separable(sepgrid(member, time) + 0 | group, list(...metadata...))
##
## This keeps the workaround local to glmmTMB, uses the existing `reTrmAddArgs`
## path, and avoids both encoded strings and formula-level sidecar metadata.
.sep_deparse <- function(x) deparse1(x, collapse = "", width.cutoff = 500L)

.sep_call_name <- function(x) {
    if (!is.call(x)) return(NULL)
    .sep_deparse(x[[1]])
}

.sep_find_call <- function(x, name) {
    ## Recursively search a call object.  We need this because the actual
    ## random-effect expression passed here is `sepgrid(member, time) + 0 | group`;
    ## the `sepgrid()` call is nested inside arithmetic and bar calls.
    if (!is.call(x)) return(NULL)
    if (identical(.sep_call_name(x), name)) return(x)
    for (i in seq_along(x)[-1]) {
        ans <- .sep_find_call(x[[i]], name)
        if (!is.null(ans)) return(ans)
    }
    NULL
}

.sep_find_calls <- function(x, name) {
    ## Return all calls with the requested head.  This is used by glmmTMB() to
    ## identify generated `sepgrid(...)` model-frame columns that must keep their
    ## full Cartesian levels even when ordinary factors are level-dropped.
    if (!is.call(x)) return(list())
    ans <- if (identical(.sep_call_name(x), name)) list(x) else list()
    for (i in seq_along(x)[-1]) {
        ans <- c(ans, .sep_find_calls(x[[i]], name))
    }
    ans
}

.sepgrid_colnames <- function(...) {
    forms <- list(...)
    calls <- unlist(lapply(forms, function(f) {
        if (!inherits(f, "formula")) return(list())
        .sep_find_calls(f[[length(f)]], "sepgrid")
    }), recursive = FALSE)
    unique(vapply(calls, .sep_deparse, character(1)))
}

.sep_is_zero <- function(x) {
    is.numeric(x) && length(x) == 1L && isTRUE(unname(x) == 0)
}

.sep_product_margin_var <- function(x) {
    ## Phase-1 product syntax deliberately accepts only no-intercept one-variable
    ## margins, e.g. `us(0 + role)`.  The returned variable is compiled to the
    ## current complete-grid representation; richer margin formulas will need a
    ## future marginal-design compiler rather than `sepgrid()`.
    if (is.call(x) && identical(.sep_call_name(x), "+") && length(x) == 3L) {
        if (.sep_is_zero(x[[2]]) && is.name(x[[3]])) return(.sep_deparse(x[[3]]))
        if (.sep_is_zero(x[[3]]) && is.name(x[[2]])) return(.sep_deparse(x[[2]]))
    }
    stop("separable() product margins currently require exactly one ",
         "no-intercept variable, e.g. us(0 + role) %x% ar1(0 + day).")
}

.sep_product_margin_spec <- function(x) {
    if (!is.call(x) || length(x) != 2L)
        stop("separable() product margins must look like us(0 + role) ",
             "or ar1(0 + day).")
    structure(.sep_product_margin_var(x[[2]]), names = .sep_deparse(x[[1]]))
}

.sepgrid_bar_call <- function(vars, group) {
    grid_call <- as.call(c(list(as.name("sepgrid")), lapply(vars, as.name)))
    as.call(list(as.name("|"),
                 as.call(list(as.name("+"), grid_call, 0)),
                 group))
}

.sep_spec_df <- function(x, what = "margin") {
    if (is.null(x)) return(NULL)
    if (is.matrix(x) || is.data.frame(x)) {
        ans <- as.data.frame(x, stringsAsFactors = FALSE)
    } else {
        ans <- data.frame(struc = unname(names(x)),
                          var = unname(x),
                          stringsAsFactors = FALSE)
    }
    if (!all(ans$struc %in% names(.sep_margin_registry))) {
        bad <- unique(ans$struc[!ans$struc %in% names(.sep_margin_registry)])
        stop("Unsupported separable() ", what, ": ", paste(bad, collapse = ", "))
    }
    ans
}

.sep_parse_spec <- function(x) {
    if (is.list(x) && !is.null(x$grid) && !is.null(x$margins)) {
        x$margins <- as.matrix(x$margins)
        return(x)
    }
    stop("Internal separable() metadata is missing or malformed.")
}

.sep_margin_registry <- list(
    homcs = list(
        code = "homcs",
        density_kind = "dense_corr",
        can_scale = TRUE,
        n_scale = function(n) 1L,
        n_corr = function(n) 1L
    ),
    us = list(
        code = "us",
        density_kind = "dense_corr",
        can_scale = TRUE,
        n_scale = function(n) as.integer(n),
        n_corr = function(n) as.integer(n * (n - 1L) / 2L)
    ),
    ar1 = list(
        code = "ar1",
        density_kind = "ar1",
        can_scale = FALSE,
        n_scale = function(n) 0L,
        n_corr = function(n) 1L
    )
)

.sep_density_kind_code <- c(dense_corr = 1L, ar1 = 2L)

.sep_dispatch_code <- c(dense_ar1 = 1L)

.sep_scale_mode_code <- c(
    ## Current implemented scale mode.  "margin" means that one separable
    ## margin supplies absolute standard deviations; all other margins are
    ## correlation-only.  Future global/product/cell modes should get new
    ## entries here instead of overloading the selected-margin integer.
    margin = 1L
)

.sep_supported_pairs <- list(
    list(
        dispatch = "dense_ar1",
        density_kinds = c("dense_corr", "ar1"),
        allowed_codes = list(dense_corr = c("homcs", "us"))
    )
)

.sep_pair_dispatch <- function(regs) {
    ## Look up the order-insensitive C++ evaluator for a pair of separable
    ## margins.  This table is intentionally small for the prototype, but it is
    ## the place future combinations should be enabled.  For example, adding
    ## dense_corr x dense_corr or ar1 x ar1 should mean adding a table entry plus
    ## the corresponding C++ dispatch branch, not changing the formula parser.
    ##
    ## The current entry means:
    ##
    ##   one dense-correlation margin, currently restricted to homcs/us
    ##   one AR(1) margin
    ##
    ## Because the lookup sorts density kinds, it accepts both product orders:
    ## `us(0 + member) %x% ar1(0 + time)` and
    ## `ar1(0 + time) %x% us(0 + member)`.
    kinds <- vapply(regs, `[[`, character(1), "density_kind")
    codes <- vapply(regs, `[[`, character(1), "code")
    for (pair in .sep_supported_pairs) {
        if (!identical(unname(sort(kinds)),
                       unname(sort(pair$density_kinds)))) {
            next
        }
        allowed <- pair$allowed_codes
        allowed_ok <- TRUE
        for (kind in names(allowed)) {
            codes_for_kind <- codes[kinds == kind]
            allowed_ok <- allowed_ok &&
                length(codes_for_kind) > 0L &&
                all(codes_for_kind %in% allowed[[kind]])
        }
        if (allowed_ok) return(pair$dispatch)
    }
    NA_character_
}

.sep_margin_label <- function(x) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    paste0(x$struc, "(", x$var, ")", collapse = " x ")
}

.sep_scale_info <- function(margins, regs, scale = NULL) {
    ## Resolve the scale layer independently of the supported covariance-pair
    ## lookup.  This is the small R-side contract that future scale modes should
    ## extend:
    ##
    ##   mode = "margin": one margin supplies absolute SD parameters
    ##   spec = zero-based selected margin index
    ##
    ## The current implementation deliberately infers margin scale only when
    ## exactly one margin is scale-capable.  If no margin can carry scale, the
    ## model needs a future `scale = global()` mode.  If multiple margins can
    ## carry scale, the user must choose explicitly.
    can_scale <- vapply(regs, `[[`, logical(1), "can_scale")
    scale_candidates <- which(can_scale)

    if (is.null(scale)) {
        scale_margin <- scale_candidates
    } else {
        scale_spec <- .sep_spec_df(scale, "scale margin")
        if (nrow(scale_spec) != 1L) {
            stop("separable() scale must be a single margin call such as ",
                 "scale = us(0 + member).")
        }
        scale_margin <- which(margins$struc == scale_spec$struc &
                              margins$var == scale_spec$var)
        if (length(scale_margin) != 1L) {
            stop("separable() scale must match one of the specified margins, ",
                 "for example scale = us(0 + member) when us(0 + member) ",
                 "is a margin.")
        }
        if (!regs[[scale_margin]]$can_scale) {
            stop("separable() scale = ", scale_spec$struc, "(",
                 scale_spec$var, ") selects a correlation-only margin. ",
                 "Use a scale-capable margin such as homcs() or us().")
        }
    }

    if (length(scale_margin) == 0L) {
        stop("separable() margins ", .sep_margin_label(margins),
             " define only a correlation product. This needs an explicit ",
             "overall scale, but scale = global() is not implemented yet.")
    }
    if (length(scale_margin) > 1L) {
        stop("separable() requires exactly one scale-carrying margin. ",
             "Specify it explicitly with scale = us(0 + member) or use one ",
             "scale-capable margin with one correlation-only margin.")
    }

    list(
        mode = "margin",
        mode_code = as.integer(.sep_scale_mode_code[["margin"]]),
        spec = as.integer(scale_margin - 1L),
        margin = scale_margin
    )
}

.sep_stop_unsupported_pair <- function(margins, regs, scale = NULL) {
    ## Unsupported separable combinations often fail because the likelihood
    ## dispatch has not been implemented yet.  However, scale ambiguity is a
    ## separate and important design issue, so diagnose that first when possible.
    can_scale <- vapply(regs, `[[`, logical(1), "can_scale")

    if (sum(can_scale) == 0L) {
        stop("separable() margins ", .sep_margin_label(margins),
             " currently define only a correlation product. A future global ",
             "scale mode is planned, for example scale = global(), but it is ",
             "not implemented yet.")
    }

    if (is.null(scale) && sum(can_scale) > 1L) {
        stop("More than one separable() margin can carry scale in ",
             .sep_margin_label(margins), ". Please specify the scale margin ",
             "explicitly, for example scale = ", margins$struc[which(can_scale)[1]],
             "(", margins$var[which(can_scale)[1]], "). This covariance pair ",
             "is also outside the current homcs/us x ar1 prototype.")
    }

    if (!is.null(scale)) {
        scale <- .sep_spec_df(scale, "scale margin")
        scale_match <- which(margins$struc == scale$struc &
                             margins$var == scale$var)
        if (length(scale_match) == 1L && !regs[[scale_match]]$can_scale) {
            stop("separable() scale = ", scale$struc, "(",
                 scale$var, ") selects a correlation-only margin. ",
                 "Use a scale-capable margin such as homcs() or us(), ",
                 "or wait for a future global scale mode.")
        }
        if (length(scale_match) == 1L) {
            stop("separable() does not yet implement the covariance pair ",
                 .sep_margin_label(margins), ". The explicit scale selector ",
                 "chooses the scale margin, but it does not enable unsupported ",
                 "density combinations.")
        }
    }

    stop("separable() currently supports homcs(0 + member) %x% ar1(0 + time) ",
         "and us(0 + member) %x% ar1(0 + time).")
}

.sep_restruc_info <- function(spec, cnms, blksize) {
    ## Centralize the R-side contract for currently supported separable terms.
    ## This keeps `getReStruc()` from knowing how many parameters each margin
    ## contributes, and gives future margins one obvious place to declare their
    ## parameter counts and density kind.
    spec <- .sep_parse_spec(spec)

    coords <- parseNumLevels(cnms)
    if (ncol(coords) != 2L)
        stop("separable() currently requires a two-dimensional sepgrid().")
    dims <- as.integer(apply(coords, 2, function(z) length(unique(z))))
    if (prod(dims) != blksize)
        stop("separable() requires a complete rectangular sepgrid().")

    margins <- .sep_spec_df(spec$margins)
    if (nrow(margins) != 2L)
        stop("separable() currently requires exactly two margin structures.")

    pair <- margins$struc
    regs <- .sep_margin_registry[pair]
    dispatch <- .sep_pair_dispatch(regs)
    if (is.na(dispatch)) {
        .sep_stop_unsupported_pair(margins, regs, spec$scale)
    }
    if (!identical(margins$var, spec$grid)) {
        stop("The separable() margin variables must match sepgrid() variables. ",
             "Use, for example, ",
             "separable(homcs(0 + member) %x% ar1(0 + time) | group).")
    }

    scale_info <- .sep_scale_info(margins, regs, spec$scale)

    ntheta <- sum(vapply(seq_along(regs), function(i) {
        nscale <- if (i == scale_info$margin) regs[[i]]$n_scale(dims[[i]]) else 0L
        nscale + regs[[i]]$n_corr(dims[[i]])
    }, integer(1)))
    density_kind <- vapply(regs, `[[`, character(1), "density_kind")

    list(
        dims = dims,
        codes = as.integer(vapply(pair, function(z) .valid_covstruct[[z]], numeric(1))),
        density_kinds = as.integer(.sep_density_kind_code[density_kind]),
        dispatch = as.integer(.sep_dispatch_code[dispatch]),
        scale_mode = scale_info$mode_code,
        scale_spec = scale_info$spec,
        ntheta = as.integer(ntheta),
        density_kind = density_kind,
        margins = margins,
        scale = spec$scale
    )
}

.sep_make_product_spec <- function(bar_expr, scale = NULL) {
    ## Compile the public product syntax
    ##
    ##   separable(us(0 + role) %x% ar1(0 + day) | group, scale = us(0 + role))
    ##
    ## into the same complete-grid representation used by the current backend.
    ## Only simple no-intercept one-variable margins are accepted for now; this
    ## keeps the door open for a later marginal-design compiler without implying
    ## that random-slope margins are already supported.
    if (!is.call(bar_expr) || !identical(.sep_call_name(bar_expr), "|") ||
        length(bar_expr) != 3L) {
        stop("separable() product syntax must look like ",
             "separable(us(0 + role) %x% ar1(0 + day) | group, ...).")
    }
    prod_expr <- bar_expr[[2]]
    if (!is.call(prod_expr) || !identical(.sep_call_name(prod_expr), "%x%") ||
        length(prod_expr) != 3L) {
        stop("separable() product syntax requires exactly two margins joined ",
             "by %x%, e.g. us(0 + role) %x% ar1(0 + day).")
    }

    margin_calls <- list(prod_expr[[2]], prod_expr[[3]])
    margins <- c(.sep_product_margin_spec(margin_calls[[1]]),
                 .sep_product_margin_spec(margin_calls[[2]]))
    if (anyDuplicated(unname(margins))) {
        stop("separable() product margins must use distinct variables.")
    }
    if (!all(names(margins) %in% names(.sep_margin_registry))) {
        bad <- unique(names(margins)[!names(margins) %in% names(.sep_margin_registry)])
        stop("Unsupported separable() margin: ", paste(bad, collapse = ", "))
    }

    scale_spec <- NULL
    if (!is.null(scale)) {
        scale_vec <- .sep_product_margin_spec(scale)
        scale_spec <- .sep_spec_df(scale_vec, "scale margin")
    }

    m <- .sep_spec_df(margins)

    list(grid = unname(margins),
         margins = m,
         scale = scale_spec,
         grid_expr = .sepgrid_bar_call(unname(margins), bar_expr[[3]]))
}

.sep_spec_call <- function(spec) {
    ## Convert structured metadata to an unevaluated `list(...)` call.  This
    ## lets `reformulas::splitForm()` carry the metadata in its existing
    ## `reTrmAddArgs` slot, attached to the exact random-effect term.
    chr_vec_call <- function(x) as.call(c(list(as.name("c")), as.list(x)))
    margins <- as.data.frame(spec$margins, stringsAsFactors = FALSE)
    call_args <- list(
        as.name("list"),
        grid = chr_vec_call(spec$grid),
        margins = as.call(list(
            as.name("data.frame"),
            struc = chr_vec_call(margins$struc),
            var = chr_vec_call(margins$var),
            stringsAsFactors = FALSE
        ))
    )
    if (!is.null(spec$scale)) {
        scale <- as.data.frame(spec$scale, stringsAsFactors = FALSE)
        call_args$scale <- as.call(list(
            as.name("data.frame"),
            struc = chr_vec_call(scale$struc),
            var = chr_vec_call(scale$var),
            stringsAsFactors = FALSE
        ))
    }
    as.call(call_args)
}

.rewrite_separable_expr <- function(x) {
    ## Walk the formula call tree and rewrite rich separable calls into the
    ## lower-level form that `reformulas::splitForm()` already understands:
    ##
    ##   separable(grid + 0 | group, list(...metadata...))
    ##
    ## The public syntax handled here is the product form:
    ##
    ##   separable(us(0 + member) %x% ar1(0 + time) | group,
    ##             scale = us(0 + member))
    ##
    ## If a formula already contains the lower-level list form, leave it alone.
    if (!is.call(x)) return(x)
    if (identical(.sep_call_name(x), "separable")) {
        args <- as.list(x[-1])
        nms <- names(args)
        if (is.null(nms)) nms <- rep("", length(args))
        scale_i <- which(nms == "scale")
        if (length(scale_i) > 1L)
            stop("separable() accepts at most one scale argument.")
        scale_arg <- if (length(scale_i)) args[[scale_i]] else NULL
        named_i <- which(nzchar(nms) & nms != "scale")
        if (length(named_i) > 0L) {
            stop("separable() only accepts named argument scale.")
        }

        unnamed_args <- args[!nzchar(nms)]
        if (length(unnamed_args) == 1L) {
            spec <- .sep_make_product_spec(unnamed_args[[1]], scale = scale_arg)
            return(as.call(list(as.name("separable"),
                                spec$grid_expr,
                                .sep_spec_call(spec))))
        }

        ## Leave the internal low-level form `separable(grid, list(...))` alone.
        if (length(unnamed_args) == 2L && is.call(unnamed_args[[2]]) &&
            identical(.sep_call_name(unnamed_args[[2]]), "list") &&
            is.null(scale_arg)) {
            return(x)
        }

        stop("separable() requires separable(",
             "margin1(0 + variable) %x% margin2(0 + variable) | group).")
    }
    for (i in seq_along(x)[-1]) x[[i]] <- .rewrite_separable_expr(x[[i]])
    x
}

rewrite_separable_formula <- function(f) {
    ## Apply the rewrite only to the RHS.  We call this early in `glmmTMB()` for
    ## conditional, zero-inflation, and dispersion formulas so all downstream
    ## machinery sees a normal random-effect special.
    if (!inherits(f, "formula")) return(f)
    f[[length(f)]] <- .rewrite_separable_expr(f[[length(f)]])
    f
}
