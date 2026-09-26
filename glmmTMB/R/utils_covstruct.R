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

## Rewrite kron(A(0 + f1) %x% B(0 + f2) | g) as the ordinary covstruct term
## kron(0 + f2:f1 | g, c("A", "B")): as in kronecker(), the last margin
## varies fastest
kron_sugar <- function(x) {
    if (identical(x[[1]], quote(kron)) && length(x) == 2) {
        flat <- function(e) {
            if (!is.call(e) || !identical(e[[1]], quote(`%x%`))) return(list(e))
            c(flat(e[[2]]), flat(e[[3]]))
        }
        margin <- function(e) {
            tt <- if (is.call(e) && length(e) == 2) terms(eval(call("~", e[[2]])))
            lab <- attr(tt, "term.labels")
            if (length(lab) != 1 || attr(tt, "intercept"))
                stop("kron() margins must look like us(0 + f)", call. = FALSE)
            c(struc = deparse(e[[1]]), f = lab)
        }
        m <- vapply(flat(x[[2]][[2]]), margin, c(struc = "", f = ""))
        lhs <- str2lang(paste("0 +", paste(rev(m["f", ]), collapse = ":")))
        return(call("kron", call("|", lhs, x[[2]][[3]]), m["struc", ]))
    }
    for (i in seq_along(x)[-1]) if (is.call(x[[i]])) x[[i]] <- kron_sugar(x[[i]])
    x
}

