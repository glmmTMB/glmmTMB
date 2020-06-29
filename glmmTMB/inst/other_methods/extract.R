
## methods for texreg
##' @param model a glmmTMB model
##' @export
##' @importFrom stats coef
##' @rawNamespace if(requireNamespace("texreg")) importFrom(texreg,extract)
##' 
extract.glmmTMB <- function(model) {
    cc <- coef(summary(model))
    cc <- cc[!vapply(cc,is.null,logical(1))]
    nn <- names(cc)
    pref <- ifelse(nn=="cond","",paste0(nn,"_"))
    names <- unlist(Map(function(x,y) paste0(x,rownames(y)),
                        pref,cc))
    co <- unlist(lapply(cc,function(x) x[,"Estimate"]))
    se <- unlist(lapply(cc,function(x) x[,"Std. Error"]))
    pval <- unlist(lapply(cc,function(x) x[,"Pr(>|z|)"]))
    ## alternately could do something like
    ##  rownames(unlist(cc)) -> gsub("\\.","_",.)
    ##  -> gsub("cond_","")
    if (!requireNamespace("texreg")) {
        stop("must have texreg package installed")
    }
    tr <- texreg::createTexreg(
        coef.names = names,
        coef = co,
        se = se,
        pvalues = pval
        ## use AIC/BIC/logLik?
        ## gof.names = gof.names,
        ## gof = gof
    )
    return(tr)
}

## setMethod("extract", signature = className("glmmTMB", "glmmTMB"),
##          definition = extract.glmmTMB)

#' @param x an expression to evaluate
#' @param stop_ex regular expression for truncating output
#' @param stop_which which occurrence of \code{stop_ex} to find
#' @importFrom utils capture.output
verbatim <- function(x, stop_ex=NULL, stop_which=1)
{
    txt <- trimws(capture.output(x))
    if (!is.null(stop_ex)) {
        txt <- txt[seq(grep(stop_ex,txt)[stop_which]-1)]
    }
    txt <- c("\\begin{verbatim}", txt, "\\end{verbatim}")
    cat(txt, sep="\n")
}

setOldClass("glmmTMB")
setMethod("extract","glmmTMB",extract.glmmTMB)


## m1 <- glmmTMB(formula = count ~ mined + (1 | site), data = Salamanders, 
##               family = poisson, ziformula = ~mined, dispformula = ~1)

## verbatim(summary(m1), "Random effects:")


