
## run model & extract elapsed time
form0 <- use ~ urban+age+livch+(urban|district)
tt <- function(fun,data,form=form0,family="binomial",debug=FALSE) {
    argList <- list(form,data)
    if (!is.null(family)) argList <- c(argList,list(family=family))
    unname(system.time(do.call(fun,argList))["elapsed"])
}

## replicate timings across models for a specified expansion
## of the data
getTimes <- function(n=1,which=c("glmmTMB","glmer"),
                     basedata=Contraception,form=form0,
                     family="binomial",
                     debug=FALSE) {
  if (n>1) {
    if (!n==round(n)) stop("only integer magnification allowed")
    ## replicate data 
    c2 <- do.call(rbind,replicate(n,basedata,simplify=FALSE))
  } else {
    ## subsample
    c2 <- basedata[sample.int(round(n*nrow(basedata))),]
  }
  ## run for all models
  res <- setNames(vapply(which,tt,data=c2,form=form,
                         family=family,numeric(1)),which)
  return(res)
}
