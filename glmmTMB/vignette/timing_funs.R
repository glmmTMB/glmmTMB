
## run model & extract elapsed time
form0 <- use ~ urban+age+livch+(urban|district)
tt <- function(fun,data,form=form0) {
  if (is.character(fun)) fun <- get(fun)
  unname(system.time(fun(form,data,"binomial"))["elapsed"])
}

## replicate timings across models for a specified expansion
## of the data
getTimes <- function(n=1,which=c("glmmTMB","glmer"),
                     basedata=Contraception) {
  if (n>1) {
    if (!n==round(n)) stop("only integer magnification allowed")
    ## replicate data 
    c2 <- do.call(rbind,replicate(n,basedata,simplify=FALSE))
  } else {
    ## subsample
    c2 <- basedata[sample.int(round(n*nrow(basedata))),]
  }
  ## run for all models
  res <- setNames(vapply(which,tt,data=c2,numeric(1)),which)
  return(res)
}
