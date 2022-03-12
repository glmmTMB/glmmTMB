loadNamespace("TMB")
library(TMB)
loadNamespace("glmmTMB")
library(glmmTMB)
print(packageVersion("Matrix"))
print(packageVersion("TMB"))
print(packageVersion("glmmTMB"))

## TMB internal check: should warn if inconsistent
TMB:::checkMatrixPackageVersion()

## the same by hand
scan(system.file("Matrix-version", package = "TMB"),
     what = character(), quiet = TRUE)
scan(system.file("Matrix-version", package = "glmmTMB"),
     what = character(), quiet = TRUE)
scan(system.file("TMB-version", package = "glmmTMB"),
     what = character(), quiet = TRUE)

oldfit <- readRDS(system.file("example_files","salamander1.rds",package="glmmTMB"))
isNullPointer <- function (x) {
    attributes(x) <- NULL
    identical(x, new("externalptr"))
}

## guts of up2date()
print(isNullPointer(oldfit$obj$env$ADFun$ptr))
obj <- oldfit$obj
oldfit$obj <- with(obj$env, TMB::MakeADFun(data, parameters,
                                           map = map,
                                           random = random,
                                           silent = FALSE, DLL = "glmmTMB"))
print(class(oldfit$obj))
oldfit2 <- readRDS(system.file("example_files","salamander1.rds",package="glmmTMB"))

print(up2date)
oldfit2 <- up2date(oldfit2)
print(class(oldfit2$obj))
