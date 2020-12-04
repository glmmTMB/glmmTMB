## machinery for generating binaries on r-hub and putting them in the right place
## run this in the head directory

## install rhub package
while (!require("rhub")) {
    remotes::install_github("r-hub/rhub")
}

## run check of platform for Windows and MacOS release versions
check(platform=c("windows-x86_64-release","macos-highsierra-release"),
      email="bbolker@gmail.com")

## extract appropriate hashes/tmpdirs from email you get back
mhash <- "3a60fd03bbd74c7f9c48e1054dc018f4"
whash <- "4cc959649991463f94a340380bf10b86"
pkg_version <- "1.0.2.9000"


## helper fun for constructing URLs
ufun <- function(pkg_ver=pkg_version,hash,ext,base_dir=FALSE) {
    url0 <- sprintf("https://artifacts.r-hub.io/glmmTMB_%s.tar.gz-%s",pkg_ver,hash)
    if (base_dir) return(url0)
    return(sprintf("%s/glmmTMB_%s.%s",url0,pkg_ver,ext))
}


src <- sprintf("glmmTMB_%s.tar.gz",pkg_version)
macbin <- sprintf("glmmTMB_%s.tgz",pkg_version)
winbin <- sprintf("glmmTMB_%s.zip",pkg_version)
download.file(ufun(hash=mhash,ext="tgz"), dest=macbin)
download.file(ufun(hash=whash,ext="zip"), dest=winbin)

## put stuff in the right place
library(drat)
insertPackage(macbin,"repos")
insertPackage(winbin,"repos")
insertPackage(src,"repos")

## to install:
## install.packages("glmmTMB", repos="https://github.com/glmmTMB/glmmTMB/tree/master/repos",
##                 type="binary")
unlink(macbin)
unlink(winbin)
unlink(src)
