## machinery for generating binaries on r-hub and putting them in the right place
## run this in the head directory

## load rhub package, installing if necessary
while (!require("rhub")) {
    remotes::install_github("r-hub/rhub")
}
library(drat)

## run check/build binaries for Windows and MacOS with the current R release
## use rhub::platforms() to see which platforms are available
releases <- c("windows-x86_64-release","macos-highsierra-release")
email <- "bbolker@gmail.com"  ## set this to your e-mail!

check(platform=releases,email=email)

## extract appropriate hashes/tmpdirs from the email you get back
## mhash <- "3a60fd03bbd74c7f9c48e1054dc018f4"
## whash <- "4cc959649991463f94a340380bf10b86"
## pkg_version <- "1.0.2.9000"


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

## system("git checkout master")
## system('R CMD build --compact-vignettes="both" glmmTMB')
system("git checkout gh-pages")

## put stuff in the right place
insertPackage(macbin,"repos")
insertPackage(winbin,"repos")
insertPackage(src,"repos")

## to test installation (if you are on MacOS or Windows)
if (FALSE) {
    install.packages("glmmTMB", repos="https://github.com/glmmTMB/glmmTMB/tree/master/repos",
                 type="binary")
}

## test installation from source:
if (FALSE) {
    install.packages("glmmTMB", repos="https://glmmTMB.github.io/glmmTMB/repos",
                     type="source")
}

## clean up
unlink(macbin)
unlink(winbin)
unlink(src)
