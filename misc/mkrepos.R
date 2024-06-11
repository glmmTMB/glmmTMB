## machinery for generating binaries on r-hub and putting them in the right place
## run this in the head directory

## load rhub package, installing if necessary
while (!require("rhub")) {
    remotes::install_github("r-hub/rhub")
}
library(drat)
library(rhub) ## redundant

## run check/build binaries for Windows and MacOS with the current R release
## use rhub::platforms() to see which platforms are available
releases <- c("windows-x86_64-release") ##  "macos-highsierra-release")
## maybe macos-highsierra-release-cran? shouldn't make much difference
##   for building binaries ... ?
email <- "bbolker@gmail.com"  ## set this to your e-mail!
oldrel <- "windows-x86_64-oldrel"

## assuming working directory = glmmTMB repo head
## use here::here() ?
check("glmmTMB", platform = releases, email=email)
check("glmmTMB", platform = oldrel, email=email)

## extract appropriate hashes/tmpdirs from the email you get back
## mhash <- "3a60fd03bbd74c7f9c48e1054dc018f4"
## whash <- "4cc959649991463f94a340380bf10b86"
## pkg_version <- "1.0.2.9000"

## mhash <- "09f2a6a3d51847e19eabcf005e99f715"
## whash <- "5b768de377c94f12a264deff0924553f"
## pkg_version <- "1.1.2"

## whash <- "8676a5b981a94267acf737ae1bd6297f"
## mhash <- "e7461c3b0baa44f0972c12fb7463a7e4"
## pkg_version <- "1.1.3"

pkg_version <- "1.1.8"


## helper fun for constructing URLs
ufun <- function(pkg_ver=pkg_version,hash,ext,base_dir=FALSE) {
    url0 <- sprintf("https://artifacts.r-hub.io/glmmTMB_%s.tar.gz-%s",pkg_ver,hash)
    if (base_dir) return(url0)
    return(sprintf("%s/glmmTMB_%s.%s", url0, pkg_ver, ext))
}

src <- sprintf("glmmTMB_%s.tar.gz", pkg_version)
macbin <- sprintf("glmmTMB_%s.tgz", pkg_version)
winbin <- sprintf("glmmTMB_%s.zip", pkg_version)
download.file(ufun(hash=mhash, ext="tgz"), dest = macbin)
download.file(ufun(hash=whash, ext="zip"), dest = winbin)

## system("git checkout master")
## system('R CMD build --compact-vignettes="both" glmmTMB')
## system("git checkout gh-pages")

## put stuff in the right place
drat::insertPackage(macbin, "docs/repos")
drat::insertPackage(winbin, "docs/repos")
drat::insertPackage(src, "docs/repos")

## to test installation (if you are on MacOS or Windows)
if (FALSE) {
    install.packages("glmmTMB", repos="https://glmmTMB.github.io/glmmTMB/repos",
                 type = "binary")
}

## test installation from source:
if (FALSE) {
    install.packages("glmmTMB", repos = "https://glmmTMB.github.io/glmmTMB/repos",
                     type = "source")
}

## clean up
unlink(macbin)
unlink(winbin)
unlink(src)

setwd("docs/repos")
rmarkdown::render("index.Rmd")
