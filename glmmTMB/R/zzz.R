.onLoad <- function(lib, pkg) {
    cat("Loading compiled code...\n")
    library.dynam("glmmTMB", pkg, lib)
}
