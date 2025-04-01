old_tb <- "~/Downloads/glmmTMB_1.1.10.tar.gz"
new_tb <- "glmmTMB_1.1.11.tar.gz"
ord_var <- "diff"

hack_date <- function(x) {
    gsub("([A-Z][a-z]{2}) +([0-9]{1,2})","\\1_\\2", x)
}
get_tb_info <- function(tb) {
    (untar(tb, list = TRUE, extras = "-v")
        |> hack_date()
        |> strsplit(" +")
        |> do.call(what=rbind)
        |> as.data.frame()
        |> setNames(c("perms", "owner", "size", "date", "time", "path"))
        |> subset(select = c(size, path))
        |> transform(size = as.numeric(size))
    )
}

mm <- merge(get_tb_info(old_tb), get_tb_info(new_tb), by = "path", all = TRUE,
            suffixes = c(".old", ".new"))
## include only changed files
mm <- (mm
    |> subset(is.na(size.old) | (!is.na(size.new) & size.new != size.old))
    |> transform(size.old = ifelse(is.na(size.old), 0, size.old))
    |> transform(diff = size.new - size.old)
)

mm <- mm[order(mm[[ord_var]], decreasing = TRUE),]

View(mm)

if (FALSE) {
    
    ## testing hack_date()
    xx <- c("-rw-r--r--  0 molliebrooks staff    4100 Mar  5 10:40 glmmTMB/DESCRIPTION",                                       
            "-rw-r--r--  0 molliebrooks staff    5657 Jan 27 17:39 glmmTMB/NAMESPACE",                                         
            "drwxr-xr-x  0 molliebrooks staff       0 Mar  5 10:36 glmmTMB/R/",                                  
            "-rw-r--r--  0 molliebrooks staff   10198 Jan 27 17:39 glmmTMB/R/Anova.R")
    xx |> hack_date()  |> strsplit(" +") |> do.call(what=rbind)

}
