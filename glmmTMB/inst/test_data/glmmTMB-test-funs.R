## remove leading/trailing white space, empty lines, multiple spaces
squash_white <- function(x) {
    x <- gsub(" +"," ",trimws(x))
    return(x[nchar(x)>0])
}
