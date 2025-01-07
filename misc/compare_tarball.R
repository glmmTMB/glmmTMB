old_tb <- "~/Downloads/glmmTMB_1.1.10.tar.gz"
new_tb <- "glmmTMB_1.1.11.tar.gz"

get_tb_info <- function(tb) {
    (untar(tb, list = TRUE, extras = "-v")
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
mm <- mm[order(mm$size.new, decreasing = TRUE),]
mm <- mm[is.na(mm$size.old) | (!is.na(mm$size.new) & mm$size.new != mm$size.old), ]
