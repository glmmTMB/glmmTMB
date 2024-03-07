library(glmmTMB)
## devtools::load_all("~/R/pkgs/glmmTMB/glmmTMB")
## dd <- read.csv2("df_glmmTMB.csv") |>
##     subset(select = c(preMDS, F_Absetzen, Betrieb)) |>
##     transform(F_Absetzen = factor(F_Absetzen),
##               Betrieb = factor(Betrieb)) |>
##     na.omit()
## attr(dd, "na.action") <- NULL


dd <- structure(list(preMDS = c(6L, 2L, 1L, 2L, 3L, 34L, 3L, 239L, 
1L, 2L, 4L, 81L, 1L, 1L, 1L, 255L, 8L, 72L, 110L, 3L, 6L, 61L, 
253L, 113L, 49L, 124L, 72L, 4L, 35L, 4206L, 3660L, 3100L, 4308L, 
5871L, 1362L, 4301L, 2673L, 204L, 216L), F_Absetzen = structure(c(1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L), levels = c("0", "1"), class = "factor"), 
    Betrieb = structure(c(2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
    2L, 2L, 2L, 3L, 3L, 5L, 8L, 8L, 8L, 8L, 8L, 8L, 8L, 8L, 8L, 
    8L, 8L, 8L, 8L, 8L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 9L
    ), levels = c("B02", "B03", "B04", "B05", "B06", "B07", "B08", 
    "B10", "B11", "B13", "B13Zucht", "B14"), class = "factor")), row.names = c(58L, 
59L, 60L, 61L, 62L, 63L, 64L, 65L, 66L, 67L, 68L, 69L, 83L, 84L, 
111L, 148L, 149L, 150L, 151L, 152L, 153L, 154L, 155L, 156L, 157L, 
158L, 159L, 160L, 161L, 174L, 175L, 176L, 177L, 178L, 179L, 180L, 
181L, 182L, 183L), class = "data.frame")

library(glmmTMB)
m1 <- glmmTMB(preMDS ~ 1 + F_Absetzen + (1 | Betrieb), data = dd,
        family = truncated_nbinom1)
diagnose(m1)

