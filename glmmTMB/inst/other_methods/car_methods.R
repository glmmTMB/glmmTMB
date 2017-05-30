## EXPERIMENTAL (not working, not yet exported)
## modified from car::Anova.default
Anova.glmmTMB <- function (mod, type = c("II", "III", 2, 3),
                           test.statistic = c("Chisq", "F"),
                           vcov. = vcov(mod)[["cond"]], singular.ok, ...) 
{
    stop("not finished yet")
    if (is.function(vcov.)) 
        vcov. <- vcov.(mod)
    type <- as.character(type)
    type <- match.arg(type)
    test.statistic <- match.arg(test.statistic)
    if (missing(singular.ok)) 
        singular.ok <- type == "2" || type == "II"
    switch(type, II = Anova.II.default(mod, vcov., test.statistic, 
        singular.ok = singular.ok), III = Anova.III.default(mod, 
        vcov., test.statistic, singular.ok = singular.ok), `2` = Anova.II.default(mod, 
        vcov., test.statistic, singular.ok = singular.ok), `3` = Anova.III.default(mod, 
        vcov., test.statistic, singular.ok = singular.ok))
}
