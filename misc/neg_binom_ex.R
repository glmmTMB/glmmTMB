## from https://github.com/glmmTMB/glmmTMB/issues/988

library(glmmTMB)
set.seed(101); dd <- data.frame(x = 1:10)
m1 <- glmmTMB(cbind(x, 10-x) ~ 1, family = betabinomial, data = dd)
sessionInfo()

dd2 <- rbind(dd, data.frame(x=-1))
update(m1, data = dd2)  ## NOT detected
update(m1, data = dd2, family = binomial)  ## gives weird/infinite answers

try(glm(cbind(x, 10-x) ~ 1, family = binomial, data = dd2))



our_binomInitialize <- function(family) {
    newtest <- substitute(
        ## added test for glmmTMB
        if (any(y<0)) {
            stop(sprintf('negative values not allowed in %s responses', FAMILY))
        }
      , list(FAMILY=family))
    b0 <- binomial()$initialize
    b0[[3]] <- newtest
    return(b0)
}

b_new <- betabinomial()
b_new$initialize <- our_binomInitialize("betabinomial")

try(update(m1, data = dd2, family = b_new))
