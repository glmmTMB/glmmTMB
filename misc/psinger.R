library(glmmTMB)

if(FALSE) {
    ## Full psinger example:
    data <- read.table("r_sample.csv",header=TRUE,sep=",")
    system.time( mnb1 <- glmmTMB(y ~ 1 + a + (1|id), data=data, family="nbinom2", verbose=TRUE) )
    ## =================== Result:
    ## outer mgc:  NaN 
    ## Error in nlminb(start = par, objective = fn, gradient = gr) : 
    ## NA/NaN gradient evaluation
    ## In addition: Warning message:
    ## In nlminb(start = par, objective = fn, gradient = gr) :
    ## NA/NaN function evaluation
    ## ===========================
}

## Reproduce same error on 'small' subset:
data <- read.table("r_sample.csv",header=TRUE,sep=",")
set.seed(1); data <- data[sample(1:nrow(data), 10000), ]
mnb1 = glmmTMB(y ~ 1 + a + (1|id), data=data, family="nbinom2", verbose=TRUE)

## Produce different kind of error:
data <- read.table("r_sample.csv",header=TRUE,sep=",")
set.seed(1); data <- data[sample(1:nrow(data), 7000), ]
mnb1 = glmmTMB(y ~ 1 + a + (1|id), data=data, family="nbinom2", verbose=TRUE)

## Error with binomial
data <- read.table("r_sample.csv",header=TRUE,sep=",")
data$y <- (data$y>0)
system.time( mb1 <- glmmTMB(y ~ 1 + a + (1|id), data=data, family="binomial", verbose=TRUE) )
