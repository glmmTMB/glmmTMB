library(mgcv); library(glmmTMB); library(reformulas)
devtools::load_all()
test1 <- function(x,z,sx=0.3,sz=0.4) { 
       x <- x*20
       (pi**sx*sz)*(1.2*exp(-(x-0.2)^2/sx^2-(z-0.3)^2/sz^2)+
       0.8*exp(-(x-0.7)^2/sx^2-(z-0.8)^2/sz^2))
}
n <- 500
set.seed(101)
x <- runif(n)/20;z <- runif(n)
y <- test1(x,z) + rnorm(n)*0.2
dd <- data.frame(x, z, y)
b1 <- gam(y~s(x,z), method = "ML", data = dd)
b1B <- glmmTMB(y ~ s(x,z), data = dd)
stopifnot(all.equal(unname(c(predict(b1))), predict(b1B)))
b2 <- gam(y~t2(x,z), data = dd, method = "ML")
b2B <- glmmTMB(y ~ t2(x,z), data = dd) ## error


f <- y ~ s(x,z)
noSpecials(sub_specials(f), delete=FALSE, , specials = c(names(.valid_covstruct), mgcv_specials))

f <- y ~ s(x,z)

f <- y ~ te(x,z)
f <- sub_specials(f, specials = c("|", "||", mgcv_specials))
f <- noSpecials(f, delete=FALSE, , specials = c(names(.valid_covstruct), mgcv_specials))


## fit is not the same ... bad starting values etc. ?
plot(predict(b2), predict(b2B))
all.equal(unname(c(predict(b2))), predict(b2B))

