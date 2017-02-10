stopifnot(require("testthat"),
          require("glmmTMB"))

context("Testing dispersion")

sim1=function(nfac=40, nt=100, facsd=.1, tsd=.15, mu=0, residsd=1)
{
	dat=expand.grid(fac=factor(letters[1:nfac]), t= 1:nt)
	n=nrow(dat)
	dat$REfac=rnorm(nfac, sd= facsd)[dat$fac]
	dat$REt=rnorm(nt, sd= tsd)[dat$t]

	dat$x=rnorm(n, mean=mu, sd=residsd) + dat$REfac + dat$REt
	return(dat)
}

set.seed(101)

d1=sim1(mu=100, residsd =10)
d2=sim1(mu=200, residsd =5)

d1=transform(d1, 
	fac=paste0(fac, 1),
	disp="ten")
	
d2=transform(d2, 
	fac=paste0(fac, 2),
	disp="five")

## global assignment for testthat
dat <<- rbind(d1, d2)
m0 <<- glmmTMB(x~disp+(1|fac), dispformula=~disp, dat)

test_that("disp calc", {
    expect_equal(unname(fixef(m0)$disp), c(log(10^2), log(5^2)-log(10^2)), tol=1e-2)
})

dat2 <<- rbind(head(d1, 50), head(d2, 50)) #smaller for faster fitting when not checking estimates
nbm0 <<- glmmTMB(round(x)~disp+(1|fac), ziformula=~0, dispformula=~disp, dat2, family=nbinom1, se=FALSE)
pm0 <<- update(nbm0, family=poisson)
nbm1 <<- update(pm0, family=nbinom1)
test_that("update maintains dispformula in call", {
	expect_equal(getCall(nbm0), getCall(nbm1))
})