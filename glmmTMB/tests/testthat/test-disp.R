stopifnot(require("testthat"),
          require("glmmTMB"))

context("Testing dispersion")

sim1=function(nfac=20, nt=1000, facsd=.1, tsd=.15, mu=0, residsd=1)
{
	dat=expand.grid(fac=factor(letters[1:nfac]), t= 1:nt)
	n=nrow(dat)
	dat$REfac=rnorm(nfac, sd= facsd)[dat$fac]
	dat$REt=rnorm(nt, sd= tsd)[dat$t]

	dat$x=rnorm(n, mean=mu, sd=residsd) + dat$REfac + dat$REt
	return(dat)
}

set.seed(101)

d1=sim1(mu=1, residsd =10)
d2=sim1(mu=10, residsd =5)

d1=transform(d1, 
	fac=paste0(fac, 1),
	disp="ten")
	
d2=transform(d2, 
	fac=paste0(fac, 2),
	disp="five")

dat=rbind(d1, d2)

m0=glmmTMB(x~disp+(1|t)+(1|fac), dispformula=~disp, dat)
expect_equal(unname(fixef(m0)$disp), c(log(10^2), log(5^2)-log(10^2)), tol=1e-2)
