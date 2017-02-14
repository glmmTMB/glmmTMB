library(glmmTMB)
sim1=function(nfac=50, nt=100, facsd=.1, tsd=.15, mu=0, residsd=1)
{
  dat=expand.grid(fac=1:nfac, t= 1:nt)
  n=nrow(dat)
  dat$REfac=rnorm(nfac, sd= facsd)[dat$fac]
  dat$REt=rnorm(nt, sd= tsd)[dat$t]
  dat$x=rnorm(n, mean=mu, sd=residsd) + dat$REfac + dat$REt
  return(dat)
}
set.seed(101)
d1 = sim1(mu=100, residsd =10)
d2 = sim1(mu=200, residsd =5)
d1 = transform(d1,
	fac=paste0(fac, 1),
	disp="ten")
d2 = transform(d2,
	fac=paste0(fac, 2),
	disp="five")
dat = rbind(d1, d2)
m0 = glmmTMB(x~disp+(1|t)+(1|fac), dispformula=~disp, dat)

#increasing nfac to 100 allows the model to converge
set.seed(101)
d3 = sim1(nfac=100, mu=100, residsd =10)
d4 = sim1(nfac=100, mu=200, residsd =5)
d3 = transform(d3,
	fac=paste0(fac, 1),
	disp="ten")
d4 = transform(d4,
	fac=paste0(fac, 2),
	disp="five")
dat34 = rbind(d3, d4)
m34 = glmmTMB(x~disp+(1|t)+(1|fac), dispformula=~disp, dat34)
