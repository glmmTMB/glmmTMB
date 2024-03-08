# glmmTMB

<!-- badges: start -->
[![cran version](http://www.r-pkg.org/badges/version/glmmTMB)](https://cran.r-project.org/package=glmmTMB)
[![downloads](http://cranlogs.r-pkg.org/badges/glmmTMB)](http://cranlogs.r-pkg.org/badges/glmmTMB)
[![total downloads](http://cranlogs.r-pkg.org/badges/grand-total/glmmTMB)](http://cranlogs.r-pkg.org/badges/grand-total/glmmTMB)
[![R-CMD-check](https://github.com/glmmTMB/glmmTMB/workflows/R-CMD-check/badge.svg)](https://github.com/glmmTMB/glmmTMB/actions)
<!-- badges: end -->

**Note** Most user-level information has migrated to the [GitHub pages site](https://glmmTMB.github.io/glmmTMB/); please check there.

`glmmTMB` is an R package for fitting generalized linear mixed models (GLMMs) and extensions, built on [Template Model Builder](https://github.com/kaskr/adcomp), which is in turn built on [CppAD](https://www.coin-or.org/CppAD/) and [Eigen](eigen.tuxfamily.org/). It handles a wide range of statistical distributions (Gaussian, Poisson, binomial, negative binomial, Beta ...) as well as model extensions such as zero-inflation, heteroscedasticity, and autocorrelation. Fixed and random effects models can be specified for the conditional and zero-inflated components of the model, as well as fixed effects models for the dispersion parameter.

## Where to ask questions

- for questions about mixed models and data analysis:
    - the [r-sig-mixed-models mailing list](https://stat.ethz.ch/mailman/listinfo/r-sig-mixed-models) (please [subscribe to the list](https://stat.ethz.ch/mailman/listinfo/r-sig-mixed-models) before posting) or 
	- [Cross Validated](https://stats.stackexchange.com)  
The mailing list has more people who know about using mixed models in R; CrossValidated has a larger overall audience and a nicer question/answer format (allows pictures, code formatting, etc.).
- for questions about `glmmTMB` usage:
    - `r-sig-mixed-models@r-project.org` or 
	- [Stack Overflow](https://stackoverflow.com)  
	(similar pros/cons as the previous point)
- for bug, infelicity, and wishlist reporting: the [glmmTMB issues list](https://github.com/glmmTMB/glmmTMB/issues) 
- for TMB-specific questions: the [TMB users forum](https://groups.google.com/forum/#!forum/tmb-users) 
- maintainer e-mail only for urgent/private communications

Please do **not** cross-post, i.e. ask the same question in more than one forum, unless it's suggested that you have posted in the wrong place, or unless you receive total silence in one forum. In the latter case it may be better to send a reminder/"bump" message to the original forum; in any case you should mention in your new message where/when you've previously asked the question.

