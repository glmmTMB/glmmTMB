R gets from user 
--------------
* formula
* ziformula
* data
* family (link)

TMB gets from R
--------------
TODO: create translation tables for integers below in C++ that gets sent to R

* data
	* linkInt (log        = 0,
	           logit      = 1,
	           probit     = 2,
	           inverse    = 3,
	           cloglog    = 4,
	           identity   = 5 )
	* familyInt (gaussian           =0,
	             binomial           =100,
	             betabinomial       =101,
	             beta               =200,
	             gamma              =300,
	             poisson            =400,
	             truncatedpoisson   =401,
	             nbinom1            =500, 
	             nbinom2            =501, 
	             truncatednbinom    =502,
	             t                  =600,
	             tweedie            =700 )
	* covInt (0=diag, 1=unstructured, 2=compound symmetry)
	* MATRIX X
	* MATRIX Z
	* MATRIX Xzi
	* MATRIX Zzi
	* MATRIX Xd (for dispersion model)
	* VECTOR yobs
	* IVECTOR blockreps
	* IVECTOR blocksize
	* IVECTOR blocknumtheta
	* IVECTOR covcode
	*
* initial values for all parameters (size is the important characteristic)
	* VECTOR beta
	* VECTOR betazi
	* VECTOR betad (for dispersion model)
	* VECTOR b
	* VECTOR bzi
	* VECTOR theta (to build variance covariance 
	* VECTOR 

TMB sends back to R (TODO: be more specific with names and types here)
------------------
* theta
* variance and covariance estimates
	* beta = list(est, se)
	* logSigmaSq = list( list(est, se))
	