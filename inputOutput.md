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
	             truncated_poisson  =401,
	             nbinom1            =500, 
	             nbinom2            =501, 
	             truncated_nbinom1  =502,
	             truncated_nbinom2  =503,
	             t                  =600,
	             tweedie            =700 )
	* MATRIX X
	* MATRIX Z
	* MATRIX Xzi
	* MATRIX Zzi
	* MATRIX Xd (for dispersion model)
	* VECTOR yobs
	* VECTOR offset
	* IVECTOR blockReps
	* IVECTOR blockSize
	* IVECTOR blockNumTheta
	* IVECTOR covCode (0=diagonal, 1=unstructured, 2=compound symmetry, 3=ar1)
	* IVECTOR blockRepszi
	* IVECTOR blockSizezi
	* IVECTOR blockNumThetazi
	* IVECTOR covCodezi
	
* initial values for all parameters (size is the important characteristic)
	* VECTOR beta
	* VECTOR betazi
	* VECTOR betad (for dispersion model)
	* VECTOR b
	* VECTOR bzi
	* VECTOR theta (to build variance covariance of the RE of conditional model)
	* VECTOR thetazi (to build variance covariance of the RE of zero-inflation model)

TMB sends back to R (TODO: be more specific with names and types here)
------------------
* theta
* variance and covariance estimates
	* beta = list(est, se)
	* logSigmaSq = list( list(est, se))
	