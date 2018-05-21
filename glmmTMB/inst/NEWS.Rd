\newcommand{\PR}{\Sexpr[results=rd]{tools:::Rd_expr_PR(#1)}}

\name{NEWS}
\title{glmmTMB News}
\encoding{UTF-8}

\section{CHANGES IN VERSION 0.2.1}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item likelihood profiles (via \code{profile}) and likelihood
      profile confidence intervals (via \code{confint(profile(.))})
      can now be computed;
      \code{confint(fitted,method="profile")} and
      \code{confint(fitted,method="uniroot")} (find CIs by using
      a root-finding algorithm on the likelihood profile)
      \item offsets are now allowed in the zero-inflation and dispersion
      formulas as well as in the main (conditional-mean) formula
      (if \code{offset} is specified as a separate argument, it applies
      only to the conditional mean)
      \item zero-truncated generalized Poisson \code{family=truncated_genpois}
      \item zero-truncated Conway-Maxwell-Poisson
      \code{family=truncated_compois}
      \item \code{predict} now allows \code{type} ("link", "response",
      "conditional", "zprob", "zlink")
    }
   }
   \subsection{BUG FIXES}{
     \itemize{
       \item built-in \code{betar()} family for Beta regression fixed
       (and name changed to \code{beta_family()}) (GH #278)
       \item fixed segfault in predict method when response is specified
       as two columns (GH #289)
       \item fixed summary-printing bug when some random effects have
       covariance terms and others don't (GH #291)
       \item fix bugs in binomial residuals and prediction (GH #307)
     }
   }
   \subsection{USER-VISIBLE CHANGES}{
     \itemize{
       \item in \code{predict.glmmTMB},
       the \code{zitype} argument has been rolled into the new
       \code{type} argument: \strong{default prediction type is now
	 "link" instead of "response", in order to match \code{glm()} default}
     }
   }
 }
 
