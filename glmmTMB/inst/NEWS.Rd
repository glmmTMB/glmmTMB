\newcommand{\PR}{\Sexpr[results=rd]{tools:::Rd_expr_PR(#1)}}

\name{NEWS}
\title{glmmTMB News}
\encoding{UTF-8}
\section{CHANGES IN VERSION 1.0.2.1}{
	\itemize{
      \item fix documentation links for CRAN checks
	}
}
\section{CHANGES IN VERSION 1.0.2 (unreleased)}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item the \code{refit()} function is now re-exported (i.e., you no
      longer need to load \code{lme4} to use it)
      \item a \code{modelparm.glmmTMB} method is now provided (so that
      \code{multcomp::glht} should work out of the box with
      \code{glmmTMB} objects)
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item new \code{sparseX} argument to specify sparse fixed-effect
      model matrices for one or more components
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{summary} and model printing now work if
      \code{control=glmmTMBControl(optimizer=optim)} is used (GH #589)
      \item structured covariance models now work in zero-inflation and
      dispersion components (GH #579)
      \item documentaion of formula for variance in beta family (GH #595)
    }
  } % bug fixes
}
\section{CHANGES IN VERSION 1.0.1}{
	\itemize{
      \item updated for R-devel changes (R 4.0.0 will set stringsAsFactors=FALSE by default)
	}
}
\section{CHANGES IN VERSION 1.0.0}{

The 1.0.0 release does not introduce any major changes or incompatibilities,
but signifies that glmmTMB is considered stable and reliable for general use.

  \subsection{NEW FEATURES}{
    \itemize{
      \item new \code{map} argument to \code{glmmTMB} allows for some
      parameter values to be fixed
      (see \code{?TMB::MakeADFun} for details)
      \item new \code{optimizer} and
      \code{optArgs} arguments to \code{glmmTMBControl} allow use of
      optimizers other than \code{nlminb}
      \item \code{predict} can make population-level predictions
       (i.e., setting all random effects to zero).
       See \code{?predict.glmmTMB} for details.
       \item \code{beta_family} now allows zero-inflation;
       new \code{ziGamma} family (minor modification of
       \code{stats::Gamma}) allows zero-inflation
       (i.e., Gamma-hurdle models)
    }
  } % new features
  \subsection{BUG FIXES}{
  \itemize{
      \item \code{vcov(., full=TRUE)} (and hence profiling) now work for models with \code{dispformula=~0}
      \item Documentation fix: when \code{family=genpois}, the
      index of dispersion is known as phi^2.
      \item \code{Anova} now respects the \code{component} argument (GH
      #494, from @eds-slim)
      \item \code{predict} now works when contrasts are set on factors
      in original data (GH #439, from @cvoeten)
      \item \code{bootMer} now works with models with Bernoulli
      responses (even though \code{simulate()} returns a two-column
      matrix in this case) (GH #529, @frousseu)
      \item better support for \code{emmeans} applied to zero-inflation
      or dispersion models (correct link functions) (Russ Lenth)
   }
  } % bug fixes
  \subsection{USER-VISIBLE CHANGES}{
  \itemize{
      \item \code{sigma(.)} now returns \code{NA} for models with
      non-trivial dispersion models (i.e. models with more than one
      dispersion parameter) (raised by GH #533, from @marek-tph)
      \item \code{VarCorr} no longer prints residual variances for
      models with \code{dispformula=~0}
      \item the \code{model.matrix()} and \code{terms()} methods
      for \code{glmmTMB} objects have been slightly modified
    }
  } % user-visible changes
} % version 1.0.0

\section{CHANGES IN VERSION 0.2.3}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{ranef} now returns information about conditional variances (as
      attributes of the individual random effects terms) by default;
      this information can easily be retrieved by
      \code{as.data.frame(ranef(.))}.
      \item \code{coef} method now available: as in \code{lme4}, returns
      sum of fixed + random effects for each random-effects
      level. (Conditional variances for \code{coef} \emph{not} yet available.)
      \item simulate works for models with genpois family
      \item parametric bootstrapping should work, using
      \code{\link[lme4]{bootMer}} from the \code{lme4} package as a front end.
    } % itemize
  } % new features
  \subsection{BUG FIXES}{
    \itemize{
      \item models with multiple types of RE (e.g. ar1 and us)
      may have failed previously (GH #329)
      \item \code{predict} was not handling data-dependent predictors (e.g. \code{poly}, \code{spline}, \code{scale}) correctly
      \item \code{profile} now works for models without random effects
    }
  } % bug fixes
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item The value returned from \code{simulate} for binomial models
      is now a non-standard data frame where each element contains a
      two-column matrix (as in the base-R \code{\link{simulate}} method for
      binomial GLMS).
    } % itemize
  } % user-visible
} % version 0.2.3
\section{CHANGES IN VERSION 0.2.2}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item REML is now an option (GH #352). It is typically only for
      Gaussian response variables, but can also be useful
      for some non-Gaussian response variables if used with caution
      (i.e. simulate a test case first).
    }
  }
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Because family functions are now available for all
      families that have been implemented in the underlying TMB
      code, specifying the \code{family} argument as a raw list (rather than as a family
      function, the name of a family function, or the output of such a
      function) is now deprecated.
    }
  }
}

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

