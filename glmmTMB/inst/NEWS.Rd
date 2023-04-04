\newcommand{\PR}{\Sexpr[results=rd]{tools:::Rd_expr_PR(#1)}}

\name{NEWS}
\title{glmmTMB News}
\encoding{UTF-8}

\section{CHANGES IN VERSION 1.1.7}{
  \itemize{
    \item \code{up2date()} adds a dispersion component to family
    objects where required, for compatibility when checking
    stored fits across R versions
  }
}

\section{CHANGES IN VERSION 1.1.6 (2023-03-19)}{
  \subsection{BUG FIXES}{
    \itemize{
      \item fix bug in diagnose for Tweedie, other models with 'psi'
      parameter (@nalimilan, GH #135)
    }
  }  % bug fixes
  \subsection{NEW FEATURES}{
    \itemize{
      \item added a \code{doOptim} argument to \code{fitTMB}, to return
      the constructed TMB object without fitting the parameters
      \item new (experimental) function \code{simulate_new}, to simulate
      from a formula, list of parameters, and covariate data (rather
      than from a fitted object)
      \item \code{emmeans.glmmTMB} method adds two options
      for the \code{component} argument: "response" and "cmean",
      corresponding  to \code{type = "response"} and \code{type = "conditional"} 
      in \code{predict.glmmTMB} (Russ Lenth)
      \item new covariance structure \code{homdiag} (homogeneous diagonal;
      the existing diagonal covariance structure, \code{diag}, assumes heterogeneous
      variances)					    
    }
  } % new features
 \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item The \code{emmeans} method for glmmTMB fits now returns infinite
      "df" (i.e., normal- rather than t-based CIs and tests) for non-Gaussian
      families, consistent with \code{glm} and other packages (GH #893)
    }
  } % user-visible changes
}

\section{CHANGES IN VERSION 1.1.5 (2022-11-16)}{
    \subsection{BUG FIXES}{
      \itemize{
	\item improve predict-handling of complex bases (GH #632, #845,
	#853)
	\item all standard deviations are now printed in output for
	models using \code{cs()} (GH #851)
	\item corrected conditional and response predictions for
	truncated distributions (GH #634, #860, #873)
	\item \code{ranef()} now works correctly for families
	with extra parameters (Tweedie etc.) (GH #870)
      } % itemize
    } % bug fixes
    \subsection{USER-VISIBLE CHANGES}{
      \itemize{
	\item glmmTMB has switched to using a different (newer, under
	active development) autodifferentiation library under the hood
	(from \code{CppAD} to \code{TMBad}). This is likely to lead to
	small changes in estimates, including tipping marginally stable
	computations to instability or vice versa (e.g. presence or
	absence of convergence warnings, positive-definite Hessian
	warnings, \code{NaN} values of standard errors, etc.) You can revert
	to using the older autodiff engine by commenting out the line
	\code{PKG_CPPFLAGS = -DTMBAD_FRAMEWORK} in \code{src/Makevars}
	and reinstalling the package (make sure to delete .o/.so files
	from the \code{src} directory first if installing from the package
	directory, as the Makefile doesn't recognize know that this
	change requires recompilation).
	\item glmmTMB now warns if fixed-effect model matrices are
	rank deficient (i.e., perfectly collinear predictors); this
	warning can be suppressed by setting
	\code{glmmTMBControl(rank_check="skip")}
	(\code{rank_check = "stop"} will throw an error).
	If \code{rank_check="adjust"}, glmmTMB
	will automatically drop appropriate columns from the model
	matrix  (Daniel B. Stouffer)
	\item the vector of "extra" family parameters
	(Tweedie power, Student-t df,
	etc.) has been renamed from "thetaf" to "psi"; \code{start}
	and \code{map} arguments that set this parameter will need to
	be changed. Users will need to run \code{up2date()}
	when loading stored model objects from previous versions of the
	package.
	\item \code{predict} now warns if extra (ignored) arguments are
	provided in \code{...}
      } % itemize
    } % user-visible changes
    \subsection{NEW FEATURES}{
      \itemize{
	\item Student-t response distribution is now implemented (see
	\code{\link{t_family}})
	\item ordered beta regression as in Kubinec (2022), for
	proportion data containing exact 0 and 1 values, is now
	implemented (\code{\link{ordbeta}})
	\item \code{glmmTMBControl} now has a \code{conv_check} argument
	that allows suppressing convergence warnings
	(the intended use is when these warnings are irrelevant,
	e.g. when running small examples for testing purposes)
    } % itemize
  } % new features
} % version 1.1.5

\section{CHANGES IN VERSION 1.1.4 (2022-07-12)}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item row names of \code{confint} output for random effects
      parameters have changed (new format is
      \code{Std.Dev} . (term) | (grouping variable) for standard
      deviations, \code{Cor} . (term1) . (term2) | (grouping variable)
      for correlations)
    } % itemize
  } % user-visible changes
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{predict(., "zprob")} now returns 0 and
      \code{predict(., "zlink")} returns \code{-Inf} for all observations
      for models without zero-inflation (GH #798, Brenton Wiernik) [was
      previously supposed to throw an error, but incorrectly returned
      conditional values]
      \item bug fixes and other improvements for \code{diagnose} (inverted Z-score; now handles
      models without random effects)
      \item \code{confint} now works for models with more than one
      random effect
      \item \code{confint} works better (although not completely) for
  models with mapped parameters
    } % itemize
  } % bug fixes
  \subsection{NEW FEATURES}{
    \itemize{
      \item now provides Pearson residuals for zero-inflated and
      variable-dispersion models (Brenton Wiernik)
    } % itemize
  } % new features
} % version 1.1.4

\section{CHANGES IN VERSION 1.1.3 (2022-03-13)}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item minor improvements in \code{diagnose()}
      }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item offset variables with attributes now work properly
      (previously threw an error; now stripped before being passed to
      TMB)
      \item \code{emmeans} methods now work when \code{component}
      is non-default (GH #780, @rvlenth/@marosteg)
      \item \code{vcov(., full = TRUE)} is now named for models with
      multiple variance components
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item implemented working residuals (\code{residuals(., type =
	"working")}; GH #776, @lionel68)
      \item new option \code{print_trivial} for the \code{print} method
      for fixed effects (\code{fixef} objects); contributed by @d-morrison
    }
  }
}
\section{CHANGES IN VERSION 1.1.2.3 (2021-09-20)}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Double-bar notation (\code{(x+y||g)}) is now translated to a
  diagonal-covariance term (\code{diag(x+y|g)}) rather than being split
  into separate random effects terms as in  \code{lme4}. This should not
  change modeling results, but may change their
  presentation/ordering/etc.. (This is also a bug fix, as double-bar
  notation was not working in several previous versions.)
  \item \code{glmmTMB} now issues a warning when (1) \code{$} is used
  within formulas or (2) the \code{data} argument is not specified (the
  latter warning can be suppressed by specifying
  \code{data=NULL}).
}
 }
  \subsection{NEW FEATURES}{
    \itemize{
      \item New (experimental) function \code{up2date} for updating
      stored \code{glmmTMB} fits that were created with an earlier
      version of \code{TMB} than the one used when \code{glmmTMB} was
      compiled to binary/installed from source
      \item Utility functions \code{dtruncnbinom1},
      \code{dtruncnbinom2}, \code{dtruncpoisson} for k-truncated
      count distributions
    }
  }
}

\section{CHANGES IN VERSION 1.1.2.2 (2021-09-01)}{
  \itemize{
    \item This is an administrative release (minor revisions for CRAN).
  }
}

\section{CHANGES IN VERSION 1.1.2 (2021-07-20)}{
  \subsection{BUG FIXES}{
    \itemize{
      \item resolved OpenMP thread-safety issues on Windows
      \item resolved bug that caused Tweedie models to crash on Solaris
      \item resolved problems with vignettes on Solaris (GH #721)
    }
  } % bug fixes
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item improved control of OpenMP threading for prediction,
      profiling etc.
    }
  } % user-visible changes
} % version 1.1.2

\section{CHANGES IN VERSION 1.1.0 (2021-06-24)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item reduced rank covariance for GLVMs implemented by M.McGillycuddy
      (see covstruct vignette for details)
      \item \code{diagnose} function to investigate potential causes of convergence problems
      \item improved parallel processing (GH #620 #652)
      \item truncated nbinom2 family now includes a \code{variance}
      component
    }
  } % new features
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{Anova} with \code{type="III"} now handles
      \code{component} argument correctly, more robust to trivial models
      \item fixed a typo/omission in the type-3 Anova method that made
      zi Anova break in some conditions (GH #674)
      \item fixed bugs/inconsistencies in handling of mapped parameters
      (GH #678)
      \item \code{confint} with \code{parm="beta_"} or
      \code{parm="theta_"} now work correctly with more complex models
      (e.g. including both zero inflation and random effects)
      (reported by @MKie45 on Stack Overflow)
      \item \code{confint} works for single-parameter models and those
      with a \code{dispformula} (GH #622)
      \item mapped (fixed) variables could give incorrect predictions (GH #644)
      \item \code{simulate} is more robust for truncated_nbinom1 and
      truncated_nbinom2 (GH #572)
    }
  }
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item "mapped" parameters (i.e., fixed by user rather than
      optimized) are now given variances/standard deviations of NA rather
      than 0 in \code{vcov(., include_mapped=TRUE)} and by extension in
      \code{summary}; hence Z-statistics and P-values will also be NA for
      these parameters
      \item row ordering has changed in \code{confint} output data
      frames (random effects parameters come \emph{last}, matching the
      row/column order in \code{vcov(., full=TRUE)})
      \item new \code{fast} flag for predictions decreases memory use
      and computational time (only if \code{newdata}, \code{newparams}
      not specified); default in \code{fitted()} method
      \item improved robustness of beta-binomial fits
      (results of fitting such models may change slightly from previous
      versions)
      \item consistent predictions between link and inverse-link (GH #696)
      \item improved vignette titles
      \item The \code{emm_basis} method for \code{glmmTMB} objects now accepts a user-specified covariance matrix (\code{vcov.} argument)
    }
  }
}

\section{CHANGES IN VERSION 1.0.2.1 (2020-07-02)}{
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
      \item structured covariance models now work in zero-inflation components (GH #579)
      \item documentation of formula for variance in beta family (GH #595)
    }
  } % bug fixes
}
\section{CHANGES IN VERSION 1.0.1 (2020-03-15)}{
  \itemize{
    \item updated for R-devel changes (R 4.0.0 will set stringsAsFactors=FALSE by default)
  }
}
\section{CHANGES IN VERSION 1.0.0 (2020-02-03)}{

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

\section{CHANGES IN VERSION 0.2.3 (2019-01-11)}{
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
\section{CHANGES IN VERSION 0.2.2 (2018-07-03)}{
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

\section{CHANGES IN VERSION 0.2.1 (2018-05-11)}{
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

