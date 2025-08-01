\newcommand{\PR}{\Sexpr[results=rd]{tools:::Rd_expr_PR(#1)}}

\name{NEWS}
\title{glmmTMB News}
\encoding{UTF-8}

\section{CHANGES IN VERSION 1.1.11.9000}{
  \subsection{BUG FIXES}{
    \itemize{
      \item default value for \code{allow.new.levels} in
      \code{predict.glmmTMB} is now NULL, meaning that
      new levels in grouping variables are allowed when
      predicting at the population level (i.e., \code{re.form}
      is \code{NA} or \code{~0}) - unless \code{allow.new.levels} is explicitly set to
      TRUE (GH#1189, @fisher-j)
      \item implementation of \code{hetar1} covstruct was wrong, fixed now
      (GH #1204, Daniel Sabanés Bové and Nik Krieger)
      \item prediction on rank-deficient X matrices could throw a false
      positive error (Erin McCallum, \url{https://stackoverflow.com/q/77517125/190277})
    }
  } % bug fixes
  \subsection{NEW FEATURES}{
    \itemize{
      \item new argument \code{full_corr} for \code{glmmTMBControl} to
      enable/disable computation and storage of full correlation
      matrices for random effect terms (helps with GH#995)
    }
  } % new features
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item improved formatting of random effects standard
      deviations/variances/correlation matrices
      \item mismatches in contrasts between original data and prediction
      data now throw a warning rather than an error (related to GH #1197)
    }
  } % user-visible changes
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{predict} method no longer fails on subsets of data
      involving ordered factors (resolving an issue comparing contrasts
      associated with factors in original data and \code{newdata}) (GH
      #1197, @fisher-j); we no longer check contrasts of factors that
      appear only as grouping variables
    }
  }
} % 1.1.11.9000 (devel)

\section{CHANGES IN VERSION 1.1.11 (2025-04-02)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{glmmTMB} gains a \code{subset} argument (GH #1128,
      @strengejacke)
      \item added a heterogeneous-variance AR1 (\code{hetar1})
      covariance structure (GH #1095) (experimental, still prints badly)
      \item added the Bell distribution (\code{bell()}) as in
      Castellares et al. 2018 \doi{10.1016/j.apm.2017.12.014} (Hatice
      Tül Kübra Akdur)
      \item added aggregated predictions with bias correction as in Kindt-Larsen, Glemerec, et al. 2023 \doi{10.1098/rspb.2022.2570} and Thorson & Kristensen 2016 \doi{10.1016/j.fishres.2015.11.016}
    } % itemize
  } % new features
  \subsection{BUG FIXES}{
    \itemize{
      \item likelihood profiling now works for models with mapped parameters
      \item \code{glmmTMB} no longer changes the order of terms in
      fixed-effect model matrices (GH #1122, @dongwenluo). Note that
      this bug fix will change the order of results (parameter vectors)
      in some cases.
      \item \code{vcov()} behaviour improved for models with mapped
      parameters, especially the case where parameters are set equal rather
      than fixed to starting values (e.g. \code{map = list(beta =
	factor(c(1,1)))}) (GH #1120, @DrJerryTAO)
      \item fixed newly introduced bug in Pearson residuals for
      zero-inflated models (GH #1101, @strengejacke)
      \item fix bug in 'exotic' families (those such as truncated
      distributions using the internal \code{make_family} helper
      function) that caused errors when calling \code{effects::Effect()}
      (GH #1133, @strengejacke)
      \item fixed bug blocking reduced-rank models with binomial response
      (GH #1151, @toddvogel1628)
      \item fixed minor bug with \code{confint} applied to models with
      random effects in the dispersion model
    } % itemize
  } % bug fixes
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item headline of \code{print} and \code{summary} output
      now labels the minimum of the objective function
      (correctly) as "-2*log(L)" rather than "deviance" (GH #1156, @ladin100)
    } % itemize
  } % user-visible changes

} % 1.1.11
  
\section{CHANGES IN VERSION 1.1.10  (2024-09-26)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item added random effect structure propto to fit multivariate random effects proportional to a known variance-covariance matrix. This feature is
      EXPERIMENTAL; please post any problems at \url{https://github.com/glmmTMB/glmmTMB/pull/1068}. See \code{vignette("covstruct")}
      \item added "nbinom12" family after Lindén and Mäntyniemi (2011)
      \item random effects, including smooth terms with \code{s()}, can
      now be included in the dispersion model. This feature is
      EXPERIMENTAL; please post any problems to the issues list
      \item added "skewnormal" family (@psmorris) (EXPERIMENTAL:
      some manual adjustment of starting values may be required)
      \item \code{predict()} now offers \code{type = "latent"}, which
      returns the values of the latent variables (with conditional
      standard deviations if \code{se.fit = TRUE})
      \item now works with automatic parallelization based on the
      underlying TMBad autodiff library (specify
      \code{parallel=list(..., autopar=TRUE)} in \code{glmmTMBControl}
      or set \code{options(glmmTMB.autopar = TRUE)})
    } % itemize
  } % new features
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item the name of the
      vector of fixed-effect parameters for the dispersion
      model has changed from \code{betad} to \code{betadisp};
      code that specifies parameter values (e.g. as part
      of a \code{start} argument or in \code{newparams} for
      \code{simulate_new}) may need to be adjusted. (For
      developers: There are similar name changes to other internal
      model components, e.g. \code{Xd} becomes \code{Xdisp}.)
    } % itemize
  } % user-visible changes
  \subsection{BUG FIXES}{
    \itemize{
      \item fixed bug in \code{simulate_new} for \code{family = "beta_family"}
      \item simulation now works for models fitted with the (scaled) t
	distribution (GH #1024)
      \item \code{vcov} works better for cases where \code{map} is used
      to fix sets of parameters equal to each other
      \item fix printing bug for zero-inflation covariance matrices (SO
      78393784)
      \item \code{up2date} now works for models with mapped parameters
      (GH #874)
      \item \code{ranef()} now works properly for reduced-rank models
      \item Pearson residuals now work for models with non-trivial
      dispersion components (GH #1053)
      \item bug fixes to \code{get_cor()}, \code{put_cor()} utility functions
      \item \code{mgcv} smooths with no unpenalized components
      (e.g. \code{s(..., bs  = "sos")}) now work
      \item better calculation of \code{df.residual} when
      \code{REML=TRUE} (#1039)=
    } % itemize
  } % bug fixes
  \subsection{OTHER CHANGES}{
    \itemize{
      \item now imports random effect machinery from the
      \code{reformulas} package rather than from \code{lme4}
    } % itemize
  } % other changes
} % 1.1.10
  
\section{CHANGES IN VERSION 1.1.9 (2024-03-20)}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item the underlying parameterization of the dispersion for
      Gaussian models has changed from the variance to the standard deviation scale,
      with the following user-visible consequences:
      \itemize{
	\item results of fitting Gaussian models may change slightly,
	especially for unstable fits
	\item stored model objects need to have their \code{betadisp}
	parameters halved for consistency: use \code{up2date(..., adjust_gauss_disp = TRUE)}
	\item parameter estimates for dispersion components will change,
	e.g. from \code{fixef(.)$disp} or \code{confint()} for models
	with dispersion components
	\item any operations that directly handle dispersion parameters
	(e.g. \code{offset} terms for dispersion models) should be
	specified on the SD rather than the variance scale
      }
      \item interpretation of the \emph{weights} variable for
      binomial-type GL(M)Ms has changed. Previously, the \emph{weights}
      argument was \emph{ignored} for a vector- (rather than
      matrix-valued) numeric response, if all observations were either 0 or 1.
      Now the \emph{weights} variable is multiplied by the resonse
      variable to compute the number of successes (consistently with
      \code{stats::glm(., family = "binomial")}). (This change makes
      it easier to use \emph{weights} to specify the number of trials
      per observation for \code{simulate_new()}.)
      \item \code{simulate_new} gives useful warning and error
      messages about unrecognized parameter names and length mismatches
      (length mismatches previously gave only a warning)
    } % itemize
  } % user-visible changes
  \subsection{NEW FEATURES}{
    \itemize{
      \item lognormal-hurdle models now work (i.e., zero values
      can occur in the response if \code{ziformula} is specified)
      \item better checking for illegal (negative or non-integer) values
      in response variables
      \item experimental support for priors: see \code{?priors} and
      \code{vignette("priors", package = "glmmTMB")}
    }
  } % new features
  \subsection{BUG FIXES}{
    \itemize{
      \item predictions now work when \code{weights} variables have
      attributes (GH #977)
    } % itemize
  } % bug fixes
} % 1.1.9

\section{CHANGES IN VERSION 1.1.8 (2023-10-07)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item "lognormal" family available (log-Normal, parameterized
      by the mean and standard deviation on the \emph{data scale})
      \item an experimental implementation of penalized splines
      (and related smooth terms) is available, 
      based on \code{mgcv}. See the example using \code{s} in \code{?glmmTMB}.  
      Constructive feedback welcome at
      \url{https://github.com/glmmTMB/glmmTMB/issues/928}
      \item population-level prediction with new data no longer requires
      that the variables involved in the random effects be explicitly
      specified as \code{NA} (GH #923, Russ Lenth)
      \item the \code{simulate} method now works for models fitted with
      \code{family=ordbetareg} (GH #942, Daniel Lüdecke)
      \item deviance residuals are now available for some families
      (built-in families from base R (see \code{?family})
      plus \code{nbinom1}, \code{nbinom2}); deviance
      residuals for other families may be implemented on request.
      \item setting \code{option(glmmTMB_openmp_debug=TRUE)} will
      produce debugging output about the number of OpenMP threads used
      \item \code{getME(., "b")} returns the vector of conditional modes
      (BLUPs, in the case of linear mixed models)
    } % itemize
  } % new features
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item changes related to handling rank-deficient fits:
      \itemize{
	\item default value of \code{rank_check} in
	\code{glmmTMBControl} changed to "adjust" (i.e., rank-deficient
	columns of fixed-effect model matrices are automatically dropped, with
	a message, by default)
	\item the \code{include_mapped} argument of \code{vcov} and
	\code{confint} is changed to \code{include_nonest},
	controlling both mapped parameters and those
	dropped due to rank-deficiency, and is now \code{TRUE} by default for
	\code{vcov}
	\item \code{model.matrix} now returns the fixed-effect model matrix
	actually used in fitting (including dropping columns for
	rank-deficiency)
      } %inner itemize
      \item \code{glmmTMB} now accepts single-column matrices (e.g. as
	produced by \code{scale} as response variables (GH #937, @santoshbs)
    } % itemize
  } % user-visible changes
} % version 1.1.7.9000
  
\section{CHANGES IN VERSION 1.1.7 (2023-04-05)}{
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

