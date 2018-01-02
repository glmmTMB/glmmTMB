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
      \code{confint(fitted,method="profile")} also works
     }
   }
   \subsection{BUG FIXES}{
     \itemize{
       \item built-in \code{betar()} family for Beta regression fixed
       (and name changed to \code{beta_family()})
     }
   }
   \subsection{USER-VISIBLE CHANGES}{
     \itemize{
       
     }
   }
 }
 
