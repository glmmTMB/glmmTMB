## Override RTMB method which, by mistake, lacks the '...' argument in
## RTMB version (<=2.0). This is an identical copy of the corrected
## upstream version:
setMethod("solve", signature("num", "num."),
          function(a, b, ...) {
              base::solve(a, b, ...)
          })
