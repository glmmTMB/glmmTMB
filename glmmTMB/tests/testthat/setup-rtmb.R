local_useRTMB <- function(flag, .env = parent.frame()) {
  old_use_rtmb <- glmmTMB::useRTMB()
  withr::defer(glmmTMB::useRTMB(old_use_rtmb), envir = .env)
  glmmTMB::useRTMB(flag)
}
