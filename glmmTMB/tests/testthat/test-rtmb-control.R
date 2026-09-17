context("RTMB control option")

is_rtmb_fit <- function(x) {
    !is.null(x$obj$env$rtmb_data_env)
}

test_that("glmmTMBControl use_rtmb overrides backend for one fit", {
    skip_if_not_installed("RTMB")
    data("sleepstudy", package = "lme4")
    old_use_rtmb <- glmmTMB::useRTMB()
    withr::defer(glmmTMB::useRTMB(old_use_rtmb))

    glmmTMB::useRTMB(FALSE)
    m_rtmb <- glmmTMB(Reaction ~ Days, data = sleepstudy, se = FALSE,
                      control = glmmTMBControl(use_rtmb = TRUE))
    expect_true(is_rtmb_fit(m_rtmb))
    expect_false(glmmTMB::useRTMB())

    glmmTMB::useRTMB(TRUE)
    m_tmb <- glmmTMB(Reaction ~ Days, data = sleepstudy, se = FALSE,
                     control = glmmTMBControl(use_rtmb = FALSE))
    expect_false(is_rtmb_fit(m_tmb))
    expect_true(glmmTMB::useRTMB())
})

test_that("glmmTMBControl use_rtmb NULL leaves backend unchanged", {
    skip_if_not_installed("RTMB")
    data("sleepstudy", package = "lme4")
    old_use_rtmb <- glmmTMB::useRTMB()
    withr::defer(glmmTMB::useRTMB(old_use_rtmb))

    glmmTMB::useRTMB(FALSE)
    m_tmb <- glmmTMB(Reaction ~ Days, data = sleepstudy, se = FALSE,
                     control = glmmTMBControl(use_rtmb = NULL))
    expect_false(is_rtmb_fit(m_tmb))
    expect_false(glmmTMB::useRTMB())

    glmmTMB::useRTMB(TRUE)
    m_rtmb <- glmmTMB(Reaction ~ Days, data = sleepstudy, se = FALSE,
                      control = glmmTMBControl(use_rtmb = NULL))
    expect_true(is_rtmb_fit(m_rtmb))
    expect_true(glmmTMB::useRTMB())
})

test_that("glmmTMBControl use_rtmb validates scalar logical input", {
    expect_error(glmmTMBControl(use_rtmb = c(TRUE, FALSE)), "use_rtmb")
    expect_error(glmmTMBControl(use_rtmb = NA), "use_rtmb")
    expect_error(glmmTMBControl(use_rtmb = 1), "use_rtmb")
    expect_error(glmmTMBControl(use_rtmb = "TRUE"), "use_rtmb")
})
