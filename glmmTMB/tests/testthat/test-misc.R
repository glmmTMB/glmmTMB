stopifnot(require("testthat"),
          require("glmmTMB"))

test_that("predhess handling", {
  load(system.file("test_data", "cygu.rda", package = "glmmTMB"))
  jointmodel_tmb <- suppressWarnings(glmmTMB(model_form_joint,
                 data=model_df,
                 control=glmmTMBControl(
                     profile=TRUE
                 ),
                 family = binomial(link = "logit")
                 ))
  ## lame test, but we just want to make sure not to get error
  ## Error in fitTMB(TMBStruc) : object 'g' not found
  expect_is(jointmodel_tmb, "glmmTMB")
})
