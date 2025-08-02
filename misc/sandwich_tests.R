# Few additional tests for the sandwich calculations.

owls_nb1 <- glmmTMB(SiblingNegotiation ~ FoodTreatment*SexParent +
                        (1|Nest)+offset(log(BroodSize)),
                    contrasts=list(FoodTreatment="contr.sum",
                                   SexParent="contr.sum"),
                    family = nbinom1,
                    zi = ~1, data=Owls)
clusters <- getGroups(owls_nb1)
nlevels(clusters)
# so only 27 here.

system.time(meat_res <- meatHC(owls_nb1))
# takes about 2 seconds with retape() each time.

expected_res <- matrix(
    c(149.13, -186.32, 6.06, 41.7, -186.32, 1376.14, -105.8, 
-396.93, 6.06, -105.8, 575.56, -39.35, 41.7, -396.93, -39.35,
647.42),
    nrow = 4, ncol = 4
)

expect_equal(meat_res, expected_res, tolerance = 1e-4, check.attributes = FALSE)
