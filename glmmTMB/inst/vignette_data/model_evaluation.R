library(glmmTMB)
library(MuMIn)
library(DHARMa)
owls_nb1 <- glmmTMB(SiblingNegotiation ~ FoodTreatment*SexParent +
                        (1|Nest)+offset(log(BroodSize)),
                    contrasts=list(FoodTreatment="contr.sum",
                                   SexParent="contr.sum"),
                    family = nbinom1,
                    zi = ~1, data=Owls,
                    na.action=na.fail)
owls_nb1_dredge <- MuMIn::dredge(owls_nb1)
owls_nb1_simres <- DHARMa::simulateResiduals(owls_nb1)
source(system.file("other_methods","influence_mixed.R", package="glmmTMB"))
owls_nb1_influence_time <- system.time(
    owls_nb1_influence <- influence_mixed(owls_nb1, groups="Nest", ncores=1))
save("owls_nb1",
     "owls_nb1_simres",
     "owls_nb1_dredge",
     "owls_nb1_influence",
     "owls_nb1_influence_time",
     file="model_evaluation.rda",
     version=2 ## for compatibility with R < 3.6.0
     )
