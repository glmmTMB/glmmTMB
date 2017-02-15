library(glmmTMB)
data(Salamanders)

#This is a collection of models that do not converge.
#They are probably overpaameterized.
zinb1m1 = glmmTMB(count~spp * mined + (1|site), zi=~(1|site), Salamanders, family="nbinom1")
zinb1m2 = glmmTMB(count~spp * mined + (1|site), zi=~spp, Salamanders, family="nbinom1")
zinb1m3 = glmmTMB(count~spp * mined + (1|site), zi=~spp + mined, Salamanders, family="nbinom1")
zinb1m4 = glmmTMB(count~spp * mined + (1|site), zi=~spp * mined, Salamanders, family="nbinom1")
zinb1m5 = glmmTMB(count~spp * mined + (1|site), zi=~spp * mined+ (1|site), Salamanders, family="nbinom1")

zinb2m3 = glmmTMB(count~spp * mined + (1|site), zi=~spp * mined + (1|site), Salamanders, family="nbinom2")

bbmle::AICtab(zinb1m1,
zinb1m2,
zinb1m3,
zinb1m4,
zinb1m5,
zinb2m3)

#Genpois (this distribution is known for poor convergence)
m1=glmmTMB(count~spp + mined + (1|site), zi=~spp + mined, Salamanders, family="genpois")