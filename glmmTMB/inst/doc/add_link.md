---
title: "adding link functions to glmmTMB"
---

- in `glmmTMB.cpp`
    - add your link to `valid_link`
    - add your link to the `switch` statement in `inverse_linkfun` (don't forget to `break;`)
- run `make enum-update` (this will rebuild `R/enum.R`)
- add link to the roxygen documentation in `R/family.R`; rebuild documentation (`make doc-update`)

	
- if your link will be used with binomial responses, add it to `logit_inverse_linkfun` 
- if your link will be used with neg binomial responses, add it to `log_inverse_linkfun()`
