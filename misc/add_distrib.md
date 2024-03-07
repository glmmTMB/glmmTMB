- add a value for your distribution to `valid_family` (at l. 19 of `src/glmmTMB.cpp`
   - if it is considered a new class of distributions, number by 100s
   - if it is a subcategory of an existing class (e.g. a truncated or under/overdispersed variant), number by 1s
- run `make enum-update`
- add family function to `R/family.R`
   
