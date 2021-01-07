#! /bin/bash
## run from head directory
ls glmmTMB/vignettes/*.Rnw | \
    xargs grep usepackage | \
    egrep -v '% *\\' | \
    egrep -v knit_hooks | \
    sed -e 's/^[^{]*{//' | \
    sed -e 's/}$//' | \
    sort | uniq | \
    paste -sd ',' - | \
    sed -e 's/,/","/g' | \
    sed -re 's/(^|$)/"/g'

