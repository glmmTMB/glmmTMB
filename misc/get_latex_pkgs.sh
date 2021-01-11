#! /bin/bash
## run from head directory
## this isolates all of the \usepackage{} packages into
## a quoted, comma-separated list
## missing: chicago.bst -> chicago, american.ldf -> babel-english
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

