## configuration stuff

Removed because it seemed to be doing more harm than good; perhaps `src/Makevars` is good enough for what we want? The idea was to allow use of OpenMP wherever it is available, but in the absence of a `src/Makevars` file it seems that `configure` isn't automatically run/`src/Makevars` isn't automatically regenerated (and therefore the OpenMP flags aren't added) when installing from tarball/GitHub repo?

Also, don't know exactly what we should do about `Makevars.win`

- Are there good models of packages using OpenMP that we can follow?
- Could re-read https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Configure-and-cleanup ...
