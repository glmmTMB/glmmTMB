PACKAGE=glmmTMB
VERSION=0.0.0.9000
TARBALL=${PACKAGE}_${VERSION}.tar.gz
ZIPFILE=${PACKAGE}_${VERSION}.zip

all:
	make enum-update
	make doc-update
	make build-package
	make install
	make pdf

enum-update:
	echo "## Auto generated - do not edit by hand" > glmmTMB/R/enum.R

	echo "\n.valid_link <- c(" >> glmmTMB/R/enum.R
	grep _link.*= glmmTMB/src/glmmTMB.cpp | sed s/_link//g >> glmmTMB/R/enum.R
	echo ")" >> glmmTMB/R/enum.R

	echo "\n.valid_family <- c(" >> glmmTMB/R/enum.R
	grep _family.*= glmmTMB/src/glmmTMB.cpp | sed s/_family//g >> glmmTMB/R/enum.R
	echo ")" >> glmmTMB/R/enum.R

	echo "\n.valid_covstruct <- c(" >> glmmTMB/R/enum.R
	grep _covstruct.*= glmmTMB/src/glmmTMB.cpp | sed s/_covstruct//g >> glmmTMB/R/enum.R
	echo ")" >> glmmTMB/R/enum.R

doc-update:
	echo "library(roxygen2);roxygenize(\"$(PACKAGE)\",roclets = c(\"collate\", \"rd\"))" | R --slave

build-package:
	R CMD build --resave-data=no $(PACKAGE)

install:
	make build-package
	R CMD INSTALL --preclean $(TARBALL)

## To enable quick compile, run from R:
##    library(TMB); precompile(flags="-O0 -g")
quick-install:
	make enum-update
	cd $(PACKAGE)/src; echo "library(TMB); compile('glmmTMB.cpp','-O0 -g')" | R --slave
	R CMD INSTALL $(PACKAGE)

unexport TEXINPUTS
pdf:
	rm -f $(PACKAGE).pdf
	R CMD Rd2pdf --no-preview $(PACKAGE)

check:
	R CMD check $(PACKAGE)

unlock:
	rm -rf ${R_LIBS}/00LOCK-glmmTMB

