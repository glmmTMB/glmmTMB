R=R
# -> you can do    R=R-devel  make ....

PACKAGE=glmmTMB
# get VERSION from glmmTMB/DESCRIPTION  
## ("::" = expand only  once, but doesn't work in make <= 3.81)
VERSION := $(shell sed -n '/^Version: /s///p' glmmTMB/DESCRIPTION)

TARBALL := $(PACKAGE)_$(VERSION).tar.gz
ZIPFILE := =$(PACKAGE)_$(VERSION).zip

CPP_SRC := $(PACKAGE)/src/*.cpp

all:
	make enum-update
	make doc-update
	make build-package
	make install
	make pdf

enum-update:: $(PACKAGE)/R/enum.R
$(PACKAGE)/R/enum.R: $(PACKAGE)/src/glmmTMB.cpp
	echo '## Auto generated - do not edit by hand' > $@
	echo ".valid_link <- c(" >> $@
	grep _link.*= $(PACKAGE)/src/glmmTMB.cpp | sed s/_link//g >> $@
	echo ")" >> $@

	echo ".valid_family <- c(" >> $@
	grep _family.*= $(PACKAGE)/src/glmmTMB.cpp | sed s/_family//g >> $@
	echo ")" >> $@

	echo ".valid_covstruct <- c(" >> $@
	grep _covstruct.*= $(PACKAGE)/src/glmmTMB.cpp | sed s/_covstruct//g >> $@
	echo ")" >> $@

	echo ".valid_zipredictcode <- c(" >> $@
	grep _zipredictcode.*= $(PACKAGE)/src/glmmTMB.cpp | sed s/_zipredictcode//g >> $@
	echo ")" >> $@

doc-update: $(PACKAGE)/R/*.R
	echo "library(roxygen2);roxygenize(\"$(PACKAGE)\",roclets = c(\"collate\", \"rd\", \"namespace\"))" | $(R) --slave
	@touch doc-update

## FIXME: build *all* .Rnw files in the directory
vignette-update: $(PACKAGE)/vignettes/*.Rnw
	cd $(PACKAGE)/vignettes; echo "library(knitr);knit2pdf('glmmTMB.Rnw')" | $(R) --slave
	 mv $(PACKAGE)/vignettes/glmmTMB.pdf $(PACKAGE)/inst/doc
	@touch vignette-update

namespace-update :: $(PACKAGE)/NAMESPACE
$(PACKAGE)/NAMESPACE: $(PACKAGE)/R/*.R
	echo "library(roxygen2);roxygenize(\"$(PACKAGE)\",roclets = c(\"namespace\"))" | $(R) --slave

##	sed -i -e "s/importFrom(lme4,sigma)/if(getRversion()>='3.3.0') importFrom(stats, sigma) else importFrom(lme4,sigma)/" $(PACKAGE)/NAMESPACE

build-package: $(TARBALL)
$(TARBALL): $(PACKAGE)/NAMESPACE $(CPP_SRC)
	$(R) CMD build --resave-data=no $(PACKAGE)

install: $(TARBALL)
	$(R) CMD INSTALL --preclean $<
	@touch $@

## To enable quick compile, run from R:
##    library(TMB); precompile(flags="-O0 -g")
quick-install: enum-update $(PACKAGE)/src/glmmTMB.so
	$(R) CMD INSTALL $(PACKAGE)

$(PACKAGE)/src/glmmTMB.so: $(PACKAGE)/src/glmmTMB.cpp
	cd $(PACKAGE)/src; echo "library(TMB); compile('glmmTMB.cpp','-O0 -g')" | $(R) --slave

unexport TEXINPUTS
pdf: $(PACKAGE).pdf
$(PACKAGE).pdf: $(PACKAGE)/man/*.Rd
	rm -f $(PACKAGE).pdf
	$(R) CMD Rd2pdf --no-preview $(PACKAGE)

build:
	$(R) CMD build $(PACKAGE)

check: $(TARBALL)
	$(R) CMD check $(TARBALL)

check-cran: $(TARBALL)
	$(R) CMD check --as-cran $(TARBALL)

## *NOT* using 'R --vanilla' : then cannot find testthat, TMB, etc they are installed into R's "system" library

test:
	echo "devtools::test('glmmTMB')" | $(R) --slave

quick-check: quick-install ex-test

ex-test:
	echo "library(glmmTMB); example(glmmTMB)" | $(R) --slave


unlock:
	\rm -rf `Rscript --vanilla -e 'writeLines(.Library)'`/00LOCK-glmmTMB
#               ------------------------------------------ = R's system library
#	rm -rf ${R_LIBS}/00LOCK-glmmTMB
##               ^^^^^^^ This only works if R_LIBS contains a single directory and the same that 'R CMD INSTALL' uses..

clean:
	\rm -f install doc-update
