R=R
# -> you can do    R=R-devel  make ....

PACKAGE=glmmTMB
VERSION=0.0.0.9000
TARBALL=${PACKAGE}_${VERSION}.tar.gz
ZIPFILE=${PACKAGE}_${VERSION}.zip

CPP_SRC = $(PACKAGE)/src/*.cpp

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
	echo "library(roxygen2);roxygenize(\"$(PACKAGE)\",roclets = c(\"collate\", \"rd\"))" | $(R) --slave
	@touch doc-update

namespace-update :: $(PACKAGE)/NAMESPACE
$(PACKAGE)/NAMESPACE: $(PACKAGE)/R/*.R
	echo "library(roxygen2);roxygenize(\"$(PACKAGE)\",roclets = c(\"namespace\"))" | $(R) --slave

build-package: $(TARBALL)
$(TARBALL): $(PACKAGE)/NAMESPACE $(CPP_SRC)
	$(R) CMD build --resave-data=no $(PACKAGE)

install: $(TARBALL)
	$(R) CMD INSTALL --preclean $(TARBALL)
	@touch install

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

check:
	$(R) CMD check $(PACKAGE)

## *NOT* using 'R --vanilla' : then cannot find testthat, TMB, etc they are installed into R's "system" library

quick-check: quick-install ex-test
	echo "source('glmmTMB/tests/AAAtest-all.R', echo=TRUE)" | $(R) --slave



unlock:
	rm -rf `Rscript --vanilla -e 'writeLines(.Library)'`/00LOCK-glmmTMB
#               ------------------------------------------ = R's system library
#	rm -rf ${R_LIBS}/00LOCK-glmmTMB
##               ^^^^^^^ This only works if R_LIBS contains a single directory and the same that 'R CMD INSTALL' uses..

test: ex-test
ex-test:
	echo "library(glmmTMB); example(glmmTMB)" | $(R) --slave


clean:
	\rm -f install doc-update
