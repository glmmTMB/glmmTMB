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
	grep "_link[ ]*=" $(PACKAGE)/src/glmmTMB.cpp | sed s/_link//g >> $@
	echo ")" >> $@

	echo ".valid_family <- c(" >> $@
	grep _family.*= $(PACKAGE)/src/glmmTMB.cpp | sed s/_family//g >> $@
	echo ")" >> $@

	echo ".valid_covstruct <- c(" >> $@
	grep "_covstruct *=" $(PACKAGE)/src/glmmTMB.cpp | sed s/_covstruct//g >> $@
	echo ")" >> $@

	echo ".valid_zipredictcode <- c(" >> $@
	grep _zipredictcode.*= $(PACKAGE)/src/glmmTMB.cpp | sed s/_zipredictcode//g >> $@
	echo ")" >> $@

	echo ".valid_prior <- c(" >> $@
	grep _prior.*= $(PACKAGE)/src/glmmTMB.cpp | sed s/_prior//g >> $@
	echo ")" >> $@

	echo ".valid_vprior <- c(" >> $@
	grep "_vprior *=" $(PACKAGE)/src/glmmTMB.cpp | sed s/_vprior//g >> $@
	echo ")" >> $@

doc-update: $(PACKAGE)/R/*.R
	echo "suppressWarnings(roxygen2::roxygenize(\"$(PACKAGE)\",roclets = c(\"collate\", \"rd\")))" | $(R) --slave
	@touch doc-update

## vignettes

## list of vignette inputs:
rnw_vig += glmmTMB model_evaluation 
rmd_vig += covstruct mcmc miscEx sim troubleshooting parallel hacking

docdir = $(PACKAGE)/inst/doc
vigdir = $(PACKAGE)/vignettes

docpdf :=  $(rnw_vig:%=${docdir}/%.pdf)
dochtml := $(rmd_vig:%=${docdir}/%.html)

## LaTeX/BibTeX must run in the same directory as the file ...
$(vigdir)/%.pdf: $(vigdir)/%.[Rr]nw
	cd $(vigdir); export NOT_CRAN=true; echo "knitr::knit2pdf(basename(\"$<\"))" | $(R) --slave

## ditto for vignette building
%.html: %.rmd
	cd $(vigdir); export NOT_CRAN=true; echo "rmarkdown::render(basename(\"$<\"))" | $(R) --slave

## files in doc dir => generate files in vig dir using downstream rules,
## then move them
$(docdir)/%: $(vigdir)/%
	mv $< $@

$(vigdir)/model_evaluation.html: $(vigdir)/model_evaluation.rmd texreg

vignette-update: ${docpdf} ${dochtml}

vigdatadir=glmmTMB/inst/vignette_data
vignette-data: $(vigdatadir)/mcmc.rda $(vigdatadir)/troubleshooting.rda $(vigdatadir)/model_evaluation.rda
## haven't figured out all of these rules yet
## R CMD BATCH corresponding *.R files in vigdatadir ...


####
namespace-update :: $(PACKAGE)/NAMESPACE
$(PACKAGE)/NAMESPACE: $(PACKAGE)/R/*.R
	echo "suppressWarnings(roxygen2::roxygenize(\"$(PACKAGE)\",roclets = c(\"namespace\")))" | $(R) --slave

build-package: $(TARBALL)
$(TARBALL): $(PACKAGE)/NAMESPACE $(CPP_SRC)
	$(R) CMD build --resave-data=no $(PACKAGE)

install: $(TARBALL)
	export NOT_CRAN=true; $(R) CMD INSTALL --preclean $<
	@touch $@

## To enable quick compile, run from R:
##    library(TMB); precompile(flags="-O0 -g")
quick-install: enum-update $(PACKAGE)/src/glmmTMB.so
	$(R) CMD INSTALL $(PACKAGE)

$(PACKAGE)/src/glmmTMB.so: $(PACKAGE)/src/glmmTMB.cpp
	cd $(PACKAGE)/src; echo "library(TMB); compile('glmmTMB.cpp','-O0 -g',libinit=FALSE)" | $(R) --slave

unexport TEXINPUTS
pdf: $(PACKAGE).pdf
$(PACKAGE).pdf: $(PACKAGE)/man/*.Rd
	rm -f $(PACKAGE).pdf
	$(R) CMD Rd2pdf --no-preview $(PACKAGE)

build:
	$(R) CMD build $(PACKAGE)

check: $(TARBALL)
	export _R_S3_METHOD_LOOKUP_BASEENV_AFTER_GLOBALENV=TRUE
	$(R) CMD check $(TARBALL)

check-cran: $(TARBALL)
	export _R_S3_METHOD_LOOKUP_BASEENV_AFTER_GLOBALENV=TRUE
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

revdep_check:
	Rscript -e 'library("revdepcheck"); revdep_reset("glmmTMB"); revdep_check("glmmTMB", num_workers = 2, timeout = as.difftime(60, units="mins"))'
