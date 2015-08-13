PACKAGE=glmmTMB
VERSION=0.0.0.9000
TARBALL=${PACKAGE}_${VERSION}.tar.gz
ZIPFILE=${PACKAGE}_${VERSION}.zip

all:
	make doc-update
	make build-package
	make install
	make pdf

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

