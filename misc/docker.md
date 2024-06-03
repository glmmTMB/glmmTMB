## testing with docker

Should make this into a proper Dockerfile ...

```
docker pull rocker/r2u
docker run --rm -ti rocker/r2u
## starts in shell, not R
apt -y update
apt -y upgrade
apt -y install git-all
apt -y install texlive-science texlive-latex-extra texlive-bibtex-extra texinfo
apt -y install pandoc pandoc-citeproc qpdf ghostscript fonts-inconsolata
tlmgr install inconsolata ## instead of 'apt -y install texlive-fonts-extra'
updmap-user
install2.r devtools
git clone https://github.com/glmmTMB/glmmTMB.git
r - "devtools::install_local('glmmTMB', dependencies = TRUE, upgrade = 'always')"
```

```
export BRANCH=dispRE
git clone https://github.com/glmmTMB/glmmTMB.git
cd glmmTMB
git checkout $BRANCH
git pull
cd ..
## could also do this (install all dependencies) by first installing
##  glmmTMB from CRAN ...
## more efficient to use pacman ... ?
cat <<EOF >my_install.R
pp <- tools::package_dependencies('glmmTMB', which='all')[[1]]
pp <- union(pp, c("reformulas"))  ## not in CRAN deps yet ...
sapply(pp, install.packages, character.only = TRUE)
EOF
r my_install.R
R CMD build glmmTMB/glmmTMB
export TARBALL=`ls -t *.tar.gz | head -1`
R CMD check --as-cran $TARBALL
```
