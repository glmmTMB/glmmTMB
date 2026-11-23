#!/bin/bash
## script to modify the maintainer email of current branch
## and build the tarball
##
## run from head glmmTMB/ directory (above pkg directory),
## i.e. . misc/autoreplace.sh
##
MY_EMAIL=${1:-bbolker@gmail.com}
MAINTAINER_EMAIL=${2:-mollieebrooks@gmail.com}
version=`grep 'Version' glmmTMB/DESCRIPTION | sed -e 's/Version: //'`
echo "glmmTMB version $version"
tarball="glmmTMB_${version}.tar.gz"
## substitute e-mail in DESCRIPTION file and build tarball
sed -i -e "s/$MAINTAINER_EMAIL/$MY_EMAIL/" glmmTMB/DESCRIPTION
R CMD build glmmTMB
tar zxvfO $tarball glmmTMB/DESCRIPTION | grep Maintainer
echo "tarball: $tarball"
## revert changes to DESCRIPTION file
git checkout -- glmmTMB/DESCRIPTION
