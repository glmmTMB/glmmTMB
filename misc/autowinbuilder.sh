#!/bin/bash
## script to modify the maintainer email of current branch,
## build the tarball, and upload to win-builder
## run from head directory
MY_EMAIL=bbolker@gmail.com
MAINTAINER_EMAIL=mollieebrooks@gmail.com
version=`grep 'Version' glmmTMB/DESCRIPTION | sed -e 's/Version: //'`
echo "glmmTMB version $version"
sed -i -e "s/$MAINTAINER_EMAIL/$MY_EMAIL/" glmmTMB/DESCRIPTION
R CMD build glmmTMB
tar zxvfO glmmTMB_1.1.2.tar.gz glmmTMB/DESCRIPTION | grep Maintainer
tarball="glmmTMB_${version}.tar.gz"
echo "tarball: $tarball"
## https://serverfault.com/questions/279176/ftp-uploading-in-bash-script
HOST=win-builder.r-project.org
USER=ftp
PASS=$MY_EMAIL
ftp -inv $HOST << EOT

user $USER $PASS
binary
passive
cd R-devel
put $tarball
cd ../R-release
put $tarball
bye
EOT
git checkout -- glmmTMB/DESCRIPTION
Rscript -e "foghorn::winbuilder_queue('glmmTMB')"

## NOTES
## HEAD~10 (fbb7deb "Merge pull request #712 from glmmTMB/rr_fixes": version 1.1.0) fails due to tweedie issues
## HEAD~5  (1215e38 Merge pull request #707 from glmmTMB/vignette_fix: version 1.1.2)
