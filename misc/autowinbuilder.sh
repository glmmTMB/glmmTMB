#!/bin/bash
## script to check out a commit from Github, modify the maintainer email,
## build the tarball, and upload to win-builder
## run from head directory
MYEMAIL=bbolker@gmail.com
head="HEAD~$1"
echo "checking out $head"
git checkout $head
version=`grep 'Version' glmmTMB/DESCRIPTION | sed -e 's/Version: //'`
echo "glmmTMB version $version"
sed -i -e "s/molliebrooks@gmail.com/$MYEMAIL/" glmmTMB/DESCRIPTION
R CMD build glmmTMB
tarball="glmmTMB_${version}.tar.gz"
echo "tarball: $tarball"
## https://serverfault.com/questions/279176/ftp-uploading-in-bash-script
HOST=win-builder.r-project.org
USER=ftp
PASS=$MYEMAIL
ftp -inv $HOST << EOT

user $USER $PASS
binary
passive
cd R-devel
put $tarball
bye
EOT
git checkout -- glmmTMB/DESCRIPTION
git checkout master
Rscript -e "foghorn::winbuilder_queue('glmmTMB')"

## NOTES
## HEAD~10 (fbb7deb "Merge pull request #712 from glmmTMB/rr_fixes": version 1.1.0) fails due to tweedie issues
## HEAD~5  (1215e38 Merge pull request #707 from glmmTMB/vignette_fix: version 1.1.2)
