#!/bin/bash
## script to modify the maintainer email of current branch,
## build the tarball, and upload to win-builder
## run from head directory (glmmTMB, above glmmTMB package dir),
## e.g. . misc/autowinbuilder.sh
## arg 1: which platform to test (both, release, or devel) ?
whichrel=${1:-both}
echo $whichrel
MY_EMAIL=bbolker@gmail.com
MAINTAINER_EMAIL=mollieebrooks@gmail.com
. misc/autoreplace.sh $MY_EMAIL $MAINTAINER_EMAIL
## https://serverfault.com/questions/279176/ftp-uploading-in-bash-script
HOST=win-builder.r-project.org
USER=ftp
PASS=$MY_EMAIL
if [ $whichrel == "both" ] || [ $whichrel == "devel" ]; then
echo "uploading to win-builder/devel"
ftp -inv $HOST << EOT

user $USER $PASS
binary
cd R-devel
put $tarball
bye 
EOT
fi
if [ $whichrel == "both" ] || [ $whichrel == "release" ]; then
echo "uploading to win-builder/release"
## upload to R-release
ftp -inv $HOST << EOT

user $USER $PASS
binary
cd R-release
put $tarball
bye
EOT
fi
## check status
Rscript -e "foghorn::winbuilder_queue('glmmTMB')"

## NOTES
## HEAD~10 (fbb7deb "Merge pull request #712 from glmmTMB/rr_fixes": version 1.1.0) fails due to tweedie issues
## HEAD~5  (1215e38 Merge pull request #707 from glmmTMB/vignette_fix: version 1.1.2)
