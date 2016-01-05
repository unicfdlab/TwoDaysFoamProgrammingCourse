#!/bin/bash


SCRIPT_PATH=`pwd`

if [ $# -eq 2 ]
then
    SCRIPT_PATH=$2
    cd $1
fi

CFILES=`ls -a`

for CFILE in $CFILES
do
    if [ -d $CFILE ]
    then
	if [ $CFILE != "." -a  $CFILE != ".." ]
	then
	    if [ $CFILE == ".svn" ]
	    then
		echo "removing .svn"
		rm -rf .svn
	    else
		echo "entering $CFILE"
		$SCRIPT_PATH/rmSvn.sh $CFILE  $SCRIPT_PATH
	    fi
	fi
    fi
done

if [ $# -eq 2 ]
then
    cd ../
fi

#
#END-OF-FILE
#

