#!/bin/bash

#
# remove old
#

rm -rf ./main
rm -rf ./configure
rm -rf ./autom*
rm -rf ./aclocal.m4
rm -rf ./config.status
rm -rf ./Makefile
rm -rf ./Makefile.in
rm -rf *.o
rm -rf ./.deps
rm -rf ./config.log
rm -rf ./depcomp
rm -rf ./install-sh ./missing

#
# do build
#

aclocal

automake --add-missing

autoconf

./configure

make

