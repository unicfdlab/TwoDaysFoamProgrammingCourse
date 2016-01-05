#!/bin/sh

rm -rf Build
mkdir Build
cd Build
cmake ../
make
./main


#
# END_OF_FILE
#

