#!/bin/bash

rm -rf day*.swf
PDF_DIR=../Pdf
PDF_FILES=`ls $PDF_DIR`
THIS_DIR=`pwd`
export PATH="PATH:$THIS_DIR"
for PDF in $PDF_FILES
do
    ./pdf2swf -t -B ./keyboard_viewer.swf $PDF_DIR/$PDF -o $PDF.swf
done
#
# END_OF_FILE
#
