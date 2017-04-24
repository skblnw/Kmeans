#!/bin/bash

TEMP=input.vmd
SITE=SITE_cg8_lamb0.0001.dat
OUTPUT_VMD=out
sed -n '1,/mol delrep 0 top/p' $TEMP > $OUTPUT_VMD.vmd
for nn in {0..7}; do
    list=`grep "Site_$nn" $SITE | sed 's/Site_'$nn'= //g'`
    sed -e "s:SEL:$list:g" \
        -e "s:CC:$nn:g" \
        -e "s:NN:$nn:g" \
    template_selection.txt > SELECT
    cat SELECT >> $OUTPUT_VMD.vmd
    echo "" >> $OUTPUT_VMD.vmd
done
sed -n '/mol rename top initial-noW.psf/,$p' $TEMP >> $OUTPUT_VMD.vmd

rm -f SELECT
