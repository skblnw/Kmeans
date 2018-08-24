#!/bin/bash

TEMP=input.vmd
SITE=../SITE_cg5.dat
OUTPUT_VMD=out
sed -n '1,/mol delrep 0 top/p' $TEMP > $OUTPUT_VMD.vmd
cat template_vm_draw.txt > vm_draw.tcl
for nn in {0..4}; do
    list=`grep "Site_$nn" $SITE | sed 's/Site_'$nn'= //g'`
    sed -e "s:SEL:$list:g" \
        -e "s:CC:$nn:g" \
        -e "s:NN:$nn:g" \
    template_selection.txt > SELECT
    cat SELECT >> $OUTPUT_VMD.vmd
    echo "" >> $OUTPUT_VMD.vmd
    
    sed -e "s:SEL:$list:g" \
        -e "s:COLOR:$nn:g" \
    template_command_draw.txt >> vm_draw.tcl
    echo "" >> vm_draw.tcl
done
sed -n '/mol rename top initial-noW.psf/,$p' $TEMP >> $OUTPUT_VMD.vmd

rm -f SELECT
