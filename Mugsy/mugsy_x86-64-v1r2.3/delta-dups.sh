#!/bin/sh
#Identify duplicated regions in a pairwise delta file from NUCmer

mugsypath=$MUGSY_INSTALL
mummerpath=$MUGSY_INSTALL/MUMmer3.20/
deltafile=$1;

#Run delta-filter -b for duplications that are detected using LIS
$mummerpath/delta-filter -b $deltafile > $deltafile.b
#Capture additional dup/repeat regions by looking for overlapping alignments
#Alignments that overlap by more than half their lengths are reports as dups
$mummerpath/delta-filter -m $deltafile > $deltafile.m
$mummerpath/delta-filter -v -u 50 $deltafile.m > $deltafile.u
#Dump union of two sets to maf format
$mummerpath/delta2maf $deltafile.b 2> /dev/null | $mugsypath/fixMAFnames.pl 
#Skip first line
$mummerpath/delta2maf $deltafile.u 2> /dev/null | $mugsypath/fixMAFnames.pl | tail -n +1 
rm $deltafile.b &
rm $deltafile.m &
rm $deltafile.u &

