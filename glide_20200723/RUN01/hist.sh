#!/bin/sh
ps3="plots/hist.ps"
rm $ps3
awk ' {print ($3-$1)}' stuff/misfits.txt | pshistogram -R-50/50/0/25 -Z1 -W1 -Jx.2c/0.2c -Gblue -B10S:"ages":/5W:"Relative Frequency (%)":WeSn -F > $ps3
