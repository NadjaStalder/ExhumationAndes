#!/bin/sh

region=-R-72/-62/-36/-18.2
proj="M3c"
path="../data"

for i in {1..37}

do

echo $i

tt="plots/ttpath_$i.ps"

psbasemap -R0/80/0/1.8  -JX8c -B10S:"Time (Ma)":/0.2W:"Exhumation rate (km/Myr)":WeSn -K > $tt
psxy stuff/edot2_$i -R -J -Wthicker -O -V -K -B >> $tt
psxy stuff/unc2_$i -R0/80/0./1 -JX8c/6c -Wthicker -B10S:"Time (Ma)":/0.2W:"Normalised Variance":WeSn -Y10c -O -V -K >> $tt

psbasemap $region -J$proj -B1S:"X":/1W:"Y":WeSn -X10c -Y-10c -O -K  >> $tt
makecpt -Cgray -T-1/1/0.00001 > apt1.cpt
grdimage $path/ill.grd -J -R -Capt1.cpt -O -V -K >> $tt
psxy stuff/cont.xy $region -J$proj -Sc0.5c -Wthicker -O -V -B >> $tt

done
#open -a Preview $tt
