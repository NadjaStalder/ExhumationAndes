#!/bin/sh
ps3="plots/misfits.ps"
gmtset MAP_TICK_LENGTH_PRIMARY = -0.2c
psbasemap -R0/80/0/80 -JX15c -B5S:"Predicted Age (Ma)":/5W:"Measured Age (Ma)":WeSn  -K > $ps3
psxy -R -J -O -Wthicker -K  >> $ps3 << END
0 0
80 80
END

awk ' $4==1 {print $3,$1,$2}' stuff/misfits.txt | psxy -J -R -Gorange -Ey/0.1 -Sc0.25c -Wthick,black -O  -K >> $ps3
awk ' $4==3 {print $3,$1,$2}' stuff/misfits.txt | psxy -J -R -Gred -Ey/0.1 -Sc0.25c -Wthick,black -O -K >> $ps3
awk ' $4==2 {print $3,$1,$2}' stuff/misfits.txt | psxy -J -R -Ggreen -Ey/0.1 -Sc0.25c -Wthick,black -O -K >> $ps3
awk ' $4==4 {print $3,$1,$2}' stuff/misfits.txt | psxy -J -R -Gblue -Ey/0.1 -Sc0.25c -Wthick,black -O >> $ps3

open -a Preview $ps3

ps4="plots/AERAFT.ps"
rm $ps4

psbasemap -R0/80/0/5000 -JX15c/10c -B10S:"Age (Ma)":/300W:"elevation":WeSn  -K > $ps4
awk ' $4==1 {print $1,$5,$2}' stuff/misfits.txt | psxy -J -R -Gorange -Ex/0.1 -Sc0.25c -Wthick,black -O  -K >> $ps4
awk ' $4==1 {print $3,$5}' stuff/misfits.txt | psxy -J -R -Gblue -Sd0.25c -Wthick,black -O  >> $ps4

ps5="plots/AERAHE.ps"
rm $ps5
psbasemap -R0/80/0/5000 -JX15c/10c -B10S:"Age (Ma)":/300W:"elevation":WeSn  -K > $ps5
awk ' $4==3 {print $1,$5,$2}' stuff/misfits.txt | psxy -J -R -Gyellow -Ex/0.1 -Sc0.25c -Wthick,black -O  -K >> $ps5
awk ' $4==3 {print $3,$5}' stuff/misfits.txt | psxy -J -R -Ggreen -Sd0.25c -Wthick,black -O  >> $ps5

open -a Preview $ps4
open -a Preview $ps5
