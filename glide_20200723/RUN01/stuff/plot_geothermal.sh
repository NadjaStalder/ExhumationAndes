#!/bin/sh

#rm plots/*

#script to make figures and movie of output maps

#Setting some initial things
region=-R-72/-62/-36.5/-18
res="0.01"
proj="M5c"
#when changing proj change pstext too
ps="tt_paths.ps"
#colr1="glacier"
colr1="rainbow_N"
colr2="rainbow_N"
#colr2="seis_new"
#path is the
path="../../data"
makecpt -Cpolar -T0/0.01/0.0001 -Mgray > apt.cpt
makecpt -C$colr2 -Dwhite/gray -T40/200/0.01 > apt2.cpt
makecpt -C$colr2 -Dwhite/gray -T15/75/0.01 > apt3.cpt


gmtset FONT_TITLE=10p
gmtset FONT_LABEL=12p
gmtset PS_MEDIA=A4
gmtset MAP_TICK_LENGTH = 0.15c
# gmtset COLOR_BACKGROUND=gray

ps2="heatflow_0020.ps"

awk '{print $1,$2,$3}' heat_flux.txt | nearneighbor -GHF.grd $region -I$res -N4/0 -S15k -E-9999 -V

# Calculate geothermal gradient from heat flow (divided by 2.3)
grdmath HF.grd 2.6 DIV = GT.grd
awk '{print $1,$2,$3}' ../unc/0020 | nearneighbor -G../unc/0020.grd $region -I$res -N4/0 -S15k -E1 -V

#Mask out regions, colour them grey, where the uncertainty is greater than....
grdclip ../unc/0020.grd -Gclip.grd -Sa.995/-100000
grdclip clip.grd -Gnew_clip.grd -Sa-1/10000.
grdinfo new_clip.grd

grdmath HF.grd new_clip.grd MIN = Heatflow.grd
grdmath GT.grd new_clip.grd MIN = Geothermal.grd
grdinfo HF.grd
grdinfo GT.grd
grdinfo Heatflow.grd
grdinfo Geothermal.grd
grdmath ../unc/0020.grd new_clip.grd MIN = ../unc/new.grd

grdinfo ../unc/new.grd

#Plot Heatflow - create basemap,grdimage, sample locations and scale
psbasemap $region -J$proj -B3S:"X":/2W:"Y":WesN -K -P -Y6.5c > $ps2
grdimage Heatflow.grd -J$proj $region -Capt2.cpt -I$path/ill.grd -K -O >> $ps2
pscoast -J$proj -R -Di -I1 -S255/255/255 -W1p/0/0/0 -O -K >> $ps2
psxy ../../data/Boundaries_WD_GMT_EAST.txt -R -J -O -P -K -Wthick -Sf0.25 -V >> $ps2
psxy ../../data/Boundaries_WD_GMT_WEST.txt -R -J -O -P -K -Wthick -Sf0.25 -V >> $ps2
#awk '{print $1,$2}' ../../data/data20200424_18-36_noBlocks.txt | psxy -Gblack -J$proj -R -Wblack -Sc0.1 -O -K >> $ps2
psscale -D2.5/-0.5/5.3/0.5h -Baf5:"Heat flow (mW/m/K)": -I -Capt2.cpt -O -K >> $ps2

#plot Geothermal gradient
grdimage Geothermal.grd -J$proj -R -B3S:"X":/2E:"Y":wEsN -Capt3.cpt -X6c -I$path/ill.grd -O -K >> $ps2
pscoast -J$proj -R -Di -I1 -S255/255/255 -W1p/0/0/0 -O -K >> $ps2
#awk '$4<'$t1'&& $4>'$t2' {print $1,$2}' ../data/data20200424_18-36_noBlocks.txt | psxy -Gblack -J$proj -R -Wblack -Sc0.1 -O -K >> $ps2
psxy ../../data/Boundaries_WD_GMT_EAST.txt -R -J -O -P -K -Wthick -Sf0.25 -V >> $ps2
psxy ../../data/Boundaries_WD_GMT_WEST.txt -R -J -O -P -K -Wthick -Sf0.25 -V >> $ps2
psscale -D2.5/-0.5/5.3/0.5h -Baf5:"Geothermal Gradient (C/km)": -I -Capt3.cpt -O -K >> $ps2

makecpt -Crainbow -T300/600/0.5 > apt.cpt
makecpt -Cgray -T-1/1/0.00001 > apt1.cpt


#cleaning up the folders
rm apt*
#rm *.cpt
rm *.grd.ps

#Converting the ps files to jpg for movie
ps2raster -Tj plots/*ps -P -V -A
