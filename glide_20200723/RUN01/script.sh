#!/bin/sh

#rm plots/*

#script to make figures and movie of output maps

#Setting some initial things
region=-R-72/-62/-36.5/-18
res="0.01"
proj="M5c"
#when changing proj change pstext too
ps="tt_paths.ps"
colr1="glacier"
colr2="seis_Nadja"
#path is the
path="../data"
makecpt -Cpolar -T0/0.01/0.0001 -Mgray > apt.cpt
makecpt -C$colr1 -Dwhite/red -T0/2.0/0.001 > apt2.cpt
#makecpt -Cviridis -Dwhite/red -T0/1.8/0.001 > apt2.cpt
makecpt -C$colr2 -Dwhite/grey -T0./1.0/0.001 > apt3.cpt


#Light grey version of Matthew
#makecpt -Cpolar -T0/0.01/0.0001 -Mgray > apt.cpt
#makecpt -C$colr1 -M -T0./1.5/0.001 > apt2.cpt
#makecpt -Cseis -T0./1./0.001 -Mgray -V > apt3.cpt

gmtset FONT_TITLE=10p
gmtset FONT_LABEL=12p
gmtset PS_MEDIA=A4
gmtset MAP_TICK_LENGTH = 0.15c
# gmtset COLOR_BACKGROUND=gray

#blockmedian $path/DEM20190108.xyz $region -I$res -Vl > $path/cut_block.xyz
#surface $path/cut_block.xyz -G$path/cut.grd $region -I$res -V -T0.27
#grdgradient $path/cut.grd -A0 -Nt0.75 -G$path/ill.grd -V

#The loop is for every file in maps
for i in *

do

if [ ${i: -4} == ".txt" ]
then
  echo ""$i" file ends in .txt"
  continue
fi


ps2=plots/"$i.ps"

#Interpolate surfaces
#awk '{print $1,$2,$3}' $i | surface -G$i.grd $region -I$res -T0.4 -V -S0.1
#awk '{print $1,$2,$3}' unc/$i | surface  -Gunc/x.grd $region -I$res -T0.4 -V -S0.1
#awk '{print $1,$2,$4}' unc/$i | surface  -Gunc/$i.grd $region -I$res -T0.4 -V -S0.1

awk '{print $1,$2,$3}' $i | nearneighbor -G$i.grd $region -I$res -N4/0 -S15k -E-9999 -V
awk '{print $1,$2,$3}' unc/$i | nearneighbor -Gunc/$i.grd $region -I$res -N4/0 -S15k -E1 -V
#awk '{print $1,$2,$4}' unc/$i | nearneighbor -Gunc/$i.grd $region -I$res -N4/0 -S0.1 -E1 -V

#Mask out regions, colour them grey, where the uncertainty is greater than....
grdclip unc/$i.grd -Gclip.grd -Sa.999/-100000
grdclip clip.grd -Gnew_clip.grd -Sa-1/10000.
grdinfo new_clip.grd
#gmtset COLOR_FOREGROUND=114/114/114
#gmtset COLOR_NAN=128/128/128

grdmath $i.grd new_clip.grd MIN = new.grd
grdinfo $i.grd
grdinfo new.grd
#grdclip unc/$i.grd -Gclip.grd -Sa.95/1.1
grdmath unc/$i.grd new_clip.grd MIN = unc/new.grd

grdinfo unc/new.grd

t1=`awk '{print $7}' $i.txt`
t2=`awk '{print $9}' $i.txt`

#create basemap,grdimage, sample locations and scale
psbasemap $region -J$proj -B3S:"X":/2W:"Y":WesN -K -P -Y6.5c > $ps2
grdimage new.grd -J$proj $region -Capt2.cpt -I$path/ill.grd -K -O >> $ps2
pscoast -J$proj -R -Di -I1 -S255/255/255 -W1p/0/0/0 -O -K >> $ps2
psxy ../data/Boundaries_WD_GMT_EAST.txt -R -J -O -P -K -Wthick -Sf0.25 -V >> $ps2
psxy ../data/Boundaries_WD_GMT_WEST.txt -R -J -O -P -K -Wthick -Sf0.25 -V >> $ps2
awk '$4<'$t1'&& $4>'$t2' {print $1,$2}' ../data/data20200424_18-36_noBlocks.txt | psxy -Gblack -J$proj -R -Wblack -Sc0.1 -O -K >> $ps2
#pstext ../data/labels.txt $region -J$proj -F+f8p,Arial-bold,white=0.5p,white -N -V -O -K >> $ps2
#pstext ../data/labels_cities.txt $region -J$proj -F+f7p,Arial-bold,white -N -V -O -K >> $ps2
#awk '{print $1, $2}' ../data/labels_cities.txt | psxy $region -J$proj -Sc0.05c -Gblack -V -O -K >> $ps2

psscale -D2.5/-0.5/5.3/0.5h -Baf0.1:"Exhumation rate [km/Ma]": -I -Capt2.cpt -O -K >> $ps2
pstext $i.txt $region -J$proj -Wblack -Gwhite -N -O -K -X3.31c -Y-10.34c >> $ps2

#Same as above but for the unc
grdimage unc/new.grd -J$proj -R -B3S:"X":/2E:"Y":wEsN -Capt3.cpt -X3.5c -Y10.34c -I$path/ill.grd -O -K >> $ps2
pscoast -J$proj -R -Di -I1 -S255/255/255 -W1p/0/0/0 -O -K >> $ps2
#psxy ../data/block_boundariesGMT.txt -R -J -O -P -K -Wthick -Sf0.25 -V  >> $ps2
psxy ../data/Boundaries_WD_GMT_EAST.txt -R -J -O -P -K -Wthick -Sf0.25 -V >> $ps2
psxy ../data/Boundaries_WD_GMT_WEST.txt -R -J -O -P -K -Wthick -Sf0.25 -V >> $ps2
#pstext ../data/labels.txt $region -J$proj -F+f8p,Arial-bold,white=0.5p,white -N -V -O -K >> $ps2
#pstext ../data/labels_cities.txt $region -J$proj -F+f7p,Arial-bold,white -N -V -O -K >> $ps2
#awk '{print $1, $2}' ../data/labels_cities.txt | psxy $region -J$proj -Sc0.1c -Gblack -V -O -K >> $ps2
#awk '{print $1, $2}' ../data/test_label.txt | psxy $region -J$proj -Sc0.2c -Gblack -V -O -K >> $ps2
awk '$4<'$t1'&& $4>'$t2' {print $1,$2}' ../data/data20200424_18-36_noBlocks.txt | psxy -Gblack -J$proj -R -Wblack -Sc0.1 -O -K >> $ps2

psscale -D2.5/-0.5/5.3/0.5h -Baf0.1:"Reduced variance": -I -Capt3.cpt -O -K >> $ps2
pstext $i.txt $region -J$proj -Wblack -Gwhite -N -O -X3.31c -Y-10.34c >> $ps2

rm $i.grd
rm new.grd
rm unc/$i.grd
rm unc/new.grd
rm clip.grd
done

makecpt -Crainbow -T300/600/0.5 > apt.cpt
makecpt -Cgray -T-1/1/0.00001 > apt1.cpt

tt="plots/ttpath.ps"
psbasemap -R0/50/0/2 -JX8c -B5S:"Time (Ma)":/0.25W:"Exhumation rate (km/Myr)":e -K > $tt
psxy stuff/edot2 -R -J -M -Wthicker -O -V -K -B >> $tt
psxy stuff/unc2 -R0/50/0./1 -JX8c/6c -M -Wthicker -B5S:"Time (Ma)":/0.1W:"Normalised Variance":e -Y10c -O -V -K >> $tt

grdimage $path/ill.grd -J$proj $region -B1SWne -Capt1.cpt -X10c -Y-10c -O -V -K >> $tt
psxy stuff/cont.xy $region -J$proj -Sc0.5c -Wthicker -M -O -V -B >> $tt

ps3="plots/misfits.ps"
gmtset TICK_LENGTH = -0.2c
psbasemap -R0/80/0/80 -JX15c -B5S:"Predicted Age (Ma)":/5W:"Measured Age (Ma)":nSeW -V -K > $ps3
psxy -R -J -O -Wthicker -K -V >> $ps3 << END
0 0
200 200
END

awk ' $4==1 {print $3,$1,$2}' stuff/misfits.txt | psxy -J -R -Gorange -Ey/0.1/thicker -Sc0.25c -Wthick,black -O -V -K >> $ps3
awk ' $4==3 {print $3,$1,$2}' stuff/misfits.txt | psxy -J -R -Gred -Ey/0.1/thicker -Sc0.25c -Wthick,black -O -V -K >> $ps3
awk ' $4==2 {print $3,$1,$2}' stuff/misfits.txt | psxy -J -R -Ggreen -Ey/0.1/thicker -Sc0.25c -Wthick,black -O -V -K >> $ps3
awk ' $4==4 {print $3,$1,$2}' stuff/misfits.txt | psxy -J -R -Gblue -Ey/0.1/thicker -Sc0.25c -Wthick,black -O -V -K >> $ps3


#cleaning up the folders
rm apt*
rm plots/apt*
rm plots/*.cpt
rm plots/*.grd.ps
rm plots/*.cpt.ps
rm plots/plots.ps
rm plots/*.txt.ps
rm plots/unc.ps
rm plots/visualize2.sh.ps
rm *.eps

#Converting the ps files to jpg for movie
ps2raster -Tj plots/*ps -P -V -A
#ps2raster -Tf plots/*ps -P -V -A
