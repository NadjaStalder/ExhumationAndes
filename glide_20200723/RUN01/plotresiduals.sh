ps3="plots/residuals.ps"
rm $ps3
#gmtset TICK_LENGTH = -0.2c
#gmtset ANOT_FONT_SIZE=12
#gmtset LABEL_FONT_SIZE=14

psbasemap -R0/12/-0.8/0.8 -JX4c/7c -B4S:"closure depth (km)":/.2W:"ln(zc) - g(m)":WeSn -V -P -K >> $ps3
awk '{print $2,$3}' resi_log.txt | psxy -J -R -Gblue -Sc0.25c -Wthick,black -O -V -P -K >> $ps3

psbasemap -R0/100/0/400 -JX6c/6c -Y9 -B20S:"iteration":/100W:"2S":WeSn -O -V -P -K >> $ps3
awk '{print $1,$2}' residualsLOG.txt | psxy -J -R -Wthick,black -O -P -V -K >> $ps3

awk '{print $3}' resi_log.txt | pshistogram -Z1 -R-1.2/1.2/0/50 -W0.1 -Gblue -L -P -Jx2.5c/0.1c -Y10 -B0.4S:"ln(zc) - g(m)":/5W:"Relative Frequency (%)":WeSn -O -V >> $ps3

open -a Preview plots/residuals.ps
