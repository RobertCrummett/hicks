# gmt begin tfa pdf
# gmt xyz2grd mag_grid_data.csv -Gtmi.grd+n -R367150.1875/397243.2188/4145707.75/4178161.25 -I60.3067/65.0371+e
# gmt grd2cpt tmi.grd -Cgeosoft.cpt -E -L51300/51700 -Z
# gmt grdimage tmi.grd -I -JX12c/12c -Ba
# gmt grdcontour tmi.grd -C10 -L-51400/51650 -t80
# gmt end show
# 
# gmt begin raw pdf
# gmt makecpt -Cgeosoft.cpt -T51300/51700
# gmt basemap -R367150.1875/397243.2188/4145707.75/4178161.25 -JX12c/12c -Ba
# gmt plot combined_data_rad_mag.csv -Sc0.05c -C -i0,1,5
# gmt end show

gmt begin tmirbf pdf
gmt xyz2grd tmi_rbf.csv -Gtmirbf.grd+n -R367150.1875/397243.2188/4145707.75/4178161.25 -I60.3067/65.0371+e
gmt grd2cpt tmirbf.grd -Cgeosoft.cpt -E -L51300/51700 -Z
bs_Y_min
gmt grdimage tmirbf.grd -I -JX12c/12c -Ba
gmt grdcontour tmirbf.grd -C10 -L-51300/51700 -t80
gmt end show
