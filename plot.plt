reset
set xlabel "Number of outer iterations"
set ylabel "Relative residual"
set logscale y
set grid
set multiplot
set key spacing 0 height 1
plot "reshis.dat" title "BA-NRSOR" w l
set nomultiplot
pause -1
