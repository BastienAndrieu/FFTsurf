 set terminal png size    16.0000000     cm,   12.0000000     cm
 set output '2var.png'
set logscale x
set format x "10^{%T}"
set logscale y
set format y "10^{%T}"
set grid xtics ytics 
set key on box lw 1.0 opaque reverse vertical Left inside left top
set xrange [ 0.100000000000000E+01: 0.100000000000000E+04]
set yrange [ 0.100000001168610E-06: 0.100000001490116E+00]
set xlabel  'N'
set ylabel  'time [s]'
 plot '2var.dat' using 1:2  title 'ifcht2' with linespoints ls 1 ,\
 '' using 1:3  title 'globalw' with linespoints ls 2 ,\
