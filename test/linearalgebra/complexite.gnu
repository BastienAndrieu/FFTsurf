set terminal png size 16 cm, 14 cm
set output "complexite.png"

set datafile separator ','


set logscale x
set logscale y
set format xy "10^{%L}"
set ytics 1.0e1

set xrange [ 1 : 1e3 ]
set yrange [ 1e-7 : 1e0 ]

set grid xtics mxtics ytics

set xlabel  'n'
set ylabel  'time [s]'

set key on box lw 1.0 opaque reverse vertical Left inside left top font ",9" height -0


plot 2e1 <= x && x <= 2e2 ? x*x*x*x*6e-9 : 1/0 title "n^4" with lines lc rgb '#000000' lt 4 lw 2.0 dt "-." ,\
     2e1 <= x && x <= 2e2 ? x*x*x*6e-9 : 1/0 title "n^3" with lines lc rgb '#000000' lt 4 lw 2.0 dt "-" ,\
     for [i=2:*] 'complexite.csv' using 1:i with linespoints ls i title columnhead 

#plot 2e1 <= x && x <= 2e2 ? x*x*x*6e-7 : 1/0 title "n^3" with lines lc rgb '#000000' lt 4 lw 2.0 dt "-." ,\
#     2e1 <= x && x <= 2e2 ? x*x*6e-7 : 1/0 title "n^2" with lines lc rgb '#000000' lt 4 lw 2.0 dt "-" ,\