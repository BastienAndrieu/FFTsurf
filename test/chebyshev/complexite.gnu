set terminal png size  16 cm, 12cm
set output 'complexite.png'

set logscale x
set format x "10^{%T}"
set logscale y
set format y "10^{%T}"
set grid xtics ytics 
set key on box lw 1.0 opaque reverse vertical Left inside left top font ",8"
set xrange [1e0:1e3]
set yrange [1e-7:1e-2]
set xlabel  'N'
set ylabel  'time [s]'

plot 'complexite.dat' using 1:2  title 'diff matrix' with linespoints ls 1 ,\
                   '' using 1:3  title 'matmul1' with linespoints ls 2 ,\
                   '' using 1:4  title 'chebdiff' with linespoints ls 3 ,\
                   '' using 1:5  title 'ifcht' with linespoints ls 4 ,\
1.e2 <= x && x <= 5.e2 ? x*x*1e-8 : 1/0 title "N^2" with lines lc rgb '#000000' lt 4 lw 2.0 dt "-" ,\
1.e2 <= x && x <= 5.e2 ? x*4e-9 : 1/0 title "N" with lines lc rgb '#000000' lt 4 lw 2.0 dt "."
#1.e0 <= x && x <= 1.e3 ? x*log(x)*7e-10 : 1/0 title "N logN" with lines lc rgb '#000000' lt 4 lw 2.0 dt "."#, \