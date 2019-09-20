set terminal png size  16 cm, 18cm
set grid xtics ytics
set xlabel  '#test'

set output 'obb2.png'
set multiplot           # engage multiplot mode

#set output 'obb2_complexite.png'
set size 1.0, 0.4
set ylabel  'time [s]'
plot 'obb2.dat' using 0:1 title 'Chebyshev' with linespoints ls 1 ,\
             '' using 0:2 title 'Bernstein' with linespoints ls 2


#set output 'obb2_volume.png'
set size 1.0, 0.4
set origin 0.0, 0.5
set ylabel  'volume'
plot 'obb2.dat' using 0:3 title 'Chebyshev' with linespoints ls 1 ,\
             '' using 0:4 title 'Bernstein' with linespoints ls 2


unset multiplot
set terminal png size  24 cm, 12cm
set output 'obb2_ratio.png'
set size 1, 1
set origin 0, 0
#set logscale y
set yrange [0:2]
set ylabel  'Bernstein/Chebyshev'
set key outside

stats 'obb2.dat' using ($2)/($1):($4)/($3) nooutput

plot 'obb2.dat' using ($0+1):($2)/($1) title 'time' with linespoints ls 1 ,\
             '' using ($0+1):($4)/($3) title 'volume' with linespoints ls 2, \
	     '' using ($0+1):(STATS_median_x) w l lt 3 t 'med time', \
	     '' using ($0+1):(STATS_median_y) w l lt 4 t 'med volume'
#	     '' using ($0+1):(STATS_mean_x) w l lt 1 t 'avg time', \
#	     '' using ($0+1):(STATS_mean_y) w l lt 2 t 'avg volume', \