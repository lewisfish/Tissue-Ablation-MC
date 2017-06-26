set terminal pngcairo
set output 'speedup-heat.png'
set xrange [0.7:8.2]
set yrange [0.7:6.2]
unset key
set xlabel '# of processes'
set ylabel 'Speed up'
set title 'Speed-up of 3D finite heat transport'
plot 'speedup.dat' u 1:6 w lp

set output 'eff-heat.png'
set xrange [0.7:8.1]
set yrange [:1.05]
set ylabel 'Efficency'
set title 'Efficency of 3D finite heat transport'

plot 'speedup.dat' u 1:7 w lp


set output 'time-heat.png'
set xrange [0.7:8.1]
set yrange [0:90]

unset key
set ylabel 'Time/s'
set title 'Timings of 3D finite heat transport'

plot 'speedup.dat' u 1:5 w lp