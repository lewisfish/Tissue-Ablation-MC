set terminal pngcairo size 1920,1080 font ",20"
set output 'time takes for various grid sizes'
f(x) = a * exp(b*x)
a = 0.05; b = 0.05
fit f(x) 'runtimes.dat' via a, b
#set xrange[0:200]
set xlabel '# of grid points'
set ylabel 'time/s'
#set logscale y
set title "Time taken for N grid points"
p  a*exp(b*x) lw 5,'runtimes.dat' lw 5 w p