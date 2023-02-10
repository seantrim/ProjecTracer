#gnuplot: effect of the real to imaginary eigenvalue ratio

ratio=5.0

reset
set term qt 0
f(x)=(1.0-(x+1.0)**2.0)**0.5

m=-1/ratio
g(x)=m*x

set xlabel "Real"
set ylabel "Imaginary"

set size ratio 0.5
set xrange [-2:0]
set yrange [0:1]
set samples 1000
plot f(x), \
g(x)
