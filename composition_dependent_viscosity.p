#gnuplot: test the analytical viscosity filter equation

reset
f(x)=0.5+0.5*tanh(2.5*(x-0.5))/tanh(1.25)
g(x)=x
set term qt 0
set xrange [0:1]
set xlabel "Original C Value"
set ylabel "Filtered C Value"
set key top left
plot g(x) title "Unfiltered", \
f(x) title "Filtered"
