#gnuplot

psize=0.5    #point size in plots
nframes=1   #number of snapshot files
aspect=1.0   #aspect ratio
array=3      #field to be plotted: 3=T,4=psi,5=u,6=w,8=vis,11=C1,12=C2 etc.
term_name='qt'

reset
set term term_name 0 size 600,600
unset border
unset key
unset xtics
unset ytics
set size ratio 1.0/aspect
do for [i=0:nframes] {
n=10000+i
fname='../T'.n
plot fname u 1:2:(column(array)) w p pt 5 ps psize palette
pause mouse key
}
