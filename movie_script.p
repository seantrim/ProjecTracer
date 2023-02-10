#gnuplot
reset
set term gif animate delay 10 size 640,600
comp=2                        ## number of enriched compositions
psize=0.5                     ## point size in plots
aspect=2.0                    ## model aspect ratio
nframes=system("ls -1 movie_files | wc -l")-2
set output "movie.gif"
set cbrange [0:1]
unset key
unset xtics
unset ytics
unset border
set size ratio 1.0/aspect
do for [i=0:nframes] {
n=10000+i
fname='f'.n
if (comp==0) {
plot "<paste movie_files/movie_grid movie_files/".fname u 1:2:3 w p pt 5 ps psize palette}
if (comp==1) {
set multiplot layout 2,1 rowsfirst
plot "<paste movie_files/movie_grid movie_files/".fname u 1:2:3 w p pt 5 ps psize palette
plot "<paste movie_files/movie_grid movie_files/".fname u 1:2:4 w p pt 5 ps psize palette
}
if (comp==2) {
set multiplot layout 3,1 rowsfirst
plot "<paste movie_files/movie_grid movie_files/".fname u 1:2:3 w p pt 5 ps psize palette
plot "<paste movie_files/movie_grid movie_files/".fname u 1:2:4 w p pt 5 ps psize palette
plot "<paste movie_files/movie_grid movie_files/".fname u 1:2:5 w p pt 5 ps psize palette
}
}
