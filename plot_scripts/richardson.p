#gnuplot: graphical analysis of the Richardson method (non-stationary iterative method used for multigrid smoother)
reset
######## Test
RK1(x,n)=(1.0+ht(n)*x) ##x=lambda
RK2(x,n)=(1.0+ht(n)*x+((ht(n)*x)**2.0)/2.0)
RK3(x,n)=(1.0+ht(n)*x+((ht(n)*x)**2.0)/2.0+((ht(n)*x)**3.0)/6.0)
RK4(x,n)=(1.0+ht(n)*x+((ht(n)*x)**2.0)/2.0+((ht(n)*x)**3.0)/6.0+((ht(n)*x)**4.0)/24.0)
flh(n,x)=ht(n)*x ##lambda*h in the stability diagram
fzz(s,n,x,order)=omega0(eps)+omega1(eps,order)*flh(n,x) ##complex z second order RKC
#T(s,zz)=(s==0?1:(s==1?zz:2*zz*T(s-1,zz)-T(s-2,zz))) ##Chebyshev Polynomial of order s as a function of complex zz -- recursive
#Tz(s,zz)=(s==0?0:(s==1?1:2*T(s-1,zz)+2*zz*Tz(s-1,zz)-Tz(s-2,zz))) ##First derivative of Chebyshev Polynomial of order s as a function of complex zz -- recursive
T(s,zz)=abs(zz)>1.0 ? 0.5*(zz+(zz**2-1.0)**0.5)**s+0.5*(zz-(zz**2-1.0)**0.5)**s : cos(s*acos(zz)) ##Chebyshev polynomial -- explicit version
U(s,zz)=((zz+(zz**2-1.0)**0.5)**(s+1)-(zz-(zz**2-1.0)**0.5)**(s+1))/(2*(zz**2-1.0)**0.5)
Tz(s,zz)=s*U(s-1,zz) ##explicit version
Tzz(s,zz)=(s==0?0:(s==1?0:4*Tz(s-1,zz)+2*zz*Tzz(s-1,zz)-Tzz(s-2,zz))) ##Second derivative of Chebyshev Polynomial of order s as a function of complex zz
###--------Damped RKC
eps=1000.0 ##damping factor
stages=35
order=1
omega0(eps)=1+eps/stages**2
omega1(eps,order)=(order==1?T(stages,omega0(eps))/Tz(stages,omega0(eps)):Tz(stages,omega0(eps))/Tzz(stages,omega0(eps))) ##for second order RKC
B_R=-2.0*omega0(eps)/omega1(eps,1) ##For damped RKC1
h2=0.95*(-B_R) ##courant # times -B_R
h1=0.9*(-B_R)
Nt=1
if (Nt>1) {
ht(n)=(n<=Nt?(h2-h1)*(n-1)/(1-Nt)+h2:1/0)} ##linear time step profile
Tomega0=T(stages,omega0(eps))
if (Nt==1) {ht(n)=h2}
Pdamp(s,x,n)=T(s,fzz(s,n,x,1))/Tomega0
###---------End damped RKC
test(x,n)=Pdamp(stages,x,n)   ###RKC1_2(x,n) ##select method to be tested
######## End Test

######## For Original Richardson Iterations
B_R_RK1=-2.0 ##stability limit for RK1
lb=-1.0
la=-0.5
N=70
######## End Original Richardson Iterations

TT(x)=abs(x)>1.0 ? 0.5*(x+(x**2-1.0)**0.5)**N+0.5*(x-(x**2-1.0)**0.5)**N : cos(N*acos(x)) ##Chebyshev polynomial

P(x)=TT((2*x-la-lb)/(la-lb))/TT((-la-lb)/(la-lb))   ##original Richardson method
h(x)=(x<=N?2.0/(-lb-la+(lb-la)*cos((2.0*x-1)*pi/2.0/N)):1/0)

P_test(x,n)=(n>1?test(x,n)*P_test(x,n-1):test(x,1))

##############################################################for graphics
set title sprintf("N=%1.1f,lb=%1.1f,la=%1.1f",N,lb,la)
set samples 200
set term qt 0
set xtics add (-2.0/4.0,-1.0/16.0)
set xrange [-1:-1.0/16.0]
set yrange [-0.5:0.5]
set ylabel "Amplitude"
set xlabel "Spectral Radius Fraction"
set arrow from -1.0/16.0,-1.0 to -1.0/16.0,1.0 nohead
set arrow from -2.0/4.0,-1.0 to -2.0/4.0,1.0 nohead
set arrow from -1.0,-1.0 to -1.0,1.0 nohead
plot P(x), \
P_test(x,Nt)

reset
unset arrow
set title sprintf("N=%1.1f,lb=%1.1f,la=%1.1f",N,lb,la)
set samples 200
set term qt 1
set xrange [-1.0/16.0:0]
set yrange [-0.5:0.5]
set ylabel "Amplitude"
set xlabel "Spectral Radius Fraction"
plot P(x), \
P_test(x,Nt)

reset
unset arrow
set title sprintf("N=%1.1f,lb=%1.1f,la=%1.1f",N,lb,la)
set samples 200
set term qt 2
set xrange [-0.005:0]
set yrange [-0.5:0.5]
set ylabel "Amplitude"
set xlabel "Spectral Radius Fraction"
plot P(x), \
P_test(x,Nt)

reset
unset arrow
set title sprintf("N=%1.1f,lb=%1.1f,la=%1.1f",N,lb,la)
set samples 200
set term qt 3
set xrange [-1e-9:0]
set yrange [-0.5:0.5]
set ylabel "Amplitude"
set xlabel "Spectral Radius Fraction"
plot P(x), \
P_test(x,Nt)
##############################################################end for graphics


##################for stats of functions
!rm trash.dat fit.log
set samples 1000
set table "trash.dat"
set xrange [-1:0]
plot P(x), \
P_test(x,Nt)
unset table

print " "
print "***********************Upper Frequency Range***********************"
f(x)=P_avg_upper
fit [-1.0:-1.0/16.0] f(x) "trash.dat" i 0 u 1:(abs($2)) via P_avg_upper  ##range desgined for fine grid smoothing for the 4th order biharmonic equation (need coarse grid to resolve restricted residual)
f(x)=P_test_avg_upper
fit [-1.0:-1.0/16.0] f(x) "trash.dat" i 1 u 1:(abs($2)) via P_test_avg_upper
print " "
print "***********************Lower Frequency Range***********************"
f(x)=P_avg_lower
fit [-1.0/16.0:0.0] f(x) "trash.dat" i 0 u 1:(abs($2)) via P_avg_lower  ##range desgined for fine grid smoothing for the 4th order biharmonic equation (need coarse grid to resolve restricted residual)
f(x)=P_test_avg_lower
fit [-1.0/16.0:0.0] f(x) "trash.dat" i 1 u 1:(abs($2)) via P_test_avg_lower
print " "
print "***********************Total Frequency Range***********************"
f(x)=P_avg_total
fit [-1.0:0.0] f(x) "trash.dat" i 0 u 1:(abs($2)) via P_avg_total  ##range desgined for fine grid smoothing for the 4th order biharmonic equation (need coarse grid to resolve restricted residual)
f(x)=P_test_avg_total
fit [-1.0:0.0] f(x) "trash.dat" i 1 u 1:(abs($2)) via P_test_avg_total
################end stats of functions

#set term qt 3
#set xrange [1:(N>=Nt?N:Nt)]
#unset arrow
#set yrange [0:1]
#set ylabel "Courant Factor"
#set xlabel "Time Step Number"
#plot h(x)/(-B_R_RK1), \
#ht(x)/(-B_R)

reset
set term qt 4
set samples 200
set xrange [-1:0]
set yrange [-1.0:1.0]
plot test(x,1) w l, \
test(x,Nt) w l

set print "-"
print " "
print "RKC1 B_R=",B_R
print "************************Summary**************************"
print "      P                  P_test"
print "total ",P_avg_total,P_test_avg_total
print "lower ",P_avg_lower,P_test_avg_lower
print "upper ",P_avg_upper,P_test_avg_upper
