#gnuplot: plot stability curves for RK and RKC methods

##############################INPUT PARAMETERS
alpha=0.5                         ###time step scale factor for first stage of RK2 and RKC2
tol=0.02                          ###tolerance for stability function =1
eps=1000.0                      ##damping factor
order=1
stages=35                              ##number of RKC stages
Imax=35.0   ##max value of imaginary axis
eig_ratio=4.0 ##expected value of eig_ratio
stage_steps=0 ##plot individual time steps within all RKC2 stages??
sig_test(u,v)=Pdamp(stages,u,v) ##select stability polynomial
##############################END INPUT PARAMETERS 

############################################################### FUNCTION DEFINITIONS
j={0,1} ###imaginary unit
#T(s,zz)=(s==0?1:(s==1?zz:2*zz*T(s-1,zz)-T(s-2,zz))) ##Chebyshev Polynomial of order s as a function of complex zz -- recursive version
T(s,zz)=abs(zz)>1.0 ? 0.5*(zz+(zz**2-1.0)**0.5)**s+0.5*(zz-(zz**2-1.0)**0.5)**s : cos(s*acos(zz)) ##Chebyshev polynomial -- explicit version
U(s,zz)=((zz+(zz**2-1.0)**0.5)**(s+1)-(zz-(zz**2-1.0)**0.5)**(s+1))/(2*(zz**2-1.0)**0.5)
#Tz(s,zz)=(s==0?0:(s==1?1:2*T(s-1,zz)+2*zz*Tz(s-1,zz)-Tz(s-2,zz))) ##First derivative of Chebyshev Polynomial of order s as a function of complex zz
Tz(s,zz)=s*U(s-1,zz) ##explicit version
Tzz(s,zz)=(s==0?0:(s==1?0:4*Tz(s-1,zz)+2*zz*Tzz(s-1,zz)-Tzz(s-2,zz))) ##Second derivative of Chebyshev Polynomial of order s as a function of complex zz
flh(u,v)=fu(u,v)+j*fv(u,v) ##lambda*h in the stability diagram

################## RKC STABILITY FUNCTIONS
omega0=1+eps/stages**2
omega1_f(eps,order)=(order==1?T(stages,omega0)/Tz(stages,omega0):Tz(stages,omega0)/Tzz(stages,omega0)) ##for second order RKC
omega1=(order==1?T(stages,omega0)/Tz(stages,omega0):Tz(stages,omega0)/Tzz(stages,omega0)) ##for second order RKC
fzz(s,u,v,order)=omega0+omega1*flh(u,v) ##complex z second order RKC
Tomega0=T(stages,omega0)
Tzomega0=Tz(stages,omega0)
if (order==2) {
Tzzomega0=Tzz(stages,omega0)
temp=(Tzzomega0/Tzomega0**2)
Bdamp(s,u,v)=1+temp*(T(s,fzz(s,u,v,2))-Tomega0)} ##damped second order stability polynomial
Pdamp(s,u,v)=T(s,fzz(s,u,v,1))/Tomega0
B_R=-2.0*omega0/omega1
Rmin=B_R  ##min value of real axis
################## END RKC STABILITY FUNCTIONS

############### RK STABILITY FUNCTIONS
sigRK1(u,v)=1+flh(u,v)
sigRK2(u,v)=(1+flh(u,v)+0.5*flh(u,v)**2)
sigRK3(u,v)=(1+flh(u,v)+0.5*flh(u,v)**2+(flh(u,v)**3)/6)
sigRK4(u,v)=(1+flh(u,v)+0.5*flh(u,v)**2+(flh(u,v)**3)/6+(flh(u,v)**4)/24)
############### END RK STABILITY FUNCTIONS
############################################################### END FUNCTION DEFINITIONS


####***************************************************************############################# Compute RKC parameters
if (order==2) {
b0=0            ####free parameter
b1(alpha,eps)=alpha/omega1_f(eps,2) ####free parameter
b(s,alpha,eps)=(s==0?b0:(s==1?b1(alpha,eps):Tzz(s,omega0)/Tz(s,omega0)**2))
a(s,alpha,eps)=1-b(s,alpha,eps)*T(s,omega0)

mewt1(alpha,eps)=b(1,alpha,eps)*omega1_f(eps,2)
mewt(s,alpha,eps)=(s==1?mewt1(alpha,eps):2*b(s,alpha,eps)*omega1_f(eps,2)/b(s-1,alpha,eps))

gt(s,alpha,eps)=(-a(s-1,alpha,eps)*mewt(s,alpha,eps))
set print "-"
print " "
print "alpha=",alpha
print "epsilon=",eps
print " "
print "s ","mewt              ","gt"
print 1,mewt(1,alpha,eps)
do for [s=2:stages] {
print s,mewt(s,alpha,eps),gt(s,alpha,eps)}
}
print "RKC1 B_R=",B_R
####***************************************************************############################# END Compute RKC parameters

b=Imax ##35.0
a=abs(B_R/2.0)
fe_x(v)=a*cos(v)+B_R/2.0
fe_y(v)=b*sin(v)
reset
set term qt 0 font "Times,14" enhanced size 800,500
set parametric
set key inside bottom center
set view map
#set isosamples 700,300
set isosamples 500,250
set xlabel "Real Axis"
set ylabel "Imaginary Axis"
set xrange [Rmin:0.1]
set yrange [0:Imax]
set urange [-2*pi:0]    ###must capture at least all 360 degrees range (for tracer stability region) and enough of the negative real axis (for RK stability curves)
set vrange [0:2*pi]
fu(u,v)=-Rmin*u/(2*pi) ##rescale parameter range for RK and RKC stability curves (for better curve sampling)
fv(u,v)=Imax*v/(2*pi)
set arrow from 0.0,0.0 to B_R,(-B_R/eig_ratio) nohead
splot fu(u,v),fv(u,v),(((abs(sig_test(u,v))<1.0+tol)&(abs(sig_test(u,v))>1.0-tol))?1:1/0) w p pt 5 ps 0.5 lc 1 title "actual", \
fe_x(v),fe_y(v),1.0 w p pt 5 ps 0.5 lc 2 title "ellipse"


if (stage_steps==1) {
do for [s=2:stages] {
reset
set term qt sprintf("%1.2f",s) font "Times,14" enhanced size 1300,500
unset key
set view map
set xlabel "alpha"
set ylabel "epsilon"
set xrange [0.25:1.0] ##alpha
set yrange [1:30]     ##epsilon
set isosamples 100
set multiplot layout 1,2 rowsfirst
set label 1 "" center at first alpha,eps point pt 6 ps 2 front tc rgb "black"
set title "mewt"
splot mewt(s,x,y) w p pt 5 ps 0.5 palette
set title "gt"
splot gt(s,x,y) w p pt 5 ps 0.5 palette 
unset multiplot
}
}
