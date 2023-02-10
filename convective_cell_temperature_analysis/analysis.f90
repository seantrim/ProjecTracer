program analysis_of_convective_cell_temperature
integer*4 :: i,k
real*8 :: RaT,aratio,x,z,dx,dz
real*8 :: lateral_avg

!!!!!inputs
aratio=2.d0
RaT=3.d8
nx=2000
nz=1000
!!!!!end inputs

open(unit=101,file="geotherm.dat")
dx=aratio/real(nx,8)
dz=1.d0/real(nz,8)
z=0.d0
do k=0,nz
 write(101,*) z,lateral_avg(z,nx,dx,RaT,aratio)
 z=z+dz
 if (k.eq.(nz-1)) z=1.d0
end do

end program analysis_of_convective_cell_temperature

real*8 function convective_cell(x,z,RaT,aratio) !!analytic expression from the van Keken et al. (1997) benchmark
implicit none
integer*4 :: i,k
real*8 :: Tu,Tl,Tr,Ts,Q,u0,v0,x,z,perturb,aratio,pii,RaT
if (z.eq.0.d0) then
 convective_cell=1.d0
elseif (z.eq.1.d0) then
 convective_cell=0.d0
else 
 pii=3.1415926535897932d0
 u0=(aratio**(7.d0/3.d0))*((RaT/2.d0/(pii**0.5d0))**(2.d0/3.d0))*&
 &(1.d0+aratio**4.d0)**(-2.d0/3.d0)
 v0=u0
 Q=2.d0*(aratio/pii/u0)**0.5d0
   Tu=0.5d0*derf((1-z)*((u0/x)**0.5d0)/2.d0)
   Tl=1.d0-0.5d0*derf(z*((u0/(aratio-x))**0.5d0)/2.d0)
   Tr=0.5d0+(Q/2.d0/pii**0.5d0)*((v0/(z+1.d0))**0.5d0)*&
   &dexp(-(v0*x**2.d0)/(4.d0*z+4.d0))
   Ts=0.5d0-(Q/2.d0/pii**0.5d0)*((v0/(2.d0-z))**0.5d0)*&
   &dexp(-(v0*(aratio-x)**2.d0)/(8.d0-4.d0*z))
   convective_cell=Tu+Tl+Tr+Ts-1.5d0
end if
end

real*8 function lateral_avg(z,nx,dx,RaT,aratio)
integer*4 :: i,nx
real*8 :: x,z,convective_cell,integral,x1,x2,dx,RaT,aratio
intergral=0.d0
x1=0.d0
do i=0,nx-1
 x2=x1+dx
 integral=integral+0.5d0*dx*(convective_cell(x1,z,RaT,aratio)+convective_cell(x2,z,RaT,aratio))
 x1=x2
end do
integral=integral/aratio
lateral_avg=integral
end
