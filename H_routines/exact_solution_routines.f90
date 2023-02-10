real*8 function vRMS(lambda,t)
implicit none

!!input
real*8 :: lambda !!aspect ratio
real*8 :: t !!time

!!internal variables
real*8, parameter :: pii=3.1415926535897932d0
real*8 :: f_func !!external function f(t) -- time dependence of the stream function
vRMS=pii*sqrt(lambda**2+1.d0)/2.d0/lambda*abs(f_func(t))
end function vRMS

subroutine compute_entrainment(tmin,tmax,nt,zR,lambda,k,zI,RaT,RaC,nx,nz,fname)
implicit none

!!input
real*8 :: zR !!reference z value
real*8 :: tmin,tmax !!lower and upper limit for time
integer*4 :: nt !!number of data points in the time series
integer*4 :: nx,nz !!# of mesh points in the x and z directions
real*8 :: xmin,xmax,zmin,zmax !!limits of integration
real*8 :: lambda !!aspect ratio
real*8 :: zI,k
real*8 :: RaT,RaC !!Rayleigh numbers
character*256 :: fname !!file name

!!output is in file named using the input variable fname

!!internal variables
real*8 :: C_array(1:nx,1:nz)
real*8 :: t,dt
real*8 :: integral,entrainment
integer*4 :: kt

open(unit=666,file=fname)

dt=(tmax-tmin)/real(nt-1,8)

do kt=1,nt
 t=max(dt*real(kt-1,8),0.d0) !!negative t values are not allowed
 call compute_array('C',t,lambda,k,zI,RaT,RaC,nx,nz,C_array)
 call volume_integral(0.d0,lambda,zR,1.d0,lambda,nx,nz,C_array,integral)
 entrainment=integral/lambda/zI
 write(666,'(2(g16.9))') t,entrainment
end do
close(unit=666)
end subroutine compute_entrainment

subroutine volume_integral(xmin,xmax,zmin,zmax,lambda,nx,nz,array,integral)
!!computes volume integral of values stored in array over x-z space
!!Limits of integration given by xmin, xmax, zmin, and zmax.
implicit none

!!input
integer*4 :: nx,nz !!# of mesh points in the x and z directions
real*8 :: xmin,xmax,zmin,zmax !!limits of integration
real*8 :: lambda !!aspect ratio
real*8 :: array(1:nx,1:nz) !!array containing the data

!!output
real*8 :: integral

!!internal variables
integer*4 :: iint,kint
integer*4 :: kint_min,kint_max
real*8 :: dz
real*8 :: integral_temp(1:nz)

dz=1.d0/real(nz-1,8)

!integral=0.d0
kint_min=nint(zmin/dz,4)+1
kint_max=nint(zmax/dz,4)+1

do kint=kint_min,kint_max
 call integral_1D(xmin,xmax,lambda,nx,array(1:nx,kint),integral_temp(kint)) !!create array of integrals over x
end do

call integral_1D(zmin,zmax,1.d0,nz,integral_temp,integral)
end subroutine volume_integral

subroutine integral_1D(xmin,xmax,length,nx,array_1D,integral)
!!compute integral of values stored in array_1D over a single spatial dimension over the interval [0,xmax]
implicit none

!!input
integer*4 :: nx !!# of mesh points
real*8 :: length !!length of spatial axis
real*8 :: xmin,xmax !!lower and upper limits of integration
real*8 :: array_1D(1:nx)

!!output
real*8 :: integral

!!internal variables
integer*4 :: iint,iint_min,iint_max
real*8 :: dx

dx=length/real(nx-1,8)

integral=0.d0
iint_min=nint(xmin/dx,4)+1
iint_max=nint(xmax/dx,4)+1
do iint=iint_min,iint_max
 integral=integral+dx*(array_1D(iint)+array_1D(iint+1))/2.d0
end do
end subroutine integral_1D

subroutine create_datafile(lambda,nx,nz,array,fname)
!!create datafile for array in column format: x z array
implicit none

!!input
character*256 :: fname !!output file name
integer*4 :: nx,nz !!# of mesh points in the x and z directions
real*8 :: lambda !!aspect ratio
real*8 :: array(1:nx,1:nz) !!array containing the data

!!output is stored in file named fname

!!internal variables
integer*4 :: iint,kint
real*8 :: dx,dz
real*8 :: x,z

open(unit=666,file=fname)

dx=lambda/real(nx-1,8); dz=1.d0/real(nz-1,8)

do kint=1,nz
 z=max(min(dz*real(kint-1,8),1.d0),0.d0)
 do iint=1,nx
  x=max(min(dx*real(iint-1,8),lambda),0.d0)
  write(666,'(3(g17.9))') x,z,array(iint,kint)
 end do
end do

close(666)
end subroutine create_datafile

subroutine compute_array(option,t,lambda,k,zI,RaT,RaC,nx,nz,array)
!!compute array of composition values for a time t
implicit none

!!input
character*1 :: option !!specify function to be computed -- valid options are: "C", "T", "D", and "H".
real*8 :: t !!time
real*8 :: lambda,k,zI !!aspect ratio, sharpness parameter, and z value of interface at t=0
real*8 :: RaT,RaC !!Rayleigh numbers
integer*4 :: nx,nz !!# of mesh points in the x and z directions

!!output
real*8 :: array(1:nx,1:nz)

!!internal variables
integer*4, parameter :: n_ghost=10
integer*4 :: iint,kint
real*8 :: dx,dz
real*8 :: x,z
real*8 :: C,D,T_func,H_func !!external functions
complex*16 :: H

dx=lambda/real(nx-1,8); dz=1.d0/real(nz-1,8)

if (option.eq.'C') then
 do kint=1,nz
  z=max(min(dz*real(kint-1,8),1.d0),0.d0)
  do iint=1,nx
   x=max(min(dx*real(iint-1,8),lambda),0.d0)
   array(iint,kint)=C(x,z,t,lambda,k,zI)
  end do
 end do
elseif (option.eq.'T') then
 do kint=1,nz
  z=max(min(dz*real(kint-1,8),1.d0),0.d0)
  do iint=1,nx
   x=max(min(dx*real(iint-1,8),lambda),0.d0)
   array(iint,kint)=T_func(x,z,t,lambda,k,zI,RaT,RaC)
  end do
 end do
elseif (option.eq.'D') then
 do kint=1,nz
  z=max(min(dz*real(kint-1,8),1.d0),0.d0)
  do iint=1,nx
   x=max(min(dx*real(iint-1,8),lambda),0.d0)
   array(iint,kint)=D(x,z,lambda)
  end do
 end do
elseif (option.eq.'H') then
 do kint=1,nz
  z=max(min(dz*real(kint-1,8),1.d0),0.d0)
  do iint=1,nx
   x=max(min(dx*real(iint-1,8),lambda),0.d0)
   call compute_H_func(x,z,t,lambda,k,zI,RaT,RaC,H_func)
   array(iint,kint)=H_func
  end do
 end do
else
 write(*,*) "Error in compute_array: unrecognized option."
end if
end subroutine compute_array

real*8 function T_func(x,z,t,lambda,k,zI,RaT,RaC)
!!temperature function
!!assumes x belongs to (-lambda/2,3/2*lambda) and z belongs to (-1,2)
implicit none

!!inputs
real*8 :: x,z,t  !!position and time
real*8 :: lambda !!aspect ratio
real*8 :: k,zI   !!sharpness parameter and initial layer height
real*8 :: RaT,RaC !!Rayleigh numbers

!!internal variables
real*8, parameter :: pii=3.1415926535897932d0
real*8 :: z0,C,f_func

T_func=(-(pii**3*(lambda**2+1.d0)**2/lambda**3)*cos(pii*x/lambda)*sin(pii*z)*f_func(t)+RaC*C(x,z,t,lambda,k,zI)+&
&(RaT-RaC)*(1.d0-z))/RaT

end function T_func

real*8 function C(x,z,t,lambda,k,zI)
!!compute composition
!!assumes x belongs to (-lambda/2,3/2*lambda) and z belongs to (-1,2)
implicit none

!!inputs
real*8 :: x,z,t  !!position and time
real*8 :: lambda !!aspect ratio
real*8 :: k,zI   !!sharpness parameter and initial layer height

!!internal variables
real*8 :: z0

call compute_z0(x,z,t,lambda,z0)
C=1.d0/(1.d0+exp(-2.d0*k*(zI-z0)))
end function C

subroutine compute_z0(x,z,t,lambda,z0)
!!Compute z0: the initial z value for a fluid parcel at position (x,z) at time t.
!!assumes x belongs to (-lambda/2,3/2*lambda) and z belongs to (-1,2)
implicit none

!!inputs
real*8 :: x,z,t  !!position and time
real*8 :: lambda !!aspect ratio

!!output
real*8 :: z0

!!internal variables
real*8, parameter :: pii=3.1415926535897932d0
real*8 :: Q,arg,eQ,bigZ0 !!variables for sidewalls
real*8 :: phi,m !!elliptic amplitude and parameter
real*8 :: D,S,f_integral,arccot !!external functions
real*8 :: x_,z_
complex*16 :: F,E !!elliptic integrals of the first and second kinds
complex*16 :: u   !!argument for the Jacobi elliptic functions
complex*16 :: sn,cn,dn !!Jacobi elliptic function values

if ((z.eq.0.d0).or.(z.eq.1.d0).or.((x.eq.(lambda/2.d0)).and.(z.eq.0.5d0))) then
 z0=z
elseif ((x.eq.0.d0).or.(x.eq.lambda)) then
 arg=pii*z
 Q=log(abs(1.d0/sin(arg)+1.d0/tan(arg)))-pii/lambda*S(x,lambda)*f_integral(t)
 eQ=exp(Q)
 bigZ0=(eQ-1.d0/eQ)/2.d0
 if (bigZ0.ge.0.d0) then
  z0=arccot(bigZ0)/pii
 else
  z0=1.d0+arccot(bigZ0)/pii
 end if
else
 if (x.lt.0.d0) then !!beyond left sidewall
  x_=-x
 elseif (x.gt.lambda) then !!beyond right sidewall
  x_=2.d0*lambda-x
 else !!domain interior
  x_=x
 end if
 if (z.lt.0.d0) then !!beyond bottom
  z_=-z
 elseif (z.gt.1.d0) then !!beyond top
  z_=2.d0-z  
 else !!domain interior
  z_=z
 end if
 phi=pii*z_; m=1.d0/D(x_,z_,lambda)**2.d0
 call incomplete_elliptic_integrals(phi,m,F,E)
 u=F
 if (t.ne.0.d0) u%IM=u%IM-S(x_,lambda)*(pii**2.d0*D(x_,z_,lambda)/lambda)*f_integral(t)
 call Jacobi_elliptic_functions(u,m,sn,cn,dn)
 if (abs(cn%RE).gt.1.d0) then
  write(*,*) "compute_z0: large cn magnitude detected -- cn=",cn%RE
  stop
 end if
 z0=acos(cn%RE)/pii


! phi=pii*z; m=1.d0/D(x,z,lambda)**2.d0
! call incomplete_elliptic_integrals(phi,m,F,E)
! u=F
! if (t.ne.0.d0) u%IM=u%IM-S(x,lambda)*(pii**2.d0*D(x,z,lambda)/lambda)*f_integral(t)
! call Jacobi_elliptic_functions(u,m,sn,cn,dn)
! if (abs(cn%RE).gt.1.d0) then
!  write(*,*) "compute_z0: large cn magnitude detected -- cn=",cn%RE
!  stop
! end if
! z0=acos(cn%RE)/pii
end if
end subroutine compute_z0

real*8 function D(x,z,lambda)
!!characteristic orbital value
implicit none

!!inputs
real*8 :: x,z,lambda

!!internal variables
real*8, parameter :: pii=3.1415926535897932d0

D=abs(sin(pii*x/lambda)*sin(pii*z))
end function D

real*8 function S(x,lambda)
!!transformed unit step function
!!assumes -lambda/2 <= x <= 3/2*lambda
implicit none
real*8 :: x,lambda
if (x.le.lambda/2.d0) then
 S=1.d0
else
 S=-1.d0
end if
end function S

real*8 function arccot(x)
!!inverse cotangent function -- range is (-pi/2,pi/2]-{0} for computational convenience
!!assuming x is real
implicit none

!!input
real*8 :: x

!!internal variables
real*8, parameter :: pii=3.1415926535897932d0

if (x.lt.0.d0) then
 arccot=-pii/2.d0-atan(x)
else
 arccot=pii/2.d0-atan(x)
end if
end function arccot

