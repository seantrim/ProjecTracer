!!!define special initial conditions here
integer*4 function composition(x,z)
 use basics
implicit none
real*8 :: x,z
if (initial_comp.eq.0) then !!two phase flow: flat dense bottom layer
 if (z.gt.dlayer) then
  composition=0
 else
  composition=1
 end if
elseif (initial_comp.eq.1) then !!two phase flow: flat dense bottom layer with cosine perturbation (van Keken et al., 1997)
 if (z.gt.(0.02d0*dcos(pii*x/aspect)+dlayer)) then
  composition=0
 else
  composition=1
 end if
else !!!!!!!!!!!!!!!!!!!WORK AREA
 if ((z.gt.(1.d0-0.052d0)).and.(x.le.0.15d0)) then !!continents
  composition=3
 elseif ((z.gt.(1.d0-0.017d0)).and.(x.gt.0.15d0)) then !!oceans
  composition=2
 elseif (z.le.0.138d0) then !!basal layer
  composition=1
 else
  composition=0
 end if
end if
end

real*8 function composition_continuous(x,z)
 use basics
 use arrays
implicit none
real*8 :: x,z
if (initial_comp.eq.0) then !!two phase flow: flat dense bottom layer with extremely sharp interface
 if (z.gt.dlayer) then
  composition_continuous=0.d0
 else
  composition_continuous=1.d0
 end if
elseif (initial_comp.eq.-1) then !!custom
 composition_continuous=1.d0/(1.d0+exp(-2.d0*klayer*(dlayer-z))) !!exact solution paper - two layer flow with specific interface gradient
end if
end function composition_continuous

real*8 function temperature(x,z)
 use basics
 use arrays
implicit none
integer*4 :: i,k
real*8 :: convective_cell
real*8 :: x,z,perturb,aratio,xx
real*8 :: f0
if (initial_temp.eq.0) then !!conductive with sine perturbation
 temperature=(1.d0-z)+0.1d0*dsin(2.d0*pii*x/aspect)
elseif (initial_temp.eq.1) then !!!!!!!convective cell from van Keken et al. (1997)
 aratio=aspect
 temperature=convective_cell(x,z,aratio)
elseif (initial_temp.eq.2) then !!conductive with sine*cos style perturbation
 temperature=(1.d0-z)+0.01d0*dcos(pii*x/aspect)*dsin(pii*z/aspect)
elseif (initial_temp.eq.3) then !!plume experiment
 perturb=0.5d0*(dsin(pii*x/aspect))**4.d0
 temperature=((1.d0-z)**3.d0)*(1.d0+perturb)
 if (temperature.gt.1.d0) temperature=1.d0
 if (temperature.lt.0.d0) temperature=0.d0
elseif (initial_temp.eq.4) then !!convective cells with an aspect ratio of 2
 aratio=2.d0
 if (x.le.2.d0) then
  xx=x-0.d0
  temperature=convective_cell(xx,z,aratio)
 elseif (x.le.4.d0) then
  xx=2.d0-(x-2.d0)
  temperature=convective_cell(xx,z,aratio)
 elseif (x.le.6.d0) then
  xx=x-4.d0
  temperature=convective_cell(xx,z,aratio)
 elseif (x.le.8.d0) then
  xx=2.d0-(x-6.d0)
  temperature=convective_cell(xx,z,aratio)
 elseif (x.le.10.d0) then
  xx=x-8.d0
  temperature=convective_cell(xx,z,aratio)
 end if
else!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!WORK AREA
! if ((z.gt.0.9d0)) then
!  temperature=0.d0       !!!cold stagnant lid
!!  temperature=-8.d0*z+8.d0 !!!linear increase with depth
!! elseif ((x.lt.0.1d0).and.(z.gt.0.8)) then
!!  temperature=0.0d0
! elseif ((x.gt.(aspect-0.1d0)).and.(z.gt.0.6)) then !!downwelling on the right
!  temperature=0.0d0
!! elseif (((aspect/2.d0-0.025d0).lt.x).and.(x.lt.(aspect/2.d0+0.025d0)).and.(z.lt.0.2d0)) then
!!  temperature=1.d0  !!hot anomaly to initiate convection beneath lid
! else
!  temperature=0.8d0 !!temperature beneath lid 
! end if
!! temperature=0.5d0*(1.d0+dtanh((0.2d0+0.04d0*dcos(pii*x/aspect)-z)/0.04d0)) !!thermal RT problem
! temperature=1.d0-z !!conductive; for analytic solution case (internal heating will initiate convection)
 f0=0.d0 !!first sample problem from exact solution paper: f(t=0)
 temperature=(-pii**3*(aspect**2+1.d0)**2/aspect**3*cos(pii*x/aspect)*sin(pii*z)*f0&
           &+RaC(1)/(1.d0+exp(-2.d0*klayer*(dlayer-z)))+(RaT-RaC(1))*(1.d0-z))/RaT !!exact solution paper with sinusoidal perturbation
end if
end


real*8 function convective_cell(x,z,aratio) !!analytic expression from the van Keken et al. (1997) benchmark
 use basics
implicit none
integer*4 :: i,k
real*8 :: Tu,Tl,Tr,Ts,Q,u0,v0,x,z,perturb,aratio
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
end
