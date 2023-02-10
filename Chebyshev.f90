subroutine Chebyshev_coefficients(x) !!for RKC2
 use basics
 use arrays
integer*4 :: j
real*8 :: x

TCheb(0)=1.d0
TCheb(1)=x
TCheb_x(0)=0.d0
TCheb_x(1)=1.d0
TCheb_xx(0)=0.d0
TCheb_xx(1)=0.d0
do j=2,nstage
 TCheb(j)=2.d0*x*TCheb(j-1)-TCheb(j-2)
 TCheb_x(j)=2.d0*TCheb(j-1)+2.d0*x*TCheb_x(j-1)-TCheb_x(j-2)
 TCheb_xx(j)=4.d0*TCheb_x(j-1)+2.d0*x*TCheb_xx(j-1)-TCheb_xx(j-2)
end do
end

subroutine Chebyshev_coefficients_stream(x) !!for RKC1
 use basics
 use arrays
integer*4 :: j
real*8 :: x

TChebS(0)=1.d0
TChebS(1)=x
TCheb_xS(0)=0.d0
TCheb_xS(1)=1.d0
do j=2,nstageS
 TChebS(j)=2.d0*x*TChebS(j-1)-TChebS(j-2)
 TCheb_xS(j)=2.d0*TChebS(j-1)+2.d0*x*TCheb_xS(j-1)-TCheb_xS(j-2)
end do
end


subroutine Chebyshev_prep !!for RKC2
 use basics
 use arrays
implicit none
integer*4 :: j
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!startup
 omega0=1.d0+damping/real(nstage**2,8)
 call Chebyshev_coefficients(omega0)
 omega1=TCheb_x(nstage)/TCheb_xx(nstage)
 do j=2,nstage
  bCheb(j)=TCheb_xx(j)/TCheb_x(j)**2.d0
 end do
! bCheb(0)=bCheb(2)
! bCheb(1)=bCheb(2) !!!!!!!!!!!!b0 and b1 are free parameters -- more than one option used in the literature
 bCheb(0)=-bCheb(2)
 bCheb(1)=aRK2/omega1
 aCheb(0)=1.d0-bCheb(0)
 aCheb(1)=1.d0-bCheb(1)*omega0
 do j=2,nstage
  aCheb(j)=1.d0-bCheb(j)*TCheb(j)
 end do
 mewtCheb(1)=bCheb(1)*omega1
 do j=2,nstage
  mewCheb(j)=2.d0*bCheb(j)*omega0/bCheb(j-1)
  vCheb(j)=-bCheb(j)/bCheb(j-2)
  cCheb(j)=1.d0-mewCheb(j)-vCheb(j)
 end do
 do j=2,nstage
  mewtCheb(j)=2.d0*bCheb(j)*omega1/bCheb(j-1)
  gtCheb(j)=-aCheb(j-1)*mewtCheb(j)  
 end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end startup

write(*,*) "omega0,omega1=",omega0,omega1
write(*,*) "mewtCheb(1)=",mewtCheb(1)
write(*,*) "mewtCheb=",mewtCheb(2:nstage)
write(*,*) "gtCheb=",gtCheb(2:nstage)
write(*,*) "A21=",mewtCheb(1)
write(*,*) "A31=",mewCheb(2)*mewtCheb(1)+gtCheb(2)
write(*,*) "A32=",mewtCheb(2)
write(*,*) "B1=",mewCheb(3)*mewCheb(2)*mewtCheb(1)+mewCheb(3)*gtCheb(2)+vCheb(3)*mewtCheb(1)+gtCheb(3)
write(*,*) "B2=",mewCheb(3)*mewtCheb(2)
write(*,*) "B3=",mewtCheb(3)
end

subroutine Chebyshev_prep_stream !!for RKC1
 use basics
 use arrays
implicit none
integer*4 :: j
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!startup
 omega0S=1.d0+dampingS/real(nstageS**2,8)
 call Chebyshev_coefficients_stream(omega0S)
 omega1S=TChebS(nstageS)/TCheb_xS(nstageS)
 do j=0,nstageS
  bChebS(j)=1.d0/TChebS(j)
 end do
 mewtChebS(1)=bChebS(1)*omega1S
 do j=2,nstageS
  mewChebS(j)=2.d0*bChebS(j)*omega0S/bChebS(j-1)
  vChebS(j)=-bChebS(j)/bChebS(j-2)
  cChebS(j)=1.d0-mewChebS(j)-vChebS(j)
 end do
 do j=2,nstageS
  mewtChebS(j)=2.d0*bChebS(j)*omega1S/bChebS(j-1)
 end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end startup

RKC1_stability=-2.d0*omega0S/omega1S
write(*,*) "RKC1 Fine Grid Smoother:"
!write(*,*) "omega0S,omega1S=",omega0S,omega1S
!write(*,*) "mewtChebS=",mewtChebS(1:nstageS)
write(*,*) "RKC1_stability=",RKC1_stability
end


subroutine energy_time_Chebyshev
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,j,k,ID
real*8 :: energy_space,temp,tempr,temps,tempT,xr_interp,zs_interp
real*8 :: T0(-span:nx+span,-span:nz+span)
real*8, allocatable :: C0(:,:,:)
real*8 :: time1,time2,dto6,ta1,ta2,dtt,dtstage,dtstage2,tstage
real*8 :: u0(1:ntr),w0(1:ntr),delsqT0(1:ntr)

if (tstep.eq.0) then
 call Chebyshev_prep
end if

!$ time1=omp_get_wtime()
 call stability_RK4
if (equalize.eq.1) then
 if (original.eq.0) then
  call equalize_tracers !!updated equalize (better parallel scaling with variable viscosity constraints)
 else
  call equalize_tracers_original !!original equalize (some parts are serial -- without variable viscosity constraints)
 end if
end if
if (comp.eq.1) then
 allocate(C0(1:ntype,-span:nx+span,-span:nz+span))
 C0=Cnew
end if
 T0=Tratio
 Ttr0=Ttr
 rtr0=rtr
 str0=str

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Begin 1st stage
 if (RaT.ne.0.d0) then
  tstage=time
  call load_tracer_space_array(tstage)
  call compute_derivatives(tracer_space_array,DERr_tracer)
 end if
  call compute_derivatives(u,DERr_u)
  call compute_derivatives(w,DERr_w)
!$ ta1=omp_get_wtime()
  dtstage=mewtCheb(1)*dt
!$OMP PARALLEL DO PRIVATE(tempr,temps,tempT,xr_interp,zs_interp)
  do ID=1,ntr
   call interpolate(rtr0(ID),str0(ID),DERr_u,u,tempr)
   if (gridX.ne.0) then
    call interpolate_xr(rtr0(ID),xr_interp)
    tempr=tempr/xr_interp
   end if
   call interpolate(rtr0(ID),str0(ID),DERr_w,w,temps)
   if (gridZ.ne.0) then
    call interpolate_zs(str0(ID),zs_interp)
    temps=temps/zs_interp
   end if
   call interpolate(rtr0(ID),str0(ID),DERr_tracer,tracer_space_array,tempT)
    u0(ID)=tempr
    w0(ID)=temps
    delsqT0(ID)=tempT

   rtr1(ID)=rtr0(ID)+dtstage*tempr    !!!stuff for all RK methods
   str1(ID)=str0(ID)+dtstage*temps
   Ttr1(ID)=Ttr0(ID)+dtstage*tempT

   rtr2(ID)=rtr0(ID)    !!!stuff for all RK methods
   str2(ID)=str0(ID)
   Ttr2(ID)=Ttr0(ID)
  end do
!$OMP END PARALLEL DO
!$ ta2=omp_get_wtime()
!$ tadvect_tracer=tadvect_tracer+ta2-ta1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! End 1st Stage

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!remaining stages
do j=2,nstage
 if (comp.eq.1) Cnew=C0 !!reset initial composition values for second RK4 stage in case of empty cells
 call tracers_to_corners(rtr1,str1,Ttr1,T0)
 tstage=time+dtstage !!dtstage from previous stage
 call load_tracer_space_array(tstage)
 call compute_derivatives(tracer_space_array,DERr_tracer)
 if (RK_update.eq.1) then
  call update_velocities
  call compute_derivatives(u,DERr_u)
  call compute_derivatives(w,DERr_w)
 end if

!$ ta1=omp_get_wtime()
  dtstage=mewtCheb(j)*dt
  dtstage2=gtCheb(j)*dt
!$OMP PARALLEL DO PRIVATE(tempr,temps,tempT,xr_interp,zs_interp)
  do ID=1,ntr
   call interpolate(rtr1(ID),str1(ID),DERr_u,u,tempr)
   if (gridX.ne.0) then
    call interpolate_xr(rtr1(ID),xr_interp)
    tempr=tempr/xr_interp
   end if
   call interpolate(rtr1(ID),str1(ID),DERr_w,w,temps)
   if (gridZ.ne.0) then
    call interpolate_zs(str1(ID),zs_interp)
    temps=temps/zs_interp
   end if
    call interpolate(rtr1(ID),str1(ID),DERr_tracer,tracer_space_array,tempT)

    rtr(ID)=cCheb(j)*rtr0(ID)+mewCheb(j)*rtr1(ID)+vCheb(j)*rtr2(ID)+dtstage*tempr+dtstage2*u0(ID)
    str(ID)=cCheb(j)*str0(ID)+mewCheb(j)*str1(ID)+vCheb(j)*str2(ID)+dtstage*temps+dtstage2*w0(ID)
    Ttr(ID)=cCheb(j)*Ttr0(ID)+mewCheb(j)*Ttr1(ID)+vCheb(j)*Ttr2(ID)+dtstage*tempT+dtstage2*delsqT0(ID)

    rtr2(ID)=rtr1(ID) !!set up recursion for next stage if needed
    str2(ID)=str1(ID)
    Ttr2(ID)=Ttr1(ID)

    rtr1(ID)=rtr(ID) !!set up recursion for next stage if needed
    str1(ID)=str(ID)
    Ttr1(ID)=Ttr(ID)
  end do
!$OMP END PARALLEL DO
!$ ta2=omp_get_wtime()
!$ tadvect_tracer=tadvect_tracer+ta2-ta1
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end remaining stages

!!!!!!wrap up: convert tracers to corner values
 if (comp.eq.1) Cnew=C0 !!reset initial composition values for final RK4 stage in case of empty cells
 call tracers_to_corners(rtr,str,Ttr,T0)
 T=Tratio
!$ time2=omp_get_wtime()
!$ tenergy=tenergy+time2-time1

end

