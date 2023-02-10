subroutine viscosity_contrast(level,diff,T0)
!$ use OMP_LIB
 use basics
implicit none
real*8 :: a1,a2,a3,a4,diff,a,T0(-span:nx+span,-span:nz+span)
integer*4 :: i,k,level,nxx,nzz
if (iterate.eq.0) then
 nxx=nx
 nzz=nz
elseif (iterate.eq.1) then
 nxx=nx_grid(level)
 nzz=nz_grid(level)
end if

diff=0.d0
!$OMP PARALLEL DO PRIVATE(i,k,a1,a2,a3,a4,a) REDUCTION(max:diff)
do i=0,nxx-1
 do k=1,nzz-1
  a1=max(T0(i,k-1),T0(i+1,k-1),T0(i+1,k),T0(i,k+1),T0(i+1,k+1))
  a2=min(T0(i,k-1),T0(i+1,k-1),T0(i+1,k),T0(i,k+1),T0(i+1,k+1))
  a3=max(a1,T0(i,k))/min(a1,T0(i,k))
  a4=max(a2,T0(i,k))/min(a2,T0(i,k))
  a=max(a3,a4)
  if (a.gt.diff) diff=a
 end do
end do
!$OMP END PARALLEL DO
end


subroutine compute_smoothing_parameters
 use basics
 use arrays
implicit none
integer*4 :: i
real*8 :: A(1:ntype)
!!!!!!!!!!!!!!!!!!!!!!!!!!Temperature Smoothing
if ((ivis.eq.0).or.(RaT.eq.0.d0)) then
 delta_T=2.d0 !!no smoothing required
 delta_T_init=2.d0
elseif (ivis.ge.1) then
 if (visT.ne.1.d0) then
  delta_T=dlog(delta_visT)/dlog(visT)
  delta_T_init=dlog(delta_visT_init)/dlog(visT)
 else
  delta_T=2.d0
  delta_T_init=2.d0
 end if
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!End Temperature Smoothing

!!!!!!!!!!!!!!!!!!!!!!!!!!Composition Smoothing
if (comp.eq.1) then
 allocate(delta_C(1:ntype))
 delta_C=2.d0 !!default -- no smoothing
 do i=1,ntype
  if (conduct_factor(i).ne.1.d0) then
   delta_C(i)=dabs(dlog(delta_visC(i))/dlog(1.d0/conduct_factor(i))) 
  end if
 end do
 write(*,*) "conductivity delta_C=",delta_C
 if (ivis.ge.2) then
  do i=1,ntype
   if (visC(i).ne.1.d0) then
    delta_C(i)=dabs(dlog(delta_visC(i))/dlog(visC(i))) !!!!!vectorize
   end if
  end do
 end if
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!End Composition Smoothing
write(*,*) "T Smoothing Parameter:", delta_T
if (comp.eq.1) write(*,*) "C Smoothing Parameters:",delta_C
end

subroutine smoothness(diff,T0)
!$ use OMP_LIB
 use basics
implicit none
real*8 :: a1,a2,a3,a4,a5,diff,a,T0(-span:nx+span,-span:nz+span)
integer*4 :: i,k
diff=0.d0
!$OMP PARALLEL DO PRIVATE(i,k,a1,a2,a3,a4,a5,a) REDUCTION(max:diff)
do i=0,nx-1
 do k=1,nz-1
  a1=dabs(T0(i,k-1)-T0(i,k))
  a2=dabs(T0(i+1,k-1)-T0(i,k))
  a3=dabs(T0(i+1,k)-T0(i,k))
  a4=dabs(T0(i,k+1)-T0(i,k))
  a5=dabs(T0(i+1,k+1)-T0(i,k))
  a=max(a1,a2,a3,a4,a5)
  if (a.gt.diff) diff=a
 end do
end do
!$OMP END PARALLEL DO
end

subroutine smoothness_residual(level,diff,r0)
!$ use OMP_LIB
 use basics
implicit none
real*8 :: a1,a2,a3,a4,a5,diff,a,r0(1:nx-1,1:nz-1)
integer*4 :: i,k,level
diff=0.d0
!$OMP PARALLEL DO PRIVATE(i,k,a1,a2,a3,a4,a5,a) REDUCTION(max:diff)
do i=1,nx_grid(level)-2
 do k=2,nz_grid(level)-2
  a1=dabs(r0(i,k-1)-r0(i,k))
  a2=dabs(r0(i+1,k-1)-r0(i,k))
  a3=dabs(r0(i+1,k)-r0(i,k))
  a4=dabs(r0(i,k+1)-r0(i,k))
  a5=dabs(r0(i+1,k+1)-r0(i,k))
  a=max(a1,a2,a3,a4,a5)
  if (a.gt.diff) diff=a
 end do
end do
!$OMP END PARALLEL DO
diff=diff/real(nx_grid(level)-2,8)/real(nz_grid(level)-2,8) !!scale since r0 array has not been normalized
end


real*8 function smoother_space(i,k,T0) !!T0 represents the composition field
 use basics
 use arrays
implicit none
integer*4 :: i,k
real*8 :: Tr,Ts,Trr,Tss
real*8 :: Tx,Tz,Txx,Tzz
real*8 :: T0(-span:nx+span,-span:nz+span)
Tr=dot_product(D1(-span1:span1),T0(i-span1:i+span1,k))/dr
Ts=dot_product(D1(-span1:span1),T0(i,k-span1:k+span1))/ds
Trr=dot_product(D2(-span1:span1),T0(i-span1:i+span1,k))/dr2
Tss=dot_product(D2(-span1:span1),T0(i,k-span1:k+span1))/ds2

Tx=Tr/xr(i)
Tz=Ts/zs(k)
Txx=(Trr-Tx*xrr(i))/xr(i)**2.d0
Tzz=(Tss-Tz*zss(k))/zs(k)**2.d0
smoother_space=Txx+Tzz
end

real*8 function vis_space(i,k,T0,level) !!T0 represents the composition field
 use basics
 use arrays
implicit none
integer*4 :: i,k,level,i1,k1
real*8 :: Tr,Ts,Trr,Tss
real*8 :: Tx,Tz,Txx,Tzz
real*8 :: T0(-span:nx+span,-span:nz+span)

if (gridX.eq.0) then
 Txx=dot_product(D2_so(-1:1),T0(i-1:i+1,k))/dr2_grid(level)
else
 i1=i*2**(level-1) !!fine grid coordinates for grid metrics (may need to restrict the metrics eventually)
 Tr=dot_product(D1_so(-1:1),T0(i-1:i+1,k))/dr_grid(level)     !!!!second order FD stencils
 Trr=dot_product(D2_so(-1:1),T0(i-1:i+1,k))/dr2_grid(level)
 Tx=Tr/xr(i1)
 Txx=(Trr-Tx*xrr(i1))/xr(i1)**2.d0
end if

if (gridZ.eq.0) then
 Tzz=dot_product(D2_so(-1:1),T0(i,k-1:k+1))/ds2_grid(level)
else
 k1=k*2**(level-1)
 Ts=dot_product(D1_so(-1:1),T0(i,k-1:k+1))/ds_grid(level)
 Tss=dot_product(D2_so(-1:1),T0(i,k-1:k+1))/ds2_grid(level)
 Tz=Ts/zs(k1)
 Tzz=(Tss-Tz*zss(k1))/zs(k1)**2.d0
end if

vis_space=Txx+Tzz
end

real*8 function stream_smoother_space_uniform(i,k,T0,level) !!T0 represents the composition field
 use basics
 use arrays
implicit none
integer*4 :: i,k,level
real*8 :: T0(-span:nx+span,-span:nz+span)
 stream_smoother_space_uniform=(sum(D2(-span1:span1)*T0(i-span1:i+span1,k)+&
 &D2(-span1:span1)*T0(i,k-span1:k+span1)))/dr2_grid(level)
end

real*8 function stream_smoother_space(i,k,T0,level) !!T0 represents the composition field
 use basics
 use arrays
implicit none
integer*4 :: i,k,level,i1,k1
real*8 :: Tr,Ts,Trr,Tss
real*8 :: Tx,Tz,Txx,Tzz
real*8 :: T0(-span:nx+span,-span:nz+span)
if ((gridX.eq.0).and.(gridZ.eq.0)) then
 stream_smoother_space=(sum(D2(-span1:span1)*T0(i-span1:i+span1,k)+&
 &D2(-span1:span1)*T0(i,k-span1:k+span1)))/dr2_grid(level)
else
 if (gridX.eq.0) then
  Txx=dot_product(D2(-span1:span1),T0(i-span1:i+span1,k))/dr2_grid(level)
 else
  i1=i*2**(level-1) !!fine grid coordinates for grid metrics (may need to restrict the metrics eventually)
  Tr=dot_product(D1(-span1:span1),T0(i-span1:i+span1,k))/dr_grid(level)
  Trr=dot_product(D2(-span1:span1),T0(i-span1:i+span1,k))/dr2_grid(level)
  Tx=Tr/xr(i1)
  Txx=(Trr-Tx*xrr(i1))/xr(i1)**2.d0
 end if

 if (gridZ.eq.0) then
  Tzz=dot_product(D2(-span1:span1),T0(i,k-span1:k+span1))/ds2_grid(level)
 else
  k1=k*2**(level-1)
  Ts=dot_product(D1(-span1:span1),T0(i,k-span1:k+span1))/ds_grid(level)
  Tss=dot_product(D2(-span1:span1),T0(i,k-span1:k+span1))/ds2_grid(level)
  Tz=Ts/zs(k1)
  Tzz=(Tss-Tz*zss(k1))/zs(k1)**2.d0
 end if
 stream_smoother_space=Txx+Tzz
end if
end

subroutine smoother_strain(level,interpolate) !!!!!!!smooth strain rate to avoid HUGE gradients in non-Newtonian viscosity
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k,n,level,interpolate
integer*4 :: nn,ntemp,Niter
real*8 :: stream_smoother_space,stream_smoother_space_uniform,temp,c
real*8 :: T0(-span:nx+span,-span:nz+span)
real*8 :: la,lb
if (interpolate.eq.0) then
 if (level.eq.1) then
  Niter=6
  lb=-5.d0
  la=-2.d0
 else
  Niter=3
  lb=-6.d0
  la=-4.d0
 end if
end if

if (interpolate.eq.1) then
 Niter=vis_iter   !!!inspired by smoother_vis: for interpolation of stream function to coarser grids during multigrid: only retain the lower half of the frequency spectrum
 lb=-2.d0
 la=-1.d0
end if

if ((gridX.eq.0).and.(gridZ.eq.0)) then !!!!!!!!uniform grid: for speed
 do nn=1,Niter-1 !!for Richardson's method
  ntemp=Niter-nn+1 !!!start with smaller time steps to kill off high frequencies first
  c=0.99d0*(2.d0/(-lb-la+(lb-la)*dcos(real(2*ntemp-1,8)*pii/real(2*Niter,8))))

  !$OMP PARALLEL PRIVATE(i,k)
  !$OMP DO
  do k=-span1,nz_grid(level)+span1
   do i=-span1,nx_grid(level)+span1
    T0(i,k)=strain(i,k,level) !!feed in fine stream function field to be smoothed
   end do
  end do
  !$OMP END DO
  !$OMP DO
  do k=0,nz_grid(level)
   do i=0,nx_grid(level)
    strain(i,k,level)=T0(i,k)+c*dt_stream_smoother(i,k,level)*((sum(D2(-span1:span1)*T0(i-span1:i+span1,k)+&
    &D2(-span1:span1)*T0(i,k-span1:k+span1)))/dr2_grid(level))
   end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  call strainBCs(level)
 end do

  ntemp=1 !!!!!compute yield to strain ratio on last iteration
  c=0.99d0*(2.d0/(-lb-la+(lb-la)*dcos(real(2*ntemp-1,8)*pii/real(2*Niter,8))))

  if (plastic.ge.0) then
   !$OMP PARALLEL PRIVATE(i,k)
   !$OMP DO
   do k=-span1,nz_grid(level)+span1
    do i=-span1,nx_grid(level)+span1
     T0(i,k)=strain(i,k,level) !!feed in fine stream function field to be smoothed
    end do
   end do
   !$OMP END DO
   !$OMP DO
   do k=0,nz_grid(level)
    do i=0,nx_grid(level)
     strain(i,k,level)=T0(i,k)+c*dt_stream_smoother(i,k,level)*((sum(D2(-span1:span1)*T0(i-span1:i+span1,k)+&
    &D2(-span1:span1)*T0(i,k-span1:k+span1)))/dr2_grid(level))
     yield_strain(i,k,level)=visNN_static(i,k,level)/(strain(i,k,level)+eps) !!for plastic yielding -- avoid blow up
    end do
   end do
   !$OMP END DO
   !$OMP END PARALLEL
  else
   !$OMP PARALLEL PRIVATE(i,k)
   !$OMP DO
   do k=-span1,nz_grid(level)+span1
    do i=-span1,nx_grid(level)+span1
     T0(i,k)=strain(i,k,level) !!feed in fine stream function field to be smoothed
    end do
   end do
   !$OMP END DO
   !$OMP DO
   do k=0,nz_grid(level)
    do i=0,nx_grid(level)
     strain(i,k,level)=T0(i,k)+c*dt_stream_smoother(i,k,level)*((sum(D2(-span1:span1)*T0(i-span1:i+span1,k)+&
    &D2(-span1:span1)*T0(i,k-span1:k+span1)))/dr2_grid(level))
     yield_strain(i,k,level)=(Cpla*strain(i,k,level)**mpla+visNN_static(i,k,level))/(strain(i,k,level)+eps) !!for plastic yielding
    end do
   end do
   !$OMP END DO
   !$OMP END PARALLEL
  end if
else !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! grid refinement
 do nn=1,Niter-1 !!for Richardson's method
  ntemp=Niter-nn+1 !!!start with smaller time steps to kill off high frequencies first
  c=0.99d0*(2.d0/(-lb-la+(lb-la)*dcos(real(2*ntemp-1,8)*pii/real(2*Niter,8))))
 
  T0=strain(-span:nx+span,-span:nz+span,level) !!feed in fine stream function field to be smoothed
  !$OMP PARALLEL DO PRIVATE(i,k)
  do k=0,nz_grid(level)
   do i=0,nx_grid(level)
    strain(i,k,level)=T0(i,k)+c*dt_stream_smoother(i,k,level)*stream_smoother_space(i,k,T0,level)
   end do
  end do
  !$OMP END PARALLEL DO
  call strainBCs(level)
 end do

  ntemp=1 !!!!!compute yield to strain ratio on last iteration
  c=0.99d0*(2.d0/(-lb-la+(lb-la)*dcos(real(2*ntemp-1,8)*pii/real(2*Niter,8))))

  T0=strain(-span:nx+span,-span:nz+span,level) !!feed in fine stream function field to be smoothed
  if (plastic.gt.0) then
   !$OMP PARALLEL DO PRIVATE(i,k)
   do k=0,nz_grid(level)
    do i=0,nx_grid(level)
     strain(i,k,level)=T0(i,k)+c*dt_stream_smoother(i,k,level)*stream_smoother_space(i,k,T0,level)
     yield_strain(i,k,level)=visNN_static(i,k,level)/(strain(i,k,level)+eps) !!for plastic yielding
    end do
   end do
   !$OMP END PARALLEL DO
  else
   !$OMP PARALLEL DO PRIVATE(i,k)
   do k=0,nz_grid(level)
    do i=0,nx_grid(level)
     strain(i,k,level)=T0(i,k)+c*dt_stream_smoother(i,k,level)*stream_smoother_space(i,k,T0,level)
     yield_strain(i,k,level)=(Cpla*strain(i,k,level)**mpla+visNN_static(i,k,level))/(strain(i,k,level)+eps) !!for plastic yielding
    end do
   end do
   !$OMP END PARALLEL DO
  end if
end if
end

subroutine strainBCs(level)
 use basics
 use arrays
implicit none
integer*4 :: i,k,level
do k=0,nz_grid(level)
 do i=1,span
  strain(-i,k,level)=strain(i,k,level)
  strain(nx_grid(level)+i,k,level)=strain(nx_grid(level)-i,k,level)
 end do
end do

do k=1,span
 do i=-span,nx_grid(level)+span
  strain(i,-k,level)=strain(i,k,level)
 end do
end do
do k=1,span
 do i=-span,nx_grid(level)+span
  strain(i,nz_grid(level)+k,level)=strain(i,nz_grid(level)-k,level)
 end do
end do
end

subroutine smoother_stream(level,interpolate) !!!!!!!smooth stream function in preparation for restriction during multigrid
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k,n,level,interpolate
integer*4 :: nn,ntemp,Niter
real*8 :: stream_smoother_space,c
real*8 :: T0(-span:nx+span,-span:nz+span)
real*8 :: la,lb
if (interpolate.eq.0) then
 Niter=6   !!!remove unresolved frequencies
 lb=-5.d0
 la=-2.d0
end if

if (interpolate.eq.1) then
 Niter=vis_iter   !!!inspired by smoother_vis: for interpolation of stream function to coarser grids during multigrid: only retain the lower half of the frequency spectrum
 lb=-2.d0
 la=-1.d0
end if

do nn=1,Niter !!for Richardson's method

 ntemp=Niter-nn+1 !!!start with smaller time steps to kill off high frequencies first
 c=2.d0/(-lb-la+(lb-la)*dcos(real(2*ntemp-1,8)*pii/real(2*Niter,8)))
 c=0.99d0*c !!to ensure stability

 T0=error(-span:nx+span,-span:nz+span,level) !!feed in fine stream function field to be smoothed
!$OMP PARALLEL DO PRIVATE(i,k)
 do k=1,nz_grid(level)-1
  do i=1,nx_grid(level)-1  
   stream_smooth(i,k,level)=T0(i,k)+c*dt_stream_smoother(i,k,level)*stream_smoother_space(i,k,T0,level)
  end do
 end do
!$OMP END PARALLEL DO
 if (nn.lt.Niter) call enforceBCs_stream_smoother(level,stream_smooth(:,:,level))
end do
 call enforceBCs_stream(level,stream_smooth(:,:,level)) !!ensure that the BCs are appropriate for the biharmonic equation
end

subroutine enforceBCs_stream_smoother(level,S0) !!!!update this
 use basics
implicit none
integer*4 :: i,k,level
real*8 :: S0(-span:nx+span,-span:nz+span)
S0(0,-span:nz_grid(level)+span)=0.d0   !!impermeable boundaries
S0(nx_grid(level),-span:nz_grid(level)+span)=0.d0
S0(-span:nx_grid(level)+span,0)=0.d0
S0(-span:nx_grid(level)+span,nz_grid(level))=0.d0
do k=0,nz_grid(level)
 do i=1,span  !!!antisymmetry for stream=0 on boundaries
  S0(-i,k)=-S0(i,k) 
 end do
end do
do k=0,nz_grid(level)
 do i=1,span  !!!antisymmetry for stream=0 on boundaries
  S0(nx_grid(level)+i,k)=-S0(nx_grid(level)-i,k)
 end do
end do

do k=1,span
 do i=-span,nx_grid(level)+span
  S0(i,-k)=-S0(i,k)
 end do
end do
do k=1,span
 do i=-span,nx_grid(level)+span
  S0(i,nz_grid(level)+k)=-S0(i,nz_grid(level)-k)
 end do
end do

end


subroutine smoother_vis(level,interpolate) !!smoother for the viscosity field
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,j,k,n,level,interpolate
integer*4 :: nn,ntemp,Niter,trigger,option
real*8 :: vis_space,temp,diff,c,ratio_min
real*8 :: t11,t22,t1,t2,time1,time2,la,lb
!$ time1=omp_get_wtime()

!$ t1=omp_get_wtime() !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! A
 call viscosityBCs(level)

if (interpolate.eq.0) then
  option=1
  Niter=3
  lb=-3.d0
  la=-2.d0
elseif (interpolate.eq.1) then          !!!for interpolation of viscosity to coarser grids during multigrid: only retain the lower half of the frequency spectrum
 option=0
 Niter=vis_iter
 lb=-2.d0
 la=-1.d0
end if
 Nvis_store=0
 nn_store=1
!$ t2=omp_get_wtime()
!$ t_svA=t_svA+(t2-t1) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! end A

n=0
trigger=0
do
do nn=nn_store,Niter !!for Richardson's method

!$ t1=omp_get_wtime() !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! B
!$ t11=omp_get_wtime() !!!!!!!!!!!!!!!!BA
  call viscosity_gradients(level,Nvis_store)
!$ t22=omp_get_wtime()
!$ t_svBA=t_svBA+(t22-t11) !!!!!!!!!!end BA
 if (interpolate.eq.0) then
  if ((Nvis_store.ge.1).or.(ivis.eq.0)) then
!$ t11=omp_get_wtime() !!!!!!!!!!!!!!!!BB
   call test_eigenvalues(level)
!$ t22=omp_get_wtime()
!$ t_svBB=t_svBB+(t22-t11) !!!!!!!!!!end BB
!$ t11=omp_get_wtime() !!!!!!!!!!!!!!!!BC
   ratio_min=1.d100
   !$OMP PARALLEL DO PRIVATE(i,k) REDUCTION(min:ratio_min)
   do k=1,nz_grid(level)-1
    do i=1,nx_grid(level)-1
     if (ratio(i,k,level).lt.ratio_min) ratio_min=ratio(i,k,level)
    end do
   end do
   !$OMP END PARALLEL DO
!$ t22=omp_get_wtime()
!$ t_svBC=t_svBC+(t22-t11) !!!!!!!!!!end BC
!  write(*,*) "smoother_vis:",fp,level,n,nn,ratio_min
  end if
 end if
!$ t2=omp_get_wtime()
!$ t_svB=t_svB+(t2-t1) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! end B

 if (option.eq.1) then
   
   if (((Nvis_store.ge.1).and.(ratio_min.gt.eig_ratio_min)).or.(ivis.eq.0).or.(n.gt.20)) then !!!!!!original -- require at least one full filter pass
    eig_fine=ratio_min
    n_eta=n
!    write(*,*) "smoother_vis:",fp,level,n,nn,ratio_min,eig_ratio_min(level)
    trigger=1
    exit
   end if

 elseif (option.eq.0) then
  if (interpolate.eq.0) then !!!prep for smoother iterations
   if ((ratio_min.gt.eig_ratio).or.(ivis.eq.0)) then
!    write(*,*) "smoother_vis:",level,n,nn,ratio_min,interpolate
    trigger=1
    nn_store=nn
    exit
   end if
  elseif (interpolate.ge.1) then !!!prep for restriction of viscosity to coarser grid
    if ((Nvis_store.ge.1).or.(ivis.eq.0)) then
!     write(*,*) "smoother_vis:",level,n,nn,ratio_min,interpolate
     trigger=1
     exit
    end if
  end if
 end if

 ntemp=Niter-nn+1 !!!start with smaller time steps to kill off high frequencies first
 c=2.d0/(-lb-la+(lb-la)*dcos(real(2*ntemp-1,8)*pii/real(2*Niter,8)))
 c=0.99d0*c !!!to ensure stability (really for dt_vis array)

!$ t1=omp_get_wtime() !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! C
!!!!!!!!!!interior points
!$OMP PARALLEL DO PRIVATE(i,k)
 do k=0,nz_grid(level)
  do i=0,nx_grid(level)  !!!!depends on assumed BCs for this smoother (will not affect physical gradients)
   vis_grid(i,k,level)=vis_grid(i,k,level)+c*vis_smooth(i,k,level)!!!!!dt_vis(i,k,level)*vis_xxpzz(i,k,level)
  end do
 end do
!$OMP END PARALLEL DO
!!!!!!!!!end interior points
!$ t2=omp_get_wtime()
!$ t_svC=t_svC+(t2-t1) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! end C

!$ t1=omp_get_wtime() !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! D
call viscosityBCs(level)
!$ t2=omp_get_wtime()
!$ t_svD=t_svD+(t2-t1) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! end D
end do

 if (trigger.eq.1) exit
 Nvis_store=Nvis_store+1
 nn_store=1
 if (option.eq.1) then !!!Once high unresolved frequencies are damped, start damping resolved frequencies to increase the eigenvalue ratio
  nn_store=1
 Niter=6
 lb=-2.d0
 la=-1.d0
 end if
n=n+1
end do
!$ time2=omp_get_wtime()
!$ t_sv=t_sv+(time2-time1)
end

subroutine smoother_time !!smoother for the C field
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k,ID,n,nn,tt,ntemp,Niter,trigger
real*8 :: smoother_space,temp,diff,c,dts
real*8 :: T0(-span:nx+span,-span:nz+span)
real*8 :: time1,time2,la,lb
 Niter=10  !!!first pass: filter out unresolved frequencies and the upper half of the resolved spectrum to smooth C while supressing spurious oscillations
 lb=-5.d0
 la=-1.d0

 Cvis=Cnew
 Cbuoy=Cnew
do tt=1,ntype
nn=0
trigger=0
do
do n=1,Niter
 call smoothness(diff,Cvis(tt,:,:))
 if (diff.lt.delta_C(tt)) then
  trigger=1
  exit
 end if

 if (nn.eq.0) then !!!filter out unresolved frequencies on first pass
  ntemp=Niter-n+1 !!!start with smaller time steps to kill off high frequencies first
  c=2.d0/(-lb-la+(lb-la)*dcos(real(2*ntemp-1,8)*pii/real(2*Niter,8)))
 else
  c=1.0d0  !!continue with more efficient smoothing once highest frequencies are damped
 end if

T0=Cvis(tt,-span:nx+span,-span:nz+span) !!feed in sharp C field to be smoothed
!$OMP PARALLEL DO PRIVATE(i,k,temp,dts)
do i=0,nx
 do k=0,nz
  temp=smoother_space(i,k,T0)
  dts=dts_array(i,k)*courantC
  Cvis(tt,i,k)=T0(i,k)+c*dts*temp
  if (Cvis(tt,i,k).lt.0.d0) Cvis(tt,i,k)=0.d0
  if (Cvis(tt,i,k).gt.1.d0) Cvis(tt,i,k)=1.d0
 end do
end do
!$OMP END PARALLEL DO
 call enforceBCs_comp(Cvis(tt,:,:))
end do
 if (trigger.eq.1) exit
 Niter=1
 nn=nn+1
end do
 n_comp(tt)=nn
! write(*,*) "C smoother:",nn,n,diff,tt
end do
end

subroutine smoother_time_T(Tin) !!smoother for the temperature field
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k,ID,n,nsmall,Nsmoother,ntemp
real*8 :: smoother_space,temp,diff,dts
real*8 :: T0(-span:nx+span,-span:nz+span)
real*8 :: time1,time2,Tin(-span:nx+span,-span:nz+span)
 Tvis=Tin
 Tbuoy=Tin
if ((RaT.eq.0).or.(ivis.eq.0)) return
n=0
do
 call smoothness(diff,Tvis)
 if (diff.lt.delta_T) then
  n_temp=n
!  write(*,*) "T smoother:",n,diff,tstep
  exit
 end if

T0=Tvis !!feed in sharp C field to be smoothed
!$OMP PARALLEL DO PRIVATE(i,k,temp,dts)
do k=0,nz
 do i=0,nx
  temp=smoother_space(i,k,T0)
  dts=dts_array(i,k)*courantT
  Tvis(i,k)=T0(i,k)+dts*temp
 end do
end do
!$OMP END PARALLEL DO
 call enforceBCs(Tvis)
n=n+1
end do
end

subroutine smoother_initial_T !!smoother for the initial temperature field
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k,ID,n,nsmall,Nsmoother,ntemp
real*8 :: smoother_space,temp,diff,dts
real*8 :: T0(-span:nx+span,-span:nz+span),Tdumb(-span_interp:nx+span_interp,-span_interp:nz+span_interp)
real*8 :: time1,time2,Tin(-span:nx+span,-span:nz+span)
n=0
do
 call smoothness(diff,T)
 if (diff.lt.delta_T_init) then
  write(*,*) "initial T smoother:",n,diff,tstep
  exit
 end if

T0=T !!feed in sharp C field to be smoothed
!$OMP PARALLEL DO PRIVATE(i,k,temp,dts)
do i=0,nx
 do k=0,nz
  temp=smoother_space(i,k,T0)
  dts=dts_array(i,k)*courantT
  T(i,k)=T0(i,k)+dts*temp
 end do
end do
!$OMP END PARALLEL DO
 call enforceBCs(T)
n=n+1
end do

if (tracer.eq.1) then
 call extendT(T) !!for interpolation of T to Ttr during equalizing
 call compute_derivatives(Textend(-span_interp:nx+span_interp,-span_interp:nz+span_interp),DERr_Textend)
 Tdumb=Textend(-span_interp:nx+span_interp,-span_interp:nz+span_interp)
!$OMP PARALLEL DO PRIVATE(ID)
 do ID=1,ntr !!!interpolate temperature values to repositioned tracers
  call interpolate(rtr(ID),str(ID),DERr_Textend,Tdumb,Ttr(ID))
 end do
!$OMP END PARALLEL DO
 call tracers_to_corners(rtr,str,Ttr,T)
end if
end


subroutine stability_RK4_smoother !!actually, this is for RK1 (switched from RK4 to enable Richardson iterations)
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k
real*8 :: eigA_smoother,eigD_smoother,dtA,dtD,eR,eI
Ts_max=Tr_base/ds     !!scale according to the denominators of finite difference formulas
Tss_max=-Trr_base/ds2 !!even # of derivative operators=real, odd #=imaginary. For real values keep signs (sign=(imaginary unit)^(# of derivatives)), for imaginary eigenvalues take absolute values
Tr_max=Tr_base/dr
Trr_max=-Trr_base/dr2

!$OMP PARALLEL DO PRIVATE(i,k,eR,eI)
do k=0,nz
 do i=0,nx !!using local time steps for smoother to account for grid spacing
  eR=eigD_smoother(i,k)
  eI=eigA_smoother(i,k)
  dts_array(i,k)=2.d0*eR/(eR**2.d0+eI**2.d0)
 end do
end do
!$OMP END PARALLEL DO
end

real*8 function eigA_smoother(i,k)
 use basics
 use arrays
implicit none
integer*4 :: i,k
!!!!!!mathematical advection from grid spacing
  eigA_smoother=dabs(xrr(i)/(xr(i)**3.d0))*Tr_max+dabs(zss(k)/(zs(k)**3.d0))*Ts_max
end


real*8 function eigD_smoother(i,k)
 use basics
 use arrays
implicit none
integer*4 :: i,k
  eigD_smoother=dabs(Trr_max/(xr(i)**2.d0))+dabs(Tss_max/(zs(k)**2.d0))
end


subroutine stability_vis_smoother(level)
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k,level
real*8 :: eigA_grid,eigD_grid,eR,eI
Ts_max_so=Tr_base_so/ds_grid(level)     !!scale according to the denominators of (second order) finite difference formulas
Tss_max_so=-Trr_base_so/ds2_grid(level) !!even # of derivative operators=real, odd #=imaginary. For real values keep signs (sign=(imaginary unit)^(# of derivatives)), for imaginary eigenvalues take absolute values
Tr_max_so=Tr_base_so/dr_grid(level)
Trr_max_so=-Trr_base_so/dr2_grid(level)

!$OMP PARALLEL DO PRIVATE(i,k,eR,eI)
do k=0,nz_grid(level)
 do i=0,nx_grid(level) !!using local time steps for smoother to account for grid spacing
  eR=eigD_grid(level,i,k)
  eI=eigA_grid(level,i,k)
  dt_vis(i,k,level)=2.d0*eR/(eR**2.d0+eI**2.d0)
 end do
end do
!$OMP END PARALLEL DO
end

real*8 function eigA_grid(level,i,k)
 use basics
 use arrays
implicit none
integer*4 :: i,k,level,i1,k1
  i1=i*2**(level-1) !!fine grid coordinates for grid metrics (may need to restrict the metrics eventually)
  k1=k*2**(level-1)
!!!!!!mathematical advection from grid spacing
  eigA_grid=dabs(xrr(i1)/(xr(i1)**3.d0))*Tr_max_so+dabs(zss(k1)/(zs(k1)**3.d0))*Ts_max_so
end


real*8 function eigD_grid(level,i,k)
 use basics
 use arrays
implicit none
integer*4 :: i,k,level,i1,k1
  i1=i*2**(level-1) !!fine grid coordinates for grid metrics (may need to restrict the metrics eventually)
  k1=k*2**(level-1)
  eigD_grid=dabs(Trr_max_so/(xr(i1)**2.d0))+dabs(Tss_max_so/(zs(k1)**2.d0))
end

subroutine stability_stream_smoother(level)
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k,level
real*8 :: eigA_stream_smooth,eigD_stream_smooth,eR,eI
Ts_max=Tr_base/ds_grid(level)     !!scale according to the denominators of finite difference formulas
Tss_max=-Trr_base/ds2_grid(level) !!even # of derivative operators=real, odd #=imaginary. For real values keep signs (sign=(imaginary unit)^(# of derivatives)), for imaginary eigenvalues take absolute values
Tr_max=Tr_base/dr_grid(level)
Trr_max=-Trr_base/dr2_grid(level)

!$OMP PARALLEL DO PRIVATE(i,k,eR,eI)
do k=0,nz_grid(level)
 do i=0,nx_grid(level) !!using local time steps for smoother to account for grid spacing
  eR=eigD_stream_smooth(level,i,k)
  eI=eigA_stream_smooth(level,i,k)
  dt_stream_smoother(i,k,level)=2.d0*eR/(eR**2.d0+eI**2.d0)
 end do
end do
!$OMP END PARALLEL DO
end

real*8 function eigA_stream_smooth(level,i,k)
 use basics
 use arrays
implicit none
integer*4 :: i,k,level,i1,k1
  i1=i*2**(level-1) !!fine grid coordinates for grid metrics (may need to restrict the metrics eventually)
  k1=k*2**(level-1)
!!!!!!mathematical advection from grid spacing
  eigA_stream_smooth=dabs(xrr(i1)/(xr(i1)**3.d0))*Tr_max+dabs(zss(k1)/(zs(k1)**3.d0))*Ts_max
end


real*8 function eigD_stream_smooth(level,i,k)
 use basics
 use arrays
implicit none
integer*4 :: i,k,level,i1,k1
  i1=i*2**(level-1) !!fine grid coordinates for grid metrics (may need to restrict the metrics eventually)
  k1=k*2**(level-1)
  eigD_stream_smooth=dabs(Trr_max/(xr(i1)**2.d0))+dabs(Tss_max/(zs(k1)**2.d0))
end

