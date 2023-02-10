subroutine optimize_grid(nx,nz,grids,n_points_min)
implicit none
integer*4 :: n_in,grids,nx,nz,nxbig
integer*4 :: n,n_points,n_coarse,n_temp,n_coarse_temp,n_points_temp
integer*4 :: n_min,n_max,n_points_max,n_coarse_max,n_points_min,n_coarse_min
n_points_max=min(nx,nz)   !!!max # of grid levels on coarsest grid
n_coarse_min=1
n_coarse_max=grids-1   !!!max # of coarse grids
n_min=n_points_min*2**n_coarse_min
n_max=n_points_max*2**n_coarse_max 

if (nz.le.nx) then !!# of coarse grids constrained by smaller of nx and nz
 nxbig=1
 n_in=nz
else
 nxbig=0
 n_in=nx
end if
if (n_min.gt.n_in) then
 write(*,*) "Grid Cannot Be Optimized With Current Settings"
 write(*,*) "Current settings allow a minimum grid size of",n_min,"."
end if
if (n_max.lt.n_in) then
 write(*,*) "Grid Cannot Be Optimized With Current Settings"
 write(*,*) "Current settings allow a maximum grid size of",n_max,"."
end if
n_temp=n_max !!start with upper limit
do n_points=n_points_min,n_points_max!!!range for # of grid levels on coarsest grid
 do n_coarse=n_coarse_min,n_coarse_max!!!range for total # of coarse grids 
  n=n_points*2**n_coarse !!# of grid levels on fine grid
  if ((n.ge.n_in).and.(n.lt.n_temp)) then
   n_temp=n
   n_coarse_temp=n_coarse
   n_points_temp=n_points
  end if
 end do
end do
!write(*,*) n_temp,n_points_temp,n_coarse_temp
grids=n_coarse_temp+1
if (nxbig.eq.0) then
 nx=n_temp
else
 nz=n_temp
end if

!!!given the # of coarse grids, figure out the resolution for the larger of nx and nz
if (nz.le.nx) then !!# of coarse grids constrained by smaller of nx and nz
 n_in=nx
else
 n_in=nz
end if
n_points_temp=nint(real(n_in,8)/real(2**n_coarse_temp,8),4)
n_temp=n_points_temp*2**n_coarse_temp
!write(*,*) n_temp,n_points_temp,n_coarse_temp
if (nxbig.eq.0) then
 nz=n_temp
else
 nx=n_temp
end if
write(*,*) "# of multigrid levels=",grids
end

subroutine print_field(level,field,fname)
 use basics
 use arrays
implicit none
integer*4 :: i,k,level
real*8 :: field(1:nx-1,1:nz-1)
character*6 :: fname
open(unit=999,file=fname)
do i=1,nx_grid(level)-1
 do k=1,nz_grid(level)-1
   write(999,*) xg(i),zg(k),field(i,k)
 end do
end do
 close(999)
end

subroutine iterative_solver
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: level,dir,p,q
real*8 :: time1,time2
!$ time1=omp_get_wtime()
level=1
error(:,:,1)=SF(:,:)
 call compute_RHS(1)
 call compute_reference_residual(1)

if (non_Newtonian.eq.0) then !!if Newtonian rheology, these quantities do not need to be computed as often
  stream_tstep=1 !!compute dt_stream in test_eigenvalues subroutine
  diff_adv=1     !!time steps depend on diffusive and advective terms
  diffusion=0    !!do not compute dt_stream_diff
  static=1
  call compute_viscosity(level)
end if

  dir=0
  p=1
  q=1
  do
   call iterate_stream(level,Nmulti,dir,p) !! fine grid viscosity is updated in this subroutine for non-Newtonian flows
   write(*,'(i5,3(g20.8))') q,residual_mag(level),res_ref(level),sum(vis_grid(:,:,1))
   if ((residual_mag(1)/res_ref(1)).le.tolerance) exit
   q=q+1
   if (q.gt.max_cycles) then
    write(*,*) "Convergence Problem: Too Many Iterative Cycles Needed"
    call test_snapshot('a',fp,0,nx_grid(1),0,nz_grid(1),vis_grid(0:nx_grid(1),0:nz_grid(1),1)) !!!debug snapshots
    call test_snapshot('c',fp,0,nx_grid(1),0,nz_grid(1),strain(0:nx_grid(1),0:nz_grid(1),1))
    call print_restart
    call snapshots
    stop
   end if
  end do

SF(:,:)=error(:,:,1)
write(*,'(2(i5),2(g20.8))') tstep,q,residual_mag(1),res_ref(1)
!$ time2=omp_get_wtime()
!$ tmultigrid=tmultigrid+(time2-time1)
end

subroutine multigrid
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: level,dir,multiply,INFO,q,p,ii,lvl,i,rank,ip,kp
real*8 :: ROWCND,COLCND,AMAX,diff,res_store,time1,time2,t1,t2
real*8 :: res_coarse(1:ngrid_coarse),res_coarse_mag,res_coarse_ref,rc_store
real*8 :: error_store(-span:nx+span,-span:nz+span),rpart(0:Nthreads-1),rpartB(0:Nthreads-1)
!$ time1=omp_get_wtime()
error(:,:,1)=SF(:,:)
 call compute_RHS(1)
 call compute_reference_residual(1)

if (non_Newtonian.eq.0) then !!if Newtonian rheology, these quantities do not need to be computed as often
 do level=1,grids
  stream_tstep=1 !!compute dt_stream/dt_stream_diff in test_eigenvalues subroutine
  diff_adv=1     !!compute dt_stream based on diffusive and advective terms
  if (splitting.eq.0) then !!compute dt_stream_diff??
   diffusion=0
  else
   diffusion=1
  end if 
  static=1
  call compute_viscosity(level)
 end do
 if ((time.eq.0.d0).or.(ivis.gt.0)) then
  call Compute_Matrix_coarse
!$ t1=omp_get_wtime()
  call DGBTRF(ngrid_coarse,ngrid_coarse,kl,ku,Matrix_coarse,2*kl+ku+1,IPIV,INFO) !!compute LU factors
!$ t2=omp_get_wtime()
!$ t_lu=t_lu+(t2-t1)
 end if
fp=0 !!!for compute_viscosity subroutine
static=1 !!for compute_viscosity subroutine
end if

q=1
do
qcycle=q
 do p=1,nsteps
  level=multi(p)
  qlevel=level
  if (multi(p+1).gt.multi(p)) then
   dir=0
  else
   dir=1 
  end if
  if (multi(p-1).lt.multi(p)) then
   call compute_RHS(level)
   if ((p.le.grids).and.(q.eq.1)) then
    call compute_reference_residual(level)
   end if
  end if
  if (level.lt.grids) then
   if (p.eq.1) then
    if (q.eq.1) then
     res_store=res_ref(p)
    else
     res_store=residual_mag(p)
    end if
   end if
   call iterate_stream(level,Nmulti,dir,p) !! fine grid viscosity is updated in this subroutine for non-Newtonian flows
   if (p.eq.1) then
    !!!!!!!!!!!!!!!!!!!! For non-Newtonian rheology
    if (non_Newtonian.eq.1) then
     do lvl=2,grids-1 !!compute viscosity and viscosity dependent local time steps for each COARSE multigrid level
      call compute_viscosity(lvl)
     end do
    end if
    !!!!!!!!!!!!!!!!!!!! End for non-Newtoninan rheology
   end if
   if (level.eq.1) then
    write(*,'(i5,3(g20.8))') q,residual_mag(level),res_ref(level),sum(vis_grid(:,:,1))
   end if
  else
   !!!!!!!!!!!!!!!set up linear system
   if (non_Newtonian.eq.0) then
    call Compute_coarse_vector
     !$ t1=omp_get_wtime()
    call DGBTRS('N',ngrid_coarse,kl,ku,1,Matrix_coarse,2*kl+ku+1,IPIV,B_coarse,ngrid_coarse,INFO)
     !$ t2=omp_get_wtime()
     !$ t_fs=t_fs+(t2-t1)
    call Stream_grid_coarse
    call enforceBCs_stream(level,error(:,:,level))
!    call compute_residual(level,error(level,:,:),p)
   !!!!!!!!!!!!!!!end set up linear system
!   write(*,*) q,level,p,residual_mag(p),res_ref(level),residual_mag(p)/res_ref(level)
   end if
   !!!!!!!!!!!!!!!*******************************************************************!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!for non-linear multigrid
   if (non_Newtonian.eq.1) then
    do fp=0,fp_max
     if (fp.eq.0) then
      call smoother_stream(level-1,1) !!!smooth fine stream function before restricting
      call restrict(level,stream_smooth(1:nx-1,1:nz-1,level-1),stream_restrict(1:nx-1,1:nz-1,level))!!!restrict stream function of next finest grid
      call enforceBCs_stream(level,stream_restrict(:,:,level))
      error(:,:,level)=stream_restrict(:,:,level) !!!for input to smoother_stream subroutine
      call smoother_stream(level,0) !!!filter unresolved frequencies on the coarse grid
      stream_restrict(:,:,level)=stream_smooth(:,:,level) !!put filter result into stream_restrict array
      error(:,:,level)=stream_restrict(:,:,level) !!!use (smoothed) restricted stream function to compute coarse viscosity      
      if (qcycle.eq.1) then 
       static=1  !!!compute static portion of the coarse grid viscosity
      else
       static=0
      end if
      call stream_derivatives(level,error(:,:,level))
      call compute_viscosity(grids)
      call Compute_Matrix_coarse      !!!LAPACK row scaling built in
      !$ t1=omp_get_wtime()
      call fake_matrix_multiply
      !$ t2=omp_get_wtime()
      !$ t_mult=t_mult+(t2-t1)
      call Compute_coarse_vector !!compute B_coarse_save (without LAPACK scaling)
     end if
     error_store(:,:)=error(:,:,level)
     !$ t1=omp_get_wtime()
     call DGBTRF(ngrid_coarse,ngrid_coarse,kl,ku,Matrix_coarse,2*kl+ku+1,IPIV,INFO) !!compute LU factors
     !$ t2=omp_get_wtime()
     !$ t_lu=t_lu+(t2-t1)
     !$OMP PARALLEL PRIVATE(i)
     !$OMP DO
     do i=1,ngrid_coarse
      B_coarse(i)=B_coarse_save(i)*RowB(i)  !!reuse saved B_coarse but apply new LAPACK scalings
     end do
     !$OMP END DO
     !$OMP END PARALLEL
     !$ t1=omp_get_wtime()
     call DGBTRS('N',ngrid_coarse,kl,ku,1,Matrix_coarse,2*kl+ku+1,IPIV,B_coarse,ngrid_coarse,INFO)
     !$ t2=omp_get_wtime()
     !$ t_fs=t_fs+(t2-t1)
     call Stream_grid_coarse !!put solution (temporarily in B_coarse) in error(:,:,level) array
     call enforceBCs_stream(level,error(:,:,level))
     if ((fp_max.eq.0).or.((residual_mag(1)/res_ref(1)).lt.1.d-2).or.(residual_mag(1).ge.res_store).or.(q.gt.3)) then !!!FP iterations only used when fine grid residual is large compared to the reference value
      fp_save=0
      exit !!!!!!!!!!!!!!!!!!!!!!!!!!!!remaining FP steps are not needed in this case
     end if

     static=0 !!do not recompute static portion of the coarse grid viscosity
     call stream_derivatives(level,error(:,:,level))
     call compute_viscosity(grids)  !!!restrict the fine grid viscosity to get the next guess for the coarse grid viscosity
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!check residual of coarse grid equation
     call Compute_Matrix_coarse    !!!compute new coarse matrix (and row scalings)
     if (fp.eq.0) then
      rc_store=1.d99
     else
      rc_store=res_coarse_mag
     end if
     !$ t1=omp_get_wtime()
     call fake_matrix_multiply
     !$ t2=omp_get_wtime()
     !$ t_mult=t_mult+(t2-t1)
     if (fp.eq.0) then !!!!!!!!!!!!!!compute coarse ref residual on first FP iteration
      rpartB=0.d0
      !$OMP PARALLEL PRIVATE(i,rank)
      rank=omp_get_thread_num()
      !$OMP DO
      do i=1,ngrid_coarse
       rpartB(rank)=rpartB(rank)+dabs(B_coarse_save(i))
      end do
      !$OMP END DO
      !$OMP END PARALLEL
      res_coarse_ref=sum(rpartB)/real(nx_grid(level)-1,8)/real(nz_grid(level)-1,8)
     end if
     rpart=0.d0
     !$OMP PARALLEL PRIVATE(i,rank,ip,kp)
     rank=omp_get_thread_num()
     !$OMP DO
     do i=1,ngrid_coarse
      call indices_coarse(i,ip,kp)
      res_coarse(i)=B_coarse_save(i)-fake(ip,kp)
      rpart(rank)=rpart(rank)+dabs(res_coarse(i))
     end do
     !$OMP END DO
     !$OMP END PARALLEL
      res_coarse_mag=sum(rpart)/real(nx_grid(level)-1,8)/real(nz_grid(level)-1,8)
     write(*,*) fp,res_coarse_mag,res_coarse_ref,sum(vis_grid(:,:,grids))
     if ((res_coarse_mag/res_coarse_ref.lt.tol_coarse).or.(fp.eq.fp_max).or.(res_coarse_mag.gt.rc_store)) then !!exit criterion for FP iterations
      if (res_coarse_mag.gt.rc_store) then
       error(:,:,level)=error_store(:,:) !!use best FP iterate available
      end if
      fp_save=fp
      exit
     end if
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end check residual of coarse grid equation
    end do
    fp=0
   end if
   !!!!!!!!!!!!!!!*******************************************************************!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end for non-linear multigrid
  end if
  if (level.eq.1) then
    if ((residual_mag(1)/res_ref(1)).le.tolerance) exit !!final viscosity and residual calculated at the end of iterate_stream subroutine
  end if
  if (multi(p+1).lt.multi(p)) then
   call prolong(level)
  end if
 end do
 if ((residual_mag(1)/res_ref(1)).le.tolerance) then
  exit
 end if
!!!!!!!!!!!!!!!!!!!!!!!end multigrid cycle
q=q+1
 if (q.gt.max_cycles) then
  write(*,*) "Multigrid Convergence Problem: Too Many Multigrid Cycles Needed"
  call test_snapshot('a',fp,0,nx_grid(1),0,nz_grid(1),vis_grid(0:nx_grid(1),0:nz_grid(1),1)) !!!debug snapshots
  call test_snapshot('b',fp,0,nx_grid(2),0,nz_grid(2),vis_grid(0:nx_grid(2),0:nz_grid(2),2))
  call test_snapshot('c',fp,0,nx_grid(1),0,nz_grid(1),strain(0:nx_grid(1),0:nz_grid(1),1))
  call test_snapshot('d',fp,0,nx_grid(2),0,nz_grid(2),strain(0:nx_grid(2),0:nz_grid(2),2))
  call print_restart
  call snapshots
  stop
 end if
end do
 qcycle=1
 SF(:,:)=error(:,:,1)
 write(*,'(2(i5),2(g20.8))') tstep,q,residual_mag(1),res_ref(1)
!$ time2=omp_get_wtime()
!$ tmultigrid=tmultigrid+(time2-time1)
end

subroutine prolong(level_coarse)
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: level_coarse,i,k,lc,lf
lc=level_coarse
lf=level_coarse-1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!non-linear multigrid correction
if (non_Newtonian.eq.1) then
!$OMP PARALLEL DO PRIVATE(i,k)
do k=1,nz_grid(lc)-1
 do i=1,nx_grid(lc)-1 !!prolong
  error(i,k,lc)=(error(i,k,lc)-stream_restrict(i,k,lc))
 end do
end do
!$OMP END PARALLEL DO
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end non-linear multigrid correction

!$OMP PARALLEL PRIVATE(i,k)
!$OMP DO 
do k=1,nz_grid(lc)-1
 do i=1,nx_grid(lc)-1 !!prolong
   error(2*i,2*k,lf)=error(2*i,2*k,lf)+error(i,k,lc)              !!bilinear
 end do
end do
!$OMP END DO
!$OMP DO
do k=1,nz_grid(lf)-1,2
 do i=2,nx_grid(lf)-2,2
   error(i,k,lf)=error(i,k,lf)+(error(i/2,(k-1)/2,lc)+error(i/2,(k+1)/2,lc))/2.d0     !!bilinear
 end do
end do
!$OMP END DO
!$OMP DO
do k=2,nz_grid(lf)-2,2
 do i=1,nx_grid(lf)-1,2
   error(i,k,lf)=error(i,k,lf)+(error((i-1)/2,k/2,lc)+error((i+1)/2,k/2,lc))/2.d0     !!bilinear
 end do
end do
!$OMP END DO
!$OMP DO
do k=1,nz_grid(lf)-1,2
 do i=1,nx_grid(lf)-1,2
   error(i,k,lf)=error(i,k,lf)+&
&(error((i-1)/2,(k-1)/2,lc)+error((i+1)/2,(k-1)/2,lc)+error((i+1)/2,(k+1)/2,lc)+error((i-1)/2,(k+1)/2,lc))/4.d0 !!bilinear
 end do
end do
!$OMP END DO
!$OMP END PARALLEL
 call enforceBCs_stream(lf,error(:,:,lf))
end

subroutine test_snapshot(f_type,f_count,i1,i2,k1,k2,array) !!for debugging purposes
implicit none
character*1 :: f_type !!e.g., "a" for file type a
character*3 :: fname !!e.g., "a1" for type a and f_count=1
integer*4 :: i,k,i1,i2,k1,k2,f_count
real*8 :: array(i1:i2,k1:k2)
write(fname,'(a,i2)') f_type,f_count
open(unit=999,file=fname)
do i=i1,i2
 do k=k1,k2
  write(999,'(i5,i5,1(g14.5E3))') i,k,array(i,k)
 end do
end do
 close(999)
end


subroutine compute_residual(level,sol,p)
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k,level,p,iq,kq,rank
real*8 ::sol(-span:nx+span,-span:nz+span),stream_space,rpart(0:Nthreads-1)
! call stream_derivatives(level,sol) !!already called - no need to repeat
rpart=0.d0
!$OMP PARALLEL PRIVATE(i,k,rank)
rank=omp_get_thread_num()
!$OMP DO
do k=1,nz_grid(level)-1
 do i=1,nx_grid(level)-1
  residual(i,k,level)=stream_space(i,k,sol,level)
  rpart(rank)=rpart(rank)+dabs(residual(i,k,level))
 end do
end do
!$OMP END DO
!$OMP END PARALLEL
residual_mag(p)=sum(rpart)/real(nx_grid(level)-2,8)/real(nz_grid(level)-2,8)

!write(*,*) level,residual_mag(p)
if ((level.eq.1).and.(residual_mag(p).gt.(1.d10*res_ref(level)))) then
 write(*,*) "Multigrid Convergence Problem: Abnormally High Residual Detected -- stopping"
 call test_snapshot('a',fp,0,nx_grid(1),0,nz_grid(1),vis_grid(0:nx_grid(1),0:nz_grid(1),1)) !!debug snapshots
 call test_snapshot('b',fp,0,nx_grid(2),0,nz_grid(2),vis_grid(0:nx_grid(2),0:nz_grid(2),2))
 call test_snapshot('c',fp,0,nx_grid(1),0,nz_grid(1),strain(0:nx_grid(1),0:nz_grid(1),1))
 call test_snapshot('d',fp,0,nx_grid(2),0,nz_grid(2),strain(0:nx_grid(2),0:nz_grid(2),2))
 call print_restart
 call snapshots
 stop
end if
end

subroutine compute_reference_residual(level)
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: level,i,k,rank
real*8 :: part(0:Nthreads-1)
!res_ref(level)=sum(dabs(RHS(1:nx_grid(level)-1,1:nz_grid(level)-1,level)))/real(nx_grid(level)-2,8)/real(nz_grid(level)-2,8)
part=0.d0
!$OMP PARALLEL PRIVATE (i,k,rank)
rank=omp_get_thread_num()
!$OMP DO
do k=1,nz_grid(level)-1
 do i=1,nx_grid(level)-1
  part(rank)=part(rank)+dabs(RHS(i,k,level))
 end do
end do
!$OMP END DO
!$OMP END PARALLEL
res_ref(level)=sum(part)/real(nx_grid(level)-2,8)/real(nz_grid(level)-2,8)
end

subroutine viscosity_gradients(level,n)
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k,level,ii,kk,n
real*8 :: visf,vis_r,vis_rr,vis_s,vis_ss,vis_rs

if (ivis.eq.0) then
   vis_x(:,:,level)=0.d0
   vis_z(:,:,level)=0.d0
   vis_xx(:,:,level)=0.d0
   vis_zz(:,:,level)=0.d0
   vis_xz(:,:,level)=0.d0
   vis_xxpzz(:,:,level)=0.d0
   vis_xxmzz(:,:,level)=0.d0
elseif ((gridX.eq.0).and.(gridZ.eq.0)) then !!uniform grid

 if (n.eq.0) then

!$OMP PARALLEL DO PRIVATE(i,k)
do k=0,nz_grid(level)
 do i=0,nx_grid(level)
   vis_xx(i,k,level)=dot_product(D2_so(-1:1),vis_grid(i-1:i+1,k,level))/dr2_grid(level)
   vis_zz(i,k,level)=dot_product(D2_so(-1:1),vis_grid(i,k-1:k+1,level))/ds2_grid(level)
   vis_xxpzz(i,k,level)=vis_xx(i,k,level)+vis_zz(i,k,level)
   vis_smooth(i,k,level)=dt_vis(i,k,level)*vis_xxpzz(i,k,level)
 end do
end do
!$OMP END PARALLEL DO
 return !!!!!!!!!!!!!!!!!!!!!only need gradients for viscosity smoother

 else !!also include other gradients needed by eigR,eigI,stream_space,stream_space_matrix

!$OMP PARALLEL DO PRIVATE(i,k)
do k=0,nz_grid(level)
 do i=0,nx_grid(level)
   vis_xz(i,k,level)=sum(Drs_so(-1:1,-1:1)*vis_grid(i-1:i+1,k-1:k+1,level))/drds_grid(level)
   vis_x(i,k,level)=dot_product(D1_so(-1:1),vis_grid(i-1:i+1,k,level))/dr_grid(level)
   vis_xx(i,k,level)=dot_product(D2_so(-1:1),vis_grid(i-1:i+1,k,level))/dr2_grid(level)
   vis_z(i,k,level)=dot_product(D1_so(-1:1),vis_grid(i,k-1:k+1,level))/ds_grid(level)
   vis_zz(i,k,level)=dot_product(D2_so(-1:1),vis_grid(i,k-1:k+1,level))/ds2_grid(level)
   vis_xxpzz(i,k,level)=vis_xx(i,k,level)+vis_zz(i,k,level)
   vis_xxmzz(i,k,level)=vis_xx(i,k,level)-vis_zz(i,k,level)
   vis_smooth(i,k,level)=dt_vis(i,k,level)*vis_xxpzz(i,k,level)
 end do
end do
!$OMP END PARALLEL DO
 end if 

else !!grid refinement -------------update for cache use efficiency

!$OMP PARALLEL DO PRIVATE(i,k,vis_r,vis_rr,vis_s,vis_ss,vis_rs,ii,kk)
do k=0,nz_grid(level)
  kk=k*2**(level-1)
 do i=0,nx_grid(level)
  ii=i*2**(level-1)
   vis_r=dot_product(D1_so(-1:1),vis_grid(i-1:i+1,k,level))/dr_grid(level)       !!second order FD stencils for viscosity gradients
   vis_rr=dot_product(D2_so(-1:1),vis_grid(i-1:i+1,k,level))/dr2_grid(level)
   vis_s=dot_product(D1_so(-1:1),vis_grid(i,k-1:k+1,level))/ds_grid(level)
   vis_ss=dot_product(D2_so(-1:1),vis_grid(i,k-1:k+1,level))/ds2_grid(level)

   vis_x(i,k,level)=vis_r/xr(ii)
   vis_z(i,k,level)=vis_s/zs(kk)
   vis_xx(i,k,level)=(vis_rr-vis_x(i,k,level)*xrr(ii))/xr(ii)**2.d0
   vis_zz(i,k,level)=(vis_ss-vis_z(i,k,level)*zss(kk))/zs(kk)**2.d0
   vis_xxpzz(i,k,level)=vis_xx(i,k,level)+vis_zz(i,k,level)
   vis_smooth(i,k,level)=dt_vis(i,k,level)*vis_xxpzz(i,k,level)
 end do
end do
!$OMP END PARALLEL DO

 if (n.eq.0) return !!!!!!!!!!!!!!!!!!!!!only need gradients for viscosity smoother

!$OMP PARALLEL DO PRIVATE(i,k,vis_r,vis_rr,vis_s,vis_ss,vis_rs,ii,kk) !!other gradients needed by stream_space etc...
do k=0,nz_grid(level)
  kk=k*2**(level-1)
 do i=0,nx_grid(level)
  ii=i*2**(level-1)
   vis_rs=sum(Drs_so(-1:1,-1:1)*vis_grid(i-1:i+1,k-1:k+1,level))/drds_grid(level)
   vis_xz(i,k,level)=vis_rs/xr(ii)/zs(kk)
   vis_xxmzz(i,k,level)=vis_xx(i,k,level)-vis_zz(i,k,level)
 end do
end do
!$OMP END PARALLEL DO
end if
end

subroutine compute_viscosity(level)
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k,level,ii,kk,q,i1,k1,interpolate,p
real*8 :: visf,t1,t2
if (level.eq.1) then

!$ t1=omp_get_wtime() !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! A
if (static.eq.1) then
!$OMP PARALLEL DO PRIVATE(i,k)
do k=0,nz_grid(level)
 do i=0,nx_grid(level)
  vis_grid_static(i,k,level)=visf(i,k) !!compute T and C dependent viscosity
 end do
end do
!$OMP END PARALLEL DO
end if
!$ t2=omp_get_wtime()
!$ t_cvA=t_cvA+(t2-t1) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end A

!$ t1=omp_get_wtime() !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! B
if (plastic.ne.0) then !!!!!put in plastic yielding
 call viscosity_loop(level)
else
 vis_grid(:,:,level)=vis_grid_static(:,:,level) !!!!!!!!!!!!!!!!! no plastic yielding
end if
!$ t2=omp_get_wtime()
!$ t_cvB=t_cvB+(t2-t1) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end B

!$ t1=omp_get_wtime() !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! C
  interpolate=0
  call smoother_vis(level,interpolate)  !!smooth computational viscosity to kill unresolved frequencies
!$ t2=omp_get_wtime()
!$ t_cvC=t_cvC+(t2-t1) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end C
else
  vis_gridf(:,:,level-1)=vis_grid(:,:,level-1) !!!store original viscosity
  vis_xf(:,:,level-1)=vis_x(:,:,level-1)
  vis_zf(:,:,level-1)=vis_z(:,:,level-1)
  vis_xxf(:,:,level-1)=vis_xx(:,:,level-1)
  vis_zzf(:,:,level-1)=vis_zz(:,:,level-1)
  vis_xzf(:,:,level-1)=vis_xz(:,:,level-1)
  vis_xxpzzf(:,:,level-1)=vis_xxpzz(:,:,level-1)
  vis_xxmzzf(:,:,level-1)=vis_xxmzz(:,:,level-1)
  vis_smoothf(:,:,level-1)=vis_smooth(:,:,level-1)

  if (static.eq.1) then !!!! if static portion of viscosity has changed then compute it 
    vis_grid(:,:,level-1)=vis_grid_static(:,:,level-1) !!!!!!!!!!!!!!!!!new
    interpolate=1
    call smoother_vis(level-1,interpolate) !!!!!smooth viscosity field before interpolating to coarse grid
    call restrict_coefficients(level,vis_grid(1:nx-1,1:nz-1,level-1),vis_grid(1:nx-1,1:nz-1,level))
    call viscosityBCs_initialize(level)
    vis_grid_static(:,:,level)=vis_grid(:,:,level) !!!store restricted static viscosity (on first multigrid cycle) ----- new
  end if

  if (non_Newtonian.ne.0) then    !!!!put in plastic yielding parts (from stream_restrict array)
   call compute_strain_rate_invariant_quick(level)
   call viscosity_loop(level)
  end if

    interpolate=0
    call smoother_vis(level,interpolate)  !!smooth computational viscosity to remove high frequencies

  vis_grid(:,:,level-1)=vis_gridf(:,:,level-1) !!restore original viscosity
  vis_x(:,:,level-1)=vis_xf(:,:,level-1)
  vis_z(:,:,level-1)=vis_zf(:,:,level-1)
  vis_xx(:,:,level-1)=vis_xxf(:,:,level-1)
  vis_zz(:,:,level-1)=vis_zzf(:,:,level-1)
  vis_xz(:,:,level-1)=vis_xzf(:,:,level-1)
  vis_xxpzz(:,:,level-1)=vis_xxpzzf(:,:,level-1)
  vis_xxmzz(:,:,level-1)=vis_xxmzzf(:,:,level-1)
  vis_smooth(:,:,level-1)=vis_smoothf(:,:,level-1)
end if
end

subroutine viscosityBCs_initialize(level) !!define inital boundary values on coarse grid levels
 use basics
 use arrays
implicit none
integer*4 :: i,k,ii,kk,level,i1,k1,tt
real*8 :: Cin(1:ntype),Bin

do k=1,nz_grid(level)-1 !!vertical boundaries
 kk=2*k !!fine grid
 vis_grid(0,k,level)=vis_grid(0,kk,level-1)
 vis_grid(nx_grid(level),k,level)=vis_grid(nx_grid(level-1),kk,level-1)
end do

do i=0,nx_grid(level) !!horizontal boundaries
 ii=2*i !!fine grid
 vis_grid(i,0,level)=vis_grid(ii,0,level-1)
 vis_grid(i,nz_grid(level),level)=vis_grid(ii,nz_grid(level-1),level-1)
end do

do i=1,span !!!side ghost points
 vis_grid(-i,0:nz_grid(level),level)=vis_grid(i,0:nz_grid(level),level)
 vis_grid(nx_grid(level)+i,0:nz_grid(level),level)=vis_grid(nx_grid(level)-i,0:nz_grid(level),level)
end do

do k=1,span !!!top/bottom ghost points
 vis_grid(-span:nx_grid(level)+span,-k,level)=vis_grid(-span:nx_grid(level)+span,k,level)
 vis_grid(-span:nx_grid(level)+span,nz_grid(level)+k,level)=vis_grid(-span:nx_grid(level)+span,nz_grid(level)-k,level)
end do
end

subroutine viscosityBCs(level)
 use basics
 use arrays
implicit none
integer*4 :: i,k,ii,kk,level

do k=-span,nz_grid(level)+span
 do i=1,span !!!side ghost points: insulating
  vis_grid(-i,k,level)=vis_grid(i,k,level)
 end do
end do
do k=-span,nz_grid(level)+span
 do i=1,span !!!side ghost points: insulating
  vis_grid(nx_grid(level)+i,k,level)=vis_grid(nx_grid(level)-i,k,level)
 end do
end do

do k=1,span !!!top/bottom ghost points
 do i=0,nx_grid(level)
  vis_grid(i,-k,level)=vis_grid(i,k,level)
  vis_grid(i,nz_grid(level)+k,level)=vis_grid(i,nz_grid(level)-k,level)
 end do
end do

end

subroutine stream_space_matrix_prep !!!!does not currently function for grid refinement
 use basics
 use arrays
implicit none
integer*4 :: i,j,k,di,dk
real*8 :: Sr,Ss,Srr,Sss,Srrr,Ssss,Srrrr,Sssss
real*8 :: Sx,Sz,Sxx,Szz,Sxxx,Szzz,Sxxxx,Szzzz
real*8 :: Srs,Srrs,Srss,Srrss
real*8 :: Sxz,Sxxzz,Sxxz,Sxzz

do dk=-span,span
 do di=-span,span
  if (dk.eq.0) then
   Sxx=D2(di)/dr2_grid(grids)
   Sxxx=D3(di)/dr3_grid(grids)
   Sxxxx=D4(di)/dr4_grid(grids)
  else
   Sxx=0.d0
   Sxxx=0.d0
   Sxxxx=0.d0
  end if

  if (di.eq.0) then
   Szz=D2(dk)/ds2_grid(grids)
   Szzz=D3(dk)/ds3_grid(grids)
   Szzzz=D4(dk)/ds4_grid(grids)
  else
   Szz=0.d0
   Szzz=0.d0
   Szzzz=0.d0
  end if

  Sxz=Drs(di,dk)/drds_grid(grids)
  Sxzz=Drss(di,dk)/drds2_grid(grids)
  Sxxz=Drrs(di,dk)/dr2ds_grid(grids)
  Sxxzz=Drrss(di,dk)/dr2ds2_grid(grids)

  SM1(di,dk)=Sxx-Szz
  SM2(di,dk)=Sxxxx+Szzzz+2.d0*Sxxzz
  SM4(di,dk)=4.d0*Sxz
  SM5(di,dk)=2.d0*(Sxxx+Sxzz)
  SM6(di,dk)=2.d0*(Szzz+Sxxz)
 end do
end do

! if (ivis.eq.0) then
!  stream_space_matrix=vis_grid(i,k,level)*(Sxxxx+Szzzz+2.d0*Sxxzz)
! else
!  stream_space_matrix=(vis_xxmzz(i,k,level)*(Sxx-Szz))+(vis_grid(i,k,level)*(Sxxxx+Szzzz+2.d0*Sxxzz))+&
!                     &(4.d0*vis_xz(i,k,level)*Sxz)+(2.d0*vis_x(i,k,level)*(Sxxx+Sxzz))+(2.d0*vis_z(i,k,level)*(Szzz+Sxxz))
! end if
end

real*8 function stream_space_matrix(i,k,level,iq,kq,Sval) !!compute coarse matrix elements
 use basics
 use arrays
implicit none
integer*4 :: i,j,k,level,i1,k1,di,dk,iq,kq,Sval
real*8 :: Sr,Ss,Srr,Sss,Srrr,Ssss,Srrrr,Sssss
real*8 :: Sx,Sz,Sxx,Szz,Sxxx,Szzz,Sxxxx,Szzzz
real*8 :: Srs,Srrs,Srss,Srrss
real*8 :: Sxz,Sxxzz,Sxxz,Sxzz

if ((gridX.eq.0).and.(gridZ.eq.0)) then !!!uniform grid
 di=iq-i
 dk=kq-k

 if (ivis.eq.0) then
  stream_space_matrix=vis_grid(i,k,level)*SM2(di,dk)
 else
  stream_space_matrix=vis_xxmzz(i,k,level)*SM1(di,dk)+vis_grid(i,k,level)*SM2(di,dk)+&
                     &vis_xz(i,k,level)*SM4(di,dk)+vis_x(i,k,level)*(SM5(di,dk))+vis_z(i,k,level)*SM6(di,dk)
 end if
else !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!grid refinement

i1=i*2**(level-1)  !!!fine grid coordinates
k1=k*2**(level-1)

 di=iq-i
 dk=kq-k

 if (k.eq.kq) then
  Sr=D1(di)/dr_grid(level)
  Srr=D2(di)/dr2_grid(level)
  Srrr=D3(di)/dr3_grid(level)
  Srrrr=D4(di)/dr4_grid(level)
 else
  Sr=0.d0
  Srr=0.d0
  Srrr=0.d0
  Srrrr=0.d0
 end if

 if (i.eq.iq) then
  Ss=D1(dk)/ds_grid(level)
  Sss=D2(dk)/ds2_grid(level)
  Ssss=D3(dk)/ds3_grid(level)
  Sssss=D4(dk)/ds4_grid(level)
 else
  Ss=0.d0
  Sss=0.d0
  Ssss=0.d0
  Sssss=0.d0
 end if

 Srs=Drs(di,dk)/drds_grid(level)
 Srss=Drss(di,dk)/drds2_grid(level)
 Srrs=Drrs(di,dk)/dr2ds_grid(level)
 Srrss=Drrss(di,dk)/dr2ds2_grid(level)

 Sx=Sr/xr(i1)
 Sz=Ss/zs(k1)
 Sxx=(Srr-Sx*xrr(i1))/xr(i1)**2.d0
 Szz=(Sss-Sz*zss(k1))/zs(k1)**2.d0
 Sxxx=(Srrr-3.d0*Sxx*xr(i1)*xrr(i1)-Sx*xrrr(i1))/xr(i1)**3.d0
 Szzz=(Ssss-3.d0*Szz*zs(k1)*zss(k1)-Sz*zsss(k1))/zs(k1)**3.d0
 Sxxxx=(Srrrr-6.d0*Sxxx*xrr(i1)*xr(i1)**2.d0-(3.d0*xrr(i1)**2.d0+4.d0*xr(i1)*xrrr(i1))*Sxx-Sx*xrrrr(i1))/xr(i1)**4.d0
 Szzzz=(Sssss-6.d0*Szzz*zss(k1)*zs(k1)**2.d0-(3.d0*zss(k1)**2.d0+4.d0*zs(k1)*zsss(k1))*Szz-Sz*zssss(k1))/zs(k1)**4.d0

 Sxz=Srs/xr(i1)/zs(k1)
 Sxxz=(Srrs-Sxz*zs(k1)*xrr(i1))/zs(k1)/xr(i1)**2.d0
 Sxzz=(Srss-Sxz*xr(i1)*zss(k1))/xr(i1)/zs(k1)**2.d0
 Sxxzz=(Srrss-Sxzz*xrr(i1)*zs(k1)**2.d0-Sxxz*zss(k1)*xr(i1)**2.d0-Sxz*xrr(i1)*zss(k1))/(xr(i1)*zs(k1))**2.d0

 if (ivis.eq.0) then
  stream_space_matrix=vis_grid(i,k,level)*(Sxxxx+Szzzz+2.d0*Sxxzz)
 else
  stream_space_matrix=(vis_xxmzz(i,k,level)*(Sxx-Szz))+(vis_grid(i,k,level)*(Sxxxx+Szzzz+2.d0*Sxxzz))+&
                     &(4.d0*vis_xz(i,k,level)*Sxz)+(2.d0*vis_x(i,k,level)*(Sxxx+Sxzz))+(2.d0*vis_z(i,k,level)*(Szzz+Sxxz))
 end if
end if

if (Sval.eq.-1) stream_space_matrix=-stream_space_matrix
end

real*8 function stream_space_full(i,k,S0,level) !!includes stream function derivatives
 use basics
 use arrays
implicit none
integer*4 :: i,k,level,i1,k1
real*8 :: Sr,Ss,Srr,Sss,Srrr,Ssss,Srrrr,Sssss
real*8 :: Sx,Sz,Sxx,Szz,Sxxx,Szzz,Sxxxx,Szzzz
real*8 :: Srs,Srrs,Srss,Srrss
real*8 :: Sxz,Sxxzz,Sxxz,Sxzz
real*8 :: S0(-span:nx+span,-span:nz+span)
real*8 :: term1,term2,term3,term4,term5,term6

i1=i*2**(level-1)  !!!fine grid coordinates
k1=k*2**(level-1)

 Sr=dot_product(D1(-span1:span1),S0(i-span1:i+span1,k))/dr_grid(level)
 Ss=dot_product(D1(-span1:span1),S0(i,k-span1:k+span1))/ds_grid(level)
 Srr=dot_product(D2(-span1:span1),S0(i-span1:i+span1,k))/dr2_grid(level)
 Sss=dot_product(D2(-span1:span1),S0(i,k-span1:k+span1))/ds2_grid(level)
 Srrr=dot_product(D3(-span:span),S0(i-span:i+span,k))/dr3_grid(level)
 Ssss=dot_product(D3(-span:span),S0(i,k-span:k+span))/ds3_grid(level)
 Srrrr=dot_product(D4(-span:span),S0(i-span:i+span,k))/dr4_grid(level)
 Sssss=dot_product(D4(-span:span),S0(i,k-span:k+span))/ds4_grid(level)

 Srs=sum(Drs(-span1:span1,-span1:span1)*S0(i-span1:i+span1,k-span1:k+span1))/drds_grid(level)
 Srss=sum(Drss(-span1:span1,-span1:span1)*S0(i-span1:i+span1,k-span1:k+span1))/drds2_grid(level)
 Srrs=sum(Drrs(-span1:span1,-span1:span1)*S0(i-span1:i+span1,k-span1:k+span1))/dr2ds_grid(level)
 Srrss=sum(Drrss(-span1:span1,-span1:span1)*S0(i-span1:i+span1,k-span1:k+span1))/dr2ds2_grid(level)

 Sx=Sr/xr(i1)
 Sz=Ss/zs(k1)
 Sxx=(Srr-Sx*xrr(i1))/xr(i1)**2.d0
 Szz=(Sss-Sz*zss(k1))/zs(k1)**2.d0
 Sxxx=(Srrr-3.d0*Sxx*xr(i1)*xrr(i1)-Sx*xrrr(i1))/xr(i1)**3.d0
 Szzz=(Ssss-3.d0*Szz*zs(k1)*zss(k1)-Sz*zsss(k1))/zs(k1)**3.d0
 Sxxxx=(Srrrr-6.d0*Sxxx*xrr(i1)*xr(i1)**2.d0-(3.d0*xrr(i1)**2.d0+4.d0*xr(i1)*xrrr(i1))*Sxx-Sx*xrrrr(i1))/xr(i1)**4.d0
 Szzzz=(Sssss-6.d0*Szzz*zss(k1)*zs(k1)**2.d0-(3.d0*zss(k1)**2.d0+4.d0*zs(k1)*zsss(k1))*Szz-Sz*zssss(k1))/zs(k1)**4.d0

 Sxz=Srs/xr(i1)/zs(k1)
 Sxxz=(Srrs-Sxz*zs(k1)*xrr(i1))/zs(k1)/xr(i1)**2.d0
 Sxzz=(Srss-Sxz*xr(i1)*zss(k1))/xr(i1)/zs(k1)**2.d0
 Sxxzz=(Srrss-Sxzz*xrr(i1)*zs(k1)**2.d0-Sxxz*zss(k1)*xr(i1)**2.d0-Sxz*xrr(i1)*zss(k1))/(xr(i1)*zs(k1))**2.d0

if (ivis.eq.0) then
 term1=0.d0
 term2=vis_grid(i,k,level)*(Sxxxx+Szzzz)
 term3=2.d0*vis_grid(i,k,level)*Sxxzz
 term4=0.d0
 term5=0.d0
 term6=0.d0
else
 term1=(vis_xx(i,k,level)-vis_zz(i,k,level))*(Sxx-Szz)
 term2=vis_grid(i,k,level)*(Sxxxx+Szzzz)
 term3=2.d0*vis_grid(i,k,level)*Sxxzz
 term4=4.d0*vis_xz(i,k,level)*Sxz
 term5=2.d0*vis_x(i,k,level)*(Sxxx+Sxzz)
 term6=2.d0*vis_z(i,k,level)*(Szzz+Sxxz)
end if

  stream_space_full=term1+term2+term3+term4+term5+term6-RHS(i,k,level) !!add terms
  stream_space_full=-stream_space_full !!!for convergence (eigenvalues must be dominated by negative real values)
end

real*8 function stream_space(i,k,S0,level)
 use basics
 use arrays
implicit none
integer*4 :: i,k,level
real*8 :: S0(-span:nx+span,-span:nz+span)
real*8 :: term1,term2,term3,term4,term5,term6
if (ivis.eq.0) then
 term2=vis_grid(i,k,level)*(SFxxxx(i,k,level)+SFzzzz(i,k,level)+2.d0*SFxxzz(i,k,level))
 stream_space=term2-RHS(i,k,level) !!combined terms 2 and 3 for speed
else
 term1=vis_xxmzz(i,k,level)*(SFxx(i,k,level)-SFzz(i,k,level))
 term2=vis_grid(i,k,level)*(SFxxxx(i,k,level)+SFzzzz(i,k,level)+2.d0*SFxxzz(i,k,level))
! term2=vis_grid(i,k,level)*(SFxxxx(i,k,level)+SFzzzz(i,k,level))
! term3=2.d0*vis_grid(i,k,level)*SFxxzz(i,k,level)
 term4=4.d0*vis_xz(i,k,level)*SFxz(i,k,level)
 term5=2.d0*vis_x(i,k,level)*(SFxxx(i,k,level)+SFxzz(i,k,level))
 term6=2.d0*vis_z(i,k,level)*(SFzzz(i,k,level)+SFxxz(i,k,level))
 stream_space=term1+term2+term4+term5+term6-RHS(i,k,level) !!combined terms 2 and 3 for speed
end if
  stream_space=-stream_space !!!for convergence (eigenvalues must be dominated by negative real values)
end

real*8 function stream_space_diff(i,k,S0,level) !!!spatial terms for diffusion multigrid smoother
 use basics
 use arrays
implicit none
integer*4 :: i,k,level
real*8 :: S0(-span:nx+span,-span:nz+span)
real*8 :: term1,term2,term3,term4,term5,term6
if (ivis.eq.0) then
 term2=vis_grid(i,k,level)*(SFxxxx(i,k,level)+SFzzzz(i,k,level)+2.d0*SFxxzz(i,k,level))
 stream_space_diff=term2-RHS(i,k,level) !!combined terms 2 and 3 for speed
else
 term1=vis_xxmzz(i,k,level)*(SFxx(i,k,level)-SFzz(i,k,level))
 term2=vis_grid(i,k,level)*(SFxxxx(i,k,level)+SFzzzz(i,k,level)+2.d0*SFxxzz(i,k,level))
 term4=4.d0*vis_xz(i,k,level)*SFxz(i,k,level)
! term5=2.d0*vis_x(i,k,level)*(SFxxx(i,k,level)+SFxzz(i,k,level))
! term6=2.d0*vis_z(i,k,level)*(SFzzz(i,k,level)+SFxxz(i,k,level))
 stream_space_diff=term1+term2+term4-RHS(i,k,level) !!combined terms 2 and 3 for speed
end if
  stream_space_diff=-stream_space_diff !!!for convergence (eigenvalues must be dominated by negative real values)
end

real*8 function stream_space_fake_matrix(i,k,S0,level)
 use basics
 use arrays
implicit none
integer*4 :: i,k,level
real*8 :: S0(-span:nx+span,-span:nz+span)
real*8 :: term1,term2,term3,term4,term5,term6
if (ivis.eq.0) then
 term2=vis_grid(i,k,level)*(SFxxxx(i,k,level)+SFzzzz(i,k,level)+2.d0*SFxxzz(i,k,level))
 stream_space_fake_matrix=term2!-RHS(i,k,level) !!combined terms 2 and 3 for speed
else
 term1=vis_xxmzz(i,k,level)*(SFxx(i,k,level)-SFzz(i,k,level))
 term2=vis_grid(i,k,level)*(SFxxxx(i,k,level)+SFzzzz(i,k,level)+2.d0*SFxxzz(i,k,level))
! term2=vis_grid(i,k,level)*(SFxxxx(i,k,level)+SFzzzz(i,k,level))
! term3=2.d0*vis_grid(i,k,level)*SFxxzz(i,k,level)
 term4=4.d0*vis_xz(i,k,level)*SFxz(i,k,level)
 term5=2.d0*vis_x(i,k,level)*(SFxxx(i,k,level)+SFxzz(i,k,level))
 term6=2.d0*vis_z(i,k,level)*(SFzzz(i,k,level)+SFxxz(i,k,level))
 stream_space_fake_matrix=term1+term2+term4+term5+term6!-RHS(i,k,level) !!combined terms 2 and 3 for speed
end if
end

subroutine fake_matrix_multiply
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k
real*8 :: stream_space_fake_matrix
fake=0.d0
!$OMP PARALLEL DO PRIVATE(i,k)
do k=1,nz_grid(grids)-1
 do i=1,nx_grid(grids)-1
  fake(i,k)=stream_space_fake_matrix(i,k,error(:,:,grids),grids)
 end do
end do
!$OMP END PARALLEL DO
end

subroutine stream_derivatives(level,S0)
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,j,k,level,i1,k1
real*8 :: Sr,Ss,Srr,Sss,Srrr,Ssss,Srrrr,Sssss
real*8 :: Sx,Sz,Sxx,Szz,Sxxx,Szzz,Sxxxx,Szzzz
real*8 :: Srs,Srrs,Srss,Srrss
real*8 :: Sxz,Sxxzz,Sxxz,Sxzz
real*8 :: S0(-span:nx+span,-span:nz+span)

if ((gridX.eq.0).and.(gridZ.eq.0)) then !!uniform grid

!$OMP PARALLEL PRIVATE(i,k)
!$OMP DO
do k=0,nz_grid(level)
 do i=0,nx_grid(level)
  SFxx(i,k,level)=dot_product(D2(-span1:span1),S0(i-span1:i+span1,k))/dr2_grid(level) 
  SFxxx(i,k,level)=dot_product(D3(-span:span),S0(i-span:i+span,k))/dr3_grid(level)
  SFxxxx(i,k,level)=dot_product(D4(-span:span),S0(i-span:i+span,k))/dr4_grid(level)
 end do
end do
!$OMP END DO

!$OMP DO
do k=0,nz_grid(level)
 do i=0,nx_grid(level)
  SFzz(i,k,level)=dot_product(D2(-span1:span1),S0(i,k-span1:k+span1))/ds2_grid(level)
  SFzzz(i,k,level)=dot_product(D3(-span:span),S0(i,k-span:k+span))/ds3_grid(level)
  SFzzzz(i,k,level)=dot_product(D4(-span:span),S0(i,k-span:k+span))/ds4_grid(level)
 end do
end do
!$OMP END DO

!$OMP DO
do k=0,nz_grid(level)
 do i=0,nx_grid(level)
  SFxz(i,k,level)=sum(Drs(-span1:span1,-span1:span1)*S0(i-span1:i+span1,k-span1:k+span1))/drds_grid(level)
  SFxxz(i,k,level)=sum(Drrs(-span1:span1,-span1:span1)*S0(i-span1:i+span1,k-span1:k+span1))/dr2ds_grid(level)
  SFxzz(i,k,level)=sum(Drss(-span1:span1,-span1:span1)*S0(i-span1:i+span1,k-span1:k+span1))/drds2_grid(level)
  SFxxzz(i,k,level)=sum(Drrss(-span1:span1,-span1:span1)*S0(i-span1:i+span1,k-span1:k+span1))/dr2ds2_grid(level)
 end do
end do
!$OMP END DO
!$OMP END PARALLEL

else !!grid refinement ----------- update to improve cache usage

!$OMP PARALLEL DO PRIVATE(i,k,Sr,Ss,Srr,Sss,Srrr,Ssss,Srrrr,Sssss,Srs,Srrs,Srss,Srrss,Sx,Sz,Sxx,Szz,Sxxx,Szzz,Sxxxx,Szzzz,Sxz&
!$OMP&                    ,Sxxzz,Sxxz,Sxzz,i1,k1)
do k=0,nz_grid(level)
k1=k*2**(level-1)
 do i=0,nx_grid(level)
 i1=i*2**(level-1)  !!!fine grid coordinates

 Sr=dot_product(D1(-span1:span1),S0(i-span1:i+span1,k))/dr_grid(level)
 Ss=dot_product(D1(-span1:span1),S0(i,k-span1:k+span1))/ds_grid(level)
 Srr=dot_product(D2(-span1:span1),S0(i-span1:i+span1,k))/dr2_grid(level)
 Sss=dot_product(D2(-span1:span1),S0(i,k-span1:k+span1))/ds2_grid(level)
 Srrr=dot_product(D3(-span:span),S0(i-span:i+span,k))/dr3_grid(level)
 Ssss=dot_product(D3(-span:span),S0(i,k-span:k+span))/ds3_grid(level)
 Srrrr=dot_product(D4(-span:span),S0(i-span:i+span,k))/dr4_grid(level)
 Sssss=dot_product(D4(-span:span),S0(i,k-span:k+span))/ds4_grid(level)

 Srs=sum(Drs(-span1:span1,-span1:span1)*S0(i-span1:i+span1,k-span1:k+span1))/drds_grid(level)
 Srss=sum(Drss(-span1:span1,-span1:span1)*S0(i-span1:i+span1,k-span1:k+span1))/drds2_grid(level)
 Srrs=sum(Drrs(-span1:span1,-span1:span1)*S0(i-span1:i+span1,k-span1:k+span1))/dr2ds_grid(level)
 Srrss=sum(Drrss(-span1:span1,-span1:span1)*S0(i-span1:i+span1,k-span1:k+span1))/dr2ds2_grid(level)

 Sx=Sr/xr(i1)
 Sz=Ss/zs(k1)
 Sxx=(Srr-Sx*xrr(i1))/xr(i1)**2.d0
 Szz=(Sss-Sz*zss(k1))/zs(k1)**2.d0
 Sxxx=(Srrr-3.d0*Sxx*xr(i1)*xrr(i1)-Sx*xrrr(i1))/xr(i1)**3.d0
 Szzz=(Ssss-3.d0*Szz*zs(k1)*zss(k1)-Sz*zsss(k1))/zs(k1)**3.d0
 Sxxxx=(Srrrr-6.d0*Sxxx*xrr(i1)*xr(i1)**2.d0-(3.d0*xrr(i1)**2.d0+4.d0*xr(i1)*xrrr(i1))*Sxx-Sx*xrrrr(i1))/xr(i1)**4.d0
 Szzzz=(Sssss-6.d0*Szzz*zss(k1)*zs(k1)**2.d0-(3.d0*zss(k1)**2.d0+4.d0*zs(k1)*zsss(k1))*Szz-Sz*zssss(k1))/zs(k1)**4.d0

 Sxz=Srs/xr(i1)/zs(k1)
 Sxxz=(Srrs-Sxz*zs(k1)*xrr(i1))/zs(k1)/xr(i1)**2.d0
 Sxzz=(Srss-Sxz*xr(i1)*zss(k1))/xr(i1)/zs(k1)**2.d0
 Sxxzz=(Srrss-Sxzz*xrr(i1)*zs(k1)**2.d0-Sxxz*zss(k1)*xr(i1)**2.d0-Sxz*xrr(i1)*zss(k1))/(xr(i1)*zs(k1))**2.d0

 SFxx(i,k,level)=Sxx  !!store derivatives
 SFzz(i,k,level)=Szz
 SFxxx(i,k,level)=Sxxx
 SFzzz(i,k,level)=Szzz
 SFxxxx(i,k,level)=Sxxxx
 SFzzzz(i,k,level)=Szzzz

 SFxz(i,k,level)=Sxz
 SFxxz(i,k,level)=Sxxz
 SFxzz(i,k,level)=Sxzz
 SFxxzz(i,k,level)=Sxxzz
 end do
end do
!$OMP END PARALLEL DO

end if
end

subroutine compute_RHS(level)
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,j,k,level
real*8 :: Tr,Tx,Cr,Cx
if (level.eq.1) then
!$OMP PARALLEL DO PRIVATE(i,j,k,Tr,Tx,Cr,Cx)
do k=1,nz_grid(level)-1
 do i=1,nx_grid(level)-1
   RHS(i,k,level)=0.d0
   if (RaT.ne.0) then
    Tr=dot_product(D1(-span1:span1),Tbuoy(i-span1:i+span1,k))/dr
    Tx=Tr/xr(i)
    RHS(i,k,level)=RaT*Tx
   end if
   if (comp.eq.1) then !!add in compositional buoyancy
    do j=1,ntype
     Cr=dot_product(D1(-span1:span1),Cbuoy(j,i-span1:i+span1,k))/dr
     Cx=Cr/xr(i)
     RHS(i,k,level)=RHS(i,k,level)-RaC(j)*Cx
    end do
   end if
 end do
end do
!$OMP END PARALLEL DO
else
 call restrict(level,residual(:,:,level-1),RHS(:,:,level))
end if
end

subroutine restrict(level,array_in,array_out)
use basics
use arrays
implicit none
integer*4 :: i,k,ii,kk,level
real*8 :: corners,sides,centre,array_in(1:nx-1,1:nz-1),array_out(1:nx-1,1:nz-1)
!$OMP PARALLEL DO PRIVATE(i,k,ii,kk,corners,sides,centre)
do k=1,nz_grid(level)-1
  kk=2*k
 do i=1,nx_grid(level)-1
  ii=2*i !!fine grid coordinates from one level up
    corners=(array_in(ii-1,kk-1)+array_in(ii-1,kk+1)+array_in(ii+1,kk+1)+array_in(ii+1,kk-1))/16.d0
    sides=(array_in(ii-1,kk)+array_in(ii,kk-1)+array_in(ii+1,kk)+array_in(ii,kk+1))/8.d0
    centre=array_in(ii,kk)/4.d0
    array_out(i,k)=(corners+sides+centre) !!bilinear
 end do
end do
!$OMP END PARALLEL DO
end

subroutine restrict_coefficients(level,array_in,array_out)
use basics
use arrays
implicit none
integer*4 :: i,k,ii,kk,level
real*8 :: corners,sides,centre,array_in(1:nx-1,1:nz-1),array_out(1:nx-1,1:nz-1)
!$OMP PARALLEL DO PRIVATE(i,k,ii,kk,corners,sides,centre)
do k=1,nz_grid(level)-1
  kk=2*k
 do i=1,nx_grid(level)-1
  ii=2*i !!fine grid coordinates from one level up
    corners=(array_in(ii-1,kk-1)+array_in(ii-1,kk+1)+array_in(ii+1,kk+1)+array_in(ii+1,kk-1))/16.d0
    sides=(array_in(ii-1,kk)+array_in(ii,kk-1)+array_in(ii+1,kk)+array_in(ii,kk+1))/8.d0
    centre=array_in(ii,kk)/4.d0
    array_out(i,k)=(corners+sides+centre)
 end do
end do
!$OMP END PARALLEL DO
end

subroutine test_eigenvalues(level)
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k,level
real*8 :: eI,eR,ctemp,ctemp_diff,ae,be
real*8 :: eigI,eigR
if (iterate.eq.0) then

!$OMP PARALLEL DO PRIVATE(i,k,eI,eR)
do k=1,nz_grid(level)-1
 do i=1,nx_grid(level)-1
   eI=eigI(i,k,level)
   eR=eigR(i,k,level)
   if (eI.ne.0.d0) then
    ratio(i,k,level)=dabs(eR/eI)
   else
    ratio(i,k,level)=1.d99
   end if
 end do
end do
!$OMP END PARALLEL DO

else

if (RKC1.eq.0) then !!RK1 fine grid smoother
 ctemp=-courant_stream*2.d0
!$OMP PARALLEL DO PRIVATE(i,k,eI,eR)
 do k=1,nz_grid(level)-1
  do i=1,nx_grid(level)-1
    eI=dabs(vis_xxpzz(i,k,level)*SI1(i,k,level))+vis_grid(i,k,level)*SI2(i,k,level)+dabs(vis_xz(i,k,level)*SI4(i,k,level))+&
      &dabs(vis_x(i,k,level)*SI5(i,k,level))+dabs(vis_z(i,k,level)*SI6(i,k,level))
    eR=vis_xxmzz(i,k,level)*SR1(i,k,level)+vis_grid(i,k,level)*SR2(i,k,level)+vis_xz(i,k,level)*SR4(i,k,level)+&
       &vis_x(i,k,level)*SR5(i,k,level)+vis_z(i,k,level)*SR6(i,k,level)
     ratio(i,k,level)=dabs(eR/eI)
    if (stream_tstep.eq.1) dt_stream(i,k,level)=ctemp*eR/(eR*eR+eI*eI)
  end do
 end do
!$OMP END PARALLEL DO
elseif (RKC1.eq.1) then !!RKC1 fine grid smoother
! ctemp=-courant_stream*(-RKC1_stability)
 ae=-RKC1_stability/2.d0 !!ellipse parameters
 be=real(nstageS,8) !!!!35.d0
 ctemp=courant_stream*(RKC1_stability/ae**2.d0)
 ctemp_diff=courant_stream*RKC1_stability
!$OMP PARALLEL DO PRIVATE(i,k,eI,eR)
 do k=1,nz_grid(level)-1
  do i=1,nx_grid(level)-1
    eI=dabs(vis_xxpzz(i,k,level)*SI1(i,k,level))+vis_grid(i,k,level)*SI2(i,k,level)+dabs(vis_xz(i,k,level)*SI4(i,k,level))+&
      &dabs(vis_x(i,k,level)*SI5(i,k,level))+dabs(vis_z(i,k,level)*SI6(i,k,level))
    eR=vis_xxmzz(i,k,level)*SR1(i,k,level)+vis_grid(i,k,level)*SR2(i,k,level)+vis_xz(i,k,level)*SR4(i,k,level)+&
       &vis_x(i,k,level)*SR5(i,k,level)+vis_z(i,k,level)*SR6(i,k,level)
     ratio(i,k,level)=dabs(eR/eI)
    if (stream_tstep.eq.1) then !!compute smoother time step sizes when needed
     if (diff_adv.eq.1) dt_stream(i,k,level)=ctemp*eR/((eR/ae)**2.d0+(eI/be)**2.d0) !!diffusion and advection terms
     if (diffusion.eq.1) dt_stream_diff(i,k,level)=ctemp_diff/eR                    !!diffusion terms only
    end if
  end do
 end do
!$OMP END PARALLEL DO
end if

end if
end

subroutine iterate_local_time_steps(level) !!now built into test_eigenvalues subroutine
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k,level
real*8 :: eI,eR
real*8 :: eigI,eigR
!$OMP PARALLEL DO PRIVATE(i,k,eI,eR)
do i=1,nx_grid(level)-1
 do k=1,nz_grid(level)-1
   eI=eigI(i,k,level)
   eR=eigR(i,k,level)
   dt_stream(i,k,level)=-courant_stream*2.d0*eR/(eR**2.d0+eI**2.d0)    
 end do
end do
!$OMP END PARALLEL DO
end

subroutine eigenvalue_prep
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k,level,i1,k1
real*8 :: Sr,Ss,Srr,Sss,Srrr,Ssss,Srrrr,Sssss
real*8 :: Sx,Sz,Sxx,Szz,Sxxx,Szzz,Sxxxx,Szzzz,Sxz,Sxxz,Sxzz,Sxxzz
real*8 :: Srs,Srrs,Srss,Srrss
real*8 :: Tr,Trr,Trrr,Trrrr,Ts,Tss,Tsss,Tssss,Trs,Trrss,Trss,Trrs

do level=1,grids

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!eigenvalues for derivatives
Ts=Tr_base/ds_grid(level)     !!scale according to the denominators of finite difference formulas
Tss=-Trr_base/ds2_grid(level) !!even # of derivative operators=real, odd #=imaginary
Tsss=Trrr_base/ds3_grid(level) !!for real values keep signs (sign=(imaginary unit)^(# of derivatives)), for imaginary eigenvalues take absolute values
Tssss=Trrrr_base/ds4_grid(level)
Tr=Tr_base/dr_grid(level)
Trr=-Trr_base/dr2_grid(level)
Trrr=Trrr_base/dr3_grid(level)
Trrrr=Trrrr_base/dr4_grid(level)
Trs=-Trs_base/drds_grid(level)
Trrss=Trrss_base/dr2ds2_grid(level)
Trss=Trrs_base/drds2_grid(level)
Trrs=Trrs_base/dr2ds_grid(level)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end eigenvalues for derivatives

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! REAL PARTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!feed in real parts of the eigenvalues (even # of derivative operators)
Sr=0.d0
Ss=0.d0
Srr=Trr
Sss=Tss
Srrr=0.d0
Ssss=0.d0
Srrrr=Trrrr
Sssss=Tssss
Srs=Trs
Srrs=0.d0
Srss=0.d0
Srrss=Trrss
!!!!!!!!!!!!!!!!!!!!

!$OMP PARALLEL DO PRIVATE(i,k,Sx,Sxx,Sxxx,Sxxxx,Sz,Szz,Szzz,Szzzz,Sxz,Sxxz,Sxzz,Sxxzz,i1,k1)
do i=1,nx_grid(level)-1
 i1=i*2**(level-1) !!fine grid coordinates for grid metrics (may need to restrict the metrics eventually)
 do k=1,nz_grid(level)-1
  k1=k*2**(level-1)
  Sx=Sr/xr(i1)
  Sz=Ss/zs(k1)
  Sxx=(Srr-Sx*xrr(i1))/xr(i1)**2.d0
  Szz=(Sss-Sz*zss(k1))/zs(k1)**2.d0
  Sxxx=(Srrr-3.d0*Sxx*xr(i1)*xrr(i1)-Sx*xrrr(i1))/xr(i1)**3.d0
  Szzz=(Ssss-3.d0*Szz*zs(k1)*zss(k1)-Sz*zsss(k1))/zs(k1)**3.d0
  Sxxxx=(Srrrr-6.d0*Sxxx*xrr(i1)*xr(i1)**2.d0-(3.d0*xrr(i1)**2.d0+4.d0*xr(i1)*xrrr(i1))*&
  &Sxx-Sx*xrrrr(i1))/xr(i1)**4.d0
  Szzzz=(Sssss-6.d0*Szzz*zss(k1)*zs(k1)**2.d0-(3.d0*zss(k1)**2.d0+4.d0*zs(k1)*zsss(k1))*&
  &Szz-Sz*zssss(k1))/zs(k1)**4.d0

  Sxz=Srs/xr(i1)/zs(k1)
  Sxxz=(Srrs-Sxz*zs(k1)*xrr(i1))/zs(k1)/xr(i1)**2.d0
  Sxzz=(Srss-Sxz*xr(i1)*zss(k1))/xr(i1)/zs(k1)**2.d0
  Sxxzz=(Srrss-Sxzz*xrr(i1)*zs(k1)**2.d0-Sxxz*zss(k1)*xr(i1)**2.d0-Sxz*xrr(i1)*zss(k1))/(xr(i1)*zs(k1))**2.d0

  SR1(i,k,level)=-(Sxx-Szz)
  SR2(i,k,level)=-((Sxxxx+Szzzz)+2.d0*Sxxzz)
!  SR2(i,k,level)=-(Sxxxx+Szzzz)
!  SR3(i,k,level)=-2.d0*Sxxzz
  SR4(i,k,level)=-4.d0*Sxz
  SR5(i,k,level)=-2.d0*(Sxxx+Sxzz)
  SR6(i,k,level)=-2.d0*(Szzz+Sxxz)
 end do
end do
!$OMP END PARALLEL DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! END REAL PARTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! IMAGINARY PARTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!feed in imaginary parts of the eigenvalues (odd # of derivative operators)
Sr=Tr
Ss=Ts
Srr=0.d0
Sss=0.d0
Srrr=Trrr
Ssss=Tsss
Srrrr=0.d0
Sssss=0.d0
Srs=0.d0
Srrs=Trrs
Srss=Trss
Srrss=0.d0
!!!!!!!!!!!!!!!!!!!!

!$OMP PARALLEL DO PRIVATE(i,k,Sx,Sxx,Sxxx,Sxxxx,Sz,Szz,Szzz,Szzzz,Sxz,Sxxz,Sxzz,Sxxzz,i1,k1)
do i=1,nx_grid(level)-1
 i1=i*2**(level-1) !!fine grid coordinates for grid metrics (may need to restrict the metrics eventually)
 do k=1,nz_grid(level)-1
  k1=k*2**(level-1)
  Sx=Sr/xr(i1)  !!!!!!!!!!!!!adding absolute values here to be conservative
  Sxx=(Srr+dabs(Sx*xrr(i1)))/xr(i1)**2.d0
  Sxxx=(Srrr+dabs(3.d0*Sxx*xr(i1)*xrr(i1))+dabs(Sx*xrrr(i1)))/dabs(xr(i1)**3.d0)
  Sxxxx=(Srrrr+dabs(6.d0*Sxxx*xrr(i1)*xr(i1)**2.d0)+(dabs(3.d0*xrr(i1)**2.d0+4.d0*xr(i1)*xrrr(i1))*&
  &dabs(Sxx))+dabs(Sx*xrrrr(i1)))/xr(i1)**4.d0

  Sz=Ss/zs(k1)
  Szz=(Sss+dabs(Sz*zss(k1)))/zs(k1)**2.d0
  Szzz=(Ssss+dabs(3.d0*Szz*zs(k1)*zss(k1))+dabs(Sz*zsss(k1)))/dabs(zs(k1)**3.d0)
  Szzzz=(Sssss+dabs(6.d0*Szzz*zss(k1)*zs(k1)**2.d0)+(dabs(3.d0*zss(k1)**2.d0+4.d0*zs(k1)*zsss(k1))*&
  &dabs(Szz))+dabs(Sz*zssss(k1)))/zs(k1)**4.d0

  Sxz=Srs/xr(i1)/zs(k1)
  Sxxz=(Srrs+dabs(Sxz*zs(k1)*xrr(i1)))/dabs(zs(k1))/xr(i1)**2.d0
  Sxzz=(Srss+dabs(Sxz*xr(i1)*zss(k1)))/dabs(xr(i1))/zs(k1)**2.d0
  Sxxzz=(Srrss+dabs(Sxzz*xrr(i1))*zs(k1)**2.d0+dabs(Sxxz*zss(k1))*xr(i1)**2.d0+&
  &dabs(Sxz*xrr(i1)*zss(k1)))/(xr(i1)*zs(k1))**2.d0

  SI1(i,k,level)=(Sxx+Szz)
  SI2(i,k,level)=dabs(Sxxxx+Szzzz)+dabs(2.d0*Sxxzz)
!  SI2(i,k,level)=(Sxxxx+Szzzz)
!  SI3(i,k,level)=2.d0*Sxxzz
  SI4(i,k,level)=4.d0*Sxz
  SI5(i,k,level)=2.d0*(Sxxx+Sxzz)
  SI6(i,k,level)=2.d0*(Szzz+Sxxz)
 end do
end do
!$OMP END PARALLEL DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! END IMAGINARY PARTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

end do !!end loop over grid levels
end 

real*8 function eigI(i,k,level)
 use basics
 use arrays
implicit none
integer*4 :: i,k,level
real*8 :: term1,term2,term3,term4,term5,term6

! term1=(vis_xx(level,i,k)+vis_zz(level,i,k))*(SIxx(level,i,k)+SIzz(level,i,k))
! term2=vis_grid(i,k,level)*(SIxxxx(level,i,k)+SIzzzz(level,i,k))
! term3=2.d0*vis_grid(i,k,level)*SIxxzz(level,i,k)
! term4=4.d0*vis_xz(level,i,k)*SIxz(level,i,k)
! term5=2.d0*vis_x(level,i,k)*(SIxxx(level,i,k)+SIxzz(level,i,k))
! term6=2.d0*vis_z(level,i,k)*(SIzzz(level,i,k)+SIxxz(level,i,k))

!!static portions are stored in SI1...SI6 - including factors of 2 and 4 

! term1=vis_xxpzz(i,k,level)*SI1(i,k,level)
! term2=vis_grid(i,k,level)*SI2(i,k,level) !!combined term2 and term3
! term4=vis_xz(i,k,level)*SI4(i,k,level)
! term5=vis_x(i,k,level)*SI5(i,k,level)
! term6=vis_z(i,k,level)*SI6(i,k,level)

! eigI=dabs(term1)+term2+dabs(term4)+dabs(term5)+dabs(term6)
 eigI=dabs(vis_xxpzz(i,k,level)*SI1(i,k,level))+vis_grid(i,k,level)*SI2(i,k,level)+dabs(vis_xz(i,k,level)*SI4(i,k,level))+&
     &dabs(vis_x(i,k,level)*SI5(i,k,level))+dabs(vis_z(i,k,level)*SI6(i,k,level))
end


real*8 function eigR(i,k,level)
 use basics
 use arrays
implicit none
integer*4 :: i,k,level
real*8 :: term1,term2,term3,term4,term5,term6

! term1=(vis_xx(level,i,k)-vis_zz(level,i,k))*(SRxx(level,i,k)-SRzz(level,i,k))
! term2=vis_grid(i,k,level)*(SRxxxx(level,i,k)+SRzzzz(level,i,k))
! term3=2.d0*vis_grid(i,k,level)*SRxxzz(level,i,k)
! term4=4.d0*vis_xz(level,i,k)*SRxz(level,i,k)
! term5=2.d0*vis_x(level,i,k)*(SRxxx(level,i,k)+SRxzz(level,i,k))
! term6=2.d0*vis_z(level,i,k)*(SRzzz(level,i,k)+SRxxz(level,i,k))

! term1=vis_xxmzz(i,k,level)*SR1(i,k,level)
! term2=vis_grid(i,k,level)*SR2(i,k,level) !!combined term2 and term3
! term4=vis_xz(i,k,level)*SR4(i,k,level)
! term5=vis_x(i,k,level)*SR5(i,k,level)
! term6=vis_z(i,k,level)*SR6(i,k,level)

!  eigR=term1+term2+term4+term5+term6
  eigR=vis_xxmzz(i,k,level)*SR1(i,k,level)+vis_grid(i,k,level)*SR2(i,k,level)+vis_xz(i,k,level)*SR4(i,k,level)+&
      &vis_x(i,k,level)*SR5(i,k,level)+vis_z(i,k,level)*SR6(i,k,level) 
end


subroutine iterate_stream(level,Niter,dir,p)
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,j,k,ID,n,Niter,tt,level,dir,nn,iq,kq,p
real*8 :: temp,la,lb,time1,time2,t1,t2,t1A,t2A,la_user,lb_user
real*8 :: S0(-span:nx+span,-span:nz+span),S1(-span:nx+span,-span:nz+span),S2(-span:nx+span,-span:nz+span)
real*8 :: stream_space,stream_space_diff,c(1:nx-1,1:nz-1),cmin,ntemp

!!!!!!!!!!!!!!!!!!!frequency range for nonstationary iterations
if (RKC1.eq.0) then
 la_user=-1.d0
 lb_user=-2.d0
elseif (RKC1.eq.1) then
 lb_user=RKC1_stability
end if
!!!!!!!!!!!!!!!!!!!end frequency range

!$ time1=omp_get_wtime()
static=1
stream_tstep=1 !!time step (dt_stream) needs to be computed
do n=1,Niter
S0=error(:,:,level) !!initial stream function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!part A
!$ t1=omp_get_wtime()
!!!!!!!!!!!!!!!!!!!!!!!!!!update viscosity field before each fine grid smoother iteration
!$ t1A=omp_get_wtime()
     call stream_derivatives(level,S0)
!$ t2A=omp_get_wtime()
!$ t_isA1=t_isA1+(t2A-t1A)
   if (non_Newtonian.eq.1) then
!$ t1A=omp_get_wtime()
    call compute_strain_rate_invariant_quick(level)
!$ t2A=omp_get_wtime()
!$ t_isA2=t_isA2+(t2A-t1A)
    if (n.eq.2) static=0 !!!!no need to recompute static portion of viscosity after first iteration
!$ t1A=omp_get_wtime()
    call compute_viscosity(level)
!$ t2A=omp_get_wtime()
!$ t_isA3=t_isA3+(t2A-t1A)
   end if
!!!!!!!!!!!!!!!!!!!!!!!!!!end update viscosity field before each fine grid iteration
!$ t2=omp_get_wtime()
!$ t_isA=t_isA+(t2-t1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end part A


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!part B
!$ t1=omp_get_wtime()
if (RKC1.eq.0) then !!!use RK1

  ntemp=Niter-n+1 !!!start with smaller time steps to kill off high frequencies first
  la=la_user
temp=dcos(real(2*ntemp-1,8)*pii/real(2*Niter,8))
!$OMP PARALLEL DO PRIVATE(i,k,lb) !!local time steps dependending on eigenvalue ratio
  do k=1,nz_grid(level)-1
   do i=1,nx_grid(level)-1
    if (ratio(i,k,level).lt.1.d99) then
     lb=lb_user*ratio(i,k,level)**2.d0/(1.d0+ratio(i,k,level)**2.d0)
    else
     lb=lb_user
    end if
    c(i,k)=2.d0/(-lb-la+(lb-la)*temp)
   end do
  end do
!$OMP END PARALLEL DO

elseif (RKC1.eq.1) then !!!use RKC1

 ntemp=Niter-n+1 !!!start with smaller time steps to kill off high frequencies first
 if (Nmulti.gt.1) then
  temp=(c_2-c_1)*real(ntemp-1,8)/real(1-Niter,8)
 else
  temp=0.d0
 end if
  c=temp+c_2   !!!courant # relative to RK1 stability limit

end if
!$ t2=omp_get_wtime()
!$ t_isB=t_isB+(t2-t1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end part B

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!part C
!$ t1=omp_get_wtime()
if (RKC1.eq.0) then !!!use RK1

if (dir.eq.1)  c=1.0d0 !!!going up the multigrid cycle - no smoothing required, just speed
!$OMP PARALLEL DO PRIVATE(i,k)
do k=1,nz_grid(level)-1
 do i=1,nx_grid(level)-1
  error(i,k,level)=S0(i,k)+c(i,k)*dt_stream(i,k,level)*stream_space(i,k,S0,level)
 end do
end do
!$OMP END PARALLEL DO
 call enforceBCs_stream(level,error(:,:,level))

elseif (RKC1.eq.1) then !!!use RKC1

if (splitting.eq.1) then !!first pass with just diffusive terms to better diminish moderate frequency error amplitudes 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Begin 1st stage -- diffusive terms
!$OMP PARALLEL DO PRIVATE(i,k)
do k=1,nz_grid(level)-1
 do i=1,nx_grid(level)-1
  S1(i,k)=S0(i,k)+mewtChebS(1)*c(i,k)*dt_stream_diff(i,k,level)*stream_space_diff(i,k,S0,level)
  S2(i,k)=S0(i,k)
 end do
end do
!$OMP END PARALLEL DO
 call enforceBCs_stream(level,S1(:,:))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! End 1st Stage -- diffusive terms

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!remaining stages -- diffusive terms
do j=2,nstageS
 call stream_derivatives(level,S1) !!!!***********************************************ADD non-Newtonian viscosity update here
!$OMP PARALLEL DO PRIVATE(i,k)
do k=1,nz_grid(level)-1
 do i=1,nx_grid(level)-1
  error(i,k,level)=cChebS(j)*S0(i,k)+mewChebS(j)*S1(i,k)+vChebS(j)*S2(i,k)+&
                  &mewtChebS(j)*c(i,k)*dt_stream_diff(i,k,level)*stream_space_diff(i,k,S1,level)
 end do
end do
!$OMP END PARALLEL DO
 call enforceBCs_stream(level,error(:,:,level))

if (j.lt.nstageS) then
!$OMP PARALLEL DO PRIVATE(i,k)
do k=-span,nz_grid(level)+span
 do i=-span,nx_grid(level)+span
  S2(i,k)=S1(i,k)
  S1(i,k)=error(i,k,level)
 end do
end do
!$OMP END PARALLEL DO
end if
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end remaining stages -- diffusive terms
S0=error(:,:,level) !!set up initial stream function for next smoother pass with diffusive and advective terms
call stream_derivatives(level,S0) !!set up stream derivatives to start next smoother pass

end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Begin 1st stage -- diffusive and advective terms
!$OMP PARALLEL DO PRIVATE(i,k)
do k=1,nz_grid(level)-1
 do i=1,nx_grid(level)-1
  S1(i,k)=S0(i,k)+mewtChebS(1)*c(i,k)*dt_stream(i,k,level)*stream_space(i,k,S0,level)
  S2(i,k)=S0(i,k)
 end do
end do
!$OMP END PARALLEL DO
 call enforceBCs_stream(level,S1(:,:))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! End 1st Stage -- diffusive and advective terms

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!remaining stages -- diffusive and advective terms
do j=2,nstageS
 call stream_derivatives(level,S1) !!!!***********************************************ADD non-Newtonian viscosity update here
!$OMP PARALLEL DO PRIVATE(i,k)
do k=1,nz_grid(level)-1
 do i=1,nx_grid(level)-1
  error(i,k,level)=cChebS(j)*S0(i,k)+mewChebS(j)*S1(i,k)+vChebS(j)*S2(i,k)+&
                  &mewtChebS(j)*c(i,k)*dt_stream(i,k,level)*stream_space(i,k,S1,level)
 end do
end do
!$OMP END PARALLEL DO
 call enforceBCs_stream(level,error(:,:,level))

if (j.lt.nstageS) then
!$OMP PARALLEL DO PRIVATE(i,k)
do k=-span,nz_grid(level)+span
 do i=-span,nx_grid(level)+span
  S2(i,k)=S1(i,k)         
  S1(i,k)=error(i,k,level)
 end do
end do
!$OMP END PARALLEL DO
end if
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end remaining stages -- diffusive and advective terms

end if
!$ t2=omp_get_wtime()
!$ t_isC=t_isC+(t2-t1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end part C
end do

 call stream_derivatives(level,error(:,:,level))
if (non_Newtonian.eq.1) then
 call compute_strain_rate_invariant_quick(level)
 stream_tstep=0 !!dt_stream does not need to be computed
 call compute_viscosity(level)
end if
 call compute_residual(level,error(:,:,level),p)
!$ time2=omp_get_wtime()
!$ t_istream=t_istream+(time2-time1)
end

subroutine enforceBCs_stream(level,S0) !!!!update this
 use basics
implicit none
integer*4 :: i,k,level
real*8 :: S0(-span:nx+span,-span:nz+span)
S0(0,-span:nz_grid(level)+span)=0.d0   !!impermeable boundaries
S0(nx_grid(level),-span:nz_grid(level)+span)=0.d0
S0(-span:nx_grid(level)+span,0)=0.d0
S0(-span:nx_grid(level)+span,nz_grid(level))=0.d0
do k=0,nz_grid(level)
 do i=1,span  !!!antisymmetry for free-slip sidewalls
  S0(-i,k)=-S0(i,k)      !!antsymmetry: free-slip
!  S0(nx_grid(level)+i,k)=-S0(nx_grid(level)-i,k)
 end do
end do
do k=0,nz_grid(level)
 do i=1,span  !!!antisymmetry for free-slip sidewalls
  S0(nx_grid(level)+i,k)=-S0(nx_grid(level)-i,k)
 end do
end do


if (Vbc.eq.0) then !!antisymmetry for free-slip top/bottom
 do k=1,span
  do i=-span,nx_grid(level)+span
   S0(i,-k)=-S0(i,k)
!   S0(i,nz_grid(level)+k)=-S0(i,nz_grid(level)-k)
  end do
 end do
 do k=1,span
  do i=-span,nx_grid(level)+span
   S0(i,nz_grid(level)+k)=-S0(i,nz_grid(level)-k)
  end do
 end do
elseif (Vbc.eq.1) then !!symmetry for rigid top/bottom
 do k=1,span
  do i=-span,nx_grid(level)+span
   S0(i,-k)=S0(i,k)
!   S0(i,nz_grid(level)+k)=S0(i,nz_grid(level)-k)
  end do
 end do
 do k=1,span
  do i=-span,nx_grid(level)+span
!   S0(i,-k)=S0(i,k)
   S0(i,nz_grid(level)+k)=S0(i,nz_grid(level)-k)
  end do
 end do
end if
end

integer*4 function pf_coarse(ip,kp) !!input is ip,kp -- output is p value
 use basics
 use arrays
implicit none
integer*4 :: ip,kp,N
if (Database.eq.1) then
 N=nx_grid(grids)+1
 pf_coarse=N*kp+ip+1
elseif (Database.eq.2) then
 N=nz_grid(grids)+1
 pf_coarse=N*ip+kp+1
elseif (Database.eq.3) then !!update or take out
 pf_coarse=parray(ip,kp)
end if
end

subroutine indices_coarse(p,ip,kp)  !!input is p, output is ip,kp
 use basics
 use arrays
implicit none
integer*4 :: p,ip,kp,N
if (Database.eq.1) then
 N=nx_grid(grids)+1
 ip=p-N*((p-1)/N)-1 !!!!!!!!!!!!!!!!!!!!beware of modular arithmetic
 kp=(p-1)/N
elseif (Database.eq.2) then
 N=nz_grid(grids)+1
 kp=p-N*((p-1)/N)-1
 ip=(p-1)/N
elseif (Database.eq.3) then !!!!!update or remove
 ip=ipvec(p)
 kp=kpvec(p)
end if
end

subroutine Compute_Matrix_coarse
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: ip,kp,p,iq,kq,q,pf_coarse,iiq,kkq,INFO,limit(-span:span),j,temp
integer*4 :: nx_coarse,nz_coarse,kl_ku_1,xval,zval,Sval
real*8 :: S0(-span:nx+span,-span:nz+span)
real*8 :: ROWCND,COLCND,AMAX,ta1,ta2,stream_space_matrix
!$ ta1=omp_get_wtime()

nx_coarse=nx_grid(grids)
nz_coarse=nz_grid(grids)

!!!!!!!!!!!!!!!!!define the stencil shape for shorter loops
limit(-span)=0
limit(-span1:-1)=span1
limit(0)=span
limit(1:span1)=span1
limit(span)=0
!!!!!!!!!!!!!!!! end stencil

if (analysis.eq.1) then !!!figure of matrix bandwidth on startup
ku=0
kl=0
!$OMP PARALLEL DO PRIVATE(p,q,temp,j,iq,kq,ip,kp) REDUCTION(max:ku,kl)
 do p=1,ngrid_coarse
  call indices_coarse(p,ip,kp)
   do j=-span,span
    iq=j+ip
    if ((iq.gt.0).and.(iq.lt.nx_coarse)) then
     do kq=-limit(j)+kp,kp+limit(j)
       if ((kq.gt.0).and.(kq.lt.nz_coarse)) then
        q=pf_coarse(iq,kq)
        temp=q-p
        if (temp.gt.ku) ku=temp
        if (-temp.gt.kl) kl=-temp
       end if
     end do
    end if
   end do
 end do
!$OMP END PARALLEL DO
write(*,*) "Coarse Matrix:",ku,kl,ngrid_coarse
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Construct Matrix Once Bandwidth is known
if (analysis.eq.0) then
kl_ku_1=kl+ku+1
Matrix_coarse=0.d0

!$OMP PARALLEL DO PRIVATE(p,q,j,ip,kp,iq,kq,iiq,kkq,xval,zval,Sval)
do p=1,ngrid_coarse
 call indices_coarse(p,ip,kp)
  if ((ip.gt.0).and.(ip.lt.nx_coarse).and.(kp.gt.0).and.(kp.lt.nz_coarse)) then  !!!for interior values of stream function
   do j=-span,span
    iq=j+ip
    if ((iq.ge.0).and.(iq.le.nx_coarse)) then
      iiq=iq
     do kq=-limit(j)+kp,kp+limit(j)
      if ((kq.ge.0).and.(kq.le.nz_coarse)) then
       kkq=kq
       Sval=1
       q=pf_coarse(iiq,kkq)
       Matrix_coarse(kl_ku_1+p-q,q)=stream_space_matrix(ip,kp,grids,iq,kq,Sval)
      end if
     end do
    end if
   end do
  end if
 if ((ip.gt.0).and.(ip.lt.nx_coarse).and.(kp.gt.0).and.(kp.lt.nz_coarse)) then  !!!for ghost values of the streamfunction
  do j=-span,span
   iq=j+ip
   do kq=-limit(j)+kp,kp+limit(j)
    if ((iq.lt.0).or.(iq.gt.nx_coarse).or.(kq.lt.0).or.(kq.gt.nz_coarse)) then
     if (iq.lt.0) then  !!free-slip left sidewall: antisymmetry
      iiq=-iq
      xval=-1
     elseif (iq.gt.nx_coarse) then!!free-slip right sidewall: antisymmetry
      iiq=2*nx_coarse-iq
      xval=-1
     else
      iiq=iq
      xval=1
     end if

     if (kq.lt.0) then !!!bottom
      kkq=-kq
      if (Vbc.eq.0) then !!free-slip: antisymmetry
       zval=-1
      elseif (Vbc.eq.1) then !!rigid: symmetry
       zval=1
      end if
     elseif (kq.gt.nz_coarse) then  !!!!top
      kkq=2*nz_coarse-kq
      if (Vbc.eq.0) then !!free-slip: antisymmetry
       zval=-1
      elseif (Vbc.eq.1) then !!rigid: symmetry
       zval=1
      end if
     else
      kkq=kq
      zval=1
     end if
     Sval=xval*zval
     q=pf_coarse(iiq,kkq)
     Matrix_coarse(kl_ku_1+p-q,q)=Matrix_coarse(kl_ku_1+p-q,q)+stream_space_matrix(ip,kp,grids,iq,kq,Sval)
    end if
   end do
  end do
 end if
end do
!$OMP END PARALLEL DO

!D!$OMP PARALLEL DO PRIVATE(kp,p) !!!!parallel is overkill????
do kp=0,nz_coarse                !!impermeable boundaries
   p=pf_coarse(0,kp)
   Matrix_coarse(kl_ku_1,p)=1.d0
   p=pf_coarse(nx_coarse,kp)
   Matrix_coarse(kl_ku_1,p)=1.d0
end do
!D!$OMP END PARALLEL DO
!D!$OMP PARALLEL DO PRIVATE(ip,p)
do ip=0,nx_coarse                !!impermeable boundaries
   p=pf_coarse(ip,0)
   Matrix_coarse(kl_ku_1,p)=1.d0
   p=pf_coarse(ip,nz_coarse)
   Matrix_coarse(kl_ku_1,p)=1.d0
end do
!D!$OMP END PARALLEL DO

 call DGBEQU(ngrid_coarse,ngrid_coarse,kl,ku,Matrix_coarse(kl+1:2*kl+ku+1,1:ngrid_coarse),kl+ku+1,RowB,ColB,ROWCND,COLCND,AMAX,INFO)
!$OMP PARALLEL DO PRIVATE(q,p)
 do q=1,ngrid_coarse  !!scale matrix
  do p=max(1,q-kl),min(ngrid_coarse,q+ku)
   Matrix_coarse(kl_ku_1+p-q,q)=RowB(p)*Matrix_coarse(kl_ku_1+p-q,q)
  end do
 end do
!$OMP END PARALLEL DO

end if
!$ ta2=omp_get_wtime()
!$ tcoarse_matrix=tcoarse_matrix+ta2-ta1
end

subroutine Compute_coarse_vector !!build RHS for linear system on coarsest grid
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: p,i,k,pf_coarse
if (non_Newtonian.eq.1) then
B_coarse_save=0.d0 !!ensure boundary values are zero to ensure stream function will be zero on model boundaries
!$OMP PARALLEL DO PRIVATE(i,k,p)
do k=1,nz_grid(grids)-1
 do i=1,nx_grid(grids)-1
  p=pf_coarse(i,k)
  B_coarse_save(p)=RHS(i,k,grids)+fake(i,k)
 end do
end do
!$OMP END PARALLEL DO
else
B_coarse=0.d0 !!ensure boundary values are zero to ensure stream function will be zero on model boundaries
!$OMP PARALLEL DO PRIVATE(i,k,p)
do k=1,nz_grid(grids)-1
 do i=1,nx_grid(grids)-1
  p=pf_coarse(i,k)
  B_coarse(p)=RHS(i,k,grids)*RowB(p)
 end do
end do
!$OMP END PARALLEL DO
end if
end

subroutine Stream_grid_coarse
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: p,ip,kp
!$OMP PARALLEL DO PRIVATE(ip,kp,p)
do p=1,ngrid_coarse
 call indices_coarse(p,ip,kp)
 error(ip,kp,grids)=B_coarse(p)
end do
!$OMP END PARALLEL DO
end
