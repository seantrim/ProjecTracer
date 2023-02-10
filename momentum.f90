subroutine compute_strain_rate_invariant
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k
real*8 :: ur,us,wr,ws
real*8 :: ux,uz,wx,wz
!$OMP PARALLEL DO PRIVATE(i,k,ur,us,wr,ws,ux,uz,wx,wz)
do k=0,nz
 do i=0,nx
  ur=dot_product(D1(-span1:span1),u(i-span1:i+span1,k))/dr
  us=dot_product(D1(-span1:span1),u(i,k-span1:k+span1))/ds
  wr=dot_product(D1(-span1:span1),w(i-span1:i+span1,k))/dr
  ws=dot_product(D1(-span1:span1),w(i,k-span1:k+span1))/ds 
  ux=ur/xr(i)
  uz=us/zs(k)
  wx=wr/xr(i)
  wz=ws/zs(k)
  strain(i,k,1)=(0.5d0*(2.d0*(ux**2.d0+wz**2.d0)+(uz+wx)**2.d0))**0.5d0
 end do
end do
!$OMP END PARALLEL DO
!!!!!!!!!!!!!!!!Boundary conditions (implied symmetry given free-slip conditions on all sides)
do i=1,span
 strain(-i,0:nz,1)=strain(i,0:nz,1)
 strain(nx+i,0:nz,1)=strain(nx-i,0:nz,1)
end do
do k=1,span
 strain(-span:nx+span,-k,1)=strain(-span:nx+span,k,1)
 strain(-span:nx+span,nz+k,1)=strain(-span:nx+span,nz-k,1)
end do
if (non_Newtonian.eq.1) call smoother_strain(1,0) !!!!!!!!!!!!!!!test: smooth high frequencies out of strain rate before feeding into viscosity
end

subroutine compute_strain_rate_invariant_coarse(level)
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k,i1,k1,level
real*8 :: ur,us,wr,ws
real*8 :: ux,uz,wx,wz
!$OMP PARALLEL DO PRIVATE(i,k,ur,us,wr,ws,ux,uz,wx,wz,i1,k1)
do k=0,nz_grid(level)
k1=k*2**(level-1)
 do i=0,nx_grid(level)
 i1=i*2**(level-1) !!fine grid coordinates for grid metrics (may need to restrict the metrics eventually)
  ur=dot_product(D1(-span1:span1),u_grid(level,i-span1:i+span1,k))/dr_grid(level)
  us=dot_product(D1(-span1:span1),u_grid(level,i,k-span1:k+span1))/ds_grid(level)
  wr=dot_product(D1(-span1:span1),w_grid(level,i-span1:i+span1,k))/dr_grid(level)
  ws=dot_product(D1(-span1:span1),w_grid(level,i,k-span1:k+span1))/ds_grid(level)
  ux=ur/xr(i1)
  uz=us/zs(k1)
  wx=wr/xr(i1)
  wz=ws/zs(k1)
  strain(i,k,level)=(0.5d0*(2.d0*(ux**2.d0+wz**2.d0)+(uz+wx)**2.d0))**0.5d0
 end do
end do
!$OMP END PARALLEL DO
!!!!!!!!!!!!!!!!Boundary conditions (implied symmetry given free-slip conditions on all sides)
do i=1,span
 strain(-i,0:nz_grid(level),level)=strain(i,0:nz_grid(level),level)
 strain(nx_grid(level)+i,0:nz_grid(level),level)=strain(nx_grid(level)-i,0:nz_grid(level),level)
end do
do k=1,span
 strain(-span:nx_grid(level)+span,-k,level)=strain(-span:nx_grid(level)+span,k,level)
 strain(-span:nx_grid(level)+span,nz_grid(level)+k,level)=strain(-span:nx_grid(level)+span,nz_grid(level)-k,level)
end do
if (non_Newtonian.eq.1) call smoother_strain(level,0) !!!!!!!!!!!!!!!test: smooth high frequencies out of strain rate before feeding into viscosity
end

subroutine compute_strain_rate_invariant_direct(level) !!!directly from stream function
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k,i1,k1,level
real*8 :: Sr,Ss,Srr,Sss,Srs
real*8 :: Sx,Sz,Sxz,Szz,Sxx
if ((gridX.eq.0).and.(gridZ.eq.0)) then

!$OMP PARALLEL DO PRIVATE(i,k,Sxz,Szz,Sxx)
do k=0,nz_grid(level)
 do i=0,nx_grid(level)
  Sxx=dot_product(D2(-span1:span1),error(i-span1:i+span1,k,level))/dr2_grid(level)
  Szz=dot_product(D2(-span1:span1),error(i,k-span1:k+span1,level))/ds2_grid(level)
  Sxz=sum(Drs(-span1:span1,-span1:span1)*error(i-span1:i+span1,k-span1:k+span1,level))/drds_grid(level)
  strain(i,k,level)=(2.d0*Sxz**2.d0+0.5d0*(Szz-Sxx)**2.d0)**0.5d0
 end do
end do
!$OMP END PARALLEL DO

else

!$OMP PARALLEL DO PRIVATE(i,k,Sr,Ss,Srr,Sss,Srs,Sx,Sz,Sxz,Szz,Sxx,i1,k1)
do k=0,nz_grid(level)
k1=k*2**(level-1)
 do i=0,nx_grid(level)
 i1=i*2**(level-1) !!fine grid coordinates for grid metrics (may need to restrict the metrics eventually)

 Sr=dot_product(D1(-span1:span1),error(i-span1:i+span1,k,level))/dr_grid(level)
 Ss=dot_product(D1(-span1:span1),error(i,k-span1:k+span1,level))/ds_grid(level)
 Srr=dot_product(D2(-span1:span1),error(i-span1:i+span1,k,level))/dr2_grid(level)
 Sss=dot_product(D2(-span1:span1),error(i,k-span1:k+span1,level))/ds2_grid(level)

 Srs=sum(Drs(-span1:span1,-span1:span1)*error(i-span1:i+span1,k-span1:k+span1,level))/drds_grid(level)

 Sx=Sr/xr(i1)
 Sz=Ss/zs(k1)
 Sxx=(Srr-Sx*xrr(i1))/xr(i1)**2.d0
 Szz=(Sss-Sz*zss(k1))/zs(k1)**2.d0

 Sxz=Srs/xr(i1)/zs(k1)

  strain(i,k,level)=(2.d0*Sxz**2.d0+0.5d0*(Szz-Sxx)**2.d0)**0.5d0
 end do
end do
!$OMP END PARALLEL DO

end if

!!!!!!!!!!!!!!!!Boundary conditions (implied symmetry given free-slip conditions on all sides)
do i=1,span
 strain(-i,0:nz_grid(level),level)=strain(i,0:nz_grid(level),level)
 strain(nx_grid(level)+i,0:nz_grid(level),level)=strain(nx_grid(level)-i,0:nz_grid(level),level)
end do
do k=1,span
 strain(-span:nx_grid(level)+span,-k,level)=strain(-span:nx_grid(level)+span,k,level)
 strain(-span:nx_grid(level)+span,nz_grid(level)+k,level)=strain(-span:nx_grid(level)+span,nz_grid(level)-k,level)
end do
if (non_Newtonian.eq.1) call smoother_strain(level,0) !!!!!!!!!!!!!!!test: smooth high frequencies out of strain rate before feeding into viscosity
end

subroutine compute_strain_rate_invariant_quick(level) !!!from stream function derivative arrays
!$ use OMP_LIB
 use basics
 use arrays
implicit none
real*8 :: t1,t2
integer*4 :: i,k,level
!$ t1=omp_get_wtime()
!$OMP PARALLEL DO PRIVATE(i,k)
do k=0,nz_grid(level)
 do i=0,nx_grid(level)
  strain(i,k,level)=(2.d0*SFxz(i,k,level)**2.d0+0.5d0*(SFzz(i,k,level)-SFxx(i,k,level))**2.d0)**0.5d0
 end do
end do
!$OMP END PARALLEL DO
!$ t2=omp_get_wtime()
!$ t_cs1=t_cs1+(t2-t1)

!$ t1=omp_get_wtime()
!!!!!!!!!!!!!!!!Boundary conditions (implied symmetry given free-slip conditions on all sides)
do i=1,span
 strain(-i,0:nz_grid(level),level)=strain(i,0:nz_grid(level),level)
 strain(nx_grid(level)+i,0:nz_grid(level),level)=strain(nx_grid(level)-i,0:nz_grid(level),level)
end do
do k=1,span
 strain(-span:nx_grid(level)+span,-k,level)=strain(-span:nx_grid(level)+span,k,level)
 strain(-span:nx_grid(level)+span,nz_grid(level)+k,level)=strain(-span:nx_grid(level)+span,nz_grid(level)-k,level)
end do
!$ t2=omp_get_wtime()
!$ t_cs2=t_cs2+(t2-t1)

!$ t1=omp_get_wtime()
if (non_Newtonian.eq.1) call smoother_strain(level,0) !!!!!!!!!!!!!!!test: smooth high frequencies out of strain rate before feeding into viscosity
!$ t2=omp_get_wtime()
!$ t_cs3=t_cs3+(t2-t1)
end


subroutine update_velocities
 use basics
implicit none
if (iterate.eq.0) then
 if (ivis.ge.1) then
  call build_SF_matrix  
 end if
 call solve_SF_equation
elseif (iterate.eq.1) then
 if (MG.eq.0) then
  call iterative_solver
 elseif (MG.eq.1) then
  call multigrid
 end if
end if
 call compute_velocities
end

subroutine Database3
 use basics
 use arrays
implicit none
integer*4 :: i,k,p
integer*4 :: i1,i2,k1,k2,Xblocks,Xb,Zblocks,Zb
Xblocks=nx/(2*span+1)+1
Zblocks=nz/(2*span+1)+1

p=1
k1=0
k2=2*span
do Zb=1,Zblocks
 i1=0
 i2=2*span
 do Xb=1,Xblocks
  do k=k1,k2
   do i=i1,i2
    parray(i,k)=p
    ipvec(p)=i
    kpvec(p)=k
    p=p+1
   end do
  end do
  i1=i1+2*span+1
  i2=min(i2+2*span+1,nx)
 end do
  k1=k1+2*span+1
  k2=min(k2+2*span+1,nz)
end do
end

subroutine indices(p,ip,kp)  !!input is p, output is ip,kp
 use basics
 use arrays
implicit none
integer*4 :: p,ip,kp,N
if (Database.eq.1) then
 N=nx+1
 ip=p-N*((p-1)/N)-1
 kp=(p-1)/N
elseif (Database.eq.2) then
 N=nz+1
 kp=p-N*((p-1)/N)-1
 ip=(p-1)/N
elseif (Database.eq.3) then
 ip=ipvec(p)
 kp=kpvec(p)
end if
end

integer*4 function pf(ip,kp) !!input is ip,kp -- output is p value
 use basics
 use arrays
implicit none
integer*4 :: ip,kp,N
if (Database.eq.1) then
 N=nx+1
 pf=N*kp+ip+1
elseif (Database.eq.2) then
 N=nz+1
 pf=N*ip+kp+1
elseif (Database.eq.3) then
 pf=parray(ip,kp)
end if
end

subroutine buoyancy
 use basics
 use arrays
implicit none
integer*4 :: p,ip,kp,tt
real*8 :: Tr,Tx
real*8, allocatable :: Cr(:),Cx(:)
if (comp.eq.1) allocate(Cr(1:ntype),Cx(1:ntype))
Bvec=0.d0
do p=1,ngrid
 call indices(p,ip,kp)
 if ((0.lt.ip).and.(ip.lt.nx).and.(0.lt.kp).and.(kp.lt.nz)) then
  Tr=dot_product(D1(-span1:span1),Tbuoy(ip-span1:ip+span1,kp))/dr
  Tx=Tr/xr(ip)
  if (comp.eq.1) then
   do tt=1,ntype
    Cr(tt)=dot_product(D1(-span1:span1),Cbuoy(tt,ip-span1:ip+span1,kp))/dr
    Cx(tt)=Cr(tt)/xr(ip)
   end do
  end if
  if (comp.eq.0) then
   Bvec(p)=RaT*Tx*RowB(p)
  else
   if (RaT.ne.0.d0) then
    Bvec(p)=(RaT*Tx-dot_product(RaC,Cx))*RowB(p)
   else
    Bvec(p)=(-dot_product(RaC,Cx))*RowB(p)
   end if
  end if
 end if
end do
end

subroutine Stream_grid
 use basics
 use arrays
implicit none
integer*4 :: p,ip,kp
do p=1,ngrid
 call indices(p,ip,kp)
  SF(ip,kp)=Bvec(p)
end do
do kp=0,nz
 do ip=1,span   !!antisymmetry
  SF(-ip,kp)=-SF(ip,kp)
  SF(nx+ip,kp)=-SF(nx-ip,kp)
 end do
end do
do kp=1,span
 do ip=-span,nx+span
  if (Vbc.eq.0) then !!free-slip
   SF(ip,-kp)=-SF(ip,kp)
   SF(ip,nz+kp)=-SF(ip,nz-kp)
  elseif (Vbc.eq.1) then !!rigid
   SF(ip,-kp)=SF(ip,kp)
   SF(ip,nz+kp)=SF(ip,nz-kp)
  end if
 end do
end do
end

subroutine build_SF_matrix
!$ use OMP_LIB
 use basics
 use arrays
implicit none
!!!p=rows of flow matrix,q=columns of flow matrix
!!!p=centre of stencil,q=stencil point
! character*1 :: TRANS,NORM
!integer*4 INFO,NRHS,IWORK(1:ngrid)
!real*8 :: ANORM,RCOND,DLANGE,WORK(1:4*ngrid),ROWCND,COLCND,AMAX
 character*1 :: TRANSB,NORMB
integer*4 INFOB,NRHSB!,IWORKB(1:ngrid)
real*8 :: ANORMB,RCONDB,DLANGB,ROWCNDB,COLCNDB,AMAXB!,WORKB(1:4*ngrid)
integer*4 :: i,k,p,q,ip,kp,iq,kq,ips,ps,kps,pf,rsym,ssym
real*8 :: Sr,Srr,Srrr,Srrrr,Ss,Sss,Ssss,Sssss,Srs,Srss,Srrs,Srrss
real*8 :: Sx,Sxx,Sxxx,Sxxxx
real*8 :: Sz,Szz,Szzz,Szzzz
real*8 :: Sxz,Sxxz,Sxzz,Sxxzz
real*8 :: visf
real*8 :: term1,term2,term3,term4,term5,term6
real*8 :: time1,time2
integer*4 :: level,interpolate
!$ time1=omp_get_wtime()
if (tstep.eq.-1) then
 subD=0  !!for # of subdiagonals
 supD=0  !!for # of superdiagonals
else
 MatrixB=0.d0
end if

if (tstep.gt.-1) then
 do i=0,nx
  do k=0,nz
   vis(i,k)=visf(i,k)
  end do
 end do
!!!!!!!!!!!!!!!!!!!!!!!!!! new filtering technique
 if (ivis.ne.0) then
  vis_grid(:,:,1)=vis(:,:)
  interpolate=0
  level=1
  call smoother_vis(level,interpolate)  !!smooth computational viscosity to kill off unresolved frequencies
  vis(:,:)=vis_grid(:,:,1)
 end if
!!!!!!!!!!!!!!!!!!!!!!!!!! end new filtering technique
end if

do p=1,ngrid !!!! outer loop over p (so that viscosity derivatives are not recomputed)
 call indices(p,ip,kp)

 if (((ip.eq.0).and.(kp.eq.0)).or.((ip.eq.nx).and.(kp.eq.0)).or.((ip.eq.0).and.(kp.eq.nz)).or.((ip.eq.nx).and.(kp.eq.nz))) then !!corners
!  Matrix(p,p)=1.d0
  if (tstep.gt.-1) MatrixB(subD+supD+1,p)=1.d0
 elseif (((Vbc.eq.1).and.(kp.eq.0)).or.((Vbc.eq.1).and.(kp.eq.nz))) then
  if (tstep.gt.-1) MatrixB(subD+supD+1,p)=1.d0
 elseif ((0.le.ip).and.(ip.le.nx).and.(0.le.kp).and.(kp.le.nz)) then !!interior grid points
  if (tstep.gt.-1) then
!   visr=dot_product(D1(-span1:span1),vis(ip-span1:ip+span1,kp))/dr
!   visrr=dot_product(D2(-span1:span1),vis(ip-span1:ip+span1,kp))/dr2
!   viss=dot_product(D1(-span1:span1),vis(ip,kp-span1:kp+span1))/ds
!   visss=dot_product(D2(-span1:span1),vis(ip,kp-span1:kp+span1))/ds2
!   visrs=sum(Drs(-span1:span1,-span1:span1)*vis(ip-span1:ip+span1,kp-span1:kp+span1))/drds

   visr=dot_product(D1_so(-1:1),vis(ip-1:ip+1,kp))/dr         !!!!!!!!!!!!use second order FD stencils for viscosity gradients
   visrr=dot_product(D2_so(-1:1),vis(ip-1:ip+1,kp))/dr2
   viss=dot_product(D1_so(-1:1),vis(ip,kp-1:kp+1))/ds
   visss=dot_product(D2_so(-1:1),vis(ip,kp-1:kp+1))/ds2
   visrs=sum(Drs_so(-1:1,-1:1)*vis(ip-1:ip+1,kp-1:kp+1))/drds


   visx=visr/xr(ip)
   visz=viss/zs(kp)
   visxx=(visrr-visx*xrr(ip))/xr(ip)**2.d0
   viszz=(visss-visz*zss(kp))/zs(kp)**2.d0
   visxz=visrs/xr(ip)/zs(kp)
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!3rd and 4th order FD terms (i.e., terms with largest stencil span)
  do i=-span,span,2*span
   rsym=0
   iq=ip+i
   if (iq.lt.0) then     !!if antisymmetry in r direction
    iq=-iq
    rsym=1
   elseif (iq.gt.nx) then
    iq=nx-(iq-nx)
    rsym=1
   end if
   k=0
    ssym=0
    kq=kp+k
    if (kq.lt.0) then     !!if antisymmetry in s direction
     kq=-kq
     ssym=1 !!antisymmetry only for free-slip: symmetry for rigid
    elseif (kq.gt.nz) then
     kq=nz-(kq-nz)
     ssym=1
    end if
    q=pf(iq,kq)
    call Matrix_element(i,k,ip,kp,p,q,rsym,ssym)
  end do

   i=0
   rsym=0
   iq=ip+i
   if (iq.lt.0) then     !!if antisymmetry in r direction
    iq=-iq
    rsym=1
   elseif (iq.gt.nx) then
    iq=nx-(iq-nx)
    rsym=1
   end if
   do k=-span,span,2*span
    ssym=0
    kq=kp+k
    if (kq.lt.0) then     !!if antisymmetry in s direction
     kq=-kq
     ssym=1
    elseif (kq.gt.nz) then
     kq=nz-(kq-nz)
     ssym=1
    end if
    q=pf(iq,kq)
    call Matrix_element(i,k,ip,kp,p,q,rsym,ssym)
   end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!End 3rd and 4th order FD terms

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1st and 2nd order FD terms
  do i=-span1,span1
   rsym=0
   iq=ip+i
   if (iq.lt.0) then     !!if antisymmetry in r direction
    iq=-iq
    rsym=1
   elseif (iq.gt.nx) then
    iq=nx-(iq-nx)
    rsym=1
   end if
   do k=-span1,span1
    ssym=0
    kq=kp+k
    if (kq.lt.0) then     !!if antisymmetry in s direction
     kq=-kq
     ssym=1
    elseif (kq.gt.nz) then
     kq=nz-(kq-nz)
     ssym=1
    end if
    q=pf(iq,kq)
    call Matrix_element(i,k,ip,kp,p,q,rsym,ssym)
   end do
  end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!End 1st and 2nd order FD terms
 end if
end do
if (tstep.eq.-1) then
 subD=-subD  !!to match definition in LAPACK
 allocate(MatrixB(1:2*subD+supD+1,1:ngrid),RowB(1:ngrid),ColB(1:ngrid),IPIVB(1:ngrid))
 write(*,*) "subD,supD,ngrid=",subD,supD,ngrid
else
 call DGBEQU(ngrid,ngrid,subD,supD,MatrixB(subD+1:2*subD+supD+1,1:ngrid),subD+supD+1,RowB,ColB,ROWCNDB,COLCNDB,AMAXB,INFOB)
! write(*,*) ROWCNDB,COLCNDB
 do q=1,ngrid  !!scale matrix
  do p=max(1,q-supD),min(ngrid,q+subD)
   MatrixB(subD+supD+1+p-q,q)=RowB(p)*MatrixB(subD+supD+1+p-q,q)
  end do
 end do
! ANORMB=DLANGB('1',ngrid,subD,supD,MatrixB(subD+1:2*subD+supD+1,1:ngrid),subD+supD+1)
 call DGBTRF(ngrid,ngrid,subD,supD,MatrixB,2*subD+supD+1,IPIVB,INFOB) 
! call DGBCON('1',ngrid,subD,supD,MatrixB,2*subD+supD+1,IPIVB,ANORMB,RCONDB,WORKB(1:3*ngrid),IWORKB(1:ngrid),INFOB)
end if
!$ time2=omp_get_wtime()
!$ tbSF=tbSF+time2-time1
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Matrix_element(i,k,ip,kp,p,q,rsym,ssym)
 use basics
 use arrays
implicit none
 character*1 :: TRANSB,NORMB
integer*4 INFOB,NRHSB!,IWORKB(1:ngrid)
real*8 :: ANORMB,RCONDB,DLANGB,ROWCNDB,COLCNDB,AMAXB!,WORKB(1:4*ngrid)
integer*4 :: i,k,p,q,ip,kp,iq,kq,ips,ps,kps,pf,rsym,ssym
real*8 :: Sr,Srr,Srrr,Srrrr,Ss,Sss,Ssss,Sssss,Srs,Srss,Srrs,Srrss
real*8 :: Sx,Sxx,Sxxx,Sxxxx
real*8 :: Sz,Szz,Szzz,Szzzz
real*8 :: Sxz,Sxxz,Sxzz,Sxxzz
real*8 :: visf
real*8 :: term1,term2,term3,term4,term5,term6

    if (i.eq.0) then
     Ss=D1(k)/ds
     Sss=D2(k)/ds2
     Ssss=D3(k)/ds3
     Sssss=D4(k)/ds4
    else
     Ss=0.d0
     Sss=0.d0
     Ssss=0.d0
     Sssss=0.d0
    end if
    if (ssym.eq.1) then
     if (Vbc.eq.0) then     !!antisymmetry for free-slip
      Ss=-Ss
      Sss=-Sss
      Ssss=-Ssss
      Sssss=-Sssss
     elseif (Vbc.eq.1) then !!symmetry for rigid BCs -- does not appear to matter
      Ss=Ss
      Sss=Sss
      Ssss=Ssss
      Sssss=Sssss
     end if
    end if
    if (k.eq.0) then
     Sr=D1(i)/dr
     Srr=D2(i)/dr2
     Srrr=D3(i)/dr3
     Srrrr=D4(i)/dr4
    else
     Sr=0.d0
     Srr=0.d0
     Srrr=0.d0
     Srrrr=0.d0
    end if
    if (rsym.eq.1) then
     Sr=-Sr
     Srr=-Srr
     Srrr=-Srrr
     Srrrr=-Srrrr     
    end if
    Sx=Sr/xr(ip)
    Sxx=(Srr-Sx*xrr(ip))/xr(ip)**2.d0
    Sxxx=(Srrr-3.d0*Sxx*xr(ip)*xrr(ip)-Sx*xrrr(ip))/xr(ip)**3.d0
    Sxxxx=(Srrrr-6.d0*Sxxx*xrr(ip)*xr(ip)**2.d0-(3.d0*xrr(ip)**2.d0+4.d0*xr(ip)*xrrr(ip))*Sxx-Sx*xrrrr(ip))/xr(ip)**4.d0

    Sz=Ss/zs(kp)
    Szz=(Sss-Sz*zss(kp))/zs(kp)**2.d0
    Szzz=(Ssss-3.d0*Szz*zs(kp)*zss(kp)-Sz*zsss(kp))/zs(kp)**3.d0
    Szzzz=(Sssss-6.d0*Szzz*zss(kp)*zs(kp)**2.d0-(3.d0*zss(kp)**2.d0+4.d0*zs(kp)*zsss(kp))*Szz-Sz*zssss(kp))/zs(kp)**4.d0

    Srs=Drs(i,k)/drds
    Srss=Drss(i,k)/drds2
    Srrs=Drrs(i,k)/dr2ds
    Srrss=Drrss(i,k)/dr2ds2
    if (Vbc.eq.0) then
     if ((rsym.eq.1).neqv.(ssym.eq.1)) then  !!!!add rigid conditions
      Srs=-Srs
      Srss=-Srss
      Srrs=-Srrs
      Srrss=-Srrss
     end if
    elseif (Vbc.eq.1) then
     if (rsym.eq.1) then
      Srs=-Srs
      Srss=-Srss
      Srrs=-Srrs
      Srrss=-Srrss
     end if
    end if
    Sxz=Srs/xr(ip)/zs(kp)
    Sxxz=(Srrs-Sxz*zs(kp)*xrr(ip))/zs(kp)/xr(ip)**2.d0
    Sxzz=(Srss-Sxz*xr(ip)*zss(kp))/xr(ip)/zs(kp)**2.d0
    Sxxzz=(Srrss-Sxzz*xrr(ip)*zs(kp)**2.d0-Sxxz*zss(kp)*xr(ip)**2.d0-Sxz*xrr(ip)*zss(kp))/(xr(ip)*zs(kp))**2.d0

    if ((ip.eq.0).or.(ip.eq.nx)) then      !!free-slip
!     Matrix(p,q)=Sxx
      if (tstep.eq.-1) then
       call bandwidth(p,q)
      else
       MatrixB(subD+supD+1+p-q,q)=MatrixB(subD+supD+1+p-q,q)+Sxx
      end if
    elseif ((kp.eq.0).or.(kp.eq.nz)) then
     if (Vbc.eq.0) then     !!free-slip
!     Matrix(p,q)=Szz
      if (tstep.eq.-1) then
       call bandwidth(p,q)
      else
       MatrixB(subD+supD+1+p-q,q)=MatrixB(subD+supD+1+p-q,q)+Szz
      end if
     elseif (Vbc.eq.1) then !!rigid
!      if (tstep.eq.-1) then
!       call bandwidth(p,q)
!      else
!       MatrixB(subD+supD+1+p-q,q)=MatrixB(subD+supD+1+p-q,q)-1.d0
!      end if
     end if
    else
!     Matrix(p,q)=(visxx-viszz)*(Sxx-Szz)+vis(ip,kp)*(Sxxxx+Szzzz)+2.d0*vis(ip,kp)*Sxxzz+4.d0*visxz*Sxz
      if (tstep.eq.-1) then
       call bandwidth(p,q)
      else
       term1=(visxx-viszz)*(Sxx-Szz)
       term2=vis(ip,kp)*(Sxxxx+Szzzz)
       term3=2.d0*vis(ip,kp)*Sxxzz
       term4=4.d0*visxz*Sxz
       term5=2.d0*visx*(Sxxx+Sxzz)
       term6=2.d0*visz*(Szzz+Sxxz)
       MatrixB(subD+supD+1+p-q,q)=MatrixB(subD+supD+1+p-q,q)+term1+term2+term3+term4+term5+term6
      end if
    end if
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine bandwidth(p,q)
 use basics
implicit none
integer*4 :: diff,p,q
diff=q-p
if (diff.lt.subD) then
 subD=diff
elseif (diff.gt.supD) then
 supD=diff
end if
end

subroutine solve_SF_equation
!$ use OMP_LIB
 use basics
 use arrays
implicit none
 character*1 :: TRANS
integer*4 INFO,NRHS
real*8 :: time1,time2
!$ time1=omp_get_wtime()
 call buoyancy
 TRANS='N'
 NRHS=1
! call DGETRS(TRANS,ngrid,NRHS,Matrix(1:ngrid,1:ngrid),ngrid,IPIV(1:ngrid),Bvec(1:ngrid),ngrid,INFO)
 call DGBTRS(TRANS,ngrid,subD,supD,NRHS,MatrixB,2*subD+supD+1,IPIVB,Bvec,ngrid,INFO)
 call Stream_grid
!$ time2=omp_get_wtime()
!$ tsSF=tsSF+time2-time1
end

subroutine compute_velocities
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k
real*8 :: SFr,SFs,SFx,SFz
if (gridX.eq.0) then

!$OMP PARALLEL DO PRIVATE(i,k)
do k=0,nz
 do i=0,nx
  w(i,k)=-dot_product(D1(-span1:span1),SF(i-span1:i+span1,k))/dr
 end do
end do
!$OMP END PARALLEL DO

else

!$OMP PARALLEL DO PRIVATE(i,k,SFr,SFx)
do k=0,nz
 do i=0,nx
  SFr=dot_product(D1(-span1:span1),SF(i-span1:i+span1,k))/dr
  SFx=SFr/xr(i)
  w(i,k)=-SFx
 end do
end do
!$OMP END PARALLEL DO

end if

if (gridZ.eq.0) then

!$OMP PARALLEL DO PRIVATE(i,k)
do k=0,nz
 do i=0,nx
  u(i,k)=dot_product(D1(-span1:span1),SF(i,k-span1:k+span1))/ds
 end do
end do
!$OMP END PARALLEL DO

else

!$OMP PARALLEL DO PRIVATE(i,k,SFs,SFz)
do k=0,nz
 do i=0,nx
  SFs=dot_product(D1(-span1:span1),SF(i,k-span1:k+span1))/ds
  SFz=SFs/zs(k)
  u(i,k)=SFz
 end do
end do
!$OMP END PARALLEL DO

end if
 call velocityBCs
end

subroutine compute_velocities_coarse(level)
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k,level,i1,k1
real*8 :: SFr,SFs,SFx,SFz
if (gridX.eq.0) then
!$OMP PARALLEL DO PRIVATE(i,k)
do k=0,nz_grid(level)
 do i=0,nx_grid(level)
  w_grid(level,i,k)=-dot_product(D1(-span1:span1),error(i-span1:i+span1,k,level))/dr_grid(level)
 end do
end do
!$OMP END PARALLEL DO
else
!$OMP PARALLEL DO PRIVATE(i,k,SFr,SFx,i1)
do k=0,nz_grid(level)
 do i=0,nx_grid(level)
 i1=i*2**(level-1) !!fine grid coordinates for grid metrics (may need to restrict the metrics eventually)
  SFr=dot_product(D1(-span1:span1),error(i-span1:i+span1,k,level))/dr_grid(level)
  SFx=SFr/xr(i1)
  w_grid(level,i,k)=-SFx
 end do
end do
!$OMP END PARALLEL DO
end if

if (gridZ.eq.0) then
!$OMP PARALLEL DO PRIVATE(i,k)
do k=0,nz_grid(level)
 do i=0,nx_grid(level)
  u_grid(level,i,k)=dot_product(D1(-span1:span1),error(i,k-span1:k+span1,level))/ds_grid(level)
 end do
end do
!$OMP END PARALLEL DO
else
!$OMP PARALLEL DO PRIVATE(i,k,SFs,SFz,k1)
do k=0,nz_grid(level)
k1=k*2**(level-1)
 do i=0,nx_grid(level)
  SFs=dot_product(D1(-span1:span1),error(i,k-span1:k+span1,level))/ds_grid(level)
  SFz=SFs/zs(k1)
  u_grid(level,i,k)=SFz
 end do
end do
!$OMP END PARALLEL DO
end if
 call velocityBCs_coarse(level)
end

