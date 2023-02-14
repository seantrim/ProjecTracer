!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Static Part of Viscosity (that does not depend on stream function/velocity/stress etc.)
real*8 function visf(i,k)  !!!function for static part of viscosity (T,C and z dependence)
 use basics
 use arrays
implicit none
integer*4 :: i,k,tt
real*8 :: Cfilter(1:ntype)
if (ivis.eq.1) then
 visf=dexp(-Tvis(i,k)*dlog(visT)+(1.d0-zg(k))*dlog(visP))
elseif ((ivis.eq.2).and.(comp.eq.1)) then !!viscosity increase with C
 if (RaT.ne.0.d0) then
  visf=dexp(-Tvis(i,k)*dlog(visT)-dot_product(Cvis(:,i,k),dlog(visC(:)))+(1.d0-zg(k))*dlog(visP))
 else
!  Cfilter(:)=0.5d0+0.5d0*dtanh(2.5d0*(Cvis(:,i,k)-0.5d0))/dtanh(1.25d0)
!  visf=dexp(-dot_product(Cfilter,dlog(visC))+(1.d0-zg(k))*dlog(visP))
  visf=dexp(-dot_product(Cvis(:,i,k),dlog(visC(:)))+(1.d0-zg(k))*dlog(visP))
 end if
elseif (ivis.eq.-1) then
 visf=dexp((1.d0-zg(k))*dlog(visP))
else
 visf=1.d0
end if
end

subroutine define_NonNewtonian_viscosity(level) !!define static portion of non-Newtonian (NN) viscosity (e.g., depth-dependent yield stress)
 use basics
 use arrays
implicit none
integer*4 :: i,k,kfine,level
do k=0,nz_grid(level)
 kfine=k*2**(level-1)
 do i=0,nx_grid(level)
  if (plastic.eq.1) then
   visNN_static(i,k,level)=(yield+yield_gradient*(1.d0-zg(kfine)))
  elseif (plastic.eq.2) then
   visNN_static(i,k,level)=(yield+yield_gradient*(1.d0-zg(kfine)))/2.d0
  elseif (plastic.eq.3) then
   visNN_static(i,k,level)=(yield+yield_gradient*(1.d0-zg(kfine)))/2.d0
  elseif (plastic.eq.-1) then
   visNN_static(i,k,level)=(yield+yield_gradient*(1.d0-zg(kfine)))
  else !!no plasticity
    visNN_static(i,k,level)=1.d99
  end if
 end do
end do
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Define effective viscosity for Non-Newtonian flows
subroutine viscosity_loop(level)
!$ use OMP_LIB
use basics
use arrays
implicit none
integer*4 :: i,k,level
real*8 :: plastic_part
if (plastic.eq.1) then !!Tosi et al. 2015
!$OMP PARALLEL DO PRIVATE(i,k)
 do k=0,nz_grid(level)
  do i=0,nx_grid(level)
     vis_grid(i,k,level)=2.d0/(1.d0/vis_grid_static(i,k,level)+1.d0/(yield_strain(i,k,level)+vis_star))
  end do
 end do
!$OMP END PARALLEL DO

elseif (plastic.eq.2) then !!Moresi and Solomatov 1998
!$OMP PARALLEL DO PRIVATE(i,k)
 do k=0,nz_grid(level)
  do i=0,nx_grid(level)
    vis_grid(i,k,level)=min(vis_grid_static(i,k,level),yield_strain(i,k,level))
  end do
 end do
!$OMP END PARALLEL DO

elseif (plastic.eq.3) then !!Smooth minimum
!$OMP PARALLEL DO PRIVATE(i,k)
 do k=0,nz_grid(level)
  do i=0,nx_grid(level)
    vis_grid(i,k,level)=(vis_grid_static(i,k,level)*dexp(alpha*vis_grid_static(i,k,level))+yield_strain(i,k,level)*&
    &dexp(alpha*yield_strain(i,k,level)))/(dexp(alpha*vis_grid_static(i,k,level))+dexp(alpha*yield_strain(i,k,level)))
  end do
 end do
!$OMP END PARALLEL DO

elseif (plastic.eq.-1) then
!$OMP PARALLEL DO PRIVATE(i,k,plastic_part)
 do k=0,nz_grid(level)
  do i=0,nx_grid(level)
!     plastic_part=vis_star+(Cpla*strain(i,k,level)**mpla+visNN_static(i,k,level))/strain(i,k,level)
     plastic_part=vis_star+yield_strain(i,k,level)
     vis_grid(i,k,level)=1.d0/(1.d0/vis_grid_static(i,k,level)+1.d0/plastic_part)
  end do
 end do
!$OMP END PARALLEL DO


else !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!no plasticity
!$OMP PARALLEL DO PRIVATE(i,k)
 do k=0,nz_grid(level)
  do i=0,nx_grid(level)
    vis_grid(i,k,level)=vis_grid_static(i,k,level)
  end do
 end do
!$OMP END PARALLEL DO
end if
end
