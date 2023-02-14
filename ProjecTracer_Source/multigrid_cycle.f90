subroutine multigrid_setup
use basics
implicit none
 call optimize_grid(nx,nz,grids,mingrid)

allocate(multi(0:51)) !!max value of nsteps assumed is 50 -- modify of needed

!!!!!!!!! Multigrid cycle definitions
if (grids.eq.2) then
 nsteps=2
 multi(1:nsteps)=(/ 1,2 /) !!standard V cycle
elseif (grids.eq.3) then
! nsteps=4
! multi(1:nsteps)=(/ 1,2,3,2 /) !!standard V cycle

 nsteps=6
 multi(1:nsteps)=(/ 1,2,3,2,3,2 /) !!standard W cycle

! nsteps=8
! multi(1:nsteps)=(/ 1,2,3,2,3,2,3,2 /) !!W cycle with extra coarse grid iterations

! nsteps=12
! multi(1:nsteps)=(/ 1,2,3,2,3,2,3,2,3,2,3,2 /) !!W cycle with extra coarse grid iterations

! nsteps=14
! multi(1:nsteps)=(/ 1,2,3,2,3,2,3,2,3,2,3,2,3,2 /) !!W cycle with extra coarse grid iterations

! nsteps=16
! multi(1:nsteps)=(/ 1,2,3,2,3,2,3,2,3,2,3,2,3,2,3,2 /) !!W cycle with extra coarse grid iterations

! nsteps=18
! multi(1:nsteps)=(/ 1,2,3,2,3,2,3,2,3,2,3,2,3,2,3,2,3,2 /) !!W cycle with extra coarse grid iterations

! nsteps=20
! multi(1:nsteps)=(/ 1,2,3,2,3,2,3,2,3,2,3,2,3,2,3,2,3,2,3,2 /) !!W cycle with extra coarse grid iterations
elseif (grids.eq.4) then
! nsteps=6
! multi(1:nsteps)=(/ 1,2,3,4,3,2 /) !!standard V cycle

! nsteps=8
! multi(1:nsteps)=(/ 1,2,3,4,3,4,3,2 /) !!cheap W

! nsteps=10
! multi(1:nsteps)=(/ 1,2,3,4,3,4,3,4,3,2 /) !!cheap W

! nsteps=12
! multi(1:nsteps)=(/ 1,2,3,4,3,4,3,4,3,4,3,2 /) !!cheap W

! nsteps=14
! multi(1:nsteps)=(/ 1,2,3,4,3,4,3,4,3,4,3,4,3,2 /) !!cheap W

! nsteps=12
! multi(1:nsteps)=(/ 1,2,3,4,3,4,3,2,3,4,3,2 /) !!standard F cycle

! nsteps=12
! multi(1:nsteps)=(/ 1,2,3,4,3,2,3,4,3,4,3,2 /) !!backwards F cycle

 nsteps=14
 multi(1:nsteps)=(/ 1,2,3,4,3,4,3,2,3,4,3,4,3,2 /) !!standard W cycle

! nsteps=16
! multi(1:nsteps)=(/ 1,2,3,4,3,2,3,4,3,4,3,2,3,4,3,2 /) !!weird W

! nsteps=18
! multi(1:nsteps)=(/ 1,2,3,4,3,4,3,4,3,2,3,4,3,4,3,4,3,2 /) !!W cycle with extra coarse grid iterations

! nsteps=26
! multi(1:nsteps)=(/ 1,2,3,4,3,4,3,4,3,2,3,4,3,4,3,4,3,2,3,4,3,4,3,4,3,2 /) !!W cycle with extra coarse grid iterations
elseif (grids.eq.5) then
! nsteps=8
! multi(1:nsteps)=(/ 1,2,3,4,5,4,3,2 /) !!standard V cycle

! nsteps=20
! multi(1:nsteps)=(/ 1,2,3,4,5,4,5,4,3,4,5,4,3,2,3,4,5,4,3,2 /) !!standard F cycle

 nsteps=26
 multi(1:nsteps)=(/ 1,2,3,4,5,4,5,4,3,4,5,4,3,2,3,4,5,4,3,4,5,4,5,4,3,2 /) !!standard W cycle

! nsteps=44
! multi(1:nsteps)=(/ 1,2,3,4,5,4,5,4,5,4,3,4,5,4,3,4,5,4,3,2,3,4,5,4,3,2,3,4,5,4,3,4,5,4,3,4,5,4,5,4,5,4,3,2 /) !!W cycle with double the coarse grid passes

! nsteps=44
! multi(1:nsteps)=(/ 1,2,3,4,5,4,5,4,3,4,5,4,3,2,3,4,5,4,3,4,5,4,5,4,5,4,3,4,5,4,3,2,3,4,5,4,3,4,5,4,5,4,3,2 /) !!two successive W cycles
elseif (grids.eq.6) then
! nsteps=10
! multi(1:nsteps)=(/ 1,2,3,4,5,6,5,4,3,2 /) !!standard V cycle

! nsteps=30
! multi(1:nsteps)=(/ 1,2,3,4,5,6,5,6,5,4,5,6,5,4,3,4,5,6,5,4,3,2,3,4,5,6,5,4,3,2 /) !!standard F cycle

 nsteps=42
 multi(1:nsteps)=(/ 1,2,3,4,5,6,5,6,5,4,5,6,5,4,3,4,5,6,5,4,3,2,3,4,5,6,5,4,3,4,5,6,5,4,5,6,5,6,5,4,3,2 /) !!standard W cycle
end if
!!!!!!!!! End Multigrid cycle definitions

multi(0)=multi(nsteps)  !!enforce periodicity
multi(nsteps+1)=multi(1)
allocate(residual_mag(1:nsteps))
end

