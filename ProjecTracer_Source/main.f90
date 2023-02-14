program ProjecTracer
!$ use OMP_LIB
 use basics
 use arrays
implicit none
real*8 :: t1,t2,tmin
integer*4 :: i,j,k,N,level,v_update,RK_update_save
real*8 :: xH,zH,tH,heating !!testing variable H
real*16 :: a_qp,b_qp

open(unit=501,file="inputs.nml")
read(nml=inputs,unit=501)
N=Nlocal
tmin=tend/2.d0          !!if local_time=0, begin time averaging at this time for stats
if (comp.eq.1) allocate(RaC(1:ntype)) !!allocate array prior to reading value from namelist
read(nml=Rayleigh_numbers,unit=501)
if (comp.eq.1) allocate(conduct_factor(1:ntype),Htype(1:ntype))
read(nml=heating_conductivity,unit=501)
if (comp.eq.1) allocate(visC(1:ntype),delta_visC(1:ntype))
read(nml=viscosity,unit=501)
close(501)

!!!quad precision test
a_qp=1.000000000000000000007q0; b_qp=2.00000000000000000000003q0
write(*,*) "Quad precision test:"
write(*,*) a_qp,b_qp,a_qp+b_qp
!!!end quad precision test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BEGIN INPUTS
!!!!!!!!!!!!! grid options
!aspect=1.d0       !!aspect ratio of Cartesian domain
!nx=200            !!# of cells in the horizontal
!nz=200            !!# of cells in the vertical
!Database=2        !!ordering of streamfunction and buoyancy vectors -- 1 (large nz), 2 (large nx), or 3 (experiment)
!order=4           !!space order: even numbers
!gridX=0           !!0=uniform,1=arctan
!gridZ=0           !!0=uniform,1=arctan,2=boundary layer
!thickness=0.1d0   !!gridZ=2: thickness of grid boundary layer
!nbound=8          !!# of cell corners in boundary layer grid
!spaceX=1.d0       !!stretching factor for arctan function
!spaceZ=1.d0       !!stretching factor for arctan and boundary layer grids
!!!!!!!!!!!!! end grid options

!!!!!!!!!!!!!velocity solver options
!iterate=0         !!how to solve for the velocities?: 0=direct method, 1=iterative method
!MG=1              !!multigrid?
!non_Newtonian=0   !!is the fluid non-Newtonian? (only works if iterate=1)
!tolerance=1.d-3   !!tolerance for multigrid -- smaller values indicate more accurate solutions
!Nmulti=1!70         !!# of iterations on each multigrid level (~70 for Newtonian flows)
!mingrid=4         !!min grid size (<nx and <nz)
!grids=2           !!max # of multigrid levels (currently between 2 and 6)
!eig_ratio=1.d0    !!min ratio of diffusive to advective eigenvalues desired on coarse grids
!vis_iter=12       !!max # of viscosity smoother iterations for multigrid
!max_cycles=1000   !!max # of multigrid cycles allowed for convergence before code stops
!tol_coarse=1.d-3  !!tolerance for fixed point iterations of coarse equation
!fp_max=100        !!max # of fixed point iterations done on the coarse equation for non-Newtonian flows
!splitting=0       !!operator splitting in multigrid smoother?
!RKC1=1            !!Time stepping method for fine grid smoother: 0=use RK1, 1=use RKC1
!nstageS=35         !!# of stages in RKC1 fine grid smoother
!dampingS=1000.d0     !!damping factor for RKC1 fine grid smoother
!c_2=0.95d0          !!courant # range for RKC1 smoother
!c_1=0.9d0
!!!!!!!!!!!!!end velocity solver options

!!!!!!!!!!!!! time stepping options
!tend=1.d-10              !!total nondimensional time 
!tsnap=tend/1.d0         !!output snapshots are taken at this time interval
!Chebyshev=0             !!0=classical RK, 1=RK Chebyshev (RKC) method (torder=2 and tracer=1 required below)
!nstage=3                !!equals # of stages used in the RKC method
!damping=100.d0          !!damping parameter for RKC method
!aRK2=0.5d0              !!alpha value for RK2 time integration and first stage of RKC
!torder=4                !!order of RK method: RK2 only if temperature tracers are used - RK4 for field method and/or tracers (2 and 4 currently available)
!movie=0                 !!create movie files? (under development...)
!tmovie_snap=tend/100.d0 !!time intervals for movie frames
!tmin=tend/2.d0          !!if local_time=0, begin time averaging at this time for stats
!N=1000                  !!if local time steps ARE used, run controlled by iteration cutoff
!Nsnap=100               !!if local time steps ARE used, snapshots controlled by iteration count
!local_time=0            !!only use with RK4 without updates (for field methods primarily)
!dt_manual=0.125d0          !!manual time step size for isothermal convection
!courant=0.99d0!2.25d0          !!factor that multiplies the time step size computed from stability analysis (0.99 for RK2, 2.25 for RKC)
!courantT=0.99d0         !!for temperature smoother
!courantC=0.99d0         !!for composition smoother
!courant_stream=0.99d0   !!for iterative solver
!RK_update=1             !!update the velocities between RK stages? (use 0 for steady state problems, 1 for time accurate runs, 2 for 2 updates during RK4)
!vel_update=1            !!frequency of velocity updates: 0=every other time step, 1=every time step
!CFL=0                   !!use velocity-dependent time steps? (for RK1 only)
!!!!!!!!!!!!! end time stepping options

!!!!!tracer options
!comp=1                  !!switch for composition: 0=off,1=on
!ntype=1                 !!# of enriched compositions:=0 for uniform composition with T tracers, >=1 for heterogenious composition
!iC_continuous=1         !!continuous C values carried by tracers? (for iC_continuous=1 equalize=0 AND ntype=1 are required)
!original=1              !!use original TRA? (short term isoviscous calculations)
!shape_function=1        !!0 for constant cell, 1 for bilinear
!buffer=10               !!buffer range for bilinear shape functions (used in TRA)
!equalize=0              !!redistribute tracers?
!tracer=1                !!temperature tracers?
!tpc=60                  !!# of tracers per cell corner
!dist_factor=1.d-4       !! <<1 controls min. dist between tracers and side boundaries during equalize
!!!!!end tracer options

!!!!!!!!!!!!!!!!!!!!Initial and Boundary Conditions
!restart=0               !!0= use analytic initial T and C, 1= use T (and possibly C -- see "add_layer" tracer option) from previous run
!initial_temp=-1         !!initial T: -1=user defined (initial_conditions.f90), 0 for conductive, 1 for convective cell (TBL scales with RaT), 2 alternate conductive
!initial_comp=-1         !!initial C: -1=user defined (initial_conditions.f90), 0=dense layer at base (dlayer), 1=sinusoidal deflection, 2=step function
!initial_temp_smooth=0   !!smooth the initial temperature field? (controlled by delta_visT_init)
!delta_visT_init=2.0d0   !!controls smoothing of the initial T field based upon viscosity dependence on T
!add_layer=0             !!add heterogenious composition to a previous isocompositional run when restarting?
!dlayer=0.5d0          !!for initial_comp=0 or 1: layer thickness
!dent=0.5d0              !!reference height for entrainment calculation
!klayer=35.d0            !!controls thickness of initial compositional interface (under development)

!! reflecting sidewalls assumed
!Tbot=0 !!Bottom thermal BC: 0=isothermal, 1=insulating
!Vbc=0  !!Top and Bottom velocity BCs: 0=free-slip,1=rigid
!!!!!!!!!!!!!!!!!!!!End Initial and Boundary Conditions

!!!!!!!!!!!!!Rayleigh Numbers and Buoyancy
!RaT=1.d5                !!Thermal Rayleigh number
!if (comp.eq.1) then
! allocate(RaC(1:ntype))
! RaC(1)=0.5d0*RaT           !!Compositional Rayleigh Number(s): >0 for increased density, <0 for decreased density
!! RaC(2)=-0.2d0*RaT           !!create as many terms as there are enriched compositions
!! RaC(3)=-0.4d0*RaT
!end if
!!!!!!!!!!!!!End Rayleigh Numbers and Buoyancy

!!!!!!!!!!Internal heating and Thermal Conductivity Options
!iheattype=1             !!=0 for constant (in each compositional component); =1 for variable (specify in energy.f90, local_time=0 only for now)
!H=0.d0                  !!internal heating rate of ambient material (for iheattype=0, i.e., constant H )
!if (comp.eq.1) then
! allocate(conduct_factor(1:ntype),Htype(1:ntype))
! conduct_factor(1)=1.d0!2.d0  !!multiplier for nondimesnional conductivity within each compositional component (C is smoothed by delta_visC)
!! conduct_factor(2)=1.d0  !!create as many terms as there are enriched compositions
!! conduct_factor(3)=0.5d0  
!
! Htype(1)=0.d0            !!set H value within components
!! Htype(2)=0.d0           !!create as many terms as there are enriched compositions
!! Htype(3)=20.d0
!end if
!!!!!!!!!!!Internal heating and Thermal Conductivity Options

!!!!!!!!!viscosity options
!plastic=0                      !!plastic yielding: -1=development,0=none,1=Tosi et al., (2015), 2=Moresi and Solomatov (1998), 3=smoothed minimum
!ivis=0                         !!viscosity: 0=isoviscous, 1=T and z dependent, 2=dependence on z, T and C
!visP=1.d0                      !!for depth dependence (>1 for an increase with depth)
!visT=1.d3!1.d6                      !!for temperature dependence (>1 for decrease)
!delta_visT=2.d0                !!max allowed viscosity contrast due to T between adjacent grid points (T used to compute viscosity field is smoothed until this criterion is satisfied)
!eig_ratio_min=2.0d0            !!min ratio of diffusive to advective eigenvalues in streamfunction eqn (viscosity will be smoothed until this criterion is also satisfied)
!if (comp.eq.1) then            !!composition dependence (>1 for decrease, <1 for increase)
! allocate(visC(1:ntype),delta_visC(1:ntype))
! visC(1)=1.d1 
! !visC(2)=1.d0  !!create as many terms as there are enriched compositions 
! !visC(3)=1.d-2
! delta_visC(1)=1.d0            !!max allowed viscosity/conductivity contrast due to C between adjacent grid points (multiplicative factor)
! !delta_visC(2)=2.d0
! !delta_visC(3)=2.d0
!end if
!!plastic yielding options:
!yield=1.d10                   !!yield stress at the surface (at zero strain-rate invariant) (plastic.ne.0)
!yield_gradient=0.d0           !!yield stress gradient (linear increase with depth) (plastic.ne.0)
!vis_star=0.d0                 !!constrains min. plastic viscosity (plastic=1)
!mpla=0.052822d0               !!exponent for strain-rate dependent yield stress (plastic=-1)
!Cpla=2.5453d0                 !!coefficient for strain-rate dependent yield stress (plastic=-1)
!alpha=-10.d0                  !!for smooth min function (plastic=3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! END INPUTS

 if (equalize.eq.1) open(unit=1024,file="equalize.dat")
 call compute_smoothing_parameters
 if (iterate.eq.1) call multigrid_setup
 call compute_basics
 call storage
 if (RKC1.eq.1) then
  call Chebyshev_prep_stream
 end if
 if (Database.eq.3) call Database3
 call FD_coefficients
 call build_grid
 !!!test for H
 !if (iheattype.eq.1) then
 ! i=11; k=11; xH=xg(i); zH=zg(k); tH=0.01d0;
 ! call compute_heating(i,k,tH,heating)
 ! write(*,*) "xH,zH,tH=",xH,zH,tH
 ! write(*,*) "heating=",heating
 !end if
 !!!end test for H
 call derivative_eigenvalues
 call stability_RK4_smoother
 if (iterate.eq.1) then
  call eigenvalue_prep !!!for stream function related eigenvalues
  if ((gridX.eq.0).and.(gridZ.eq.0)) call stream_space_matrix_prep
  do level=1,grids
   call stability_vis_smoother(level)
   if (non_Newtonian.eq.1) then
    call stability_stream_smoother(level)
    call define_NonNewtonian_viscosity(level)
   end if
  end do
 else
   grids=1
   level=1
   call stability_vis_smoother(level)
 end if
 if (restart.eq.1) then
  call smoother_time_T(T)
  if (comp.ne.0) then
   call smoother_time
  end if
 call load_conductivity_array
 end if
 if (restart.eq.0) call initial_temperature
 if ((restart.eq.0).and.((comp.eq.1).or.(tracer.eq.1))) then
  call initialize_tracers
 end if
!!!!!Initial temperature smoothing here
 if ((restart.eq.0).and.(initial_temp_smooth.eq.1)) call smoother_initial_T
 if (((restart.eq.1).and.(add_layer.eq.1)).and.(comp.eq.1)) then
  if (tracer.eq.0) then
   call initialize_tracers
  else
   call add_tracers
  end if
 end if
 tstep=-1
 if (iterate.eq.0) then
  call build_SF_matrix !!find bandwith and banded matrix size
 elseif (iterate.eq.1) then
  if (restart.eq.0) SF=0.d0 !!initial guess for stream function
  error(:,:,2:grids)=0.d0
 end if
! call tracer_snapshots !!!!time sink for large jobs
 call cpu_time(t1)
!$ t1=omp_get_wtime()
 open(unit=303,file="stats.dat")
 open(unit=304,file="viscosity_smoother.dat")
fcount=0
fmovie=0
tmovie=0.d0
Nout=Nsnap
tout=0.d0
time=0.d0
v_update=1
RK_update_save=RK_update
if (local_time.eq.0) N=1000000 !!!max iterations for time accurate run
do tstep=0,N
 if (v_update.eq.1) then !!!!!!!!!!!!!!!!!!!!!!!!!new
  if (iterate.eq.0) then
   if ((ivis.ge.1).or.(tstep.eq.0)) then
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
  if (vel_update.eq.0) then
   RK_update=0
   v_update=0 !!deactivate velocity solver on next time step
  end if
 else
  RK_update=RK_update_save
  v_update=1 !!reactivate velocity solver on next time step
 end if !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end new
 call stats
 if ((movie.eq.1).and.(time.gt.tmovie)) then
  call movie_frames
  tmovie=tmovie+tmovie_snap
  fmovie=fmovie+1
 end if
 if (local_time.eq.0) then  !!snapshots controlled by time
  if (time.ge.tout) then
   call snapshots
   call print_restart !!!!!!!!!!test: for jobs that run out of wall time
   tout=tout+tsnap
   fcount=fcount+1
   if (time.ge.tend) exit
  end if
 else                       !!snapshots controlled by iteration
  if (Nout.ge.Nsnap) then
   call snapshots
   call print_restart !!!!!!!!!!test: for jobs that run out of wall time
   Nout=1
   fcount=fcount+1
  else
   Nout=Nout+1
  end if
 end if
 if (Chebyshev.eq.1) then
  call energy_time_Chebyshev
 else
  call energy_time
 end if
 if (local_time.eq.0) then
  time=time+dt
 end if
end do
 close(303)
 close(304)
 call cpu_time(t2)
!$ t2=omp_get_wtime()
 write(*,*) "Run Time=",t2-t1,"s"
 write(*,*) "energy_time=",tenergy,"s"
 if ((tracer.eq.1).or.(comp.eq.1)) write(*,*) "advect_tracer=",tadvect_tracer,"s"
 if (iterate.eq.0) then
  write(*,*) "build_SF_matrix=",tbSF,"s"
  write(*,*) "solve_SF_equation=",tsSF,"s"
 elseif (iterate.eq.1) then
  write(*,*) "Multigrid=",tmultigrid,"s"
  write(*,*) " Compute_Matrix_coarse=",tcoarse_matrix,"s"
  write(*,*) " LU factors=",t_lu,"s"
  write(*,*) " Forward solves=",t_fs,"s"
  if (non_Newtonian.eq.1) write(*,*) " Matrix Multiply=",t_mult,"s"
  write(*,*) "iterate_stream=",t_istream,"s"
  write(*,*) " A=",t_isA,"s"
  write(*,*) "  A1=",t_isA1,"s"
  write(*,*) "  A2=",t_isA2,"s"
  write(*,*) "   A2A=",t_cs1,"s"
  write(*,*) "   A2B=",t_cs2,"s"
  write(*,*) "   A2C=",t_cs3,"s"
  write(*,*) "  A3=",t_isA3,"s"
  write(*,*) "   A3A=",t_cvA,"s"
  write(*,*) "   A3B=",t_cvB,"s"
  write(*,*) "   A3C=",t_cvC,"s"
  write(*,*) " B=",t_isB,"s"
  write(*,*) " C=",t_isC,"s"
  write(*,*) "smoother_vis=",t_sv,"s"
  write(*,*) " A=",t_svA,"s"
  write(*,*) " B=",t_svB,"s"
  write(*,*) "  BA=",t_svBA,"s"
  write(*,*) "  BB=",t_svBB,"s"
  write(*,*) "  BC=",t_svBC,"s"
  write(*,*) " C=",t_svC,"s"
  write(*,*) " D=",t_svD,"s"
 end if
 if ((tracer.eq.1).or.(comp.eq.1)) then
  write(*,*) "compute_derivatives=",tcd,"s"
  write(*,*) "tracers_to_corners=",tconvert,"s"
 end if
 if (equalize.eq.1) then
  write(*,*) "equalize=",tequalize,"s"
 end if
  write(*,*) "time steps=",tstep
 call print_restart
! call tracer_snapshots
 if ((local_time.eq.0).and.(RaT.ne.0.d0).and.(tstep.gt.10)) call time_average(tmin,tend)
end program ProjecTracer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine movie_frames
 use basics
 use arrays
implicit none
integer*4 :: i,k
character*6 :: fname
write(fname,'(a1,i5)') "f",10000+fmovie
open(unit=90210,file="movie_files/"//fname)
if (comp.eq.0) then
 do i=0,nx
  do k=0,nz
   write(90210,'(f7.3)') T(i,k)
  end do
 end do
else
 do i=0,nx
  do k=0,nz
   write(90210,'(10(f7.3))') T(i,k),Cnew(1:ntype,i,k)
  end do
 end do
end if
close(90210)
end

subroutine snapshots
 use basics
 use arrays
implicit none
real*8 :: integrate1Dx
real*8 :: Tlat,strainlat,stresslat
integer*4 :: i,k,id
 character*6 :: fname
write(fname,'(a1,i5)') 'T',10000+fcount
open(unit=101,file=fname)
if (RaT.ne.0.d0) then
 do i=0,nx
  do k=0,nz
   if ((comp.eq.0).and.(tracer.eq.0)) then
    write(101,'(8(g14.5E3))') xg(i),zg(k),T(i,k),SF(i,k),u(i,k),w(i,k),vis(i,k),Tvis(i,k)
   elseif ((comp.eq.0).and.(tracer.eq.1)) then
    write(101,'(9(g14.5E3),i6)') xg(i),zg(k),T(i,k),SF(i,k),u(i,k),w(i,k),2.d0*vis(i,k)*strain(i,k,1),vis(i,k),Tvis(i,k),&
&nempty(i,k)
   else
 write(101,'(9(g14.5E3),i6,20(g14.5E3))') xg(i),zg(k),T(i,k),SF(i,k),u(i,k),w(i,k),2.d0*vis(i,k)*strain(i,k,1),vis(i,k),Tvis(i,k)&
&,nempty(i,k),Cnew(1:ntype,i,k),Cvis(1:ntype,i,k)
   end if
  end do
  write(101,*) " "
 end do
else                !!!!isothermal convection (pure compositional convection)
 do i=0,nx
  do k=0,nz
 write(101,'(7(g14.5E3),i6,20(g14.5E3))') xg(i),zg(k),SF(i,k),u(i,k),w(i,k),2.d0*vis(i,k)*strain(i,k,1),vis(i,k),nempty(i,k),&
&Cnew(1:ntype,i,k),Cvis(1:ntype,i,k)
  end do
  write(101,*) " "
 end do
end if
 close(101)

!!!!!!!!!!!!!!!!! lateral avg profiles
write(fname,'(a1,i5)') 'P',10000+fcount
open(unit=101,file=fname)
do k=0,nz
 Tlat=integrate1Dx(T(0:nx,k))/aspect
 strainlat=integrate1Dx(strain(0:nx,k,1))/aspect
 stresslat=integrate1Dx(2.d0*vis(0:nx,k)*strain(0:nx,k,1))/aspect !!deviatoric stress invariant
 write(101,'(4(g14.5E3))') zg(k),Tlat,strainlat,stresslat
end do
 close(101)
!!!!!!!!!!!!!!!!! end lateral avg profiles
end

subroutine tracer_snapshots
 use basics
 use arrays
implicit none
real*8 :: xgrid,zgrid
integer*4 :: id
if ((comp.eq.1).or.(tracer.eq.1)) then
 if (tstep.eq.-1) then
  open(unit=666,file="tracers.dat")
 else
  open(unit=666,file="tracers.dat",position='append')
 end if
 if ((comp.eq.1).and.(tracer.eq.1)) then
    do id=1,ntr
     write(666,'(3(f9.5),i4)') xgrid(rtr(id)),zgrid(str(id)),Ttr(id),Ttype(id)
    end do
 elseif (comp.eq.1) then
    do id=1,ntr
     write(666,'(2(f9.5),i4)') xgrid(rtr(id)),zgrid(str(id)),Ttype(id)
    end do
 elseif (tracer.eq.1) then
    do id=1,ntr
     write(666,'(3(f9.5))') xgrid(rtr(id)),zgrid(str(id)),Ttr(id)
    end do
 end if
    write(666,*) " "
    write(666,*) " "
    close(666)
end if
end

subroutine print_restart
 use basics
 use arrays
implicit none
integer*4 :: i,k,ID
open(unit=101,file="Trestart",form='unformatted')
if (RaT.ne.0.d0) then
 if (comp.eq.0) then
  write(101) xg,zg,T,SF
 else
  write(101) xg,zg,T,Cnew,SF
 end if
else
 if (comp.eq.0) then
  write(101) xg,zg,SF
 else
  write(101) xg,zg,Cnew,SF
 end if
end if
 close(101)
if ((comp.eq.1).or.(tracer.eq.1)) then
 open(unit=101,file="restart_tracers",form='unformatted')
  if (tracer.eq.0) then
   if (iC_continuous.eq.0) then
    write(101) rtr,str,Ttype
   elseif (iC_continuous.eq.1) then
    write(101) rtr,str,Ctr
   end if
  elseif (comp.eq.0) then
   write(101) rtr,str,Ttr
  else
   if (iC_continuous.eq.0) then
    write(101) rtr,str,Ttype,Ttr
   elseif (iC_continuous.eq.1) then
    write(101) rtr,str,Ctr,Ttr
   end if
  end if
 close(101)
end if
end


subroutine FD_coefficients
 use basics
 use arrays
implicit none
integer*4 :: i,m,n,bigN,bigM,k
real*8 :: x0,c1,c2,c3
real*8, allocatable :: a(:),delta(:,:,:)
 bigM=max(4,Dnumber)  !!need up to 4th derivatives in SF equation - may need more for tracer velocity interpolation
 bigN=order+bigM-2 !!assuming bigM is even
! bigN=order+2
 allocate(a(0:bigN),delta(-1:bigM,0:bigN,0:bigN))
 delta=0.d0
 x0=0.d0
 a(0)=0.d0
 k=1
 do i=1,bigN,2
  a(i)=real(k,8)
  a(i+1)=-real(k,8)
  k=k+1
 end do

 delta(0,0,0)=1.d0
 c1=1.d0
 do n=1,bigN
  c2=1.d0
  do i=0,n-1
   c3=a(n)-a(i)
   c2=c2*c3
   if (n.le.bigM) delta(n,n-1,i)=0.d0
   do m=0,min(n,bigM)
    delta(m,n,i)=((a(n)-x0)*delta(m,n-1,i)-real(m,8)*delta(m-1,n-1,i))/c3
   end do
  end do
  do m=0,min(n,bigM)
   delta(m,n,n)=(c1/c2)*(real(m,8)*delta(m-1,n-1,n-1)-(a(n-1)-x0)*delta(m,n-1,n-1))
  end do
 c1=c2
 end do

 D1=0.d0
 D2=0.d0
 D1(0)=delta(1,order,0)
 D2(0)=delta(2,order,0)
 k=1
 do i=1,order,2
  D1(k)=delta(1,order,i)
  D1(-k)=delta(1,order,i+1)
  D2(k)=delta(2,order,i)
  D2(-k)=delta(2,order,i+1)
 k=k+1
 end do

!!!!!!!!for stats integrals
 D(1:Dnumber,-span_interp:span_interp)=0.d0
 do m=2,Dnumber,2
  D(m-1,0)=delta(m-1,order+m-2,0)
  D(m,0)=delta(m,order+m-2,0)
  k=1
  do i=1,order+m-2,2
   D(m-1,k)=delta(m-1,order+m-2,i)
   D(m-1,-k)=delta(m-1,order+m-2,i+1)
   D(m,k)=delta(m,order+m-2,i)
   D(m,-k)=delta(m,order+m-2,i+1)
  k=k+1
  end do
 end do
 do m=1,Dnumber,2
  D(m,0)=0.d0     !!!enforce symemetry and antisymmetry properties
  do i=1,span_interp
   D(m,i)=-D(m,-i)
   D(m+1,i)=D(m+1,-i)
  end do
 end do
!!!!!!!!

 D3=0.d0
 D4=0.d0
 D3(0)=delta(3,order+2,0)
 D4(0)=delta(4,order+2,0)
 k=1
 do i=1,order+2,2
  D3(k)=delta(3,order+2,i)
  D3(-k)=delta(3,order+2,i+1)
  D4(k)=delta(4,order+2,i)
  D4(-k)=delta(4,order+2,i+1)
 k=k+1
 end do
 

 D1(0)=0.d0     !!!enforce symemetry and antisymmetry properties
 D3(0)=0.d0
 do i=1,span
  D1(i)=-D1(-i)
  D2(i)=D2(-i)
  D3(i)=-D3(-i)
  D4(i)=D4(-i)
 end do


 Drs=0.d0
 Drrss=0.d0
 Drss=0.d0
 Drrs=0.d0
 do k=-span1,span1
  do i=-span1,span1
   Drs(i,k)=D1(k)*D1(i)
   Drss(i,k)=D2(k)*D1(i)
   Drrs(i,k)=D1(k)*D2(i)
   Drrss(i,k)=D2(k)*D2(i)
  end do
 end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!second order FD stencils for viscosity gradients
D1_so=(/-0.5d0,0.d0,0.5d0/)
D2_so=(/1.d0,-2.d0,1.d0/)

Drs_so(:,1)=(/-0.25d0,0.d0,0.25d0/)
Drs_so(:,0)=(/0.d0,0.d0,0.d0/)
Drs_so(:,-1)=(/0.25d0,0.d0,-0.25d0/)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end second order FD stencils for viscosity gradients

! write(*,'(100(g12.3,x))') D1
! write(*,'(100(g12.3,x))') D1_so
! write(*,*) " "
! write(*,'(100(g12.3,x))') D2
! write(*,'(100(g12.3,x))') D2_so
! write(*,*) " "

!do i=span,-span,-1
! write(*,'(100(g12.3,x))') Drs(i,-span:span)
!end do
! write(*,*) " "
!do i=1,-1,-1
! write(*,'(100(g12.3,x))') Drs_so(i,-1:1)
!end do
! write(*,*) " "

! write(*,'(100(g12.3,x))') D3
! write(*,'(100(g12.3,x))') D4
! write(*,*) ""
! do m=1,Dnumber
!  write(*,'(100(g12.3,x))') D(m,-span_interp:span_interp)
! end do
end
