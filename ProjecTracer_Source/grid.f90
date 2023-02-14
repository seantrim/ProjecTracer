real*8 function xgrid(r)
 use basics
implicit none
real*8 :: r
if (gridX.eq.1) then
  xgrid=aspect*(Ax*datan(spaceX*(r/aspect-0.5d0))+0.5d0)
elseif (gridX.eq.2) then
 xgrid=r
else
 xgrid=r
end if
end

real*8 function zgrid(s)
 use basics
implicit none
real*8 :: s,s1,s2,dzbound,dzin
if (gridZ.eq.1) then
 zgrid=Az*datan(spaceZ*(s-0.5d0))+0.5d0
elseif (gridZ.eq.2) then
 dzbound=thickness/real(nbound-1,8)
 dzin=(1.d0-2.d0*thickness)/real(nz-2*Nbound+2,8)
 s1=real(nbound-1,8)*ds
 s2=1.d0-s1
 if (s.lt.s1) then
  zgrid=(thickness/(dexp(spaceZ*s1)-1.d0))*(dexp(spaceZ*s)-1.d0)
 elseif (s.gt.s2) then
  zgrid=1.d0-(thickness/(dexp(spaceZ*s1)-1.d0))*(dexp(-spaceZ*(s-1.d0))-1.d0)
 else
  zgrid=thickness+(s-s1)*dzin/ds
 end if
else
 zgrid=s
end if
end

subroutine build_grid
 use basics
 use arrays
implicit none
integer*4 :: i,k,ID,tt
real*8 :: r,s,xgrid,zgrid
real*8 :: xrf,xrrf,xrrrf,xrrrrf,xrn,xrrn,xrrrn,xrrrrn
real*8 :: xr_interp,zs_interp

if (restart.eq.0) then
 do i=0,nx
  r=dr*real(i,8)
  xg(i)=xgrid(r)
 end do
 do k=0,nz
  s=ds*real(k,8)
  zg(k)=zgrid(s)
 end do
else
 open(unit=101,file="Trestart",form='unformatted')
   if (RaT.ne.0.d0) then
    if ((comp.eq.0).or.(add_layer.eq.1)) then
     read(101) xg,zg,T,SF
    else
     read(101) xg,zg,T,Cnew,SF
    end if
   else
     read(101) xg,zg,Cnew,SF
   end if
 close(101)
 call enforceBCs(T)
 if (iterate.eq.1) call enforceBCs_stream(1,SF)
 if ((comp.eq.1).and.(add_layer.eq.0)) then
  do tt=1,ntype
   call enforceBCs_comp(Cnew(1:ntype,:,:))
  end do
 end if
 if (add_layer.eq.0) then
  if ((comp.eq.1).or.(tracer.eq.1)) then
   open(unit=101,file="restart_tracers",form='unformatted')
   if (tracer.eq.0) then
    if (iC_continuous.eq.0) then
     read(101) rtr,str,Ttype
    elseif (iC_continuous.eq.1) then
     read(101) rtr,str,Ctr
    end if
   elseif (comp.eq.0) then
    read(101) rtr,str,Ttr
   else
    if (iC_continuous.eq.0) then
     read(101) rtr,str,Ttype,Ttr
    elseif (iC_continuous.eq.1) then
     read(101) rtr,str,Ctr,Ttr
    end if
   end if
   close(101)
   smooth_vis=0
   call tracers_to_corners(rtr,str,Ttr,T)
   smooth_vis=1
  end if
 elseif (add_layer.eq.1) then
  if (tracer.eq.1) then
   open(unit=101,file="restart_tracers",form='unformatted')
   read(101) rtr,str,Ttr
   close(101)
  end if
 end if
end if
 do i=1,span+span_interp      !!!symmetry for BCs
  xg(-i)=-xg(i)
  xg(nx+i)=aspect-xg(nx-i)+aspect
 end do
 do k=1,span+span_interp
  zg(-k)=-zg(k)
  zg(nz+k)=1.d0-zg(nz-k)+1.d0
 end do
 xg(0)=0.d0
 zg(0)=0.d0
 xg(nx)=aspect
 zg(nz)=1.d0

do i=0,nx !!grid metrics
 if (gridX.ne.0) then
  xrrr(i)=dot_product(D3(-span:span),xg(i-span:i+span))/dr3
  xrrrr(i)=dot_product(D4(-span:span),xg(i-span:i+span))/dr4
 else
  xrrr(i)=0.d0
  xrrrr(i)=0.d0
 end if
end do
do i=-span_interp,nx+span_interp
 if (gridX.ne.0) then
  xr(i)=dot_product(D1(-span1:span1),xg(i-span1:i+span1))/dr
  xrr(i)=dot_product(D2(-span1:span1),xg(i-span1:i+span1))/dr2
 else
  xr(i)=1.d0
  xrr(i)=0.d0
 end if
end do
do k=0,nz !!grid metrics
 if (gridZ.ne.0) then
  zsss(k)=dot_product(D3(-span:span),zg(k-span:k+span))/ds3
  zssss(k)=dot_product(D4(-span:span),zg(k-span:k+span))/ds4
 else
  zsss(k)=0.d0
  zssss(k)=0.d0
 end if
end do
do k=-span_interp,nz+span_interp
 if (gridZ.ne.0) then
  zs(k)=dot_product(D1(-span1:span1),zg(k-span1:k+span1))/ds
  zss(k)=dot_product(D2(-span1:span1),zg(k-span1:k+span1))/ds2
 else
  zs(k)=1.d0
  zss(k)=0.d0
 end if
end do
if (comp.eq.1) call compute_metric_interpolator_derivatives


!!!!!!!!!!!!!!!!!movie files
if (movie.eq.1) then
 open(unit=90210,file="movie_files/movie_grid")
 do i=0,nx
  do k=0,nz
   write(90210,'(2(f6.2))') xg(i),zg(k)
  end do
 end do
 close(90210)
end if
!!!!!!!!!!!!!!!!!end movie files
end
