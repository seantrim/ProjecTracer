module tracers
 type trace_cell
  integer*4, allocatable :: id(:)
 end type trace_cell
end module tracers

subroutine mzran(ran_real,rank)
use basics
implicit none
integer*4 :: ran_int,rank
real*8 :: ran_real
ran_int=seed_i(rank)-seed_k(rank)
if (ran_int.lt.0) ran_int=ran_int+2147483579
seed_i(rank)=seed_j(rank)
seed_j(rank)=seed_k(rank)
seed_k(rank)=ran_int
seed_n(rank)=69069*seed_n(rank)+1013904243
ran_int=ran_int+seed_n(rank)
ran_real=0.5d0+(0.2328306d-9)*real(ran_int,8)
end

subroutine shuffle_serial(a,istart,N,rank) !!for small arrays (can be called by more than one thread simultaneously)
implicit none
integer*4 :: i,N,a(istart:N),rank,r,temp,istart
real*8 :: ran_real
do i=istart,N-1
 call mzran(ran_real,rank)
 r=floor(ran_real*real(N-i,8),4)+i+1   !!random integer from i+1 to N
 temp=a(i)
 a(i)=a(r)
 a(r)=temp
end do
end

subroutine shuffle(a,N) !!for a large array "a" of size "N" ------------------------------------ Maybe use random reals instead of integers to reduce repaeted random values????
!$ use OMP_LIB
implicit none
integer*4 :: i,N,rank,a(1:N),b(1:N),r(1:N),mloc_array(1:1),mloc,Nthreads,Nshare
integer*4, allocatable :: Nmin(:),Nmax(:),mloc_rank(:)
real*8 :: ran_real
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!create array of random integers (repeats possible)
!$OMP PARALLEL PRIVATE(rank,i)
rank=omp_get_thread_num()
!$OMP DO
do i=1,N
 call mzran(ran_real,rank)
! r(i)=nint(ran_real*real(N-1,8),4)+1   !!random integer from 1 to N !!!!!!!!!!!!!!!!!!!adapt to use floor instead of nint (so boundary values are equally likely compared to interior values)
 r(i)=floor(ran_real*real(N,8),4)+1 !!!!!!!!!!!!!use this one (untested)
end do
!$OMP END DO
!$OMP END PARALLEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end create array of random integers

!!!!!!!!!!!!!!!!!!!!!!!!!figure out how to divide work between processes
Nthreads=omp_get_max_threads()
Nshare=nint(real(N,8)/real(Nthreads,8),4)
allocate(Nmin(0:Nthreads-1),Nmax(0:Nthreads-1),mloc_rank(0:Nthreads-1))
!$OMP PARALLEL DO
do rank=0,Nthreads-1
 Nmin(rank)=rank*Nshare
 Nmax(rank)=(rank+1)*Nshare-1
end do
!$OMP END PARALLEL DO
Nmin(0)=1
Nmax(Nthreads-1)=N
!!!!!!!!!!!!!!!!!!!!!!!!!end figure out how to divide work between processes

do i=1,N
!$OMP PARALLEL PRIVATE(rank,mloc_array)
 rank=omp_get_thread_num()
 mloc_array=minloc(r(Nmin(rank):Nmax(rank)))
 mloc_rank(rank)=mloc_array(1)+Nmin(rank)-1
!$OMP END PARALLEL
 mloc_array=minloc(r(mloc_rank))
 rank=mloc_array(1)-1 !!rank with smallest minvalue
 mloc=mloc_rank(rank)
 b(i)=a(mloc)
 r(mloc)=N+1 !!once random number is used exclude from the search
end do
a=b
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!original TRA
subroutine equalize_tracers_original
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: ID,p,i,k,pmax,pp,tt,tpc_max,tpc_min
integer*4, allocatable :: IDsave(:)
real*8 :: perturb_r,perturb_s,r,s,tempT,tempC,ran_i,ran_k,time1,time2
real*8 :: ran_q,minexcess
real*8 :: Tdumb(-span_interp:nx+span_interp,-span_interp:nz+span_interp)
integer*4 :: nexcess(0:nx,0:nz),nfree,nadd,ncells,q,qmax,istore(1:ngrid),kstore(1:ngrid)
real*8 :: ran_real
integer*4 :: rank,Nshare,pp_min,pp_max
tpc_max=tpc !!!!!!!!!these used to be options in main.f90
tpc_min=tpc

!$ time1=omp_get_wtime()
!!!which cells have excess tracers?? we want nexcess=0 ideally
if (comp.eq.0) then
!$OMP PARALLEL DO PRIVATE(i,k)
 do i=0,nx
  do k=0,nz
   nexcess(i,k)=0
   if (nempty(i,k).gt.0) then
    nexcess(i,k)=nempty(i,k)-tpc  !!# of excess tracers
   end if
  end do
 end do
!$OMP END PARALLEL DO
elseif (comp.eq.1) then
 if (shape_function.eq.0) then
!$OMP PARALLEL DO PRIVATE(i,k,tt)
 do i=0,nx
  do k=0,nz
!!are the tracers all of one type? Don't modify interface or near interface cells or cells near empty cells
   nexcess(i,k)=0
   do tt=0,ntype
    if ((all(count_type(tt,i-1:i+1,k-1:k+1).eq.nempty(i-1:i+1,k-1:k+1))).and.&
&(nempty(i,k).gt.0)) then
     nexcess(i,k)=nempty(i,k)-tpc  !!# of excess tracers
     exit
    end if
   end do
  end do
 end do
!$OMP END PARALLEL DO
 elseif (shape_function.eq.1) then !!compositional interface is thicker for bilinear shape functions
!$OMP PARALLEL DO PRIVATE(i,k,tt)
 do i=0,nx
  do k=0,nz
!!are the tracers all of one type? Don't modify interface or near interface cells or cells near empty cells
   nexcess(i,k)=0
   do tt=0,ntype
    if ((all(count_type(tt,i-buffer:i+buffer,k-buffer:k+buffer).eq.nempty(i-buffer:i+buffer,k-buffer:k+buffer))).and.&
&(nempty(i,k).gt.0)) then
     nexcess(i,k)=nempty(i,k)-tpc  !!# of excess tracers
     exit
    end if
   end do
  end do
 end do
!$OMP END PARALLEL DO
 end if
end if
nfree=sum(nexcess,MASK=nexcess.gt.0) !!# of freed tracers to play with

!!is equalize worth doing? If not exit
if ((ALL((nexcess.lt.(tpc_max-tpc)).and.(nexcess.gt.(tpc_min-tpc)))).or.(nfree.eq.0)) then
 return
else
 nequalize=nequalize+1
end if

if (tracer.eq.1) then
 call extendT(Tratio) !!for interpolation of T to Ttr during equalizing
 call compute_derivatives(Textend(-span_interp:nx+span_interp,-span_interp:nz+span_interp),DERr_Textend)
end if

allocate(IDsave(1:nfree))
IDsave=0
p=1
do ID=1,ntr  !!look for tracers near the edges first
 i=nint(rtr(ID)/dr,4)
 k=nint(str(ID)/ds,4)
 if (nexcess(i,k).gt.0) then
  if ((str(ID).lt.dist_s).or.(str(ID).gt.(1.d0-dist_s)).or.&
&(rtr(ID).lt.dist_r).or.(rtr(ID).gt.(aspect-dist_r))) then !!too many tracers?? These tracers will be redistributed
   IDsave(p)=ID
   nexcess(i,k)=nexcess(i,k)-1
   p=p+1
  end if
 end if
end do
do ID=1,ntr  !!now look for other extra tracers
 i=nint(rtr(ID)/dr,4)
 k=nint(str(ID)/ds,4)
 if (nexcess(i,k).gt.0) then
  if ((str(ID).ge.dist_s).and.(str(ID).le.(1.d0-dist_s)).and.&
&(rtr(ID).ge.dist_r).and.(rtr(ID).le.(aspect-dist_r))) then !!too many tracers?? These tracers will be redistributed
   IDsave(p)=ID
   nexcess(i,k)=nexcess(i,k)-1
   p=p+1
   if (p.gt.nfree) exit
  end if
 end if
end do

!!!!!!add needed tracers to lacking cells via lottery
p=1
do !!p=1,nfree
minexcess=minval(nexcess)
if (minexcess.ge.0) then
 nfree=p-1
 exit !!don't modify interface cells
end if
ncells=count(nexcess.eq.minexcess) !!# of cells that have fewest tracers
q=1
do i=0,nx   !!find which cells lack tracers
 do k=0,nz
  if (nexcess(i,k).eq.minexcess) then
   istore(q)=i
   kstore(q)=k
   q=q+1
   if (q.gt.ncells) exit
  end if
 end do
 if (q.gt.ncells) exit
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!original
rank=0
do pp=1,ncells
! call random_number(ran_q)
 call mzran(ran_q,rank)
 q=floor(real(ncells,8)*ran_q)+1                  !!random # between 1 to ncells
 i=istore(q)
 k=kstore(q)
 r=dr*(real(i,8)) !!cell corners
 s=ds*(real(k,8))
! call random_number(perturb_r)
! call random_number(perturb_s)
 call mzran(perturb_r,rank)
 call mzran(perturb_s,rank)
 perturb_r=(perturb_r-0.5d0)*dr
 perturb_s=(perturb_s-0.5d0)*ds
 ID=IDsave(p)
 rtr(ID)=r+perturb_r
 str(ID)=s+perturb_s
 if (rtr(ID).le.0.d0) rtr(ID)=dist_r
 if (rtr(ID).ge.aspect) rtr(ID)=aspect-dist_r
 if (str(ID).le.0.d0) str(ID)=dist_s
 if (str(ID).ge.1.d0) str(ID)=1.d0-dist_s
 if (comp.eq.1) then
   Ttype(ID)=maxloc(count_type(0:ntype,i,k),DIM=1)-1 !!which tracer type is in that cell? Needed to subtract 1 since index starts at 0 in count_type array
 end if
 nexcess(i,k)=nexcess(i,k)+1
 p=p+1
 if (p.gt.nfree) exit
end do
if (p.gt.nfree) exit
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end original

if (tracer.eq.1) then
Tdumb=Textend(-span_interp:nx+span_interp,-span_interp:nz+span_interp)
!$OMP PARALLEL DO PRIVATE(ID,p)
do p=1,nfree !!!interpolate temperature values to repositioned tracers
 ID=IDsave(p)
  call interpolate(rtr(ID),str(ID),DERr_Textend,Tdumb,Ttr(ID))
end do
!$OMP END PARALLEL DO
end if

!$ time2=omp_get_wtime()
smooth_vis=0 !!viscosity smoothing not required yet
call tracers_to_corners(rtr,str,Ttr,T)
smooth_vis=1
!$ tequalize=tequalize+time2-time1
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end original TRA


subroutine equalize_tracers
!$ use OMP_LIB
 use basics
 use arrays
 use tracers
implicit none
integer*4 :: ID,p,i,k,pmax,pp,tt,ii,kk
real*8 :: perturb_r,perturb_s,r,s,tempT,tempC,ran_i,ran_k,time1,time2
real*8 :: ran_q
real*8 :: Tdumb(-span_interp:nx+span_interp,-span_interp:nz+span_interp)
integer*4 :: nexcess(0:nx,0:nz),nfree,nadd,ncells,q,qmax
real*8 :: ran_real
integer*4 :: rank,counter(0:nx,0:nz),type_save,escape,minexcess(0:ntype),tpc_type(0:ntype),count_save
type (trace_cell) tracer_in_cell(0:nx,0:nz)
integer*4 :: Nshare,nx_loop(0:nx),nz_loop(0:nz),count_display(0:ntype),Nideal(0:ntype),tracer_count(0:ntype)
integer*4, allocatable :: nx_min(:),nx_max(:),donor(:,:,:),index_save(:,:),index_save2(:,:),work(:),nfree_rank(:,:)
logical :: logic_test(0:nx,0:nz)
!$ time1=omp_get_wtime()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Testing
if (comp.eq.1) then
 do tt=1,ntype
  Nideal(tt)=nint(Cinit(tt)*real(ntr,8)) 
 end do
  Nideal(0)=ntr-sum(Nideal(1:ntype))
else
  Nideal(0)=ntr
end if
if (time.eq.0.d0) write(*,*) ntr,Nideal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!End testing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Which cells have excess tracers?? we want nexcess=0 ideally
if (comp.eq.0) then
 tpc_type(0)=tpc
!$OMP PARALLEL DO PRIVATE(i,k)
 do i=0,nx
  do k=0,nz
   nexcess(i,k)=0
   if (nempty(i,k).gt.0) then
    nexcess(i,k)=nempty(i,k)-tpc  !!# of excess tracers
   end if
  end do
 end do
!$OMP END PARALLEL DO
elseif (comp.eq.1) then
 logic_test=sum(Cnew(1:ntype,0:nx,0:nz),DIM=1).eq.0.d0
 count_save=count(logic_test)
 if (count_save.gt.0) then
  tpc_type(0)=nint(real(sum(count_type(0,0:nx,0:nz),MASK=logic_test),8)/&
  &real(count_save,8))
 else
  tpc_type(0)=0
 end if
 do tt=1,ntype
  logic_test=Cnew(tt,0:nx,0:nz).eq.1.d0
  count_save=count(logic_test)
  if (count_save.gt.0) then
   tpc_type(tt)=nint(real(sum(count_type(tt,0:nx,0:nz),MASK=Cnew(tt,0:nx,0:nz).eq.1.d0),8)/real(count_save,8))
  else
   tpc_type(tt)=0
  end if
 end do
!$OMP PARALLEL DO PRIVATE(i,k,tt)
 do i=0,nx
  do k=0,nz
!!are the tracers all of one type? Don't modify interface or near interface cells or cells near empty cells
   nexcess(i,k)=0
   do tt=0,ntype
    if ((all(count_type(tt,i-1:i+1,k-1:k+1).eq.nempty(i-1:i+1,k-1:k+1))).and.&
&(nempty(i,k).gt.0)) then
!     nexcess(i,k)=nempty(i,k)-tpc  !!# of excess tracers
     nexcess(i,k)=nempty(i,k)-tpc_type(tt)  !!# of excess tracers
     exit
    end if
   end do
  end do
 end do
!$OMP END PARALLEL DO
end if
nfree=sum(nexcess,MASK=nexcess.gt.0) !!# of freed tracers to play with
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!End Which cells have excess tracers?? we want nexcess=0 ideally

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Is equalize worth doing? If not exit
if (all(nexcess(0:nx,0:nz).eq.0)) then
 return
else
 nequalize=nequalize+1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!End Is equalize worth doing?

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Initial Prep
!Nthreads=omp_get_max_threads()
Nshare=nx/Nthreads
allocate(nx_min(0:Nthreads-1),nx_max(0:Nthreads-1),index_save(0:Nthreads-1,0:ntype),nfree_rank(0:Nthreads-1,0:ntype),&
&index_save2(0:Nthreads-1,0:ntype))
do rank=0,Nthreads-1
 nx_min(rank)=rank*Nshare
 nx_max(rank)=(rank+1)*Nshare-1
end do
nx_max(Nthreads-1)=nx

!$OMP PARALLEL DO PRIVATE(i,k)
do i=0,nx
 do k=0,nz
  allocate(tracer_in_cell(i,k)%id(0:nempty(i,k)))
  tracer_in_cell(i,k)%id(0)=0
 end do
end do
!$OMP END PARALLEL DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!End Initial Prep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Organize tracer IDs by CV location (i,k)
counter=0                                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Can Parallelization be more efficient??
!$OMP PARALLEL PRIVATE(rank,ID,i,k)
rank=omp_get_thread_num()
do ID=1,ntr
 i=nint(rtr(ID)/dr,4)
 if ((nx_min(rank).le.i).and.(i.le.nx_max(rank))) then
  k=nint(str(ID)/ds,4)
  counter(i,k)=counter(i,k)+1             !!!!!!!!!!!!!!!!!!!!!maybe reuse nempty(i,k) since final version of counter(i,k) is identical -- would save memory
  tracer_in_cell(i,k)%id(counter(i,k))=ID
 end if
end do
!$OMP END PARALLEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Organize tracer IDs by CV location (i,k)

!!!!!!!!!!!!!!!!!!!!!ADD SHUFFLE HERE TO TRACER_IN_CELL ARRAY
!$OMP PARALLEL PRIVATE(i,k,rank)
rank=omp_get_thread_num()
!$OMP DO
do i=0,nx
 do k=0,nz
  call shuffle_serial(tracer_in_cell(i,k)%id(1:counter(i,k)),1,counter(i,k),rank)
 end do
end do
!$OMP END DO
!$OMP END PARALLEL
!!!!!!!!!!!!!!!!!!!!!END SHUFFLE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Prep for T interpolation
if (tracer.eq.1) then
 call extendT(Tratio) !!for interpolation of T to Ttr during equalizing
 call compute_derivatives(Textend(-span_interp:nx+span_interp,-span_interp:nz+span_interp),DERr_Textend)
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!End Prep for T interpolation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Look for donor tracers: NEWEST
allocate(donor(0:Nthreads-1,0:nfree,0:ntype))   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Save on memory here???? 0:nfree may be too big a range
donor=0  !!default values
index_save=0 !!staring value
!!!!!!!!!!!!!!!!Look near edges first
do i=0,nx,nx  !!!left and right sides
!$OMP PARALLEL PRIVATE(rank,k,q,ID)
rank=omp_get_thread_num()
!$OMP DO
do k=0,nz
 if (nexcess(i,k).gt.0) then !!only take donors from rich unrestricted CVs
  do q=1,counter(i,k)
    ID=tracer_in_cell(i,k)%id(q)
    if ((str(ID).lt.dist_s).or.(str(ID).gt.(1.d0-dist_s)).or.&
&(rtr(ID).lt.dist_r).or.(rtr(ID).gt.(aspect-dist_r))) then !!too many tracers?? These tracers will be redistributed
     index_save(rank,Ttype(ID))=index_save(rank,Ttype(ID))+1
     donor(rank,index_save(rank,Ttype(ID)),Ttype(ID))=ID
     nexcess(i,k)=nexcess(i,k)-1
     if (nexcess(i,k).eq.0) exit
    end if
  end do
 end if
end do
!$OMP END DO
!$OMP END PARALLEL
end do

do k=0,nz,nz  !!!bottom and top boundaries
!$OMP PARALLEL PRIVATE(rank,i,q,ID)
rank=omp_get_thread_num()
!$OMP DO
do i=1,nx-1   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! corners already done
 if (nexcess(i,k).gt.0) then !!only take donors from rich unrestricted CVs
  do q=1,counter(i,k)
    ID=tracer_in_cell(i,k)%id(q)
    if ((str(ID).lt.dist_s).or.(str(ID).gt.(1.d0-dist_s)).or.&
&(rtr(ID).lt.dist_r).or.(rtr(ID).gt.(aspect-dist_r))) then !!too many tracers?? These tracers will be redistributed
     index_save(rank,Ttype(ID))=index_save(rank,Ttype(ID))+1
     donor(rank,index_save(rank,Ttype(ID)),Ttype(ID))=ID
     nexcess(i,k)=nexcess(i,k)-1
     if (nexcess(i,k).eq.0) exit
    end if
  end do
 end if
end do
!$OMP END DO
!$OMP END PARALLEL
end do
!!!!!!!!!!!!!!!!End Look near edges first

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!shuffle donor tracers near the boundaries
!$OMP PARALLEL PRIVATE(rank,tt)
rank=omp_get_thread_num()
do tt=0,ntype
 call shuffle_serial(donor(rank,1:index_save(rank,tt),tt),1,index_save(rank,tt),rank)
end do
!$OMP END PARALLEL
index_save2=index_save !!index_save should be saved here for later shuffle of interior donor tracers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end shuffle tracers near the boundaries

!!!!!!!!!!!!!!!!!!!!!!!!!!wrap up edge CVs (now look for excess tracers in edge CV interiors, if any)
do i=0,nx,nx !!left and right side interiors
!$OMP PARALLEL PRIVATE(rank,k,q,ID)
rank=omp_get_thread_num()
!$OMP DO
 do k=0,nz
  if (nexcess(i,k).gt.0) then !!only take donors from rich unrestricted CVs
   do q=1,counter(i,k)
     ID=tracer_in_cell(i,k)%id(q)
     if (((str(ID).ge.dist_s).and.(str(ID).le.(1.d0-dist_s))).and.&
&((rtr(ID).ge.dist_r).and.(rtr(ID).le.(aspect-dist_r)))) then !!too many tracers?? These tracers will be redistributed
      index_save(rank,Ttype(ID))=index_save(rank,Ttype(ID))+1
      donor(rank,index_save(rank,Ttype(ID)),Ttype(ID))=ID
      nexcess(i,k)=nexcess(i,k)-1
      if (nexcess(i,k).eq.0) exit
     end if
   end do
  end if
 end do
!$OMP END DO
!$OMP END PARALLEL
end do

do k=0,nz,nz !!bottom and top boundary CV interiors
!$OMP PARALLEL PRIVATE(rank,i,q,ID)
rank=omp_get_thread_num()
!$OMP DO
 do i=1,nx-1 !!already did corners
  if (nexcess(i,k).gt.0) then !!only take donors from rich unrestricted CVs
   do q=1,counter(i,k)
     ID=tracer_in_cell(i,k)%id(q)
     if (((str(ID).ge.dist_s).and.(str(ID).le.(1.d0-dist_s))).and.&
&(rtr(ID).ge.dist_r).and.(rtr(ID).le.(aspect-dist_r))) then !!too many tracers?? These tracers will be redistributed:-------------------------- Cut down if statement (no corners)
      index_save(rank,Ttype(ID))=index_save(rank,Ttype(ID))+1
      donor(rank,index_save(rank,Ttype(ID)),Ttype(ID))=ID
      nexcess(i,k)=nexcess(i,k)-1
      if (nexcess(i,k).eq.0) exit
     end if
   end do
  end if
 end do
!$OMP END DO
!$OMP END PARALLEL
end do
!!!!!!!!!!!!!!!!!!!!!!!!!End wrap up edge cells

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Search remaining (interior) CVs for donor tracers
!$OMP PARALLEL PRIVATE(rank,i,k,q,ID)
rank=omp_get_thread_num()
!$OMP DO
do i=1,nx-1 !!left and right side interiors
 do k=1,nz-1
  if (nexcess(i,k).gt.0) then !!only take donors from rich unrestricted CVs
   do q=1,counter(i,k)
     ID=tracer_in_cell(i,k)%id(q)
     index_save(rank,Ttype(ID))=index_save(rank,Ttype(ID))+1
     donor(rank,index_save(rank,Ttype(ID)),Ttype(ID))=ID
     nexcess(i,k)=nexcess(i,k)-1
     if (nexcess(i,k).eq.0) exit
   end do
  end if
 end do
end do
!$OMP END DO
!$OMP END PARALLEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!End Search remaining (interior) CVs for donor tracers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Test to see if we got all the donor tracers
nfree_rank=0
!$OMP PARALLEL PRIVATE(rank,tt)
rank=omp_get_thread_num()
do tt=0,ntype
 nfree_rank(rank,tt)=count(donor(rank,:,tt).gt.0)
end do
!$OMP END PARALLEL
if (sum(nfree_rank).ne.nfree) then
 write(*,*) "Error: # of donor tracers found inconsistent with nfree value"
 stop
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!End test to see if we got all the donor tracers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!End Look for donor tracers: NEWEST

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!***********************Shuffle donor tracers that are not near the boundaries
!$OMP PARALLEL PRIVATE(rank,tt)
rank=omp_get_thread_num()
do tt=0,ntype
 call shuffle_serial(donor(rank,index_save2(rank,tt)+1:index_save(rank,tt),tt),index_save2(rank,tt)+1,index_save(rank,tt),rank)
end do
!$OMP END PARALLEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!***********************End Shuffle donor tracers that are not near the boundaries


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Add needed tracers to lacking cells via lottery: NEW
!$OMP PARALLEL DO !!!!!prep for randomized loops
do i=0,nx
 nx_loop(i)=i
end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
do k=0,nz
 nz_loop(k)=k
end do
!$OMP END PARALLEL DO

index_save=0 !!to keep track of which tracer IDs need to be used for T interpolation
escape=0
do tt=0,ntype
 tracer_count(tt)=sum(count_type(tt,0:nx,0:nz))
end do
if (all(tracer_count.eq.Nideal)) escape=1

do
if (escape.eq.0) then
 do tt=0,ntype
  tracer_count(tt)=sum(count_type(tt,0:nx,0:nz))
 end do
 if (all(tracer_count.eq.Nideal)) escape=1
end if
if (tpc_type(0).gt.0) then
 if (comp.eq.0) then
  minexcess(0)=minval(nexcess(0:nx,0:nz),MASK=(nexcess(0:nx,0:nz).le.0)) !!ambient composition
 else
  minexcess(0)=minval(nexcess(0:nx,0:nz),MASK=((sum(Cnew(1:ntype,0:nx,0:nz),DIM=1).eq.0.d0).and.(nexcess(0:nx,0:nz).le.0))) !!ambient composition
 end if
else
 minexcess(0)=0
end if
do p=1,ntype
 if (tpc_type(p).gt.0) then
  minexcess(p)=minval(nexcess(0:nx,0:nz),MASK=((Cnew(p,0:nx,0:nz).eq.1.d0).and.(nexcess(0:nx,0:nz).le.0))) !!enriched compositions
 else
  minexcess(p)=0
 end if
end do
if (all(minexcess.eq.0)) exit
if (sum(index_save(:,:)).ge.nfree) exit
index_save2=index_save
rank=0
 call shuffle_serial(nx_loop,0,nx,rank) !!randomize loop order
 call shuffle_serial(nz_loop,0,nz,rank) !!randomize loop order
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Find poorest CV locations
!$OMP PARALLEL PRIVATE(i,k,ii,kk,rank,r,s,perturb_r,perturb_s,ID,tt,p,type_save)
rank=omp_get_thread_num()
!$OMP DO
do ii=0,nx   !!find which cells lack tracers
 i=nx_loop(ii) !!!use randomized i values
 do kk=0,nz
  k=nz_loop(kk) !!!use randomized k values
  type_save=maxloc(count_type(0:ntype,i,k),DIM=1)-1 !!which tracer type is in that cell? Needed to subtract 1 since index starts at 0 in count_type array
  if (escape.eq.0) then
   do p=0,ntype !!!!!!!!!!!!!!!!decide which type of tracer to donate in order to get optimum balance
    if ((tracer_count(p)).gt.Nideal(p)) then
     tt=p
     exit
    end if
   end do
  else
   tt=type_save !!if ideal values are met then do not swap tracer types
  end if
  if ((nexcess(i,k).eq.minexcess(type_save)).and.(index_save(rank,tt).lt.nfree_rank(rank,tt)).and.(nexcess(i,k).lt.0)) then
   r=dr*(real(i,8)) !!cell corners
   s=ds*(real(k,8))
   call mzran(perturb_r,rank)
   call mzran(perturb_s,rank)
   perturb_r=(perturb_r-0.5d0)*dr
   perturb_s=(perturb_s-0.5d0)*ds
   index_save(rank,tt)=index_save(rank,tt)+1
   ID=donor(rank,index_save(rank,tt),tt)
   rtr(ID)=r+perturb_r
   str(ID)=s+perturb_s
   if (rtr(ID).le.0.d0) rtr(ID)=dist_r
   if (rtr(ID).ge.aspect) rtr(ID)=aspect-dist_r
   if (str(ID).le.0.d0) str(ID)=dist_s
   if (str(ID).ge.1.d0) str(ID)=1.d0-dist_s
   if (comp.eq.1) then
    Ttype(ID)=type_save !!which tracer type is in that cell? Needed to subtract 1 since index starts at 0 in count_type array
   end if
   count_type(tt,i,k)=count_type(tt,i,k)-1
   count_type(Ttype(ID),i,k)=count_type(Ttype(ID),i,k)+1
   nexcess(i,k)=nexcess(i,k)+1
  end if
 end do
end do
!$OMP END DO
!$OMP END PARALLEL
 if (all(index_save.eq.index_save2)) exit !!if no repositioning happened on the last iteration
end do

!!!!!!!!!!!!!!!!!!Use any remaining donor tracers in serial (alternating between ranks to get even usage of donors)
rank=0
do
if (tpc_type(0).gt.0) then
 if (comp.eq.0) then
  minexcess(0)=minval(nexcess(0:nx,0:nz),MASK=(nexcess(0:nx,0:nz).le.0)) !!ambient composition
 else
  minexcess(0)=minval(nexcess(0:nx,0:nz),MASK=((sum(Cnew(1:ntype,0:nx,0:nz),DIM=1).eq.0.d0).and.(nexcess(0:nx,0:nz).le.0))) !!ambient composition
 end if
else
 minexcess(0)=0
end if
do p=1,ntype
 if (tpc_type(0).gt.0) then
  minexcess(p)=minval(nexcess(0:nx,0:nz),MASK=((Cnew(p,0:nx,0:nz).eq.1.d0).and.(nexcess(0:nx,0:nz).le.0))) !!enriched compositions
 else
  minexcess(p)=0
 end if
end do
if (all(minexcess.eq.0)) exit
if (sum(index_save(:,:)).ge.nfree) exit
index_save2=index_save
rank=0
 call shuffle_serial(nx_loop,0,nx,rank) !!randomize loop order
 call shuffle_serial(nz_loop,0,nz,rank) !!randomize loop order
do ii=0,nx   !!find which cells lack tracers
 i=nx_loop(ii) !!!use randomized i values
 do kk=0,nz
  k=nz_loop(kk) !!!use randomized k values
  tt=maxloc(count_type(0:ntype,i,k),DIM=1)-1 !!which tracer type is in that cell? Needed to subtract 1 since index starts at 0 in count_type array
  if ((nexcess(i,k).eq.minexcess(tt)).and.(sum(index_save(:,tt)).lt.sum(nfree_rank(:,tt))).and.(nexcess(i,k).lt.0)) then
   do
    rank=rank+1
    if (rank.eq.Nthreads) rank=0
    if (index_save(rank,tt).lt.nfree_rank(rank,tt)) exit
   end do
   r=dr*(real(i,8)) !!cell corners
   s=ds*(real(k,8))
   call mzran(perturb_r,rank)
   call mzran(perturb_s,rank)
   perturb_r=(perturb_r-0.5d0)*dr
   perturb_s=(perturb_s-0.5d0)*ds
   index_save(rank,tt)=index_save(rank,tt)+1
   ID=donor(rank,index_save(rank,tt),tt)
   rtr(ID)=r+perturb_r
   str(ID)=s+perturb_s
   if (rtr(ID).le.0.d0) rtr(ID)=dist_r
   if (rtr(ID).ge.aspect) rtr(ID)=aspect-dist_r
   if (str(ID).le.0.d0) str(ID)=dist_s
   if (str(ID).ge.1.d0) str(ID)=1.d0-dist_s
   if (comp.eq.1) then
    Ttype(ID)=maxloc(count_type(0:ntype,i,k),DIM=1)-1 !!which tracer type is in that cell? Needed to subtract 1 since index starts at 0 in count_type array
   end if
   nexcess(i,k)=nexcess(i,k)+1
  end if
 end do
end do
 if (all(index_save.eq.index_save2)) exit !!if no repositioning happened on the last iteration
end do
!!!!!!!!!!!!!!!!!!End use any remaining donor tracers in serial 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!End Find poorest CV locations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Interpolate T values of repositioned tracers
if (tracer.eq.1) then
Tdumb=Textend(-span_interp:nx+span_interp,-span_interp:nz+span_interp)
!$OMP PARALLEL PRIVATE(ID,p,rank,tt)
rank=omp_get_thread_num()
do tt=0,ntype
 do p=1,index_save(rank,tt) !!!interpolate temperature values to repositioned tracers
  ID=donor(rank,p,tt)
  call interpolate(rtr(ID),str(ID),DERr_Textend,Tdumb,Ttr(ID))
 end do
end do
!$OMP END PARALLEL
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!End Interpolate T values of repositioned tracers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!End Add needed tracers to lacking cells via lottery: NEW

if (comp.eq.0) then
 minexcess(0)=minval(nexcess(0:nx,0:nz),MASK=(nexcess(0:nx,0:nz).le.0)) !!ambient composition
else
 minexcess(0)=minval(nexcess(0:nx,0:nz),MASK=((sum(Cnew(1:ntype,0:nx,0:nz),DIM=1).eq.0.d0).and.(nexcess(0:nx,0:nz).le.0))) !!ambient composition
 do p=1,ntype
  minexcess(p)=minval(nexcess(0:nx,0:nz),MASK=((Cnew(p,0:nx,0:nz).eq.1.d0).and.(nexcess(0:nx,0:nz).le.0))) !!enriched compositions
 end do
end if

do tt=0,ntype
 count_display(tt)=sum(count_type(tt,0:nx,0:nz))
end do
if (comp.eq.0) then
 write(1024,'(g20.8,15(i10))') time,tpc,minexcess,sum(nexcess(0:nx,0:nz),MASK=nexcess(0:nx,0:nz).le.0),&
 &sum(index_save),nfree,count_display
else
 write(1024,'(g20.8,15(i10))') time,tpc_type,minexcess,sum(nexcess(0:nx,0:nz),MASK=nexcess(0:nx,0:nz).le.0),&
 &sum(index_save),nfree,count_display
end if
!$ time2=omp_get_wtime()
smooth_vis=0 !!viscosity smoothing not required yet
call tracers_to_corners(rtr,str,Ttr,T)
smooth_vis=1
!$ tequalize=tequalize+time2-time1
end

subroutine entrainment_setup
 use basics
implicit none
integer*4 :: k
real*8 :: scell,zgrid
  do k=0,nz
   scell=ds*(real(k,8))
    if (zgrid(scell).le.dent) kent=k
  end do
end

subroutine initialize_tracers
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k,p,ID,composition,rank,nproc
integer*4, allocatable :: nx_start(:),nx_finish(:),ID_start(:),ID_finish(:)
real*8 :: rcell,scell,perturb_r,perturb_s,xtr,ztr
real*8 :: xgrid,zgrid,temperature,composition_continuous
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!define ranges of parallel loops
nproc=omp_get_max_threads()
allocate(nx_start(0:nproc-1),nx_finish(0:nproc-1),ID_start(0:nproc-1),ID_finish(0:nproc-1))

nx_start(0)=0
do rank=0,nproc-2
 nx_finish(rank)=nx_start(rank)+((nx+1)/nproc-1)
 nx_start(rank+1)=nx_finish(rank)+1
end do
nx_finish(nproc-1)=nx

ID_start(0)=1
do rank=0,nproc-2
 ID_finish(rank)=((nx_finish(rank)-nx_start(rank)+1)*(nz+1))*tpc+ID_start(rank)-1
 ID_start(rank+1)=ID_finish(rank)+1
end do
ID_finish(nproc-1)=((nx_finish(nproc-1)-nx_start(nproc-1)+1)*(nz+1))*tpc+ID_start(nproc-1)-1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end define ranges of parallel loops
 rank=0
!$OMP PARALLEL PRIVATE(i,k,rank,rcell,scell,perturb_r,perturb_s,xtr,ztr,ID)
rank=omp_get_thread_num()
ID=ID_start(rank)
 do i=nx_start(rank),nx_finish(rank) !!go through each cell corner
 rcell=dr*(real(i,8))
  do k=0,nz
   scell=ds*(real(k,8))
   do p=1,tpc
    call mzran(perturb_r,rank)
    call mzran(perturb_s,rank)
    perturb_r=(perturb_r-0.5d0)*dr
    perturb_s=(perturb_s-0.5d0)*ds
    rtr(ID)=rcell+perturb_r   
    str(ID)=scell+perturb_s
    if (rtr(ID).le.0.d0) rtr(ID)=0.d0 !dist_r
    if (rtr(ID).ge.aspect) rtr(ID)=aspect !aspect-dist_r
    if (str(ID).le.0.d0) str(ID)=0.d0 !dist_s
    if (str(ID).ge.1.d0) str(ID)=1.d0 !1.d0-dist_s
    xtr=xgrid(rtr(ID))       !!!!!maybe interpolate these?? minor detail... 
    ztr=zgrid(str(ID))
    if (comp.eq.1) then
     if (iC_continuous.eq.0) then
      Ttype(ID)=composition(xtr,ztr)
     elseif (iC_continuous.eq.1) then
      Ctr(ID)=composition_continuous(xtr,ztr)
     end if
    else
     Ttype(ID)=0  !!!uniform composition
    end if
    if (tracer.eq.1) then
     Ttr(ID)=temperature(xtr,ztr)
    end if
    ID=ID+1
   end do
  end do
 end do
!$OMP END PARALLEL
 call tracers_to_corners(rtr,str,Ttr,T)!!!!!!!!!!convert tracers to initial C field here

!!!!!!!!!!entrainment setup
 if (comp.eq.1) then
  do k=0,nz
   scell=ds*(real(k,8))
   if (zgrid(scell).le.dent) kent=k
  end do
 end if
!!!!!!!!!!end entrainment setup
end

subroutine add_tracers
 use basics
 use arrays
implicit none
integer*4 :: i,k,id,nseed,tt,composition
integer*4, allocatable :: seed(:)
real*8 :: scell,xtr,ztr
real*8 :: xgrid,zgrid,temperature
do id=1,ntr
 xtr=xgrid(rtr(id))       !!!!!maybe interpolate these?? minor detail... 
 ztr=zgrid(str(id))
 Ttype(id)=composition(xtr,ztr)
end do
 if (comp.eq.1) then
  Cnew=0.d0
  do k=0,nz
   scell=ds*(real(k,8))
   if (zgrid(scell).le.dent) kent=k
   do i=0,nx
    do tt=1,ntype
     if (composition(xg(i),zg(k)).eq.tt) Cnew(tt,i,k)=1.d0
    end do
   end do
  end do
 end if
 call tracers_to_corners(rtr,str,Ttr,T)!!!!!!!!!!convert tracers to initial C field here
end

real*8 function shape_fn(i,k,r,s) !!OG shape function (just bilinear, no if statements)
 use basics
implicit none
real*8 :: r,s
integer*4 :: i,k
shape_fn=max(0.d0,1.d0-dabs(r-real(i,8)*dr)/dr)*max(0.d0,1.d0-dabs(s-real(k,8)*ds)/ds)
end

real*8 function shape_func(i,k,r,s) !!updated shape function (choice of constant cell or bilinear)
!!Note that shape functions are defined in computation space (r,s)
!!Also note that the constant cell option has a weighting of unity on all sides of a dual grid cell.
!!Thus, the constant cell definition is slightly more inclusive than in the tracers_to_corners_OG routine.
 use basics
implicit none

!!inputs
real*8 :: r,s    !!Lagrangian tracer position 
integer*4 :: i,k !!grid indices for interpolated value

!!internal variables
real*8 :: ri,sk !!grid position of interpolated value in computation space

ri=real(i,8)*dr; sk=real(k,8)*ds !!grid position of interpolated value (r,s) space
if (shape_function.eq.0) then
 if (((ri-dro2).le.r).and.(r.le.(ri+dro2)).and.((sk-dso2).le.s).and.(s.le.(sk+dso2))) then
  shape_func=1.d0 !!if tracer is within the dual grid cell surrounding corner point (i,k)
 else
  shape_func=0.d0
 end if
elseif (shape_function.eq.1) then
 shape_func=max(0.d0,1.d0-abs(r-real(i,8)*dr)/dr)*max(0.d0,1.d0-abs(s-real(k,8)*ds)/ds)
 shape_func=max(0.d0,1.d0-abs(r-ri)/dr)*max(0.d0,1.d0-abs(s-sk)/ds)
end if
end

subroutine tracers_to_corners(rtr_,str_,Ttr_,T0)
!!wrapper function for interpolating tracer T and C values to cell corners (i.e., Eulerian T and C)
!!chooses between method for discrete tracer types (Lagrangian C values are either 0 or 1) ...
!!      OR the method that allows Lagrangian C values to be fractional
 use basics
implicit none
real*8 :: rtr_(1:ntr),str_(1:ntr),Ttr_(1:ntr)
real*8 :: T0(-span:nx+span,-span:nz+span)
if (iC_continuous.eq.1) then
 call tracers_to_corners_iC_continuous(rtr_,str_,Ttr_,T0)
elseif (iC_continuous.eq.0) then
 call tracers_to_corners_OG(rtr_,str_,Ttr_,T0)
else
 write(*,*) "Error in tracers_to_corners: invalid iC_continuous value."
 write(*,*) "iC_continuous=",iC_continuous
 stop
end if
end subroutine tracers_to_corners

subroutine tracers_to_corners_iC_continuous(rtr_,str_,Ttr_,T0) !!T0 input values are for empty cells
!!both T and C tracer values vary continuously
!!If comp=1, this routine assumes ntype=1 (for now).
!!Stores Eulerian composition and temperature in Cnew and Tratio arrays (to be used by other routines)
!$ use OMP_LIB
 use basics
 use arrays
implicit none

!!inputs (these arrays are not modified by this routine)
real*8 :: rtr_(1:ntr),str_(1:ntr),Ttr_(1:ntr) !!tracer positions and Lagrangian T values
real*8 :: T0(-span:nx+span,-span:nz+span)     !!initial T values on Eulerian grid

!!external functions
real*8 :: shape_func !!shape function

!!internal variables
integer*4 :: tt !!index for compositional type (should be 1 for now)
integer*4 :: i,k !!indices for Eulerian grid
integer*4 :: rank,Nshare !!parallel computation variables
integer*4 :: id,id_min,id_max !!tracer ID variables
integer*4, allocatable :: nempty_work(:,:,:) !!parallel work array for nempty array
real*8 :: r,s !!position in computation space
real*8 :: time1,time2 !!timing variables
real*8 :: shape_sum(-span:nx+span,-span:nz+span) !!contains sum of shape function weightings used in TxS
real*8 :: shape_sum_work(0:Nthreads-1,-span:nx+span,-span:nz+span) !!parallel work array for portions of shape_sum
real*8, allocatable :: TxS(:,:),CxS(:,:) !!contains products of Lagrangian values and shape function weightings
real*8, allocatable :: TxS_work(:,:,:),CxS_work(:,:,:) !!parallel work array for TxS

!!!!******Main Body Starts********
!$ time1=omp_get_wtime()
if (tracer.eq.1) Tratio=T0 !!for empty cells -- initial values will be retained
!intial values for C are stored in Cnew array (done before calling this routine)

shape_sum_work=0.d0 !!initialize work aray
nempty=0 !!initialize tracer count array
allocate(nempty_work(0:Nthreads-1,-span:nx+span,-span:nz+span)); nempty_work=0
if (tracer.eq.1) then
 allocate(TxS(-span:nx+span,-span:nz+span)); TxS=0.d0
 allocate(TxS_work(0:Nthreads-1,-span:nx+span,-span:nz+span)); TxS_work=0.d0
end if
if (comp.eq.1) then
 allocate(CxS(-span:nx+span,-span:nz+span)); CxS=0.d0
 allocate(CxS_work(0:Nthreads-1,-span:nx+span,-span:nz+span)); CxS_work=0.d0
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!********************compute TxS and shape_sum in parallel******************
Nshare=ntr/Nthreads !!# of tracers per OpenMP thread
rank=0 !!this ensures that the serial case works as well

!$OMP PARALLEL PRIVATE(id,id_min,id_max,r,s,i,k,rank)
!$ rank=omp_get_thread_num() !!rank is assigned for each OpenMP thread

id_min=(1+rank*Nshare)              !!!determine tracer ID ranges for each rank
if (rank.eq.(Nthreads-1)) then      !!rank ranges from 0 to Nthreads-1
 id_max=max((rank+1)*Nshare,ntr)    !!ensures ntr is not exceeded in the last rank
else
 id_max=(rank+1)*Nshare
end if

do id=id_min,id_max
 r=rtr_(id)
 s=str_(id)
 i=nint(r/dr,4)
 k=nint(s/ds,4)
 if ((r.lt.0.d0).or.(r.gt.aspect).or.(s.lt.0.d0).or.(s.gt.1.d0)) then !!validation for tracer position
  write(*,*) "Tracer advected outside of domain - stopping"
  stop
 end if
 nempty_work(rank,i,k)=nempty_work(rank,i,k)+1 !!count tracers in dual grid cell surrounding (i,k)
 if (tracer.eq.1) then !!the amount of neighbours needed may depend on the shape function
  TxS_work(rank,i-1,k-1)=TxS_work(rank,i-1,k-1)+Ttr_(id)*shape_func(i-1,k-1,r,s)
  TxS_work(rank,i-1,k)=TxS_work(rank,i-1,k)+Ttr_(id)*shape_func(i-1,k,r,s)
  TxS_work(rank,i-1,k+1)=TxS_work(rank,i-1,k+1)+Ttr_(id)*shape_func(i-1,k+1,r,s)
  TxS_work(rank,i,k-1)=TxS_work(rank,i,k-1)+Ttr_(id)*shape_func(i,k-1,r,s)
  TxS_work(rank,i,k)=TxS_work(rank,i,k)+Ttr_(id)*shape_func(i,k,r,s)
  TxS_work(rank,i,k+1)=TxS_work(rank,i,k+1)+Ttr_(id)*shape_func(i,k+1,r,s)
  TxS_work(rank,i+1,k-1)=TxS_work(rank,i+1,k-1)+Ttr_(id)*shape_func(i+1,k-1,r,s)
  TxS_work(rank,i+1,k)=TxS_work(rank,i+1,k)+Ttr_(id)*shape_func(i+1,k,r,s)
  TxS_work(rank,i+1,k+1)=TxS_work(rank,i+1,k+1)+Ttr_(id)*shape_func(i+1,k+1,r,s)
 end if
 if (comp.eq.1) then
  CxS_work(rank,i-1,k-1)=CxS_work(rank,i-1,k-1)+Ctr(id)*shape_func(i-1,k-1,r,s)
  CxS_work(rank,i-1,k)=CxS_work(rank,i-1,k)+Ctr(id)*shape_func(i-1,k,r,s)
  CxS_work(rank,i-1,k+1)=CxS_work(rank,i-1,k+1)+Ctr(id)*shape_func(i-1,k+1,r,s)
  CxS_work(rank,i,k-1)=CxS_work(rank,i,k-1)+Ctr(id)*shape_func(i,k-1,r,s)
  CxS_work(rank,i,k)=CxS_work(rank,i,k)+Ctr(id)*shape_func(i,k,r,s)
  CxS_work(rank,i,k+1)=CxS_work(rank,i,k+1)+Ctr(id)*shape_func(i,k+1,r,s)
  CxS_work(rank,i+1,k-1)=CxS_work(rank,i+1,k-1)+Ctr(id)*shape_func(i+1,k-1,r,s)
  CxS_work(rank,i+1,k)=CxS_work(rank,i+1,k)+Ctr(id)*shape_func(i+1,k,r,s)
  CxS_work(rank,i+1,k+1)=CxS_work(rank,i+1,k+1)+Ctr(id)*shape_func(i+1,k+1,r,s)
 end if

 shape_sum_work(rank,i-1,k-1)=shape_sum_work(rank,i-1,k-1)+shape_func(i-1,k-1,r,s)
 shape_sum_work(rank,i-1,k)=shape_sum_work(rank,i-1,k)+shape_func(i-1,k,r,s)
 shape_sum_work(rank,i-1,k+1)=shape_sum_work(rank,i-1,k+1)+shape_func(i-1,k+1,r,s)
 shape_sum_work(rank,i,k-1)=shape_sum_work(rank,i,k-1)+shape_func(i,k-1,r,s)
 shape_sum_work(rank,i,k)=shape_sum_work(rank,i,k)+shape_func(i,k,r,s)
 shape_sum_work(rank,i,k+1)=shape_sum_work(rank,i,k+1)+shape_func(i,k+1,r,s)
 shape_sum_work(rank,i+1,k-1)=shape_sum_work(rank,i+1,k-1)+shape_func(i+1,k-1,r,s)
 shape_sum_work(rank,i+1,k)=shape_sum_work(rank,i+1,k)+shape_func(i+1,k,r,s)
 shape_sum_work(rank,i+1,k+1)=shape_sum_work(rank,i+1,k+1)+shape_func(i+1,k+1,r,s)
end do
!$OMP END PARALLEL

!!!!Sum up parallel work results
!$OMP PARALLEL DO PRIVATE(i,k)
do k=0,nz
 do i=0,nx
  nempty(i,k)=sum(nempty_work(:,i,k))
  shape_sum(i,k)=sum(shape_sum_work(:,i,k))
  if (tracer.eq.1) TxS(i,k)=sum(TxS_work(:,i,k))
  if (comp.eq.1) CxS(i,k)=sum(CxS_work(:,i,k))
 end do
end do
!$OMP END PARALLEL DO
!!!!End Sum up parallel work results
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!********************end compute TxS and shape_sum in parallel******************

!$OMP PARALLEL DO PRIVATE(i,k)  !!apply the weighted average to determine the Eulerian interpolated value
do k=0,nz
 do i=0,nx
  if (shape_sum(i,k).gt.0.d0) then !!if there are no tracers close to the cell corner retain previous Eulerian value
   if (tracer.eq.1) Tratio(i,k)=TxS(i,k)/shape_sum(i,k)
   if (comp.eq.1) Cnew(1,i,k)=CxS(i,k)/shape_sum(i,k) !!ntype=1 assumed for now
  end if
 end do
end do
!$OMP END PARALLEL DO
!!!!******Main Body Ends********

!!post processing
if (comp.eq.1) then
 do tt=1,ntype
  call enforceBCs_comp(Cnew(tt,:,:))
 end do
 if (smooth_vis.eq.1) call smoother_time
end if
if (tracer.eq.1) then
 call enforceBCs(Tratio)
 if (smooth_vis.eq.1) call smoother_time_T(Tratio)
end if
call load_conductivity_array
!$ time2=omp_get_wtime()
!$ tconvert=tconvert+time2-time1
end

!!!!******************************** OG tracers_to_corners prior to intermediate C values for tracers (i.e., no Ctr(:) array) 
subroutine tracers_to_corners_OG(rtr_,str_,Ttr_,T0) !!T0 input values are for empty cells
!$ use OMP_LIB
 use basics
 use arrays
implicit none
real*8 :: r,s,T0(-span:nx+span,-span:nz+span),T1(-span:nx+span,-span:nz+span)
real*8 :: temp,rtr_(1:ntr),str_(1:ntr),time1,time2,shape_fn
real*8 :: Ttr_(1:ntr)
integer*4 :: id,i,k,tt
integer*4 :: rank,Nshare,id_min,id_max
real*8, allocatable :: Twork(:,:,:)
integer*4, allocatable :: count_work(:,:,:,:)
real*8, allocatable :: count_work_real(:,:,:,:)
!$ time1=omp_get_wtime()
if (tracer.eq.1) Tratio=T0 !!for empty cells
 T1=0.d0
 count_type=0
 if (shape_function.eq.1) count_type_real=0.d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!parallel
Nshare=ntr/Nthreads
if (tracer.eq.1) then
 allocate(Twork(0:Nthreads-1,-span:nx+span,-span:nz+span))
 Twork=0.d0
end if
allocate(count_work(0:Nthreads-1,0:ntype,-span:nx+span,-span:nz+span))
count_work=0
if (shape_function.eq.1) then
 allocate(count_work_real(0:Nthreads-1,0:ntype,-span:nx+span,-span:nz+span))
 count_work_real=0.d0
end if
!$OMP PARALLEL PRIVATE(id,id_min,id_max,r,s,i,k,rank)
rank=0
!$ rank=omp_get_thread_num()
id_min=(1+rank*Nshare)
if (rank.eq.(Nthreads-1)) then
 id_max=max((rank+1)*Nshare,ntr)
else
 id_max=(rank+1)*Nshare
end if
do id=id_min,id_max
 r=rtr_(id)
 s=str_(id)
 i=nint(r/dr,4)
 k=nint(s/ds,4)
 if ((r.lt.0.d0).or.(r.gt.aspect).or.(s.lt.0.d0).or.(s.gt.1.d0)) then
  write(*,*) "Tracer advected outside of domain - stopping"
  stop
 end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!modify this
  if ((tracer.eq.1).and.(shape_function.eq.0)) Twork(rank,i,k)=Twork(rank,i,k)+Ttr_(id)
  count_work(rank,Ttype(id),i,k)=count_work(rank,Ttype(id),i,k)+1
 if (shape_function.eq.1) then !!bilinear shape functions
  if (tracer.eq.1) then
   Twork(rank,i-1,k-1)=Twork(rank,i-1,k-1)+Ttr_(id)*shape_fn(i-1,k-1,r,s)
   Twork(rank,i-1,k)=Twork(rank,i-1,k)+Ttr_(id)*shape_fn(i-1,k,r,s)
   Twork(rank,i-1,k+1)=Twork(rank,i-1,k+1)+Ttr_(id)*shape_fn(i-1,k+1,r,s)
   Twork(rank,i,k-1)=Twork(rank,i,k-1)+Ttr_(id)*shape_fn(i,k-1,r,s)
   Twork(rank,i,k)=Twork(rank,i,k)+Ttr_(id)*shape_fn(i,k,r,s)
   Twork(rank,i,k+1)=Twork(rank,i,k+1)+Ttr_(id)*shape_fn(i,k+1,r,s)
   Twork(rank,i+1,k-1)=Twork(rank,i+1,k-1)+Ttr_(id)*shape_fn(i+1,k-1,r,s)
   Twork(rank,i+1,k)=Twork(rank,i+1,k)+Ttr_(id)*shape_fn(i+1,k,r,s)
   Twork(rank,i+1,k+1)=Twork(rank,i+1,k+1)+Ttr_(id)*shape_fn(i+1,k+1,r,s)
  end if

  count_work_real(rank,Ttype(id),i-1,k-1)=count_work_real(rank,Ttype(id),i-1,k-1)+shape_fn(i-1,k-1,r,s)
  count_work_real(rank,Ttype(id),i-1,k)=count_work_real(rank,Ttype(id),i-1,k)+shape_fn(i-1,k,r,s)
  count_work_real(rank,Ttype(id),i-1,k+1)=count_work_real(rank,Ttype(id),i-1,k+1)+shape_fn(i-1,k+1,r,s)
  count_work_real(rank,Ttype(id),i,k-1)=count_work_real(rank,Ttype(id),i,k-1)+shape_fn(i,k-1,r,s)
  count_work_real(rank,Ttype(id),i,k)=count_work_real(rank,Ttype(id),i,k)+shape_fn(i,k,r,s)
  count_work_real(rank,Ttype(id),i,k+1)=count_work_real(rank,Ttype(id),i,k+1)+shape_fn(i,k+1,r,s)
  count_work_real(rank,Ttype(id),i+1,k-1)=count_work_real(rank,Ttype(id),i+1,k-1)+shape_fn(i+1,k-1,r,s)
  count_work_real(rank,Ttype(id),i+1,k)=count_work_real(rank,Ttype(id),i+1,k)+shape_fn(i+1,k,r,s)
  count_work_real(rank,Ttype(id),i+1,k+1)=count_work_real(rank,Ttype(id),i+1,k+1)+shape_fn(i+1,k+1,r,s)  
 end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end modify this
end do
!$OMP END PARALLEL
!$OMP PARALLEL DO PRIVATE(i,k,tt)
do k=0,nz
 do i=0,nx
  if (tracer.eq.1) T1(i,k)=sum(Twork(:,i,k))
  do tt=0,ntype
   count_type(tt,i,k)=sum(count_work(:,tt,i,k))
   if (shape_function.eq.1) count_type_real(tt,i,k)=sum(count_work_real(:,tt,i,k))
  end do
 end do
end do
!$OMP END PARALLEL DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!parallel
!write(*,*) sum(Twork),sum(count_type),sum(count_type_real) !!!!!!!!!!!!!!!!!!!!!!!!!!test values

 if (shape_function.eq.0) then
  do tt=0,ntype
   do i=1,span
    count_type(tt,-i,0:nz)=count_type(tt,i,0:nz)
    count_type(tt,nx+i,0:nz)=count_type(tt,nx-i,0:nz)
   end do
   do k=1,span
    count_type(tt,0:nx,-k)=count_type(tt,0:nx,k)
    count_type(tt,0:nx,nz+k)=count_type(tt,0:nx,nz-k)
   end do
  end do
 end if

 if (shape_function.eq.1) then
  do tt=0,ntype
   do i=1,buffer
    count_type(tt,-i,0:nz)=count_type(tt,i,0:nz)
    count_type(tt,nx+i,0:nz)=count_type(tt,nx-i,0:nz)
   end do
   do k=1,buffer
    count_type(tt,0:nx,-k)=count_type(tt,0:nx,k)
    count_type(tt,0:nx,nz+k)=count_type(tt,0:nx,nz-k)
   end do
  end do

  do tt=0,ntype
   do i=1,span
    count_type_real(tt,-i,0:nz)=count_type_real(tt,i,0:nz)
    count_type_real(tt,nx+i,0:nz)=count_type_real(tt,nx-i,0:nz)
   end do
   do k=1,span
    count_type_real(tt,0:nx,-k)=count_type_real(tt,0:nx,k)
    count_type_real(tt,0:nx,nz+k)=count_type_real(tt,0:nx,nz-k)
   end do
  end do
 end if


nempty=sum(count_type,DIM=1)
if (shape_function.eq.1) nempty_real=sum(count_type_real,DIM=1)

if (shape_function.eq.0) then
!$OMP PARALLEL DO PRIVATE(i,k,tt)
do k=0,nz
 do i=0,nx
  if (nempty(i,k).gt.0) then !!if there are no tracers close to the cell corner use previous value for C
   if (comp.eq.1) then
    do tt=1,ntype !!assume tracer type 0 has zero mass
      if (count_type(tt,i,k).eq.nempty(i,k)) then
       Cnew(tt,i,k)=1.d0
      elseif (count_type(tt,i,k).eq.0) then
       Cnew(tt,i,k)=0.d0
      else
       Cnew(tt,i,k)=real(count_type(tt,i,k),8)/real(nempty(i,k),8)
      end if
    end do
   end if
   if (tracer.eq.1) Tratio(i,k)=T1(i,k)/real(nempty(i,k),8)
  end if
 end do
end do
!$OMP END PARALLEL DO
end if

if (shape_function.eq.1) then
!$OMP PARALLEL DO PRIVATE(i,k,tt)
do k=0,nz
 do i=0,nx
  if (nempty_real(i,k).gt.0.d0) then !!if there are no tracers close to the cell corner use previous value for C
   if (comp.eq.1) then
    do tt=1,ntype !!assume tracer type 0 has zero mass
     Cnew(tt,i,k)=count_type_real(tt,i,k)/nempty_real(i,k)
    end do
   end if
   if (tracer.eq.1) Tratio(i,k)=T1(i,k)/nempty_real(i,k)
  end if
 end do
end do
!$OMP END PARALLEL DO
end if

 if (comp.eq.1) then
  do tt=1,ntype
   call enforceBCs_comp(Cnew(tt,:,:))
  end do
  if (smooth_vis.eq.1) call smoother_time
 end if
 if (tracer.eq.1) then
  call enforceBCs(Tratio)
  if (smooth_vis.eq.1) call smoother_time_T(Tratio)
 end if
 call load_conductivity_array
!$ time2=omp_get_wtime()
!$ tconvert=tconvert+time2-time1
end


subroutine enforceBCs_comp(C0)
 use basics
implicit none
integer*4 :: i
!real*8 :: C0(-span:nx+span,-span:nz+span)
real*8 :: C0(-spanT:nx+spanT,-spanT:nz+spanT) !!expanded index range for H ghost point calculations
!do i=1,span
do i=1,spanT
 C0(-i,0:nz)=C0(i,0:nz)      !!symmetry: conservation of mass
 C0(nx+i,0:nz)=C0(nx-i,0:nz)
 C0(:,-i)=C0(:,i)      
 C0(:,nz+i)=C0(:,nz-i)
end do
end 


subroutine velocityBCs
 use basics
 use arrays
implicit none
integer*4 :: i,k
do i=1,span_interp         !!symmetry:free slip
 w(-i,0:nz)=w(i,0:nz)      
 w(nx+i,0:nz)=w(nx-i,0:nz)
 if (Vbc.eq.0) then        !!symmetry:free slip
  u(0:nx,-i)=u(0:nx,i)      
  u(0:nx,nz+i)=u(0:nx,nz-i)
 elseif (Vbc.eq.1) then    !!antisymmetry for rigid boundaries
  u(0:nx,-i)=-u(0:nx,i)      
  u(0:nx,nz+i)=-u(0:nx,nz-i)
 end if
end do
do k=1,span_interp         !!antisymmetry for impermeable boundaries
 w(:,-k)=-w(:,k)           
 w(:,nz+k)=-w(:,nz-k)
 u(-k,:)=-u(k,:)       
 u(nx+k,:)=-u(nx-k,:)
end do
u(0,:)=0.d0                !!impermeable boundaries
u(nx,:)=0.d0
w(:,0)=0.d0
w(:,nz)=0.d0
if (Vbc.eq.1) then !!rigid BCs
 u(:,0)=0.d0
 u(:,nz)=0.d0
end if
end

subroutine velocityBCs_coarse(level)
 use basics
 use arrays
implicit none
integer*4 :: i,k,level
do i=1,span_interp         !!symmetry:free slip
 w_grid(level,-i,0:nz_grid(level))=w_grid(level,i,0:nz_grid(level))
 w_grid(level,nx_grid(level)+i,0:nz_grid(level))=w_grid(level,nx_grid(level)-i,0:nz_grid(level))
 if (Vbc.eq.0) then        !!symmetry:free slip
  u_grid(level,0:nx_grid(level),-i)=u_grid(level,0:nx_grid(level),i)
  u_grid(level,0:nx_grid(level),nz_grid(level)+i)=u_grid(level,0:nx_grid(level),nz_grid(level)-i)
 elseif (Vbc.eq.1) then    !!antisymmetry for rigid boundaries
  u_grid(level,0:nx_grid(level),-i)=-u_grid(level,0:nx_grid(level),i)
  u_grid(level,0:nx_grid(level),nz_grid(level)+i)=-u_grid(level,0:nx_grid(level),nz_grid(level)-i)
 end if
end do
do k=1,span_interp         !!antisymmetry for impermeable boundaries
 w_grid(level,:,-k)=-w_grid(level,:,k)
 w_grid(level,:,nz_grid(level)+k)=-w_grid(level,:,nz_grid(level)-k)
 u_grid(level,-k,:)=-u_grid(level,k,:)
 u_grid(level,nx_grid(level)+k,:)=-u_grid(level,nx_grid(level)-k,:)
end do
u_grid(level,0,:)=0.d0                !!impermeable boundaries
u_grid(level,nx_grid(level),:)=0.d0
w_grid(level,:,0)=0.d0
w_grid(level,:,nz_grid(level))=0.d0
if (Vbc.eq.1) then !!rigid BCs
 u_grid(level,:,0)=0.d0
 u_grid(level,:,nz_grid(level))=0.d0
end if
end

subroutine interpolate(r,s,DERr,f,finterp) !!!!!!!!!!!!!!!!!!!!!!!!!!!cleaned this up: went back to original version below -- double check k value in this version?
 use basics
 use arrays
implicit none
integer*4 :: i,j,jj,k,ID
real*8 :: hr,hs,fRinterp(-span_interp:span_interp+1),finterp
real*8 :: a1,r,s
real*8 :: f(-span_interp:nx+span_interp,-span_interp:nz+span_interp)
real*8 :: DERs(1:order-1)       !!holds s derivatives of fRinterp
real*8 :: DERr(1:order-1,0:nx,-span_interp:nz+span_interp)

i=nint(r/dr,4)
k=nint(s/ds,4)
hr=r-real(i,8)*dr

 do jj=-span_interp,span_interp
  fRinterp(jj)=f(i,k+jj)
 end do
 do j=1,order-1 !!j=derivative counter
  a1=hr**real(j,8)/factorial(j)
  do jj=-span_interp,span_interp
   fRinterp(jj)=fRinterp(jj)+DERr(j,i,k+jj)*a1
  end do
 end do

hs=s-real(k,8)*ds
 do jj=1,order-1
  DERs(jj)=dot_product(D(jj,-span_interp:span_interp),fRinterp(-span_interp:span_interp))/ds_power(jj)
 end do
 finterp=fRinterp(0)
 do j=1,order-1 !!j=derivative counter
  a1=hs**real(j,8)/factorial(j)
  finterp=finterp+DERs(j)*a1
 end do
end


subroutine interpolate_old(r,s,DERr,f,finterp)
 use basics
 use arrays
implicit none
integer*4 :: i,j,jj,k,ID
real*8 :: hr,hs,fRinterp(-span_interp:span_interp+1),finterp
real*8 :: a1,r,s
real*8 :: f(-span_interp:nx+span_interp,-span_interp:nz+span_interp)
real*8 :: DERs(1:order-1)       !!holds s derivatives of fRinterp
real*8 :: DERr(1:order-1,0:nx,-span_interp:nz+span_interp)

i=int(r/dr,4)
k=int(s/ds,4) 
hr=r-real(i,8)*dr

if (hr.le.dro2) then !!forward Taylor expansion
 do jj=-span_interp,span_interp+1
  fRinterp(jj)=f(i,k+jj)
 end do
 do j=1,order-1 !!j=derivative counter
  a1=hr**real(j,8)/factorial(j)
  do jj=-span_interp,span_interp+1
   fRinterp(jj)=fRinterp(jj)+DERr(j,i,k+jj)*a1
  end do
 end do
else                       !!backward Taylor expansion
 hr=dr-hr
 do jj=-span_interp,span_interp+1
  fRinterp(jj)=f(i+1,k+jj)
 end do
 do j=1,order-1 !!j=derivative counter
  a1=(-hr)**real(j,8)/factorial(j)
  do jj=-span_interp,span_interp+1
   fRinterp(jj)=fRinterp(jj)+DERr(j,i+1,k+jj)*a1
  end do
 end do
end if

hs=s-real(k,8)*ds
if (hs.le.dso2) then !!forward Taylor expansion
 do jj=1,order-1
  DERs(jj)=dot_product(D(jj,-span_interp:span_interp),fRinterp(0-span_interp:0+span_interp))/ds_power(jj)
 end do
 finterp=fRinterp(0)
 do j=1,order-1 !!j=derivative counter
  a1=hs**real(j,8)/factorial(j)
  finterp=finterp+DERs(j)*a1
 end do
else                       !!backward Taylor expansion
 hs=ds-hs
 do jj=1,order-1
  DERs(jj)=dot_product(D(jj,-span_interp:span_interp),fRinterp(1-span_interp:1+span_interp))/ds_power(jj)
 end do
 finterp=fRinterp(1)
 do j=1,order-1 !!j=derivative counter
  a1=(-hs)**real(j,8)/factorial(j)
  finterp=finterp+DERs(j)*a1
 end do
end if
end

subroutine compute_derivatives(f,DERr)
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: jj,i,k
real*8 :: f(-span_interp:nx+span_interp,-span_interp:nz+span_interp)
real*8 :: DERr(1:order-1,0:nx,-span_interp:nz+span_interp)
real*8 :: time1,time2
!$ time1=omp_get_wtime()
!$OMP PARALLEL DO PRIVATE(k,i,jj)
do k=-span_interp,nz+span_interp
 do i=0,nx
  do jj=1,order-1
   DERr(jj,i,k)=dot_product(D(jj,-span_interp:span_interp),f(i-span_interp:i+span_interp,k))/dr_power(jj)
  end do
 end do
end do
!$OMP END PARALLEL DO
!$ time2=omp_get_wtime()
!$ tcd=tcd+time2-time1
end

subroutine compute_metric_interpolator_derivatives
 use basics
 use arrays
implicit none
integer*4 :: i,k,jj
do jj=1,order-1
 do i=0,nx
  DERxr(jj,i)=dot_product(D(jj,-span_interp:span_interp),xr(i-span_interp:i+span_interp))/dr**real(jj,8)
 end do
end do
do jj=1,order-1
 do k=0,nz
  DERzs(jj,k)=dot_product(D(jj,-span_interp:span_interp),zs(k-span_interp:k+span_interp))/ds**real(jj,8)
 end do
end do
end

subroutine interpolate_xr(r,xr_interp)
 use basics
 use arrays
implicit none
integer*4 :: i,j,jj,ID
real*8 :: hr,xr_interp
real*8 :: a1,r
if (gridX.ne.0) then
i=int(r/dr,4)
hr=r-real(i,8)*dr

if (hr.le.dro2) then !!forward Taylor expansion
 xr_interp=xr(i)
 do j=1,order-1 !!j=derivative counter
  a1=hr**real(j,8)/factorial(j)
  xr_interp=xr_interp+DERxr(j,i)*a1
 end do
else                       !!backward Taylor expansion
 hr=dr-hr
 xr_interp=xr(i+1)
 do j=1,order-1 !!j=derivative counter
  a1=(-hr)**real(j,8)/factorial(j)
  xr_interp=xr_interp+DERxr(j,i+1)*a1
 end do
end if

else
 xr_interp=1.d0
end if
end

subroutine interpolate_zs(s,zs_interp)
 use basics
 use arrays
implicit none
integer*4 :: j,jj,k,ID
real*8 :: hs,zs_interp
real*8 :: a1,s
if (gridZ.ne.0) then
k=int(s/ds,4)
hs=s-real(k,8)*ds

if (hs.le.dso2) then !!forward Taylor expansion
 zs_interp=zs(k)
 do j=1,order-1 !!j=derivative counter
  a1=hs**real(j,8)/factorial(j)
  zs_interp=zs_interp+DERzs(j,k)*a1
 end do
else                       !!backward Taylor expansion
 hs=ds-hs
 zs_interp=zs(k+1)
 do j=1,order-1 !!j=derivative counter
  a1=(-hs)**real(j,8)/factorial(j)
  zs_interp=zs_interp+DERzs(j,k+1)*a1
 end do
end if

else
 zs_interp=1.d0
end if
end


