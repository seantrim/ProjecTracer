subroutine time_average(tmin,tmax)
implicit none
real*8 :: tmin,tmax,t0,tend
real*8 :: ti,Tavgi,vrmsi,Nui,NuBi,sfi,q1i,q2i,q3i,q4i
real*8 :: tf,Tavgf,vrmsf,Nuf,NuBf,sff,q1f,q2f,q3f,q4f
real*8 :: Tavgint,vrmsint,Nuint,NuBint,sfint,q1int,q2int,q3int,q4int
Tavgint=0.d0
vrmsint=0.d0
Nuint=0.d0
NuBint=0.d0
sfint=0.d0
q1int=0.d0
q2int=0.d0
q3int=0.d0
q4int=0.d0
ti=0.d0
tf=0.d0
open(unit=666,file="stats.dat")
do while (ti.le.tmin)
 read(666,*) ti,Tavgi,vrmsi,Nui,NuBi,sfi,q1i,q2i,q3i,q4i
end do
t0=ti !!actual starting time
do while (tf.lt.tmax)
 read(666,*) tf,Tavgf,vrmsf,Nuf,NuBf,sff,q1f,q2f,q3f,q4f
 Tavgint=Tavgint+0.5d0*(Tavgf+Tavgi)*(tf-ti)
 vrmsint=vrmsint+0.5d0*(vrmsf+vrmsi)*(tf-ti)
 Nuint=Nuint+0.5d0*(Nuf+Nui)*(tf-ti)
 NuBint=NuBint+0.5d0*(NuBf+NuBi)*(tf-ti)
 sfint=sfint+0.5d0*(sff+sfi)*(tf-ti)
 q1int=q1int+0.5d0*(q1f+q1i)*(tf-ti)
 q2int=q2int+0.5d0*(q2f+q2i)*(tf-ti)
 q3int=q3int+0.5d0*(q3f+q3i)*(tf-ti)
 q4int=q4int+0.5d0*(q4f+q4i)*(tf-ti)


 ti=tf
 Tavgi=Tavgf
 vrmsi=vrmsf
 Nui=Nuf
 NuBi=NuBf
 sfi=sff
 q1i=q1f
 q2i=q2f
 q3i=q3f
 q4i=q4f
end do
tend=tf !!actual end time
Tavgint=Tavgint/(tend-t0)
vrmsint=vrmsint/(tend-t0)
Nuint=Nuint/(tend-t0)
NuBint=NuBint/(tend-t0)
sfint=sfint/(tend-t0)
q1int=q1int/(tend-t0)
q2int=q2int/(tend-t0)
q3int=q3int/(tend-t0)
q4int=q4int/(tend-t0)
write(*,*) "Time Averaged Values"
write(*,'(9(f18.5))') Tavgint,vrmsint,Nuint,NuBint,sfint,q1int,q2int,q3int,q4int
end
