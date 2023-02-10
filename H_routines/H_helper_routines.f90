real*8 function csc(theta)
!!cosecant function
implicit none

!!input
real*8 :: theta

csc=1.d0/sin(theta)
end function csc

real*8 function cot(theta)
!!cotangent function -- arccot function in exact_solution_routines.f90
implicit none

!!input
real*8 :: theta

cot=1.d0/tan(theta)
end function cot

complex*16 function JacobiZeta(u,k)
!!Jacobi Zeta function
!!assumes imaginary part of k is zero
implicit none

!!input
complex*16 :: k !!elliptic modulus
complex*16 :: u !!argument (equivalent to F(phi|m))

!!internal variables
real*8 :: m !!elliptic parameter
real*8 :: phi !!Jacobi amplitude
complex*16 :: Fc,Ec,F,E !!elliptic integral values
complex*16 :: sn,cn,dn !!Jacobi elliptic function values

call compute_JacobiZeta(u,k,JacobiZeta)
end function JacobiZeta

subroutine compute_JacobiZeta(u,k,JacobiZeta)
!!Jacobi Zeta function
!!assumes imaginary part of k is zero
!!assumes the Jacobi amplitude (phi) is real
implicit none

!!input
complex*16 :: k !!elliptic modulus
complex*16 :: u !!argument (equivalent to F(phi|m))

!!output
complex*16 :: JacobiZeta

!!internal variables
real*8, parameter :: tol_IM=1.d-3 !!tolerance for cn%IM size
real*8 :: m !!elliptic parameter
real*8 :: phi !!Jacobi amplitude
complex*16 :: Fc,Ec,F,E !!elliptic integral values
complex*16 :: twoFc !!two times Fc
complex*16 :: sn,cn,dn !!Jacobi elliptic function values
real*8 :: Nm_IM,Np_IM,Nm_RE,Np_RE !!periodicity multiplier
real*8 :: fraction_Nm_RE,fraction_Np_RE,fraction_Nm_IM,fraction_Np_IM
real*8 :: delta_Nm_RE,delta_Np_RE,delta_Nm_IM,delta_Np_IM
logical :: IM_periodic,RE_periodic
logical :: int_p_RE,int_m_RE,int_p_IM,int_m_IM

m=k%RE**2
call complete_elliptic_integrals(m,Fc,Ec)
twoFc=2.d0*Fc !!store this for frequent use
call Jacobi_elliptic_functions(u,m,sn,cn,dn)

if (abs(cn%RE).gt.1.d0) then
 write(*,*) "Error in compute_JacobiZeta -- abs(cn%RE) is greater than unity"
 write(*,*) "u,k,m=",u,k,m
 write(*,*) "cn=",cn
 stop
elseif (abs(cn%IM).gt.tol_IM) then
 write(*,*) "Error in compute_JacobiZeta -- cn%IM is not zero"
 write(*,*) "u,k,m=",u,k,m
 write(*,*) "cn=",cn
 stop
end if
phi=acos(cn%RE) !!produces phi belonging to [0,pi] which works with incomplete_elliptic_integrals routine
call incomplete_elliptic_integrals(phi,m,F,E) !!assuming imaginary part of phi is negligible

!!Re[u] not in the presumed range? 
RE_periodic=.not.((((0.d0.le.u%RE).and.(u%RE.le.(twoFc%RE))).or.(((twoFc%RE).lt.u%RE).and.(u%RE.le.0.d0))))
if (RE_periodic.eqv..true.) then
 Nm_RE=(u%RE-F%RE)/(twoFc%RE)
 Np_RE=(u%RE+F%RE)/(twoFc%RE)
 fraction_Nm_RE=abs(mod(Nm_RE,1.d0)) !fractional parts of N values -- should be very close to zero or unity for integer values of N
 fraction_Np_RE=abs(mod(Np_RE,1.d0))
 delta_Nm_RE=min(fraction_Nm_RE,1.d0-fraction_Nm_RE)
 delta_Np_RE=min(fraction_Np_RE,1.d0-fraction_Np_RE)
 if (delta_Np_RE.lt.delta_Nm_RE) then
  int_p_RE=.true.; int_m_RE=.false.
 else
  int_p_RE=.false.; int_m_RE=.true.
 end if
 if (int_p_RE.eqv..true.) then
  E%RE=2.d0*Ec%RE*Np_RE-E%RE !!use periodicity of elliptic integrals
 elseif (int_m_RE.eqv..true.) then
  E%RE=2.d0*Ec%RE*Nm_RE+E%RE
 end if
end if

!!Im[u] not in the presumed range?
IM_periodic=.not.((((0.d0.le.u%IM).and.(u%IM.le.(twoFc%IM))).or.(((twoFc%IM).lt.u%IM).and.(u%IM.le.0.d0))))
if (IM_periodic.eqv..true.) then
 Nm_IM=(u%IM-F%IM)/(twoFc%IM)
 Np_IM=(u%IM+F%IM)/(twoFc%IM)
 fraction_Nm_IM=abs(mod(Nm_IM,1.d0)) !!fractional parts of N values -- should be very close to zero or unity for integer values of N
 fraction_Np_IM=abs(mod(Np_IM,1.d0))
 delta_Nm_IM=min(fraction_Nm_IM,1.d0-fraction_Nm_IM)
 delta_Np_IM=min(fraction_Np_IM,1.d0-fraction_Np_IM)
 if (delta_Np_IM.lt.delta_Nm_IM) then
  int_p_IM=.true.; int_m_IM=.false.
 else
  int_p_IM=.false.; int_m_IM=.true.
 end if
 if (int_p_IM.eqv..true.) then
  E%IM=2.d0*Ec%IM*Np_IM-E%IM
 elseif (int_m_IM.eqv..true.) then
  E%IM=2.d0*Ec%IM*Nm_IM+E%IM
 end if
end if

JacobiZeta=E-Ec*u/Fc
end subroutine compute_JacobiZeta

real*8 function F(t)
!!wrapper function for integral of f(t)
implicit none
real*8 :: t
real*8 :: f_integral
F=f_integral(t)
end function F

real*8 function dFdt(t)
!!wrapper function for f(t)
implicit none
real*8 :: t
real*8 :: f_func
dFdt=f_func(t)
end function dFdt

real*8 function d2Fdt2(t)
!!wrapper function for derivative of f(t)
implicit none
real*8 :: t
real*8 :: f_derivative
d2Fdt2=f_derivative(t)
end function d2Fdt2

real*8 function signum(x)
!!signum function -- equivalent to Maple's abs(1,x) function
!!assumes Im[x]=0 and x != 0
implicit none
complex*16 :: x
if (x%RE.lt.0.d0) then
 signum=-1.d0
elseif (x%RE.gt.0.d0) then
 signum=1.d0
else
 signum=0.d0 !!this branch is not used in practice
end if
end function signum

real*8 function Heaviside(x)
!!Heavyside function
!!To match definition used in Maple, we take Heaviside(0)=1
!!assumes Im[x]=0
implicit none
complex*16 :: x
!!if (x.le.0.d0) then
if (x%RE.lt.0.d0) then
 Heaviside=0.d0
else
 Heaviside=1.d0
end if
end function Heaviside

real*8 function Dirac(x)
!!Dirac delta function -- to avoid overflow the spike at x=0 is ignored (relevant terms go to zero as x --> 0)
implicit none
complex*16 :: x
Dirac=0.d0
end function Dirac

real*8 function Dirac_derivative(x)
!!derivative of Dirac delta function -- to avoid overflow the spike at x=0 is ignored (relevant terms go to zero as x --> 0)
implicit none
complex*16 :: x
Dirac_derivative=0.d0
end function Dirac_derivative

complex*16 function InverseJacobiAM(phi,k)
!!wrapper for trigonometric form of the incomplete elliptic integral of the first kind
!!Assumes Im{phi}=0 and Im[k]=0
implicit none

!!input
complex*16 :: k !!elliptic modulus
complex*16 :: phi !!elliptic amplitude

!!internal variables
complex*16 :: F,E
real*8 :: m !!elliptic parameter

m=k%RE**2
call incomplete_elliptic_integrals(phi%RE,m,F,E)
InverseJacobiAM=F
end function InverseJacobiAM

subroutine compute_JacobiSN_CN_DN(u,k,sn,cn,dn)
!!wrapper for sn, cn, and dn functions
!!assumes Im[k]=0
implicit none

!!input
complex*16 :: k !!elliptic modulus
complex*16 :: u !!argument

!!output 
complex*16 :: sn,cn,dn

!!internal variables
real*8 :: m !!elliptic parameter

m=k%RE**2
call Jacobi_elliptic_functions(u,m,sn,cn,dn)
end subroutine compute_JacobiSN_CN_DN

complex*16 function JacobiSN(u,k)
!!wrapper for sn function
!!assumes Im[k]=0
implicit none

!!input
complex*16 :: k !!elliptic modulus
complex*16 :: u !!argument

!!internal variables
real*8 :: m !!elliptic parameter
complex*16 :: sn,cn,dn

m=k%RE**2
call Jacobi_elliptic_functions(u,m,sn,cn,dn)
JacobiSN=sn
end function JacobiSN

complex*16 function JacobiCN(u,k)
!!wrapper for cn function
!!assumes Im(k)=0
implicit none

!!input
complex*16 :: k !!elliptic modulus
complex*16 :: u !!argument

!!internal variables
real*8 :: m !!elliptic parameter
complex*16 :: sn,cn,dn

m=k%RE**2
call Jacobi_elliptic_functions(u,m,sn,cn,dn)
JacobiCN=cn
end function JacobiCN

complex*16 function JacobiDN(u,k)
!!wrapper for dn function
!!assumes Im[k]=0
implicit none

!!input
complex*16 :: k !!elliptic modulus
complex*16 :: u !!argument

!!internal variables
real*8 :: m !!elliptic parameter
complex*16 :: sn,cn,dn

m=k%RE**2
call Jacobi_elliptic_functions(u,m,sn,cn,dn)
JacobiDN=dn
end function JacobiDN

complex*16 function EllipticK(k)
!!wrapper for complete elliptic integral of the first kind
!!assumes Im[k]=0
implicit none

!!input
complex*16 :: k !!elliptic modulus

!!internal variables
real*8 :: m !!elliptic parameter
complex*16 :: Fc,Ec

m=k%RE**2
call complete_elliptic_integrals(m,Fc,Ec)
EllipticK=Fc
end function EllipticK

complex*16 function EllipticE(k)
!!wrapper for complete elliptic integral of the second kind
!!assumes Im[k]=0
implicit none

!!input
complex*16 :: k !!elliptic modulus

!!internal variables
real*8 :: m !!elliptic parameter
complex*16 :: Fc,Ec

m=k%RE**2
call complete_elliptic_integrals(m,Fc,Ec)
EllipticE=Ec
end function EllipticE

subroutine compute_EllipticK_EllipticE(k,EllipticK,EllipticE)
!!wrapper for complete elliptic integral of the first and second kinds
!!assumes Im[k]=0
implicit none

!!input
complex*16 :: k !!elliptic modulus

!!output
complex*16 :: EllipticK,EllipticE

!!internal variables
real*8 :: m !!elliptic parameter
complex*16 :: Fc,Ec

m=k%RE**2
call complete_elliptic_integrals(m,Fc,Ec)
EllipticK=Fc; EllipticE=Ec
end subroutine compute_EllipticK_EllipticE
