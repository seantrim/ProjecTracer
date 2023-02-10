real*8 function f_func(t)
!!user-defined function for f
implicit none

!!inputs
real*8 :: t

!!internal variables
real*8, parameter :: pii=3.1415926535897932d0

real*8, parameter :: a=100.d0 !!first sample problem
real*8, parameter :: b=100.d0

!real*8, parameter :: a=52.969966897973720d0 !a=600.d0/(pii*sqrt(13.d0)) !!second sample problem
!real*8, parameter :: b=100.d0
!real*8, parameter :: c=50.d0
!real*8, parameter :: d=397.27475173480286d0 !d=4500.d0/(pii*sqrt(13.d0))

f_func=a*sin(pii*b*t) !!first sample problem
!f_func=a*sin(pii*b*t)*exp(-c*t)+d !!second sample problem
end function f_func

real*8 function f_derivative(t)
!!user-defined function for derivative of f
implicit none

!!inputs
real*8 :: t

!!internal variables
real*8, parameter :: pii=3.1415926535897932d0

real*8, parameter :: a=100.d0 !!first sample problem
real*8, parameter :: b=100.d0

!real*8, parameter :: a=52.969966897973720d0 !a=600.d0/(pii*sqrt(13.d0)) !!second sample problem
!real*8, parameter :: b=100.d0
!real*8, parameter :: c=50.d0
!real*8, parameter :: d=397.27475173480286d0 !d=4500.d0/(pii*sqrt(13.d0))

f_derivative=pii*a*b*cos(pii*b*t) !!first sample problem
!f_derivative=(pii*b*cos(pii*b*t) - c*sin(pii*b*t))*a*exp(-c*t) !!second sample problem
end function f_derivative

real*8 function f_integral(t)
!!user-defined function for integral of f
implicit none

!!inputs
real*8 :: t

!!internal variables
real*8, parameter :: pii=3.1415926535897932d0

real*8, parameter :: a=100.d0 !!first sample problem
real*8, parameter :: b=100.d0

!real*8, parameter :: a=52.969966897973720d0 !a=600.d0/(pii*sqrt(13.d0)) !!second sample problem
!real*8, parameter :: b=100.d0
!real*8, parameter :: c=50.d0
!real*8, parameter :: d=397.27475173480286d0 !d=4500.d0/(pii*sqrt(13.d0))

f_integral=(a/(pii*b))*(1.d0-cos(pii*b*t)) !!first example from "Sample results" section using f(t)=a*sin(pi*b*t)
!f_integral=pii*a*b/((pii*b)**2+c**2)+d*t-(pii*b*cos(pii*b*t) + c*sin(pii*b*t))*a*exp(-c*t)/((pii*b)**2+c**2) !!!second sample problem
end function f_integral
