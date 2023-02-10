subroutine compute_H_func(x,z,t,lambda,k,zI,RaT,RaC,H_func)
!!Compute H
!!Assumes -lambda/2<x<3/2*lambda and -1<z<2
!!Allows evaluation for domain interior, boundaries, and beyond the boundaries (for ghost points)
implicit none

!!inputs
real*8 :: x,z,t,lambda,k,zI,RaT,RaC

!!output
real*8 :: H_func

!!internal variables
real*8, parameter :: tol=1.d-8
real*8 :: xmid
real*8 :: xL,xR,x_
complex*16 :: H,HL,HR

!!external functions
real*8 :: H_horizontal_boundaries 


!!input validation
if ((x.le.-lambda/2.d0).or.(x.ge.3.d0/2.d0*lambda).or.(z.le.-1.d0).or.(z.ge.2.d0)) then
 write(*,*) "Error in compute_H_func: requested position not supported."
 write(*,*) "Requested x,z=",x,z
 write(*,*) "lambda=",lambda
 write(*,*) "Supported x range is -lambda/2<x<3/2*lambda"
 write(*,*) "Supported z range is -1<z<2."
 stop
end if 

xmid=lambda/2.d0

if (x.lt.0.d0) then !!symmetry about left sidewall
 x_=-x
elseif (x.gt.lambda) then !!symmetry about right sidewall
 x_=2.d0*lambda-x
else !!within problem domain (interior + boundaries)
 x_=x
end if

if ((z.eq.0.d0).or.(z.eq.1.d0)) then !!at top/bottom boundary
 H_func=H_horizontal_boundaries(z,RaC,RaT,k,zI)
elseif ((x_.eq.0.d0).or.(x_.eq.lambda)) then !!at sidewalls
 call compute_H_sidewalls(x_,z,t,lambda,k,zI,RaT,RaC,H_func)
elseif (((xmid-tol).le.x_).and.(x_.le.(xmid+tol))) then !!mid line
 xL=xmid-2.d0*tol; xR=xmid+2.d0*tol
 if (z.lt.0.d0) then
  call compute_H_below(xL,z,t,lambda,k,zI,RaT,RaC,HL)
  call compute_H_below(xR,z,t,lambda,k,zI,RaT,RaC,HR)
 elseif (z.gt.1.d0) then
  call compute_H_above(xL,z,t,lambda,k,zI,RaT,RaC,HL)
  call compute_H_above(xR,z,t,lambda,k,zI,RaT,RaC,HR)
 else
  call compute_H_domain(xL,z,t,lambda,k,zI,RaT,RaC,HL)
  call compute_H_domain(xR,z,t,lambda,k,zI,RaT,RaC,HR)
 end if
 H_func=(HL%RE+HR%RE)/2.d0
elseif (z.lt.0.d0) then
 call compute_H_below(x_,z,t,lambda,k,zI,RaT,RaC,H)
 H_func=H%RE
elseif (z.gt.1.d0) then
 call compute_H_above(x_,z,t,lambda,k,zI,RaT,RaC,H)
 H_func=H%RE
else !!elsewhere (not at the domain boundaries)
 call compute_H_domain(x_,z,t,lambda,k,zI,RaT,RaC,H)
 H_func=H%RE
end if

end subroutine compute_H_func

real*8 function H_horizontal_boundaries(z,RaC,RaT,k,zI)
implicit none

!!input
real*8 :: z
real*8 :: RaC,RaT
real*8 :: k,zI

!!internal variables
real*8 :: t1,t2,t3,t4,t5

t1=z-zI
t2=k*t1
t3=exp(2.d0*t2)
t4=t3*t3
t5=t3+1.d0

H_horizontal_boundaries=4.d0*(RaC*k**2*t3/t5**2 &
&- 2.d0*RaC*k**2*t4/t5**3)/RaT
end function H_horizontal_boundaries

subroutine compute_H_sidewalls(x,z,t,lambda,k,zI,RaT,RaC,H)
!!for H at the sidewalls
implicit none

!!inputs
real*8 :: x,z,t,lambda,k,zI,RaT,RaC

!!output
real*8 :: H

!!external functions
real*8 :: csc,cot,arccot
real*8 :: f_integral,f_func,f_derivative

!!internal variables
real*8, parameter :: pii=3.1415926535897932d0
real*8 :: lam
real*8 :: H_ge,H_lt
real*8 :: Q,bigZ0
real*8 :: t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20
real*8 :: t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33

lam=lambda

t9=pii*x/lam
t10=pii*z

t11=cos(t10)
t12=cos(t9)
t13=sin(t10)
t14=sin(t9)
t15=cot(t10)
t16=csc(t10)

t2=t15+t16
t1=log(abs(t2))

t17=f_func(t)
t18=f_integral(t)
t19=f_derivative(t)
t20=lam**2
t21=(t20 + 1)**2

t25=-pii*t12*t18/lam + t1
t3=exp(t25)
t4=1.d0/t3

t24=zI - arccot(-0.5d0*t4 + 0.5d0*t3)/pii
t22=k*t24
t23=k*(t24 - 1)

t5=exp(-2*t22)
t7=exp(-4*t22)
t6=exp(-2*t23)
t8=exp(-4*t23)

t26=(t15**2 + 1)
t27=(t4 - t3)**2
t28=(pii*t15*t16 + t26*pii)
t29=(pii**2*t15**2*t16 + 2*t26*pii**2*t15 + t26*pii**2*t16)
t30=(pii**2*t4*t18*t14/t20 + pii**2*t3*t18*t14/t20)**2*RaC*k
t31=(t5 + 1)
t32=(t6 + 1)
t33=(t27 + 4)

Q=t25

bigZ0=-(exp(-Q)-exp(Q))/2.d0
!write(*,*) "CHS z,bigZ0=",z,bigZ0

if (bigZ0.ge.0.d0) then
 H_ge=(t21*pii**4*t11*t12*t17/lam**3 - RaC + RaT + 4*RaC*k*(t28*t4/(t2) + t28*t3/(t2))*t5/(t33*pii*t31**2))*pii*t12*t17*t13/(RaT*&
     &lam) - (t21*pii**5*t12*t17*t13/lam**3 + 8*RaC*k*(t28*t4/(t2) + t28*t3/(t2))**2*(t4 - t3)*t5/(t33**2*pii*t31**2) - 4*RaC*k*(&
     &2*t28**2*t4/(t2)**2 - t29*t4/(t2) - t29*t3/(t2))*t5/(t33*pii*t31**2) - 16*RaC*k**2*(t28*t4/(t2) + t28*t3/(t2))**2*t5/(t33**&
     &2*pii**2*t31**2) + 32*RaC*k**2*(t28*t4/(t2) + t28*t3/(t2))**2*t7/(t33**2*pii**2*t31**3))/RaT - (t21*pii**5*t12*t17*t13/lam*&
     &*5 + 8*t30*(t4 - t3)*t5/(t33**2*pii*t31**2) - 16*t30**2*t5/(t33**2*pii**2*t31**2) - 4*(pii**4*t4*t18**2*t14**2/lam**4 - pii&
     &**4*t3*t18**2*t14**2/lam**4 - pii**3*t12*t4*t18/lam**3 - pii**3*t12*t3*t18/lam**3)*RaC*k*t5/(t33*pii*t31**2) + 32*t30**2*t7&
     &/(t33**2*pii**2*t31**3))/RaT - (t21*pii**3*t12*t13*t19/lam**3 + 4*(pii*t12*t4*t17/lam + pii*t12*t3*t17/lam)*RaC*k*t5/(t33*p&
     &ii*t31**2))/RaT
 H=H_ge
elseif (bigZ0.lt.0.d0) then
 H_lt=(t21*pii**4*t11*t12*t17/lam**3 - RaC + RaT + 4*RaC*k*(t28*t4/(t2) + t28*t3/(t2))*t6/(t33*pii*t32**2))*pii*t12*t17*t13/(RaT*&
     &lam) - (t21*pii**5*t12*t17*t13/lam**3 + 8*RaC*k*(t28*t4/(t2) + t28*t3/(t2))**2*(t4 - t3)*t6/(t33**2*pii*t32**2) - 4*RaC*k*(&
     &2*t28**2*t4/(t2)**2 - t29*t4/(t2) - t29*t3/(t2))*t6/(t33*pii*t32**2) - 16*RaC*k**2*(t28*t4/(t2) + t28*t3/(t2))**2*t6/(t33**&
     &2*pii**2*t32**2) + 32*RaC*k**2*(t28*t4/(t2) + t28*t3/(t2))**2*t8/(t33**2*pii**2*t32**3))/RaT - (t21*pii**5*t12*t17*t13/lam*&
     &*5 + 8*t30*(t4 - t3)*t6/(t33**2*pii*t32**2) - 16*t30**2*t6/(t33**2*pii**2*t32**2) - 4*(pii**4*t4*t18**2*t14**2/lam**4 - pii&
     &**4*t3*t18**2*t14**2/lam**4 - pii**3*t12*t4*t18/lam**3 - pii**3*t12*t3*t18/lam**3)*RaC*k*t6/(t33*pii*t32**2) + 32*t30**2*t8&
     &/(t33**2*pii**2*t32**3))/RaT - (t21*pii**3*t12*t13*t19/lam**3 + 4*(pii*t12*t4*t17/lam + pii*t12*t3*t17/lam)*RaC*k*t6/(t33*p&
     &ii*t32**2))/RaT
 H=H_lt
end if
end subroutine compute_H_sidewalls

subroutine compute_H_below(x,z,t,lambda,k,zI,RaT,RaC,H)
!!Based on output from Maple's CodeGeneration[Fortran] command using optimize=true (not optimize=tryhard)
!!Edits made to Maple's output to allow routine to be compiled
!!Assumes x belongs to (-lambda/2,3*lambda/2) - [0] - [lambda] and -2<z<0
!!Note: not intended use at domain boundaries
implicit none

!!inputs
real*8 :: x,z,t,lambda,k,zI,RaT,RaC

!!output
complex*16 :: H

!!external functions
real*8 :: Heaviside,Dirac,Dirac_derivative,signum
real*8 :: F,dFdt,d2Fdt2
complex*16 :: JacobiSN,JacobiCN,JacobiDN
complex*16 :: InverseJacobiAM,EllipticK,EllipticE,JacobiZeta

!!internal functions
complex*16 :: t1,t2,t3,t5,t8,t11,t12,t13,t14,t15,t16,t17,t21,t22,t23,t24,t25,t26,t28,t29,t31,t33,t34,t35,t37,t38,t39,t44,t45,t46
complex*16 :: t47,t52,t54,t55,t56,t57,t58,t59,t61,t66,t68,t69,t71,t72,t73,t78,t79,t80,t82,t83,t84,t85,t86,t87,t88,t89,t90,t92,t94
complex*16 :: t95,t96,t101,t102,t103,t104,t105,t106,t107,t108,t111,t112,t113,t116,t120,t121,t123,t124,t126,t127,t129,t130,t131
complex*16 :: t132,t133,t134,t135,t137,t139,t140,t141,t143,t144,t146,t147,t150,t159,t166,t167,t168,t169,t172,t174,t177,t178,t180
complex*16 :: t181,t184,t194,t198,t203,t204,t205,t206,t207,t208,t209,t210,t214,t215,t218,t219,t220,t226,t230,t235,t237,t238,t239
complex*16 :: t244,t245,t246,t252,t254,t262,t264,t270,t271,t282,t283,t287,t289,t295,t298,t300,t301,t302,t308,t313,t315,t319,t333
complex*16 :: t341,t346,t355,t357,t365,t369,t374,t378,t382,t386,t393,t399,t404,t407,t408,t410,t415,t417,t419,t420,t431,t436,t437
complex*16 :: t452,t459,t461,t470,t480,t482,t484,t488,t489,t497,t498,t505,t506,t510,t513,t516,t517,t526,t527,t531,t534,t536,t547
complex*16 :: t554,t562,t590,t597,t617,t619,t652

      t1 = (0.31415926535897932D1 ** 2)
      t2 = 0.31415926535897932D1 * t1
      t3 = (lambda ** 2)
      t5 = (t3 + 1) ** 2
      t8 = 0.1D1 / lambda / t3
      t11 = 0.1D1 / lambda
      t12 = t11 * x * 0.31415926535897932D1
      t13 = cos(t12)
      t14 = (0.31415926535897932D1 * z)
      t15 = sin(t14)
      t16 = t15 * t13
      t17 = d2Fdt2(t) !diff(diff(F(t), t), t)
      t21 = 0.1D1 / 0.31415926535897932D1
      t22 = sin(t12)
      t23 = t15 * t22
      t24 = abs(t23)
      t25 = 0.1D1 / t24
      t26 = InverseJacobiAM(t14, t25)
      t28 = (-x + lambda / 0.2D1)
      t29 = Heaviside(t28)
      t31 = 2 * t29 - 1
      t33 = t11 * t24
      t34 = F(t)
      t35 = t34 * t33
      t37 = t26 + (0, 1) * t35 * t1 * t31
      call compute_JacobiSN_CN_DN(t37,t25,t55,t38,t54) !!SJT
      !t38 = JacobiCN(t37, t25)
      t39 = 0.31415926535897932D1 / 2 + (0, 1) * log(sqrt(1 - t38 ** 2&
     &) + (0, 1) * t38)
      t44 = exp(-2 * (-t39 * t21 + zI) * k)
      t45 = 1 + t44
      t46 = t45 ** 2
      t47 = 0.1D1 / t46
      t52 = dFdt(t) !diff(F(t), t)
      !t54 = JacobiDN(t37, t25)
      !t55 = JacobiSN(t37, t25)
      t56 = t55 * t54
      t57 = t38 ** 2
      t58 = -t57 + 1
      t59 = sqrt(t58)
      t61 = t44 / t59
      t66 = 0.1D1 / RaT
      t68 = 0.31415926535897932D1 * t22
      t69 = cos(t14)
      t71 = t1 ** 2
      t72 = t5 * t71
      t73 = t3 ** 2
      t78 = t47 * RaC
      t79 = k * t78
      t80 = t24 ** 2
      t82 = 0.1D1 / t24 / t80
      t83 = 0.31415926535897932D1 * t82
      t84 = t11 * t83
      t85 = signum(t23) !abs(1, t23)
      t86 = JacobiZeta(t26, t25)
      call compute_EllipticK_EllipticE(t25,t89,t87) !!SJT
      !t87 = EllipticE(t25)
      t88 = t26 * t87
      !t89 = EllipticK(t25)
      t90 = 0.1D1 / t89
      t92 = t90 * t88 + t86
      t94 = 0.1D1 / t80
      t95 = -1 + t94
      t96 = 0.1D1 / t95
      t101 = 2 * t14
      t102 = sin(t101)
      t103 = t96 * t102
      t104 = t15 ** 2
      t105 = t104 * t94
      t106 = 1 - t105
      t107 = sqrt(t106)
      t108 = 0.1D1 / t107
      t111 = -t96 * t80 * t92 / 2 - t80 * t26 / 2 + t108 * t103 / 4
      t112 = t111 * t85
      t113 = t112 * t16
      t116 = Dirac(t28)
      t120 = t2 * t31
      t121 = 0.1D1 / t3
      t123 = t34 * t85
      t124 = t123 * t16
      t126 = -2 * t113 * t84 + (0, -2) * t35 * t1 * t116 + (0, 1) * t124&
     & * t121 * t120
      t127 = t54 * t126
      t129 = 0.31415926535897932D1 * t94
      t130 = t11 * t129
      t131 = t55 ** 2
      t132 = t25 * t131
      t133 = -t95
      t134 = 0.1D1 / t133
      t135 = t38 * t134
      t137 = t37 * t24
      t139 = t134 * t24
      t140 = JacobiZeta(t37, t25)
      t141 = t87 * t37
      t143 = t90 * t141 + t140
      t144 = t143 * t139
      t146 = -t135 * t132 - t137 * t56 + t144 * t56
      t147 = t146 * t85
      t150 = -t147 * t16 * t130 - t55 * t127
      t159 = t11 * 0.31415926535897932D1
      t166 = t108 * 0.31415926535897932D1
      t167 = t22 * t82
      t168 = 0.31415926535897932D1 * t167
      t169 = t85 * t69
      t172 = 2 * t111 * t169 * t168
      t174 = t34 * t11
      t177 = t166 - t172 + (0, 1) * t174 * t169 * t22 * t120
      t178 = t54 * t177
      t180 = t22 * t94
      t181 = 0.31415926535897932D1 * t180
      t184 = -t146 * t169 * t181 - t55 * t178
      t194 = t5 * 0.31415926535897932D1 * t71
      t198 = t52 * t16
      t203 = (k ** 2)
      t204 = t203 / t45 / t46 * RaC
      t205 = 0.1D1 / t1
      t206 = t150 ** 2
      t207 = t206 * t205
      t208 = 0.1D1 / t58
      t209 = t44 ** 2
      t210 = t209 * t208
      t214 = t80 ** 2
      t215 = 0.1D1 / t214
      t218 = t13 ** 2
      t219 = t104 * t218
      t220 = t85 ** 2
      t226 = t121 * t1 * t82
      t230 = 2*Dirac(t23) !signum(1, t23)
      t235 = t13 * t11
      t237 = t85 * t15
      t238 = t90 * t87
      t239 = 1 - t105 - t238
      t244 = 0.31415926535897932D1 * t25
      t245 = t235 * t244
      t246 = t87 ** 2
      t252 = t89 ** 2
      t254 = t94 * t15
      t262 = -t246 * t26 + t87 * (-t94 + 2 - t105) * t89 * t26 + (-t254& 
     &* t107 * t69 + t86 * t95 + (t26 * t95 + t86) * t106) * t252
      t264 = 0.1D1 / t252
      t270 = t24 * t89
      t271 = t24 * t87 - t270
      t282 = t94 * t264
      t283 = 0.31415926535897932D1 * t282
      t287 = t24 * t134 * t87 - t270
      t289 = t287 * t237 * t235
      t295 = t24 * t92
      t298 = t237 * t235
      t300 = t25 * t92
      t301 = t95 ** 2
      t302 = 0.1D1 / t301
      t308 = t24 * t26
      t313 = t82 * t108 * t302 * t102
      t315 = t85 * t16 * t159
      t319 = 0.1D1 / t107 / t106
      t333 = Dirac_derivative(t28) !Dirac(1, t28)
      t341 = (0, -1) * t31
      t346 = t71 * t31
      t355 = t55 * t38
      t357 = t38 * t25
      t365 = t143 * t134 * t55 * t357 - t134 * t54 * t132 - t37 * t55 *& 
     &t357
      t369 = -t365 * t85 * t16 * t130 - t355 * t94 * t126
      t374 = t134 * t25
      t378 = t54 * t38
      t382 = -t55 * t131 * t374 - t143 * t378 * t139 + t378 * t137 + t55&
     & * t374
      t386 = -t382 * t85 * t16 * t130 + t54 * t38 * t126
      t393 = t121 * t1 * t94
      t399 = t25 * t55
      t404 = t135 * t94 * t131
      t407 = t133 ** 2
      t408 = 0.1D1 / t407
      t410 = t38 * t408 * t215 * t131
      t415 = t55 * t369
      t417 = t386 * t54
      t419 = t159 * t56
      t420 = t37 * t85
      t431 = t408 * t94
      t436 = t54 ** 2
      t437 = t436 - t238
      t452 = -t246 * t37 + t87 * (-t94 + t436 + 1) * t89 * t37 + (-t94 *&
     & t55 * t378 + t140 * t95 + (t37 * t95 + t140) * t436) * t252
      t459 = t94 * t37
      t461 = t271 * t85
      t470 = -2 * t386 * t135 * t399 + t315 * t404 + 2 * t315 * t410 - t&
     &150 * t134 * t132 - t137 * t415 - t137 * t417 - t420 * t16 * t419 &
     &- t126 * t24 * t56 + t144 * t415 + t144 * t417 + t143 * t134 * t85&
     & * t16 * t419 - 2 * t315 * t143 * t431 * t56 + (-t90 * t461 * t16 &
     &* t159 * t459 - t264 * t96 * t452 * t237 * t245 + t90 * t87 * t126&
     & + t289 * t283 * t141 + t437 * t126) * t139 * t56
      t480 = t21 * k * t78
      t482 = 0.1D1 / t59 / t58
      t484 = t38 * t44
      t488 = t203 * t78
      t489 = t44 * t208
      t497 = t184 ** 2
      t498 = t497 * t205
      t505 = t69 * 0.31415926535897932D1
      t506 = t85 * t505
      t510 = 2 * t506 * t22 * t104 * t82 - 2 * t505 * t254
      t513 = t22 ** 2
      t516 = t69 ** 2
      t517 = t220 * t516
      t526 = t1 * t513 * t82
      t527 = t230 * t516
      t531 = t166 - t172
      t534 = t505 * t22 * t25
      t536 = t264 * t96
      t547 = t287 * t169 * t68
      t554 = t169 * t68
      t562 = cos(t101)
      t590 = -t365 * t169 * t181 - t355 * t94 * t177
      t597 = -t382 * t169 * t181 + t54 * t38 * t177
      t617 = t55 * t590
      t619 = t597 * t54
      t652 = -2 * t597 * t135 * t399 + t554 * t404 + 2 * t554 * t410 - t&
     &184 * t134 * t132 - t137 * t617 - t137 * t619 - t420 * t505 * t22 &
     &* t56 - t177 * t24 * t56 + t144 * t617 + t144 * t619 + t143 * t134&
     & * t169 * t68 * t56 - 2 * t506 * t22 * t143 * t431 * t56 + (-t90 *&
     & t271 * t169 * t68 * t459 - t536 * t452 * t85 * t534 + t547 * t282&
     & * t141 + t90 * t87 * t177 + t437 * t177) * t139 * t56
      H = t66 * (-t17 * t16 * t8 * t5 * t2 + (0, -2) * t61 * t56 * &
     &t52 * t33 * t31 * 0.31415926535897932D1 * k * t47 * RaC) + t66 * (&
     &t52 * t23 / t73 * t72 + 2 * t61 * t150 * t21 * t79) * t52 * t69 * &
     &t68 - t66 * (-t52 * t69 * t13 * t8 * t72 + 2 * t61 * t184 * t21 * &
     &t79 + RaC - RaT) * t52 * t15 * t13 * t159 - t66 * (t198 / lambda /&
     & t73 * t194 + 8 * t210 * t207 * t204 + 2 * t61 * (-t55 * t54 * (6 &
     &* t111 * t220 * t219 * t121 * t1 * t215 + 2 * t112 * t23 * t226 - &
     &2 * t111 * t230 * t219 * t226 - 2 * (-t96 * t80 * (-2 * t90 * t112&
     & * t16 * t159 * t82 * t87 - t90 * t26 * t271 * t237 * t235 * t129 &
     &- 2 * t239 * t111 * t237 * t235 * t83 - t264 * t96 * t262 * t237 *&
     & t245 + t289 * t283 * t88) / 2 - t298 * 0.31415926535897932D1 * t9&
     &6 * t295 - t298 * 0.31415926535897932D1 * t302 * t300 + t113 * t11&
     & * t244 - t298 * 0.31415926535897932D1 * t308 + t315 * t313 / 2 - &
     &t85 * t235 * 0.31415926535897932D1 * t15 * t104 * t82 * t319 * t10&
     &3 / 4) * t85 * t16 * t84 + (0, 2) * t35 * t1 * t333 + (0, -4) * t1&
     &24 * t121 * t2 * t116 + t123 * t23 * t8 * t71 * t341 + (0, 1) * t3&
     &4 * t230 * t219 * t8 * t346) - t55 * t369 * t126 - t386 * t127 + 2&
     & * t146 * t220 * t219 * t226 + t147 * t23 * t393 - t146 * t230 * t&
     &219 * t393 - t470 * t85 * t16 * t130) * t21 * t79 + 2 * t484 * t48&
     &2 * t206 * t480 - 4 * t489 * t207 * t488) - t66 * (t198 * t8 * t19&
     &4 + 8 * t210 * t498 * t204 + 2 * t61 * (-t55 * t54 * (-t510 * t319&
     & * 0.31415926535897932D1 / 2 + 6 * t111 * t517 * t1 * t513 * t215 &
     &+ 2 * t111 * t237 * t1 * t167 - 2 * t111 * t527 * t526 - 2 * (-t96&
     & * t80 * (-t90 * t26 * t461 * t505 * t180 - t536 * t262 * t85 * t5&
     &34 + t547 * t282 * t88 + t90 * t531 * t87 + t239 * t531) / 2 - t55&
     &4 * t96 * t295 - t554 * t302 * t300 - t80 * t531 / 2 - t506 * t22 &
     &* t308 + t108 * t96 * t562 * 0.31415926535897932D1 / 2 + t554 * t3&
     &13 / 2 - t510 * t319 * t103 / 8) * t169 * t168 + t174 * t237 * t22&
     & * t71 * t341 + (0, 1) * t174 * t527 * t513 * t346) - t55 * t590 *&
     & t177 - t597 * t178 + 2 * t146 * t517 * t526 + t146 * t237 * t1 * &
     &t180 - t146 * t527 * t1 * t513 * t94 - t652 * t169 * t181) * t21 *&
     & t79 + 2 * t484 * t482 * t497 * t480 - 4 * t489 * t498 * t488)
end subroutine compute_H_below

subroutine compute_H_above(x,z,t,lambda,k,zI,RaT,RaC,H)
!!Based on output from Maple's CodeGeneration[Fortran] command using optimize=true (not optimize=tryhard)
!!Edits made to Maple's output to allow routine to be compiled
!!Assumes x belongs to (-lambda/2,3*lambda/2) - [0] - [lambda] and 1<z<2
!!Note: not intended use at domain boundaries
implicit none

!!inputs
real*8 :: x,z,t,lambda,k,zI,RaT,RaC

!!output
complex*16 :: H

!!external functions
real*8 :: Heaviside,Dirac,Dirac_derivative,signum
real*8 :: F,dFdt,d2Fdt2
complex*16 :: JacobiSN,JacobiCN,JacobiDN
complex*16 :: InverseJacobiAM,EllipticK,EllipticE,JacobiZeta

!!internal functions
complex*16 :: t1,t2,t3,t5,t8,t11,t12,t13,t14,t15,t16,t17,t21,t23,t24,t25,t26,t27,t28,t29,t31,t32,t34,t35,t37,t38,t39,t41,t42,t43
complex*16 :: t48,t49,t50,t51,t56,t58,t59,t60,t61,t62,t63,t65,t70,t72,t73,t75,t76,t77,t83,t84,t85,t87,t88,t89,t90,t91,t92,t93,t94
complex*16 :: t95,t96,t98,t100,t101,t102,t107,t108,t109,t110,t111,t112,t113,t114,t117,t118,t119,t122,t126,t127,t129,t130,t132,t133
complex*16 :: t135,t136,t137,t138,t139,t140,t141,t143,t145,t146,t147,t149,t150,t152,t153,t156,t165,t172,t173,t174,t175,t176,t179
complex*16 :: t182,t185,t186,t188,t189,t192,t202,t206,t211,t212,t213,t214,t215,t216,t217,t218,t222,t223,t226,t227,t228,t234,t238
complex*16 :: t243,t245,t246,t247,t252,t253,t254,t260,t262,t270,t272,t278,t279,t290,t291,t295,t297,t303,t306,t308,t309,t310,t316
complex*16 :: t321,t323,t327,t341,t348,t361,t363,t371,t375,t380,t384,t388,t392,t399,t405,t410,t413,t414,t416,t421,t423,t425,t426
complex*16 :: t437,t442,t443,t458,t465,t467,t476,t486,t488,t490,t494,t495,t503,t504,t511,t512,t516,t519,t522,t523,t532,t533,t537
complex*16 :: t540,t542,t553,t560,t568,t596,t603,t623,t625,t658

      t1 = (0.31415926535897932D1 ** 2)
      t2 = 0.31415926535897932D1 * t1
      t3 = (lambda ** 2)
      t5 = (t3 + 1) ** 2
      t8 = 0.1D1 / lambda / t3
      t11 = 0.1D1 / lambda
      t12 = t11 * x * 0.31415926535897932D1
      t13 = cos(t12)
      t14 = (0.31415926535897932D1 * z)
      t15 = sin(t14)
      t16 = t15 * t13
      t17 = d2Fdt2(t) !diff(diff(F(t), t), t)
      t21 = 0.1D1 / 0.31415926535897932D1
      t23 = (0.31415926535897932D1 * (0.2D1 - z))
      t24 = sin(t12)
      t25 = sin(t23)
      t26 = t25 * t24
      t27 = abs(t26)
      t28 = 0.1D1 / t27
      t29 = InverseJacobiAM(t23, t28)
      t31 = (-x + lambda / 0.2D1)
      t32 = Heaviside(t31)
      t34 = 2 * t32 - 1
      t35 = (0, -1) * t34
      t37 = t11 * t27
      t38 = F(t)
      t39 = t38 * t37
      t41 = t39 * t1 * t35 + t29
      call compute_JacobiSN_CN_DN(t41,t28,t59,t42,t58) !!SJT
      !t42 = JacobiCN(t41, t28)
      t43 = 0.31415926535897932D1 / 2 + (0, 1) * log(sqrt(1 - t42 ** 2&
     &) + (0, 1) * t42)
      t48 = exp(-2 * (-t43 * t21 + zI) * k)
      t49 = 1 + t48
      t50 = t49 ** 2
      t51 = 0.1D1 / t50
      t56 = dFdt(t) !diff(F(t), t)
      !t58 = JacobiDN(t41, t28)
      !t59 = JacobiSN(t41, t28)
      t60 = t59 * t58
      t61 = t42 ** 2
      t62 = -t61 + 1
      t63 = sqrt(t62)
      t65 = t48 / t63
      t70 = 0.1D1 / RaT
      t72 = 0.31415926535897932D1 * t24
      t73 = cos(t14)
      t75 = t1 ** 2
      t76 = t5 * t75
      t77 = t3 ** 2
      t83 = t51 * RaC
      t84 = k * t83
      t85 = t27 ** 2
      t87 = 0.1D1 / t27 / t85
      t88 = 0.31415926535897932D1 * t87
      t89 = t11 * t88
      t90 = t25 * t13
      t91 = signum(t26) !abs(1, t26)
      t92 = JacobiZeta(t29, t28)
      call compute_EllipticK_EllipticE(t28,t95,t93) !!SJT
      !t93 = EllipticE(t28)
      t94 = t29 * t93
      !t95 = EllipticK(t28)
      t96 = 0.1D1 / t95
      t98 = t96 * t94 + t92
      t100 = 0.1D1 / t85
      t101 = -1 + t100
      t102 = 0.1D1 / t101
      t107 = 2 * t23
      t108 = sin(t107)
      t109 = t102 * t108
      t110 = t25 ** 2
      t111 = t110 * t100
      t112 = 1 - t111
      t113 = sqrt(t112)
      t114 = 0.1D1 / t113
      t117 = -t102 * t85 * t98 / 2 - t85 * t29 / 2 + t114 * t109 / 4
      t118 = t117 * t91
      t119 = t118 * t90
      t122 = Dirac(t31)
      t126 = 0.1D1 / t3
      t127 = t126 * t2
      t129 = t38 * t91
      t130 = t129 * t90
      t132 = -2 * t119 * t89 + (0, 2) * t39 * t1 * t122 + t130 * t127 *& 
     &t35
      t133 = t58 * t132
      t135 = 0.31415926535897932D1 * t100
      t136 = t11 * t135
      t137 = t59 ** 2
      t138 = t28 * t137
      t139 = -t101
      t140 = 0.1D1 / t139
      t141 = t42 * t140
      t143 = t41 * t27
      t145 = t140 * t27
      t146 = JacobiZeta(t41, t28)
      t147 = t93 * t41
      t149 = t96 * t147 + t146
      t150 = t149 * t145
      t152 = -t141 * t138 - t143 * t60 + t150 * t60
      t153 = t152 * t91
      t156 = -t153 * t90 * t136 - t59 * t133
      t165 = t11 * 0.31415926535897932D1
      t172 = t114 * 0.31415926535897932D1
      t173 = t24 * t87
      t174 = 0.31415926535897932D1 * t173
      t175 = cos(t23)
      t176 = t91 * t175
      t179 = 2 * t117 * t176 * t174
      t182 = t38 * t11
      t185 = -t172 + t179 + (0, 1) * t182 * t176 * t24 * t2 * t34
      t186 = t58 * t185
      t188 = t24 * t100
      t189 = 0.31415926535897932D1 * t188
      t192 = t152 * t176 * t189 - t59 * t186
      t202 = t5 * 0.31415926535897932D1 * t75
      t206 = t56 * t16
      t211 = (k ** 2)
      t212 = t211 / t49 / t50 * RaC
      t213 = 0.1D1 / t1
      t214 = t156 ** 2
      t215 = t214 * t213
      t216 = 0.1D1 / t62
      t217 = t48 ** 2
      t218 = t217 * t216
      t222 = t85 ** 2
      t223 = 0.1D1 / t222
      t226 = t13 ** 2
      t227 = t110 * t226
      t228 = t91 ** 2
      t234 = t126 * t1 * t87
      t238 = 2*Dirac(t26) !signum(1, t26)
      t243 = t13 * t11
      t245 = t91 * t25
      t246 = t96 * t93
      t247 = 1 - t111 - t246
      t252 = 0.31415926535897932D1 * t28
      t253 = t243 * t252
      t254 = t93 ** 2
      t260 = t95 ** 2
      t262 = t100 * t25
      t270 = -t254 * t29 + t93 * (-t100 + 2 - t111) * t95 * t29 + (-t262&
     & * t113 * t175 + t92 * t101 + (t29 * t101 + t92) * t112) * t260
      t272 = 0.1D1 / t260
      t278 = t27 * t95
      t279 = t27 * t93 - t278
      t290 = t100 * t272
      t291 = 0.31415926535897932D1 * t290
      t295 = t27 * t140 * t93 - t278
      t297 = t295 * t245 * t243
      t303 = t27 * t98
      t306 = t245 * t243
      t308 = t28 * t98
      t309 = t101 ** 2
      t310 = 0.1D1 / t309
      t316 = t27 * t29
      t321 = t87 * t114 * t310 * t108
      t323 = t91 * t90 * t165
      t327 = 0.1D1 / t113 / t112
      t341 = Dirac_derivative(t31) !Dirac(1, t31)
      t348 = t75 * t34
      t361 = t59 * t42
      t363 = t42 * t28
      t371 = t149 * t140 * t59 * t363 - t140 * t58 * t138 - t41 * t59 *& 
     &t363
      t375 = -t371 * t91 * t90 * t136 - t361 * t100 * t132
      t380 = t140 * t28
      t384 = t58 * t42
      t388 = -t59 * t137 * t380 - t149 * t384 * t145 + t384 * t143 + t59&
     & * t380
      t392 = -t388 * t91 * t90 * t136 + t58 * t42 * t132
      t399 = t126 * t1 * t100
      t405 = t28 * t59
      t410 = t141 * t100 * t137
      t413 = t139 ** 2
      t414 = 0.1D1 / t413
      t416 = t42 * t414 * t223 * t137
      t421 = t59 * t375
      t423 = t392 * t58
      t425 = t165 * t60
      t426 = t41 * t91
      t437 = t414 * t100
      t442 = t58 ** 2
      t443 = t442 - t246
      t458 = -t254 * t41 + t93 * (-t100 + t442 + 1) * t95 * t41 + (-t100&
     & * t59 * t384 + t146 * t101 + (t41 * t101 + t146) * t442) * t260
      t465 = t100 * t41
      t467 = t279 * t91
      t476 = -2 * t392 * t141 * t405 + t323 * t410 + 2 * t323 * t416 - t&
     &156 * t140 * t138 - t143 * t421 - t143 * t423 - t426 * t90 * t425 &
     &- t132 * t27 * t60 + t150 * t421 + t150 * t423 + t149 * t140 * t91&
     & * t90 * t425 - 2 * t323 * t149 * t437 * t60 + (-t272 * t102 * t45&
     &8 * t245 * t253 - t96 * t467 * t90 * t165 * t465 + t96 * t93 * t13&
     &2 + t297 * t291 * t147 + t443 * t132) * t145 * t60
      t486 = t21 * k * t83
      t488 = 0.1D1 / t63 / t62
      t490 = t42 * t48
      t494 = t211 * t83
      t495 = t48 * t216
      t503 = t192 ** 2
      t504 = t503 * t213
      t511 = t175 * 0.31415926535897932D1
      t512 = t91 * t511
      t516 = -2 * t512 * t24 * t110 * t87 + 2 * t511 * t262
      t519 = t24 ** 2
      t522 = t175 ** 2
      t523 = t228 * t522
      t532 = t1 * t519 * t87
      t533 = t238 * t522
      t537 = -t172 + t179
      t540 = t511 * t24 * t28
      t542 = t272 * t102
      t553 = t295 * t176 * t72
      t560 = t176 * t72
      t568 = cos(t107)
      t596 = -t361 * t100 * t185 + t371 * t176 * t189
      t603 = t388 * t176 * t189 + t58 * t42 * t185
      t623 = t59 * t596
      t625 = t603 * t58
      t658 = -2 * t603 * t141 * t405 - t560 * t410 - 2 * t560 * t416 - t&
     &192 * t140 * t138 - t143 * t623 - t143 * t625 + t426 * t511 * t24 &
     &* t60 - t185 * t27 * t60 + t150 * t623 + t150 * t625 - t149 * t140&
     & * t176 * t72 * t60 + 2 * t512 * t24 * t149 * t437 * t60 + (t96 * &
     &t279 * t176 * t72 * t465 + t542 * t458 * t91 * t540 - t553 * t290 &
     &* t147 + t96 * t93 * t185 + t443 * t185) * t145 * t60
      H = t70 * (-t17 * t16 * t8 * t5 * t2 + (0, 2) * t65 * t60 * t&
     &56 * t37 * t34 * 0.31415926535897932D1 * k * t51 * RaC) + t70 * (t&
     &56 * t15 * t24 / t77 * t76 + 2 * t65 * t156 * t21 * t84) * t56 * t&
     &73 * t72 - t70 * (-t56 * t73 * t13 * t8 * t76 + 2 * t65 * t192 * t&
     &21 * t84 + RaC - RaT) * t56 * t15 * t13 * t165 - t70 * (t206 / lam&
     &bda / t77 * t202 + 8 * t218 * t215 * t212 + 2 * t65 * (-t59 * t58 &
     &* (6 * t117 * t228 * t227 * t126 * t1 * t223 + 2 * t118 * t26 * t2&
     &34 - 2 * t117 * t238 * t227 * t234 - 2 * (-t102 * t85 * (-2 * t96 &
     &* t118 * t90 * t165 * t87 * t93 - t96 * t29 * t279 * t245 * t243 *&
     & t135 - t272 * t102 * t270 * t245 * t253 - 2 * t247 * t117 * t245 &
     &* t243 * t88 + t297 * t291 * t94) / 2 - t306 * 0.31415926535897932&
     &D1 * t102 * t303 - t306 * 0.31415926535897932D1 * t310 * t308 + t1&
     &19 * t11 * t252 - t306 * 0.31415926535897932D1 * t316 + t323 * t32&
     &1 / 2 - t91 * t243 * 0.31415926535897932D1 * t25 * t110 * t87 * t3&
     &27 * t109 / 4) * t91 * t90 * t89 + (0, -2) * t39 * t1 * t341 + (0,&
     & 4) * t130 * t127 * t122 + (0, 1) * t129 * t26 * t8 * t348 + t38 *&
     & t238 * t227 * t8 * t75 * t35) - t59 * t375 * t132 - t392 * t133 +&
     & 2 * t152 * t228 * t227 * t234 + t153 * t26 * t399 - t152 * t238 *&
     & t227 * t399 - t476 * t91 * t90 * t136) * t21 * t84 + 2 * t490 * t&
     &488 * t214 * t486 - 4 * t495 * t215 * t494) - t70 * (t206 * t8 * t&
     &202 + 8 * t218 * t504 * t212 + 2 * t65 * (-t59 * t58 * (t516 * t32&
     &7 * 0.31415926535897932D1 / 2 + 6 * t117 * t523 * t1 * t519 * t223&
     & + 2 * t117 * t245 * t1 * t173 - 2 * t117 * t533 * t532 + 2 * (-t1&
     &02 * t85 * (t96 * t29 * t467 * t511 * t188 + t542 * t270 * t91 * t&
     &540 - t553 * t290 * t94 + t96 * t537 * t93 + t247 * t537) / 2 + t5&
     &60 * t102 * t303 + t560 * t310 * t308 - t85 * t537 / 2 + t512 * t2&
     &4 * t316 - t114 * t102 * t568 * 0.31415926535897932D1 / 2 - t560 *&
     & t321 / 2 - t516 * t327 * t109 / 8) * t176 * t174 + (0, 1) * t182 &
     &* t245 * t24 * t348 + t182 * t533 * t519 * t75 * t35) - t59 * t596&
     & * t185 - t603 * t186 + 2 * t152 * t523 * t532 + t152 * t245 * t1 &
     &* t188 - t152 * t533 * t1 * t519 * t100 + t658 * t176 * t189) * t2&
     &1 * t84 + 2 * t490 * t488 * t503 * t486 - 4 * t495 * t504 * t494)
end subroutine compute_H_above

subroutine compute_H_domain(x,z,t,lambda,k,zI,RaT,RaC,H)
!!Based on output from Maple's CodeGeneration[Fortran] command using optimize=true (not optimize=tryhard)
!!Edits made to Maple's output to allow routine to be compiled
!!Assumes x belongs to (-lambda/2,3*lambda/2) - [0] - [lambda] and 0<z<1
!!Note: not intended use at domain boundaries
implicit none

!!inputs
real*8 :: x,z,t,lambda,k,zI,RaT,RaC

!!output
complex*16 :: H

!!external functions
real*8 :: Heaviside,Dirac,Dirac_derivative,signum
real*8 :: F,dFdt,d2Fdt2
complex*16 :: JacobiSN,JacobiCN,JacobiDN
complex*16 :: InverseJacobiAM,EllipticK,EllipticE,JacobiZeta

!!internal functions
complex*16 :: t1,t2,t3,t5,t8,t11,t12,t13,t14,t15
complex*16 :: t16,t17,t21,t22,t23,t24,t25,t26,t28,t29,t31,t32,t34,t35,t36,t38,t39,t40,t45,t46,t47,t48,t53,t55,t56,t57,t58,t59,t60
complex*16 :: t62,t67,t69,t70,t72,t73,t74,t79,t80,t81,t83,t84,t85,t86,t87,t88,t89,t90,t91,t93,t95,t96,t97,t102,t103,t104,t105
complex*16 :: t106,t107,t108,t109,t112,t113,t114,t117,t121,t122,t124,t125,t127,t128,t130,t131,t132,t133,t134,t135,t136,t138,t140
complex*16 :: t141,t142,t144,t145,t147,t148,t151,t160,t167,t168,t169,t170,t173,t176,t179,t180,t182,t183,t186,t196,t200,t205,t206
complex*16 :: t207,t208,t209,t210,t211,t212,t216,t217,t220,t221,t222,t228,t232,t237,t239,t240,t241,t246,t247,t248,t254,t256,t264
complex*16 :: t266,t272,t273,t284,t285,t289,t291,t297,t300,t302,t303,t304,t310,t315,t317,t321,t335,t342,t355,t357,t365,t369,t374
complex*16 :: t378,t382,t386,t393,t399,t404,t407,t408,t410,t415,t417,t419,t420,t431,t436,t437,t452,t459,t461,t470,t480,t482,t484
complex*16 :: t488,t489,t497,t498,t505,t506,t510,t513,t516,t517,t526,t527,t531,t534,t536,t547,t554,t562,t590,t597,t617,t619,t652

      t1 = (0.31415926535897932D1 ** 2)
      t2 = 0.31415926535897932D1 * t1
      t3 = (lambda ** 2)
      t5 = (t3 + 1) ** 2
      t8 = 0.1D1 / lambda / t3
      t11 = 0.1D1 / lambda
      t12 = t11 * x * 0.31415926535897932D1
      t13 = cos(t12)
      t14 = (0.31415926535897932D1 * z)
      t15 = sin(t14)
      t16 = t15 * t13
      t17 = d2Fdt2(t) !diff(diff(F(t), t), t)
      t21 = 0.1D1 / 0.31415926535897932D1
      t22 = sin(t12)
      t23 = t15 * t22
      t24 = abs(t23)  !!D function
      t25 = 0.1D1 / t24 !!1/D
      t26 = InverseJacobiAM(t14, t25)
      t28 = (-x + lambda / 0.2D1)
      t29 = Heaviside(t28)
      t31 = 2 * t29 - 1
      t32 = (0, -1) * t31
      t34 = t11 * t24
      t35 = F(t)
      t36 = t35 * t34
      t38 = t36 * t1 * t32 + t26 !!Eq. 27
      call compute_JacobiSN_CN_DN(t38,t25,t56,t39,t55) !!SJT
      !t39 = JacobiCN(t38, t25)
      t40 = 0.31415926535897932D1 / 2 + (0, 1) * log(sqrt(1 - t39 ** 2&
     &) + (0, 1) * t39)
      t45 = exp(-2 * (-t40 * t21 + zI) * k)
      t46 = 1 + t45
      t47 = t46 ** 2
      t48 = 0.1D1 / t47
      t53 = dFdt(t) !diff(F(t), t)
      !t55 = JacobiDN(t38, t25)
      !t56 = JacobiSN(t38, t25)
      t57 = t56 * t55
      t58 = t39 ** 2
      t59 = -t58 + 1
      t60 = sqrt(t59)
      t62 = t45 / t60
      t67 = 0.1D1 / RaT
      t69 = 0.31415926535897932D1 * t22
      t70 = cos(t14)
      t72 = t1 ** 2
      t73 = t5 * t72
      t74 = t3 ** 2
      t79 = t48 * RaC
      t80 = k * t79
      t81 = t24 ** 2
      t83 = 0.1D1 / t24 / t81
      t84 = 0.31415926535897932D1 * t83
      t85 = t11 * t84
      t86 = signum(t23) !abs(1, t23)
      t87 = JacobiZeta(t26, t25)
      call compute_EllipticK_EllipticE(t25,t90,t88) !!SJT
      !t88 = EllipticE(t25)
      t89 = t26 * t88
      !t90 = EllipticK(t25)
      t91 = 0.1D1 / t90
      t93 = t91 * t89 + t87
      t95 = 0.1D1 / t81
      t96 = -1 + t95
      t97 = 0.1D1 / t96
      t102 = 2 * t14
      t103 = sin(t102)
      t104 = t97 * t103
      t105 = t15 ** 2
      t106 = t105 * t95
      t107 = 1 - t106
      t108 = sqrt(t107)
      t109 = 0.1D1 / t108
      t112 = -t97 * t81 * t93 / 2 - t81 * t26 / 2 + t109 * t104 / 4
      t113 = t112 * t86
      t114 = t113 * t16
      t117 = Dirac(t28)
      t121 = 0.1D1 / t3
      t122 = t121 * t2
      t124 = t35 * t86
      t125 = t124 * t16
      t127 = -2 * t114 * t85 + (0, 2) * t36 * t1 * t117 + t125 * t122 *& 
     &t32
      t128 = t55 * t127
      t130 = 0.31415926535897932D1 * t95
      t131 = t11 * t130
      t132 = t56 ** 2
      t133 = t25 * t132
      t134 = -t96
      t135 = 0.1D1 / t134
      t136 = t39 * t135
      t138 = t38 * t24
      t140 = t135 * t24
      t141 = JacobiZeta(t38, t25)
      t142 = t88 * t38
      t144 = t91 * t142 + t141
      t145 = t144 * t140
      t147 = -t136 * t133 - t138 * t57 + t145 * t57
      t148 = t147 * t86
      t151 = -t148 * t16 * t131 - t56 * t128
      t160 = t11 * 0.31415926535897932D1
      t167 = t109 * 0.31415926535897932D1
      t168 = t22 * t83
      t169 = 0.31415926535897932D1 * t168
      t170 = t86 * t70
      t173 = 2 * t112 * t170 * t169
      t176 = t35 * t11
      t179 = t176 * t170 * t22 * t2 * t32 + t167 - t173
      t180 = t55 * t179
      t182 = t22 * t95
      t183 = 0.31415926535897932D1 * t182
      t186 = -t147 * t170 * t183 - t56 * t180
      t196 = t5 * 0.31415926535897932D1 * t72
      t200 = t53 * t16
      t205 = (k ** 2)
      t206 = t205 / t46 / t47 * RaC
      t207 = 0.1D1 / t1
      t208 = t151 ** 2
      t209 = t208 * t207
      t210 = 0.1D1 / t59
      t211 = t45 ** 2
      t212 = t211 * t210
      t216 = t81 ** 2
      t217 = 0.1D1 / t216
      t220 = t13 ** 2
      t221 = t105 * t220
      t222 = t86 ** 2
      t228 = t121 * t1 * t83
      t232 = 2*Dirac(t23) !signum(1, t23)
      t237 = t13 * t11
      t239 = t86 * t15
      t240 = t91 * t88
      t241 = 1 - t106 - t240
      t246 = 0.31415926535897932D1 * t25
      t247 = t237 * t246
      t248 = t88 ** 2
      t254 = t90 ** 2
      t256 = t95 * t15
      t264 = -t248 * t26 + t88 * (-t95 + 2 - t106) * t90 * t26 + (-t256& 
     &* t108 * t70 + t87 * t96 + (t26 * t96 + t87) * t107) * t254
      t266 = 0.1D1 / t254
      t272 = t24 * t90
      t273 = t24 * t88 - t272
      t284 = t95 * t266
      t285 = 0.31415926535897932D1 * t284
      t289 = t24 * t135 * t88 - t272
      t291 = t289 * t239 * t237
      t297 = t24 * t93
      t300 = t239 * t237
      t302 = t25 * t93
      t303 = t96 ** 2
      t304 = 0.1D1 / t303
      t310 = t24 * t26
      t315 = t83 * t109 * t304 * t103
      t317 = t86 * t16 * t160
      t321 = 0.1D1 / t108 / t107
      t335 = Dirac_derivative(t28) !Dirac(1, t28)
      t342 = t72 * t31
      t355 = t56 * t39
      t357 = t39 * t25
      t365 = t144 * t135 * t56 * t357 - t135 * t55 * t133 - t38 * t56 *& 
     &t357
      t369 = -t365 * t86 * t16 * t131 - t355 * t95 * t127
      t374 = t135 * t25
      t378 = t55 * t39
      t382 = -t56 * t132 * t374 - t144 * t378 * t140 + t378 * t138 + t56&
     & * t374
      t386 = -t382 * t86 * t16 * t131 + t55 * t39 * t127
      t393 = t121 * t1 * t95
      t399 = t25 * t56
      t404 = t136 * t95 * t132
      t407 = t134 ** 2
      t408 = 0.1D1 / t407
      t410 = t39 * t408 * t217 * t132
      t415 = t56 * t369
      t417 = t386 * t55
      t419 = t160 * t57
      t420 = t38 * t86
      t431 = t408 * t95
      t436 = t55 ** 2
      t437 = t436 - t240
      t452 = -t248 * t38 + t88 * (-t95 + t436 + 1) * t90 * t38 + (-t95 *&
     & t56 * t378 + t141 * t96 + (t38 * t96 + t141) * t436) * t254
      t459 = t95 * t38
      t461 = t273 * t86
      t470 = -2 * t386 * t136 * t399 + t317 * t404 + 2 * t317 * t410 - t&
     &151 * t135 * t133 - t138 * t415 - t138 * t417 - t420 * t16 * t419 &
     &- t127 * t24 * t57 + t145 * t415 + t145 * t417 + t144 * t135 * t86&
     & * t16 * t419 - 2 * t317 * t144 * t431 * t57 + (-t91 * t461 * t16 &
     &* t160 * t459 - t266 * t97 * t452 * t239 * t247 + t91 * t88 * t127&
     & + t291 * t285 * t142 + t437 * t127) * t140 * t57
      t480 = t21 * k * t79
      t482 = 0.1D1 / t60 / t59
      t484 = t39 * t45
      t488 = t205 * t79
      t489 = t45 * t210
      t497 = t186 ** 2
      t498 = t497 * t207
      t505 = t70 * 0.31415926535897932D1
      t506 = t86 * t505
      t510 = 2 * t506 * t22 * t105 * t83 - 2 * t505 * t256
      t513 = t22 ** 2
      t516 = t70 ** 2
      t517 = t222 * t516
      t526 = t1 * t513 * t83
      t527 = t232 * t516
      t531 = t167 - t173
      t534 = t505 * t22 * t25
      t536 = t266 * t97
      t547 = t289 * t170 * t69
      t554 = t170 * t69
      t562 = cos(t102)
      t590 = -t365 * t170 * t183 - t355 * t95 * t179
      t597 = -t382 * t170 * t183 + t55 * t39 * t179
      t617 = t56 * t590
      t619 = t597 * t55
      t652 = -2 * t597 * t136 * t399 + t554 * t404 + 2 * t554 * t410 - t&
     &186 * t135 * t133 - t138 * t617 - t138 * t619 - t420 * t505 * t22 &
     &* t57 - t179 * t24 * t57 + t145 * t617 + t145 * t619 + t144 * t135&
     & * t170 * t69 * t57 - 2 * t506 * t22 * t144 * t431 * t57 + (-t91 *&
     & t273 * t170 * t69 * t459 - t536 * t452 * t86 * t534 + t547 * t284&
     & * t142 + t91 * t88 * t179 + t437 * t179) * t140 * t57
      H = t67 * (-t17 * t16 * t8 * t5 * t2 + (0, 2) * t62 * t57 * t53 * &
     &t34 * t31 * 0.31415926535897932D1 * k * t48 * RaC) + t67 * (t53 * &
     &t23 / t74 * t73 + 2 * t62 * t151 * t21 * t80) * t53 * t70 * t69 - &
     &t67 * (-t53 * t70 * t13 * t8 * t73 + 2 * t62 * t186 * t21 * t80 + &
     &RaC - RaT) * t53 * t15 * t13 * t160 - t67 * (t200 / lambda / t74 *&
     & t196 + 8 * t212 * t209 * t206 + 2 * t62 * (-t56 * t55 * (6 * t112&
     & * t222 * t221 * t121 * t1 * t217 + 2 * t113 * t23 * t228 - 2 * t1&
     &12 * t232 * t221 * t228 - 2 * (-t97 * t81 * (-2 * t91 * t113 * t16&
     & * t160 * t83 * t88 - t91 * t26 * t273 * t239 * t237 * t130 - 2 * &
     &t241 * t112 * t239 * t237 * t84 - t266 * t97 * t264 * t239 * t247 &
     &+ t291 * t285 * t89) / 2 - t300 * 0.31415926535897932D1 * t97 * t2&
     &97 - t300 * 0.31415926535897932D1 * t304 * t302 + t114 * t11 * t24&
     &6 - t300 * 0.31415926535897932D1 * t310 + t317 * t315 / 2 - t86 * &
     &t237 * 0.31415926535897932D1 * t15 * t105 * t83 * t321 * t104 / 4)&
     & * t86 * t16 * t85 + (0, -2) * t36 * t1 * t335 + (0, 4) * t125 * t&
     &122 * t117 + (0, 1) * t124 * t23 * t8 * t342 + t35 * t232 * t221 *&
     & t8 * t72 * t32) - t56 * t369 * t127 - t386 * t128 + 2 * t147 * t2&
     &22 * t221 * t228 + t148 * t23 * t393 - t147 * t232 * t221 * t393 -&
     & t470 * t86 * t16 * t131) * t21 * t80 + 2 * t484 * t482 * t208 * t&
     &480 - 4 * t489 * t209 * t488) - t67 * (t200 * t8 * t196 + 8 * t212&
     & * t498 * t206 + 2 * t62 * (-t56 * t55 * (-t510 * t321 * 0.3141592&
     &6535897932D1 / 2 + 6 * t112 * t517 * t1 * t513 * t217 + 2 * t112 *&
     & t239 * t1 * t168 - 2 * t112 * t527 * t526 - 2 * (-t97 * t81 * (-t&
     &91 * t26 * t461 * t505 * t182 - t536 * t264 * t86 * t534 + t547 * &
     &t284 * t89 + t91 * t531 * t88 + t241 * t531) / 2 - t554 * t97 * t2&
     &97 - t554 * t304 * t302 - t81 * t531 / 2 - t506 * t22 * t310 + t10&
     &9 * t97 * t562 * 0.31415926535897932D1 / 2 + t554 * t315 / 2 - t51&
     &0 * t321 * t104 / 8) * t170 * t169 + (0, 1) * t176 * t239 * t22 * &
     &t342 + t176 * t527 * t513 * t72 * t32) - t56 * t590 * t179 - t597 &
     &* t180 + 2 * t147 * t517 * t526 + t147 * t239 * t1 * t182 - t147 *&
     & t527 * t1 * t513 * t95 - t652 * t170 * t183) * t21 * t80 + 2 * t4&
     &84 * t482 * t497 * t480 - 4 * t489 * t498 * t488)
end subroutine compute_H_domain
