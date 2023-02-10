!!!!SJT: This is an edited version of the file xgscd.txt 
!!!!SJT: deleted original driver program to allow linkage of subroutines to other code
!!!!SJT: deleted sample data at end of file so that the file could be compiled
!!!!SJT: converted to free form -- comment character "!" used throughout -- "&" used at left and right sides for line continuation
!!!!SJT: added implicit none statements to each routine
!!!!SJT: disabled save statement in variable declaration for thread safety
!!!!SJT: modified computations of m from mc (quad precision) in an attempt to reduce truncation error when mc is close to unity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gscd(u,mc_qp,s,c,d) !!SJT
!subroutine gscd(u,mc,s,c,d) !!SJT: original
!
!	Double precision subroutine to compute three Jacobian elliptic functions simultaneously
!
!   For general argument: -infty < u < infty
!
!     Reference: T. Fukushima, (2012) Numer. Math. DOI 10.1007/s00211-012-0498-0
!       "Precise and Fast Computation of Jacobian Elliptic Functions by
!        Conditional Duplication"
!
!     Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
!
!     Used subprograms: scd2, elk
!
!     Inputs: u = argument, mc = 1-m, 0 < mc <= 1
!
!     Output: s = sn(u|m), c=cn(u|m), d=dn(u|m)
!
implicit none !!SJT
real*16 :: mc_qp !!SJT
real*8 u,mc,s,c,d
real*8 m,kc,ux,k,kh,kh3,kh5,kh7,k2,k3,k4,sx,cx,dx
real*8 elk
!
!m=1.d0-mc !!SJT: original
m=real(1.q0-mc_qp,8); mc=real(mc_qp,8) !!SJT: quad precision used to reduce the impact of cancellation errors
kc=sqrt(mc)
ux=abs(u)
if(ux.lt.0.785d0) then
    !call scd2(ux,mc,s,c,d) !!SJT: OG
    call scd2(ux,mc_qp,s,c,d) !!SJT
else
    !k=elk(mc) !!SJT: OG
    k=elk(mc_qp) !!SJT
    kh=k*0.5d0; kh3=k*1.5d0; kh5=k*2.5d0; kh7=k*3.5d0;
    k2=k*2.d0; k3=k*3.d0; k4=k*4.d0
    ux=ux-k4*dble(int(ux/k4))
    if(ux.lt.kh) then
        !call scd2(ux,mc,s,c,d) !!SJT: OG
        call scd2(ux,mc_qp,s,c,d) !!SJT
    elseif(ux.lt.k) then
        ux=k-ux
        !call scd2(ux,mc,s,c,d) !!SJT: OG
        call scd2(ux,mc_qp,s,c,d) !!SJT
        sx=c/d; c=kc*s/d; s=sx; d=kc/d
    elseif(ux.lt.kh3) then
        ux=ux-k
        !call scd2(ux,mc,s,c,d) !!SJT: OG
        call scd2(ux,mc_qp,s,c,d) !!SJT
        sx=c/d; c=-kc*s/d; s=sx; d=kc/d
    elseif(ux.lt.k2) then
        ux=k2-ux
        !call scd2(ux,mc,s,c,d) !!SJT: OG
        call scd2(ux,mc_qp,s,c,d) !!SJT
        c=-c
    elseif(ux.lt.kh5) then
        ux=ux-k2
        !call scd2(ux,mc,s,c,d) !!SJT: OG
        call scd2(ux,mc_qp,s,c,d) !!SJT
        s=-s; c=-c
    elseif(ux.lt.k3) then
        ux=k3-ux
        !call scd2(ux,mc,s,c,d) !!SJT: OG
        call scd2(ux,mc_qp,s,c,d) !!SJT
        sx=-c/d; c=-kc*s/d; s=sx; d=kc/d
    elseif(ux.lt.kh7) then
        ux=ux-k3
        !call scd2(ux,mc,s,c,d) !!SJT: OG
        call scd2(ux,mc_qp,s,c,d) !!SJT
        sx=-c/d; c=kc*s/d; s=sx; d=kc/d
    else
        ux=k4-ux
        !call scd2(ux,mc,s,c,d) !!SJT: OG
        call scd2(ux,mc_qp,s,c,d) !!SJT
        s=-s
    endif
endif
if(u.lt.0.d0) s=-s
return
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rgscd(u,mc,s,c,d)
!
!	Single precision subroutine to compute three Jacobian elliptic functions simultaneously
!
!   For general argument: -infty < u < infty
!
!     Reference: T. Fukushima, (2012) Numer. Math. DOI 10.1007/s00211-012-0498-0
!       "Precise and Fast Computation of Jacobian Elliptic Functions by
!        Conditional Duplication"
!
!     Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
!
!     Used subprograms: rscd2, relk
!
!     Inputs: u = argument, mc = 1-m, 0 < mc <= 1
!
!     Output: s = sn(u|m), c=cn(u|m), d=dn(u|m)
!
implicit none !!SJT
real*4 u,mc,s,c,d
real*4 m,kc,ux,k,kh,kh3,k2,kh5,k3,kh7,k4,sx,cx,dx
real*4 relk
!
m=1.0-mc
kc=sqrt(mc)
ux=abs(u)
if(ux.lt.0.785) then
    call rscd2(ux,mc,s,c,d)
else
    k=relk(mc)
    kh=k*0.5; kh3=k*1.5; kh5=k*2.5; kh7=k*3.5;
    k2=k*2.0; k3=k*3.0; k4=k*4.0
    ux=ux-k4*int(ux/k4)
    if(ux.lt.kh) then
        call rscd2(ux,mc,s,c,d)
    elseif(ux.lt.k) then
        ux=k-ux
        call rscd2(ux,mc,s,c,d)
        sx=c/d; c=kc*s/d; s=sx; d=kc/d
    elseif(ux.lt.kh3) then
        ux=ux-k
        call rscd2(ux,mc,s,c,d)
        sx=c/d; c=-kc*s/d; s=sx; d=kc/d
    elseif(ux.lt.k2) then
        ux=k2-ux
        call rscd2(ux,mc,s,c,d)
        c=-c
    elseif(ux.lt.kh5) then
        ux=ux-k2
        call rscd2(ux,mc,s,c,d)
        s=-s; c=-c
    elseif(ux.lt.k3) then
        ux=k3-ux
        call rscd2(ux,mc,s,c,d)
        sx=-c/d; c=-kc*s/d; s=sx; d=kc/d
    elseif(ux.lt.kh7) then
        ux=ux-k3
        call rscd2(ux,mc,s,c,d)
        sx=-c/d; c=kc*s/d; s=sx; d=kc/d
    else
        ux=k4-ux
        call rscd2(ux,mc,s,c,d)
        s=-s
    endif
endif
if(u.lt.0.0) s=-s
return
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine scd2(u,mc_qp,s,c,d) !!SJT
!    subroutine scd2(u,mc,s,c,d) !!SJT: original
!
!	Double precision subroutine to compute three Jacobian elliptic functions simultaneously
!
!   For limited argument: 0 <= u < K/2
!
!     Reference: T. Fukushima, (2012) Numer. Math. DOI 10.1007/s00211-012-0498-0
!       "Precise and Fast Computation of Jacobian Elliptic Functions by
!        Conditional Duplication"
!
!     Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
!
!     Inputs: u = argument, mc = 1-m, 0 < mc <= 1
!
!     Output: s = sn(u|m), c=cn(u|m), d=dn(u|m)
!
    implicit none !!SJT
    real*8 :: x !!SJT
    real*8 u,mc,s,c,d
    real*8 B10,B11,B20,B21,B22,m,uA,uT,u0,v,a,b,y,z,my,mc2,m2,xz,w
    integer n,j,i
    parameter (B10=1.d0/24.d0,B11=1.d0/6.d0,B20=1.d0/720.d0)
    parameter (B21=11.d0/180.d0,B22=1.d0/45.d0)
    real*16 :: mc_qp !!SJT
    !m=1.d0-mc; uA=1.76269d0+mc*1.16357d0; uT=5.217d-3-m*2.143d-3; u0=u !!SJT: original
    m=real(1.q0-mc_qp,8); mc=real(mc_qp,8) !!SJT
    uA=1.76269d0+mc*1.16357d0; uT=5.217d-3-m*2.143d-3; u0=u !!SJT
    do n=0,20
        if(u0.lt.uT) goto 1
        u0=u0*0.5d0
    enddo
    write(*,*) "(scd2) Too large input argument: u=", u
1 continue
    v=u0*u0; a=1.d0; b=v*(0.5d0-v*(B10+m*B11-v*(B20+m*(B21+m*B22))))
    if(u.lt.uA) then
        do j=1,n
            y=b*(a*2.d0-b); z=a*a; my=m*y; b=(y*2.d0)*(z-my); a=z*z-my*y
        enddo
    else
        do j=1,n
            y=b*(a*2.d0-b); z=a*a; my=m*y
            if(z.lt.my*2.d0) goto 2
            b=(y*2.d0)*(z-my); a=z*z-my*y
        enddo
    endif
    b=b/a; y=b*(2.d0-b); c=1.d0-b; s=sqrt(y); d=sqrt(1.d0-m*y)
     return
2 continue
    c=a-b; mc2=mc*2.d0; m2=m*2.d0
    do i=j,n
        x=c*c; z=a*a; w=m*x*x-mc*z*z; xz=x*z; c=mc2*xz+w; a=m2*xz-w
    enddo
    c=c/a; x=c*c; s=sqrt(1.d0-x); d=sqrt(mc+m*x)
    return; end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rscd2(u,mc,s,c,d)
!
!	Single precision subroutine to compute three Jacobian elliptic functions simultaneously
!
!   For limited argument: 0 <= u < K/2
!
!     Reference: T. Fukushima, (2012) Numer. Math. DOI 10.1007/s00211-012-0498-0
!       "Precise and Fast Computation of Jacobian Elliptic Functions by
!        Conditional Duplication"
!
!     Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
!
!     Inputs: u = argument, mc = 1-m, 0 < mc <= 1
!
!     Output: s = sn(u|m), c=cn(u|m), d=dn(u|m)
!
    implicit none !!SJT
    real*4 :: x !!SJT
    real*4 u,mc,s,c,d
    real*4 m,uA,uT,u0,v,a,b,y,z,my,mc2,m2,xz,w
    integer n,j,i
    m=1.0-mc; uA=1.76269+mc*1.16357; uT=9.207e-4-m*4.604e-4; u0=u
    do n=0,20
        if(u0.lt.uT) goto 1
        u0=u0*0.5
    enddo
    write(*,*) "(rscd2) Too large input argument: u=", u
1 continue
    v=u0*u0; a=1.0; b=v*0.5
    if(u.lt.uA) then
        do j=1,n
            y=b*(a*2.0-b); z=a*a; my=m*y; b=(y*2.0)*(z-my); a=z*z-my*y
        enddo
    else
        do j=1,n
            y=b*(a*2.0-b); z=a*a; my=m*y
            if(z.lt.my*2.0) goto 2
            b=(y*2.0)*(z-my); a=z*z-my*y
        enddo
    endif
    b=b/a; y=b*(2.0-b); c=1.0-b; s=sqrt(y); d=sqrt(1.0-m*y); return
2 continue
    c=a-b; mc2=mc*2.0; m2=m*2.0
    do i=j,n
        x=c*c; z=a*a; w=m*x*x-mc*z*z; xz=x*z; c=mc2*xz+w; a=m2*xz-w
    enddo
    c=c/a; x=c*c; s=sqrt(1.0-x); d=sqrt(mc+m*x); return
    end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real*8 function elk(mc_qp) !!SJT
!      real*8 function elk(mc) !!SJT: original
!c
!c	Double precision complete elliptic integral of the first kind
!c
!c     Reference: T. Fukushima, (2009) Celest. Mech. Dyn. Astron. 105, 305-328
!c        "Fast Computation of Complete Elliptic Integrlals and Jacobian
!c         Elliptic Functions"
!c
!c     Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
!c
!c     Inputs: mc   = complementary parameter 0 <= mc   <= 1
!c
!c     Output: elk
!c
        implicit none !!SJT
	real*8 mc
	real*8 mcold,PIHALF,PIINV,elkold,TINY,m,mx,P,Q
	real*8 kkc,nome
        !real*8 :: kc !!SJT
        real*16 :: mc_qp !!SJT
!c
	real*8 D1,D2,D3,D4,D5,D6,D7,D8,D9,D10,D11,D12,D13,D14
	parameter (D1=1.d0/16.d0,D2=1.d0/32.d0,D3=21.d0/1024.d0)
	parameter (D4=31.d0/2048.d0,D5=6257.d0/524288.d0)
	parameter (D6=10293.d0/1048576.d0,D7=279025.d0/33554432.d0)
	parameter (D8=483127.d0/67108864.d0)
	parameter (D9=435506703.d0/68719476736.d0)
	parameter (D10=776957575.d0/137438953472.d0)
	parameter (D11=22417045555.d0/4398046511104.d0)
	parameter (D12=40784671953.d0/8796093022208.d0)
	parameter (D13=9569130097211.d0/2251799813685248.d0)
	parameter (D14=17652604545791.d0/4503599627370496.d0)
!c
	!logical first/.TRUE./ !!SJT: note that logical variable first is .true. at the start of each routine call
	!save first,mcold,PIHALF,PIINV,elkold,TINY !!SJT: disabled for thread safety
        logical first !!SJT: initialization cannot be in declaration statement for thread safety purposes
        first=.true. !!SJT: initialization outside of declaration for thread safety
!c
	if(first) then
		first=.FALSE.
		mcold=1.d0
		PIHALF=atan(1.d0)*2.d0
		PIINV=0.5d0/PIHALF
		elkold=PIHALF
		TINY=1.d-99
	endif
	!m=1.d0-mc !!SJT: original
        m=real(1.q0-mc_qp,8); mc=real(mc_qp,8) !!SJT
        !kc=sqrt(mc); m=(1.d0+kc)*(1.d0-kc) !!!!SJT: reduce truncation error
	if(abs(m).lt.1.d-16) then
		elk=PIHALF
	elseif(abs(mc-mcold).lt.1.11d-16*mc) then
		elk=elkold
	elseif(mc.lt.TINY) then
	  elk=1.3862943611198906d0-0.5d0*log(TINY)
	elseif(mc.lt.1.11d-16) then
	  elk=1.3862943611198906d0-0.5d0*log(mc)
	elseif(mc.lt.0.1d0) then
		nome=mc*(D1+mc*(D2+mc*(D3+mc*(D4+mc*(D5+mc*(D6&
     &		+mc*(D7+mc*(D8+mc*(D9+mc*(D10+mc*(D11+mc*(D12&
     &		+mc*(D13+mc*D14)))))))))))))
		mx=mc-0.05d0
!c
!c	K'
!c
		kkc=1.591003453790792180d0+mx*(&
     &		0.416000743991786912d0+mx*(&
     &		0.245791514264103415d0+mx*(&
     &		0.179481482914906162d0+mx*(&
     &		0.144556057087555150d0+mx*(&
     &		0.123200993312427711d0+mx*(&
     &		0.108938811574293531d0+mx*(&
     &		0.098853409871592910d0+mx*(&
     &		0.091439629201749751d0+mx*(&
     &		0.085842591595413900d0+mx*(&
     &		0.081541118718303215d0))))))))))
!c
		elk=-kkc*PIINV*log(nome)
	elseif(m.le.0.1d0) then
		mx=m-0.05d0
		elk=1.591003453790792180d0+mx*(&
     &		0.416000743991786912d0+mx*(&
     &		0.245791514264103415d0+mx*(&
     &		0.179481482914906162d0+mx*(&
     &		0.144556057087555150d0+mx*(&
     &		0.123200993312427711d0+mx*(&
     &		0.108938811574293531d0+mx*(&
     &		0.098853409871592910d0+mx*(&
     &		0.091439629201749751d0+mx*(&
     &		0.085842591595413900d0+mx*(&
     &		0.081541118718303215d0))))))))))
	elseif(m.le.0.2d0) then
		mx=m-0.15d0
		elk=1.635256732264579992d0+mx*(&
     &		0.471190626148732291d0+mx*(&
     &		0.309728410831499587d0+mx*(&
     &		0.252208311773135699d0+mx*(&
     &		0.226725623219684650d0+mx*(&
     &		0.215774446729585976d0+mx*(&
     &		0.213108771877348910d0+mx*(&
     &		0.216029124605188282d0+mx*(&
     &		0.223255831633057896d0+mx*(&
     &		0.234180501294209925d0+mx*(&
     &		0.248557682972264071d0+mx*(&
     &		0.266363809892617521d0+mx*(&
     &		0.287728452156114668d0))))))))))))
	elseif(m.le.0.3d0) then
		mx=m-0.25d0
		elk=1.685750354812596043d0+mx*(&
     &		0.541731848613280329d0+mx*(&
     &		0.401524438390690257d0+mx*(&
     &		0.369642473420889090d0+mx*(&
     &		0.376060715354583645d0+mx*(&
     &		0.405235887085125919d0+mx*(&
     &		0.453294381753999079d0+mx*(&
     &		0.520518947651184205d0+mx*(&
     &		0.609426039204995055d0+mx*(&
     &		0.724263522282908870d0+mx*(&
     &		0.871013847709812357d0+mx*(&
     &		1.057652872753547036d0)))))))))))
	elseif(m.le.0.4d0) then
		mx=m-0.35d0
		elk=1.744350597225613243d0+mx*(&
     &		0.634864275371935304d0+mx*(&
     &		0.539842564164445538d0+mx*(&
     &		0.571892705193787391d0+mx*(&
     &		0.670295136265406100d0+mx*(&
     &		0.832586590010977199d0+mx*(&
     &		1.073857448247933265d0+mx*(&
     &		1.422091460675497751d0+mx*(&
     &		1.920387183402304829d0+mx*(&
     &		2.632552548331654201d0+mx*(&
     &		3.652109747319039160d0+mx*(&
     &		5.115867135558865806d0+mx*(&
     &		7.224080007363877411d0))))))))))))
	elseif(m.le.0.5d0) then
		mx=m-0.45d0
		elk=1.813883936816982644d0+mx*(&
     &		0.763163245700557246d0+mx*(&
     &		0.761928605321595831d0+mx*(&
     &		0.951074653668427927d0+mx*(&
     &		1.315180671703161215d0+mx*(&
     &		1.928560693477410941d0+mx*(&
     &		2.937509342531378755d0+mx*(&
     &		4.594894405442878062d0+mx*(&
     &		7.330071221881720772d0+mx*(&
     &		11.87151259742530180d0+mx*(&
     &		19.45851374822937738d0+mx*(&
     &		32.20638657246426863d0+mx*(&
     &		53.73749198700554656d0+mx*(&
     &		90.27388602940998849d0)))))))))))))
	elseif(m.le.0.6d0) then
		mx=m-0.55d0
		elk=1.898924910271553526d0+mx*(&
     &		0.950521794618244435d0+mx*(&
     &		1.151077589959015808d0+mx*(&
     &		1.750239106986300540d0+mx*(&
     &		2.952676812636875180d0+mx*(&
     &		5.285800396121450889d0+mx*(&
     &		9.832485716659979747d0+mx*(&
     &		18.78714868327559562d0+mx*(&
     &		36.61468615273698145d0+mx*(&
     &		72.45292395127771801d0+mx*(&
     &		145.1079577347069102d0+mx*(&
     &		293.4786396308497026d0+mx*(&
     &		598.3851815055010179d0+mx*(&
     &		1228.420013075863451d0+mx*(&
     &		2536.529755382764488d0))))))))))))))
	elseif(m.le.0.7d0) then
		mx=m-0.65d0
		elk=2.007598398424376302d0+mx*(&
     &		1.248457231212347337d0+mx*(&
     &		1.926234657076479729d0+mx*(&
     &		3.751289640087587680d0+mx*(&
     &		8.119944554932045802d0+mx*(&
     &		18.66572130873555361d0+mx*(&
     &		44.60392484291437063d0+mx*(&
     &		109.5092054309498377d0+mx*(&
     &		274.2779548232413480d0+mx*(&
     &		697.5598008606326163d0+mx*(&
     &		1795.716014500247129d0+mx*(&
     &		4668.381716790389910d0+mx*(&
     &		12235.76246813664335d0+mx*(&
     &		32290.17809718320818d0+mx*(&
     &		85713.07608195964685d0+mx*(&
     &		228672.1890493117096d0+mx*(&
     &		612757.2711915852774d0))))))))))))))))
	elseif(m.le.0.8d0) then
		mx=m-0.75d0
		elk=2.156515647499643235d0+mx*(&
     &		1.791805641849463243d0+mx*(&
     &		3.826751287465713147d0+mx*(&
     &		10.38672468363797208d0+mx*(&
     &		31.40331405468070290d0+mx*(&
     &		100.9237039498695416d0+mx*(&
     &		337.3268282632272897d0+mx*(&
     &		1158.707930567827917d0+mx*(&
     &		4060.990742193632092d0+mx*(&
     &		14454.00184034344795d0+mx*(&
     &		52076.66107599404803d0+mx*(&
     &		189493.6591462156887d0+mx*(&
     &		695184.5762413896145d0+mx*(&
     &		2.567994048255284686d6+mx*(&
     &		9.541921966748386322d6+mx*(&
     &		3.563492744218076174d7+mx*(&
     &		1.336692984612040871d8+mx*(&
     &		5.033521866866284541d8+mx*(&
     &		1.901975729538660119d9+mx*(&
     &		7.208915015330103756d9)))))))))))))))))))
	elseif(m.le.0.85d0) then
		mx=m-0.825d0
		elk=2.318122621712510589d0+mx*(&
     &		2.616920150291232841d0+mx*(&
     &		7.897935075731355823d0+mx*(&
     &		30.50239715446672327d0+mx*(&
     &		131.4869365523528456d0+mx*(&
     &		602.9847637356491617d0+mx*(&
     &		2877.024617809972641d0+mx*(&
     &		14110.51991915180325d0+mx*(&
     &		70621.44088156540229d0+mx*(&
     &		358977.2665825309926d0+mx*(&
     &		1.847238263723971684d6+mx*(&
     &		9.600515416049214109d6+mx*(&
     &		5.030767708502366879d7+mx*(&
     &		2.654441886527127967d8+mx*(&
     &		1.408862325028702687d9+mx*(&
     &		7.515687935373774627d9)))))))))))))))
      else
		mx=m-0.875d0
		elk=2.473596173751343912d0+mx*(&
     &		3.727624244118099310d0+mx*(&
     &		15.60739303554930496d0+mx*(&
     &		84.12850842805887747d0+mx*(&
     &		506.9818197040613935d0+mx*(&
     &		3252.277058145123644d0+mx*(&
     &		21713.24241957434256d0+mx*(&
     &		149037.0451890932766d0+mx*(&
     &		1.043999331089990839d6+mx*(&
     &		7.427974817042038995d6+mx*(&
     &		5.350383967558661151d7+mx*(&
     &		3.892498869948708474d8+mx*(&
     &		2.855288351100810619d9+mx*(&
     &		2.109007703876684053d10+mx*(&
     &		1.566998339477902014d11+mx*(&
     &		1.170222242422439893d12+mx*(&
     &		8.777948323668937971d12+mx*(&
     &		6.610124275248495041d13+mx*(&
     &		4.994880537133887989d14+mx*(&
     &		3.785974339724029920d15)))))))))))))))))))
	endif
!c
	mcold=mc
	elkold=elk
      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real function relk(mc)
!c
!c	Single precision complete elliptic integral of the first kind
!c
!c     Reference: T. Fukushima, (2009) Celest. Mech. Dyn. Astron. 105, 305-328
!c        "Fast Computation of Complete Elliptic Integrlals and Jacobian
!c         Elliptic Functions"
!c
!c     Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
!c
!c     Inputs: mc   = complementary parameter 0 <= mc   <= 1
!c
!c     Output: elk
!c
        implicit none !!SJT
	real mc
	real m,nome,mx,kkc,P,Q
!c
	real D1,D2,D3,D4,D5,D6
	parameter (D1=1.0/16.0,D2=1.0/32.0,D3=21.0/1024.0)
	parameter (D4=31.0/2048.0,D5=6257.0/524288.0)
	parameter (D6=10293.0/1048576.0)
	real PIHALF,PIINV
	parameter (PIHALF=1.57079633,PIINV=0.318309886)
!c
	m=1.0-mc
	if(mc.lt.1.05e-8) then
	    relk=1.38629436-0.5*log(mc)
	elseif(mc.lt.0.1) then
		nome=mc*(D1+mc*(D2+mc*(D3+mc*(D4+mc*(D5+mc*D6)))))
		mx=mc-0.05
!c
!c	K'
!c
		kkc=1.59100345+mx*(&
     &		0.41600074+mx*(&
     &		0.24579151+mx*(&
     &		0.17948148+mx*(&
     &		0.14455606))))
!c
		relk=-kkc*PIINV*log(nome)
	elseif(m.le.0.1) then
		mx=m-0.05
		relk=1.59100345+mx*(&
     &		0.41600074+mx*(&
     &		0.24579151+mx*(&
     &		0.17948148+mx*(&
     &		0.14455606))))
	elseif(m.le.0.2) then
		mx=m-0.15
		relk=1.63525673+mx*(&
     &		0.47119063+mx*(&
     &		0.30972841+mx*(&
     &		0.25220831+mx*(&
     &		0.22672562))))
	elseif(m.le.0.3) then
		mx=m-0.25
		relk=1.68575035+mx*(&
     &		0.54173185+mx*(&
     &		0.40152444+mx*(&
     &		0.36964247+mx*(&
     &		0.37606072))))
	elseif(m.le.0.4) then
		mx=m-0.35
		relk=1.74435060+mx*(&
     &		0.63486428+mx*(&
     &		0.53984256+mx*(&
     &		0.57189271+mx*(&
     &		0.67029514+mx*(&
     &		0.83258659)))))
	elseif(m.le.0.5) then
		mx=m-0.45
		relk=1.81388394+mx*(&
     &		0.76316325+mx*(&
     &		0.76192861+mx*(&
     &		0.95107465+mx*(&
     &		1.31518067+mx*(&
     &		1.92856069)))))
	elseif(m.le.0.6) then
		mx=m-0.55
		relk=1.89892491+mx*(&
     &		0.95052179+mx*(&
     &		1.15107759+mx*(&
     &		1.75023911+mx*(&
     &		2.95267681+mx*(&
     &		5.28580040)))))
	elseif(m.le.0.7) then
		mx=m-0.65
		relk=2.00759840+mx*(&
     &		1.24845723+mx*(&
     &		1.92623466+mx*(&
     &		3.75128964+mx*(&
     &		8.11994455+mx*(&
     &		18.6657213+mx*(&
     &		44.6039248))))))
	elseif(m.le.0.8) then
		mx=m-0.75
		relk=2.15651565+mx*(&
     &		1.79180564+mx*(&
     &		3.82675129+mx*(&
     &		10.3867247+mx*(&
     &		31.4033141+mx*(&
     &		100.923704+mx*(&
     &		337.326828+mx*(&
     &		1158.70793)))))))
	elseif(m.le.0.85) then
		mx=m-0.825
		relk=2.31812262+mx*(&
     &		2.61692015+mx*(&
     &		7.89793508+mx*(&
     &		30.5023972+mx*(&
     &		131.486937+mx*(&
     &		602.984764+mx*(&
     &		2877.02462))))))
      else
		mx=m-0.875
		relk=2.47359617+mx*(&
     &		3.72762424+mx*(&
     &		15.6073930+mx*(&
     &		84.1285084+mx*(&
     &		506.981820+mx*(&
     &		3252.27706+mx*(&
     &		21713.2424+mx*(&
     &		149037.045)))))))
	endif
!c
      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
