module types
implicit none
integer, parameter      :: dp = selected_real_kind(P=15)
end module types

!--------------------------------------------------------------------------------
!       About this code
!--------------------------------------------------------------------------------
!
! This code computes Mie's light scattering solution by using the Debye series,
! instead of Lorenz-Mie series. 
! The basic structure of this code follows the BHMIE code by B.T.Draine:
!  https://www.astro.princeton.edu/~draine/scattering.html
! but it was modified to improve its accuracy and convergence by RT.
!
! The Debye series can "exactly" decompose the Lorenz-Mie series into
! various scattering components, such as diffraction, reflection, and 
! refraction (including (np+1) internal refraction). 
! See also Appendix A in Tazaki et al. (2021) for more detail.
!
! In this code, I adopted the algorithm developed by 
!        Shen & Wang (2010, ApOpt, 49, 2422), 
! where they developed a stable algorithm for computing highly absorbing
! sphere as well as non-absorbing sphere. 
!
!                                                  Ryo Tazaki (2020/June/19)
!
!--------------------------------------------------------------------------------
!       Input parameters
!--------------------------------------------------------------------------------
!
! x               : size parameter of a sphere
! refrel          : complex refractive index
! nang            : number of scattering angles from 0 to 90 deg.
! ical            : switch for the Debye series calculation
!                   ical = 0 : Series are summed upto the infinite term.
!                              Returned S1p, S2p contains values for 0 ~ np terms.
!                   ical = 1 : Series are terminated at np.
!                              Returned S1p,S2p contains values for 0 ~ np terms.
! np              : Terms to be computed in the debye series terms (np >= 0). 
!                       np=0 : Fraunhofer diffraction and reflected light.
!                       np=1 : Secondary refracted light 
!                       np=2 : Tertially refracted light (one internal reflection)
!                         ...
!
!--------------------------------------------------------------------------------
!       Output parameters
!--------------------------------------------------------------------------------
!
! S1(2*nang-1),   : Mie scattering amplitude at each scattering angles 
! S2(2*nang-1)  
! S1p(p,2*nang-1) : Mie scattering amplitude of p-th order scattered wave
! S2p(p,2*nang-1)   at each scattering angles. For example, p=0 contains
!                   the amplitude of Fraunhofer diffraction + surface 
!                   reflected light, p=1 does that of twice refracted light.
! Qext,Qsca,Qback : Extinction, scattering, back scattering efficiency
! gsca            : Asymmetry parameter    
!
!--------------------------------------------------------------------------------
subroutine debye_series_rt(x,refrel,nang,ical,np,S1,S2,S1p,S2p,Qext,Qsca,Qback,gsca)
use types
implicit none

!--------------------------------------------------------------------------------
!       Input variables
!--------------------------------------------------------------------------------
real(kind=dp)                               :: x       ! size parameter
complex(kind=dp)                            :: refrel  ! refractive index
integer                                     :: nang    ! number of angle mesh
integer                                     :: ical    ! Series termination switch
integer                                     :: np      ! Number of partial waves

!--------------------------------------------------------------------------------
!       Output variables
!--------------------------------------------------------------------------------

complex(kind=dp),dimension(1:2*nang-1)      :: S1,S2   ! Mie scattering amplitude
complex(kind=dp),dimension(0:np,1:2*nang-1) :: S1p,S2p ! Partial wave sca. amplitude 
real(kind=dp)                               :: qext    ! Qext
real(kind=dp)                               :: qsca    ! Qsca
real(kind=dp)                               :: gsca    ! <cos>
real(kind=dp)                               :: qback   ! Qback

!--------------------------------------------------------------------------------
!       Local variables
!--------------------------------------------------------------------------------
complex(kind=dp),parameter:: ivec=cmplx(0.0_dp,1.0_dp)
real(kind=dp),parameter   :: floorval=1.0e-30_dp
integer::j,jj,n,nn,nstop,nmx,ip,iq
real(kind=dp)::nu,d1x_ini
real(kind=dp)::xstop,ymod,en,fn,p,fnq_diff,dang,pii,theta
real(kind=dp)::realpart,imagpart
real(kind=dp),allocatable,dimension(:)::ajx
real(kind=dp),dimension(1:nang)::amu,pi0,pi1,pi,tau
complex(kind=dp)::y,d1y_ini,xi1,xi,dpcx
complex(kind=dp)::An13yy,An31yy,u0,u11,u31,u13,u33,Tn,R212n,R121n,umR212n,umR121n
complex(kind=dp)::an,bn,an1,bn1,fnq(2),fnqp(2,0:np),alpha(2),beta(2),R121n_power
complex(kind=dp),allocatable,dimension(:)::An13xx,lnAn13yy,w31,ajy
complex(kind=dp),allocatable,dimension(:,:)::d1,d3
realpart(dpcx)=(dble(dpcx))
imagpart(dpcx)=(dimag(dpcx))


! some checks...
y     = refrel*x
ymod  = abs(y)
xstop = x + 4.05_dp * x ** (1.0/3.0) + 2.0_dp
nstop = nint(xstop)
nmx   = nint(max(xstop,ymod)) + 15 

! check numbers of ray-expansion:
if(np .lt. 0) then
        print *, 'error : negative np is detected.'
        print *, '        np should be positive. STOP!'
        stop
endif

! check calculation model:
if(ical .ne. 0 .and. ical .ne. 1) then
        print *, 'error: incorrect calculation mode: ical = ',ical
        print *, '       ical should be either 0 or 1. STOP!'
        stop
endif

! check refractive index
if(real(refrel) .lt. 1.0_dp) then
        print *, 'Now, the real part of the refractive index is ',real(refrel)
        print *, ', which is less than unity. In this case, this code cannot' 
        print *, 'correctly decompose reflection/refraction components, although'
        print *, 'infinite series sum may be correctly calculated.'
        print *, 'Safety stop!'
        stop
endif

!*** require nang.ge.1 in order to calculate scattering intensities
pii  = 4.0_dp*atan(1.0_dp)
dang = 0.0_dp
if(nang.gt.1)dang=0.5_dp*pii/real(nang-1,kind=dp)
do j=1,nang
        theta=real(j-1,kind=dp)*dang
        amu(j)=cos(theta)
enddo
do j=1,nang
        pi0(j)=0.0_dp
        pi1(j)=1.0_dp
enddo
nn=2*nang-1
do j=1,nn
        S1(j)=cmplx(0.0_dp,0.0_dp)
        S2(j)=cmplx(0.0_dp,0.0_dp)
        do ip=0,np
                S1p(ip,j)=cmplx(0.0_dp,0.0_dp)
                S2p(ip,j)=cmplx(0.0_dp,0.0_dp)
        enddo
enddo

!----------------------------------------------------------------------
! The Lentz's continued fraction method
! Reference: W. J. Lentz, Appl. Opt. 15. 668. (1976)
!
! Finding initial value of the logarithmic derivative "d" at 
! the starting order of downward recurrence (n=nmx).
!----------------------------------------------------------------------
allocate(An13xx(0:nmx),lnAn13yy(0:nmx),d1(2,0:nmx),&
        d3(2,0:nmx),w31(0:nmx),ajx(1:nmx),ajy(1:nmx))

d1      = cmplx(0.0_dp,0.0_dp)
d3      = cmplx(0.0_dp,0.0_dp)
An13xx  = cmplx(0.0_dp,0.0_dp)
lnAn13yy= cmplx(0.0_dp,0.0_dp)

! Equation (5) in Lentz (1976):
nu = real(nmx,kind=dp)+0.5_dp
ajx(1) = 2.0_dp * nu / x
ajy(1) = 2.0_dp * nu / y
do j=2,nmx
        ajx(j) = (-1.0_dp) ** (j-1) * 2.0_dp*(nu+real(j-1,kind=dp)) / x
        ajy(j) = (-1.0_dp) ** (j-1) * 2.0_dp*(nu+real(j-1,kind=dp)) / y
enddo
d1x_ini = ajx(nmx)
d1y_ini = ajy(nmx)
do j=1,nmx-1
        d1x_ini = ajx(nmx-j) + 1.0_dp / d1x_ini
        d1y_ini = ajy(nmx-j) + 1.0_dp / d1y_ini
enddo
! Equation (3) in Lentz (1976):
d1(1,nmx)=-real(nmx,kind=dp)/x+d1x_ini
d1(2,nmx)=-real(nmx,kind=dp)/y+d1y_ini
! Usual downward recurrence:
do n=1,nmx
        en=real(nmx-n+1,kind=dp)
        d1(1,nmx-n) = (en/x) - (1.0_dp/(d1(1,nmx-n+1) + en/x))
        d1(2,nmx-n) = (en/y) - (1.0_dp/(d1(2,nmx-n+1) + en/y))
enddo
!calculate D(3)_n(x,or,mx) by upward recurrence.
d3(1,0)=ivec
d3(2,0)=ivec
w31 = cmplx(0.0_dp,0.0_dp)
w31(0) = -1.0/tan(y)+ivec
do n=1,nmx
        en=real(n,kind=dp)
        d3(1,n) = -(en/x) + (1.0_dp/(en/x - d3(1,n-1)))
        w31(n) = w31(n-1)/((d3(2,n-1)-en/y)*(d1(2,n-1)-en/y))
        d3(2,n) = d1(2,n) + w31(n) 
enddo
!calculate A(13)_n(x) by upward recurrence.
An13xx(0)=2.0_dp/(1.0_dp-ivec/tan(x))
do n=1,nmx
        en = real(n,kind=dp)
        An13xx(n) = An13xx(n-1)*(d1(1,n-1)-en/x)/(d3(1,n-1)-en/x)
enddo
!calculate lnA(13)_n(y) by upward recurrence.
lnAn13yy(0)=2.0_dp*aimag(y)+log(exp(-2.0_dp*aimag(y))-exp(-ivec*2.0_dp*real(y)))
do n=1,nmx
        en = real(n,kind=dp)
        lnAn13yy(n) = lnAn13yy(n-1)+log((d1(2,n-1)-en/y)/(d3(2,n-1)-en/y))
enddo

fnq_diff = 0.5_dp       ! Fraunhofer diffraction
alpha(1) = refrel
alpha(2) = 1.0_dp
beta(1)  = 1.0_dp
beta(2)  = refrel
qsca     = 0.0_dp
gsca     = 0.0_dp
p        =-1.0_dp  

!
! Start series summation.
!
! Comments on the termination order:
! Here, the series sum is truncated at the order 'nmx' rather than 'nstop'.
! This is because nstop is not enough to attain convergence for high order 
! scattered waves (large p). 
! The Debye series code written by P. Laven:
!       http://www.philiplaven.com/mieplot.htm
! seems to adopt nstop as the termination order, leading to 
! the unconverged results for high order of p.
! 
do n=1,nmx  

        en = real(n,kind=dp)
        fn = (2.0_dp*en+1.0_dp)/(en*(en+1.0_dp))

        !*** store previous values of an and bn for use
        !    in computation of g=<cos(theta)>
        if(n.gt.1)then
                an1=an
                bn1=bn
        endif

        !initialize
        fnq  = cmplx(0.0_dp,0.0_dp)
        fnqp = cmplx(0.0_dp,0.0_dp)
        
        do iq=1,2

                ! Shen & Wang (2010b) updated algorithm.
                u0 = (d1(1,n)-d3(1,n))*(d1(2,n)-d3(2,n))
                u11 = alpha(iq) * d1(1,n) - beta(iq) * d1(2,n)
                u31 = alpha(iq) * d3(1,n) - beta(iq) * d1(2,n)
                u13 = alpha(iq) * d1(1,n) - beta(iq) * d3(2,n)
                u33 = alpha(iq) * d3(1,n) - beta(iq) * d3(2,n)
               

                ! z=|z|e^{i*theta}: complex
                ! --> ln(z) = ln|z| + i * theta
                !     ln|z| = Re(ln(z))
                if(real(lnAn13yy(n)) .lt. 0.0_dp .or. aimag(refrel) .eq. 0.0_dp) then
                        An13yy = exp(lnAn13yy(n))
                        Tn = refrel * An13xx(n) * An13yy * u0
                        Tn = Tn / (u33-An13yy*u31)**2.0_dp
                        umR212n=An13xx(n)*(u13-An13yy*u11)/&
                                (u33-An13yy*u31)
                        umR121n=-An13yy*u31/(u33-An13yy*u31)
                        R212n=1.0_dp-An13xx(n)*(u13-An13yy*u11)/&
                                (u33-An13yy*u31)
                        R121n=1.0_dp+An13yy*u31/(u33-An13yy*u31)
                else
                        An31yy = exp(-lnAn13yy(n))
                        Tn = refrel*An13xx(n)*An31yy*u0
                        Tn = Tn / (An31yy*u33-u31)**2.0_dp
                        umR212n = An13xx(n)*(An31yy*u13-u11)/&
                                (An31yy*u33-u31)
                        umR121n = -u31/(An31yy*u33-u31)
                        R212n = 1.0_dp-An13xx(n)*(An31yy*u13-u11)/&
                                (An31yy*u33-u31)
                        R121n = 1.0_dp+u31/(An31yy*u33-u31)
                endif
                
                ! In this code, p=0 is defined such that 
                ! diffraction + Fresnel reflection
                fnqp(iq,0)  = fnq_diff * umR212n

                ! For an optically thick grain, R121n^{ip-1} or fnqp
                ! often causes underflows, so it is safe to set 
                ! some floor values for them.
                R121n_power = cmplx(1.0_dp,0.0_dp)
                do ip=1,np
                fnqp(iq,ip) = - 0.5_dp * Tn * R121n_power
                R121n_power = R121n_power * R121n
                if(abs(fnqp(iq,ip)) .le. floorval) fnqp(iq,ip) = cmplx(0.0_dp,0.0_dp)
                if(abs(R121n_power) .le. floorval) exit
                enddo

                if(ical .eq. 0) then
                        fnq(iq) = fnq_diff * (umR212n - Tn / umR121n)
                elseif(ical .eq. 1) then
                        do ip=0,np
                        fnq(iq) = fnq(iq) + fnqp(iq,ip)
                        enddo
                endif
        enddo

        an = fnq(1)
        bn = fnq(2)

        !*** augment sums for qsca and g=<cos(theta)>
        qsca=qsca+(2.0_dp*en+1.0_dp)*(abs(an)**2.0+abs(bn)**2.0)
        gsca=gsca+((2.0_dp*en+1.0_dp)/(en*(en+1.0_dp)))*&
        &        (realpart(an)*realpart(bn)+imagpart(an)*imagpart(bn))
        if(n.gt.1)then
            gsca=gsca+((en-1.0_dp)*(en+1.0_dp)/en)*&
        &      (realpart(an1)*realpart(an)+imagpart(an1)*imagpart(an)+&
        &      realpart(bn1)*realpart(bn)+imagpart(bn1)*imagpart(bn))
        endif

        !c*** now calculate scattering intensity pattern
        !    first do angles from 0 to 90
        do j=1,nang
                jj=2*nang-j
                pi(j)=pi1(j)
                tau(j)=en*amu(j)*pi(j)-(en+1.0_dp)*pi0(j)
                S1(j)=S1(j)+fn*(an*pi(j)+bn*tau(j))
                S2(j)=S2(j)+fn*(an*tau(j)+bn*pi(j))
                !
                if(ical .ge. 0) then
                        do ip=0,np
                                S1p(ip,j) = S1p(ip,j) + fn*(fnqp(1,ip)*pi(j)+fnqp(2,ip)*tau(j))
                                S2p(ip,j) = S2p(ip,j) + fn*(fnqp(1,ip)*tau(j)+fnqp(2,ip)*pi(j))
                        enddo
                endif
        enddo

        !*** now do angles greater than 90 using pi and tau from
        !    angles less than 90.
        !    p=1 for n=1,3,...; p=-1 for n=2,4,...
        p=-p
        do j=1,nang-1
                jj=2*nang-j
                S1(jj)=S1(jj)+fn*p*(an*pi(j)-bn*tau(j))
                S2(jj)=S2(jj)+fn*p*(bn*pi(j)-an*tau(j))
                !
                if(ical .ge. 0) then
                do ip=0,np
                S1p(ip,jj) = S1p(ip,jj) + fn*p*(fnqp(1,ip)*pi(j)-fnqp(2,ip)*tau(j))
                S2p(ip,jj) = S2p(ip,jj) + fn*p*(fnqp(2,ip)*pi(j)-fnqp(1,ip)*tau(j))
                enddo
                endif
        enddo

        !*** compute pi_n for next value of n
        !    for each angle j, compute pi_n+1
        !    from pi = pi_n , pi0 = pi_n-1
        do j=1,nang
                pi1(j)=((2.0_dp*en+1.0_dp)*amu(j)*pi(j)-(en+1.0_dp)*pi0(j))/en
                pi0(j)=pi(j)
        enddo
enddo

!*** have summed sufficient terms.
!    now compute qsca,qext,qback,and gsca
gsca=2.0_dp*gsca/qsca
qsca=(2.0_dp/(x*x))*qsca
qext=(4.0_dp/(x*x))*realpart(S1(1))
qback=(4.0_dp*(ABS(S1(2*NANG-1))/x)**2.0_dp) 

deallocate(An13xx,lnAn13yy,d1,d3,w31,ajx,ajy)

return
end subroutine debye_series_rt
