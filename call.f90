!----------------------------------------------------------------------
!       call debye_series.f90
!----------------------------------------------------------------------
program call_debye
use types
implicit none
integer,parameter:: npmax = 100
integer,parameter:: nang  = 181
complex(kind=dp)::refrel
integer:: ical,ip,np,j
real(kind=dp)::x,ext,sca,back,gs,dang,refre,refim,tauabs
real(kind=dp),dimension(2*nang-1)::ang,S11,S12,S33,S34
real(kind=dp),dimension(0:npmax,2*nang-1)::S11p,S12p,S33p,S34p
complex(kind=dp),dimension(2*nang-1)::S1,S2
complex(kind=dp),dimension(0:npmax,2*nang-1)::S1p,S2p
!----------------------------------------------------------------------
!  INPUT PARAMETER
!----------------------------------------------------------------------
x       = 1.0e2_dp      ! size parameter     
refre   = 1.330_dp      ! real part of refractive index
refim   = 1.d-1         ! imag part of refractive index
ical    = 0             ! Switch for debye series truncation order
np      = 5             ! Maximum order of the debye series
!----------------------------------------------------------------------
!  call debye_series_rt
!----------------------------------------------------------------------
refrel = cmplx(refre,refim)
dang   = 90.0_dp/real(nang-1,kind=dp)
call debye_series_rt(x,refrel,nang,ical,np,s1,s2,s1p,s2p,ext,sca,back,gs)
do j=1,2*nang-1
        ang(j)=dang*dble(j-1)
        S11(j)=0.5D0*ABS(S2(J))*ABS(S2(J))
        S11(j)=S11(j)+0.5D0*ABS(S1(J))*ABS(S1(J))
        S12(j)=0.5D0*ABS(S2(J))*ABS(S2(J))
        S12(j)=S12(j)-0.5D0*ABS(S1(J))*ABS(S1(J))       
        S33(j)=REAL(S1(J)*CONJG(S2(J)))
        S34(j)=AIMAG(S2(J)*CONJG(S1(J)))
        do ip=0,np
                S11p(ip,j)=0.5D0*ABS(S2p(ip,J))*ABS(S2p(ip,J))
                S11p(ip,j)=S11p(ip,j)+0.5D0*ABS(S1p(ip,J))*ABS(S1p(ip,J))
                S12p(ip,j)=0.5D0*ABS(S2p(ip,J))*ABS(S2p(ip,J))
                S12p(ip,j)=S12p(ip,j)-0.5D0*ABS(S1p(ip,J))*ABS(S1p(ip,J))       
                S33p(ip,j)=REAL(S1p(ip,J)*CONJG(S2p(ip,J)))
                S34p(ip,j)=AIMAG(S2p(ip,J)*CONJG(S1p(ip,J)))
       enddo
enddo
!----------------------------------------------------------------------
!  OUTPUT
!----------------------------------------------------------------------
write(*,*) "writing... output.dat"
open(10,file="output.dat",status="unknown")
write(10,2100) x," = size parameter"
write(10,2100) refre," = Re(m)"
write(10,2100) refim," = Im(m)"
write(10,2100) ext," = Qext"
write(10,2100) sca," = Qsca"
write(10,2100) back," = Qback"
write(10,2100) gs," = g"
if(ical .eq. 0) then
        write(10,2301) "Infinite"," = Truncation order of the Debye series"
elseif(ical .eq. 1) then
        write(10,2300) np," = Truncation order of the Debye series"
endif
write(10,*)
write(10,2000) "ang (deg)","S11","S11(p=0)","S11(p=1)","..."
do j=1,2*nang-1
        write(10,1000) ang(j),S11(j),(S11p(ip,j),ip=0,min(np,5))
enddo
write(10,*)
write(10,2000) "ang (deg)","P","P(p=0)","P(p=1)","..."
do j=1,2*nang-1
        write(10,1000) ang(j),-S12(j)/S11(j),(-S12p(ip,j)/S11p(ip,j),ip=0,min(np,5))
enddo
close(10)

1000 format(' ',1P10E15.5)
2000 format('#',10A15)
2100 format('#',1PE15.5,A)
2300 format('#',I15,A)
2301 format('#',A15,A)

stop
end program call_debye
