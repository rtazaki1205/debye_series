!----------------------------------------------------------------------
!       call debye_series.f90
!----------------------------------------------------------------------
program call_debye
use types
implicit none
integer::nang,ical,ip,np,j
real(kind=dp)::x,ext,sca,back,gs,dang,refre,refim
real(kind=dp),allocatable,dimension(:)  :: ang,S11,S12,POL 
real(kind=dp),allocatable,dimension(:,:):: S11p,S12p,POLp
complex(kind=dp)::refrel
complex(kind=dp),allocatable,dimension(:)  ::S1,S2
complex(kind=dp),allocatable,dimension(:,:)::S1p,S2p

!----------------------------------------------------------------------
!  input parameters
!----------------------------------------------------------------------
x       = 1.0e2_dp      ! size parameter     
refre   = 1.330_dp      ! real part of refractive index
refim   = 1.0e-3_dp     ! imag part of refractive index
nang    = 361           ! number of scat. angles from 0 to 90 degree.
ical    = 0             ! Switch for debye series truncation order
np      = 20            ! Maximum order of the debye series

!----------------------------------------------------------------------
!  call debye_series_rt
!----------------------------------------------------------------------
allocate(S1(1:2*nang-1),S2(1:2*nang-1),ang(1:2*nang-1),S11(1:2*nang-1),&
        S12(1:2*nang-1),S11p(0:np,2*nang-1),S12p(0:np,2*nang-1),&
        S1p(0:np,2*nang-1),S2p(0:np,2*nang-1),POL(1:2*nang-1),POLp(0:np,1:2*nang-1))

refrel = cmplx(refre,refim)
dang   = 90.0_dp/real(nang-1,kind=dp)
call debye_series_rt(x,refrel,nang,ical,np,S1,S2,S1p,S2p,ext,sca,back,gs)
POLp   = 0.0_dp
do j=1,2*nang-1
        ang(j) = dang*dble(j-1)
        S11(j) = 0.5_dp*ABS(S2(J))*ABS(S2(J))
        S11(j) = S11(j)+0.5_dp*ABS(S1(J))*ABS(S1(J))
        S12(j) = 0.5_dp*ABS(S2(J))*ABS(S2(J))
        S12(j) = S12(j)-0.5_dp*ABS(S1(J))*ABS(S1(J))       
        POL(j) = -S12(j)/S11(j)
        do ip=0,np
                S11p(ip,j)=0.5_dp*ABS(S2p(ip,J))*ABS(S2p(ip,J))
                S11p(ip,j)=S11p(ip,j)+0.5_dp*ABS(S1p(ip,J))*ABS(S1p(ip,J))
                S12p(ip,j)=0.5_dp*ABS(S2p(ip,J))*ABS(S2p(ip,J))
                S12p(ip,j)=S12p(ip,j)-0.5_dp*ABS(S1p(ip,J))*ABS(S1p(ip,J)) 
                if(S11p(ip,j) .ne. 0.0_dp) POLp(ip,j)=-S12p(ip,j)/S11p(ip,j)
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
        write(10,1000) ang(j),POL(j),(POLp(ip,j),ip=0,min(np,5))
enddo
close(10)

deallocate(S1,S2,ang,S11,S12,S11p,S12p,S1p,S2p,POL,POLp)

stop

1000 format(' ',1P10E15.5)
2000 format('#',10A15)
2100 format('#',1PE15.5,A)
2300 format('#',I15,A)
2301 format('#',A15,A)
end program call_debye
