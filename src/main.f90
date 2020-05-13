!==============================!
module hyzer
save

integer, parameter :: nx=8, ny=8, nn=nx*ny, nq=nn, nb=2*nn, nn2=nn*nn
!real(8), parameter :: z=4.d0 ! 2*nb/nn

integer,parameter :: nvx=52, ivmax=2**8 - 1

integer,parameter :: ntau=100

real(8) :: vv,tt,tp,mu,beta,vv2
integer :: istep,mstep,nruns,equ

integer :: xy(2,nn),xy1(0:nx-1,0:ny-1)

integer :: l,mloop,nh
integer :: st(nn),nst(0:15,0:7),ctra(1:6,0:1,0:1,0:1,0:1)
integer :: ns2iq(0:3,0:3,0:3,0:3),iq2ns(0:3,0:15)
integer :: plqt(0:3,nq),btyp(nb),phase(nn)
real(8) :: amax,wgt(0:15),awgt(0:15),dwgt(0:15)

integer :: vxoper(nvx),vxcode(0:6,0:15),vxleg(0:7,nvx)
integer :: vxnew(0:7,0:7,nvx),op(0:6,0:15)
integer :: vtyp(nvx),ivx(0:ivmax),vxi(nvx)
real(8) :: vxprb(0:7,0:7,nvx),corr(1:nn,1:nn),strFactTemp(0:nx,0:nx),strFact(0:nx,0:nx)
!corr and strFactTemp survives between msteps, strFact writes averaged
!strFact for whole simulation.

integer, allocatable :: gstring(:)
integer, allocatable :: vert(:),link(:)

integer :: nl
real(8) :: lopers,nloops

integer :: iir,jjr,kkr,nnr

end module hyzer
!==============================!

!=============================!
module bmsr

real(8),save :: avu,avk,avp,umag,sxu,ssa,sxa,rhox,rhoy,rhotx,rhoty,rhotpx,rhotpy

end module bmsr
!============================!

!============================!
program main
!============================!
use hyzer; implicit none

call read_params
write(*,*)'read params'
write(*,*)"tt: ",tt,"tp: ",tp,"vv: ",vv,"vv2: ",vv2,"mu: ",mu,"beta: ",beta

call initran
call initconf
call simulation

end program main
!=============================!

!=========================================!
subroutine simulation
!=========================================!
use hyzer; implicit none

integer :: i,j
real(8) :: beta_store

call lattice
call pvect0
call vxweight
call initvrtx

beta_store = beta
beta = 1.d0

if (istep.ne.0) then
    open(12,file='log.txt',status='unknown',access='append')
    write(12,*)'Starting equilibration.'
    close(12)
    lopers=0.d0
    nloops=0.d0
    do i=1,istep

        call mcstep(0)
        call adjstl
        if (mod(i,istep/20).eq.0) call adjnl
        if ((mod(i,3000)==0) .and. (beta < beta_store))  then
            beta = beta + 1.d0
            open(12,file='log.txt',status='unknown',access='append')
            write(12,*)'beta= ',beta
            close(12)
        endif
    enddo
    open(12,file='log.txt',status='unknown',access='append')
    write(12,*)'Completed equilibration: L = ',l
    close(12)
endif

do i=1,nruns
    open(12,file='log.txt',status='unknown',access='append')
    write(12,*)'Starting run ',i
    close(12)
    call zerodata
    do j=1,mstep
        call mcstep(1)
    enddo
    call writeres(mstep)
    call calcCorrStr(mstep)
    call equatestg
    open(12,file='log.txt',status='unknown',access='append')
    write(12,*)'Completed run ',i
    close(12)
    open(UNIT=20,FILE='conf',STATUS='unknown',access='append')
    write(20,*)"Run",i,"conf: "
    call writeconf
    close(20)
enddo
call writestg
deallocate(vert)
deallocate(link)

end subroutine simulation
!=========================!

!=====================!
subroutine mcstep (i)
!=====================!

implicit none

integer :: i
logical :: passed

if(i.eq.0)then
1   call dupdate
    call linkoper
    call updloop(passed)
    if (.not.passed) goto 1
else
2   call dupdate
    call linkoper
    call updloop(passed)
    if (.not.passed) goto 2
endif

end subroutine mcstep
!=====================!

!==========================!
subroutine writeres (nmsr)
!==========================!

use bmsr; use hyzer;

implicit none

integer :: i,j,k,nmsr,x,y,dx,dy,k1,k2
real(8) :: cosk1,sink1,cosk2,sink2,pi

umag=umag/dble(nmsr)
avu=-avu/(dble(beta)*dble(nmsr)*dble(nn))
avu=avu+amax*dble(nb/nn)
avk=-avk/(dble(beta)*dble(nmsr)*dble(nn))
avp=-avp/(dble(beta)*dble(nmsr)*dble(nn))
sxu=sxu/dble(nmsr)
ssa=ssa/dble(nmsr)
sxa=sxa/dble(nmsr)
rhox=rhox/(dble(nmsr))
rhoy=rhoy/(dble(nmsr))
rhotx=rhotx/(dble(nmsr))
rhoty=rhoty/(dble(nmsr))
rhotpx=rhotpx/(dble(nmsr))
rhotpy=rhotpy/(dble(nmsr))

open(UNIT=10,FILE='enr.dat',STATUS='unknown',ACCESS='append')
write(10,*)avu,avk,avp
close(10)

open(UNIT=10,FILE='uni.dat',STATUS='unknown',ACCESS='append')
write(10,*)sxu,umag
close(10)

open(UNIT=10,FILE='stgpipi.dat',STATUS='unknown',ACCESS='append')
write(10,*)ssa,sxa
close(10)

open(UNIT=10,FILE='rho.dat',STATUS='unknown',ACCESS='append')
write(10,*)rhox,rhoy
close(10)

open(UNIT=10,FILE='rhot.dat',STATUS='unknown',ACCESS='append')
write(10,*)rhotx,rhoty
close(10)

open(UNIT=10,FILE='rhotp.dat',STATUS='unknown',ACCESS='append')
write(10,*)rhotpx,rhotpy
close(10)

end subroutine writeres
!======================!

!======================================!
subroutine calcCorrStr(nmsr)
!======================================!
use hyzer; implicit none
integer :: i,j,k1,k2,nmsr,ix1,iy1,ix2,iy2
real(8) :: pi,cosk1,sink1,cosk2,sink2

do i=1,nn
    do j=1,nn
        corr(i,j)=corr(i,j)/(dble(nmsr))
    enddo
enddo

pi=3.141592654
do k1=0,nx
    do k2=0,nx
        do i=1,nn
            do j=1,nn
                ix1=xy(1,i)
                iy1=xy(2,i)
                ix2=xy(1,j)
                iy2=xy(2,j)
                cosk1=DCOS(dble(k1)*2.d0*pi*dble(ix1-ix2)/dble(nx))
                sink1=DSIN(dble(k1)*2.d0*pi*dble(ix1-ix2)/dble(nx))
                cosk2=DCOS(dble(k2)*2.d0*pi*dble(iy1-iy2)/dble(nx))
                sink2=DSIN(dble(k2)*2.d0*pi*dble(iy1-iy2)/dble(nx))
                strFactTemp(k1,k2)=strFactTemp(k1,k2)+(cosk1*cosk2-sink1*sink2)*corr(i,j)
            enddo
        enddo
    enddo
enddo

end subroutine calcCorrStr
!=======================================!

!======================================!
subroutine equatestg
!======================================!
use bmsr; use hyzer;

implicit none

integer :: k1,k2

do k1=0,nx
    do k2=0,nx
        strFact(k1,k2) = strFact(k1,k2)+strFactTemp(k1,k2)/dble(nruns*nn*nn)
    enddo
enddo

end subroutine equatestg
!======================================!

!======================================!
subroutine writestg
!======================================!
use hyzer; implicit none
integer :: k1,k2

open(UNIT=10,FILE='stgfull.dat',STATUS='unknown',ACCESS='append')
do k1=0,nx
    do k2=0,nx
        write(10,*)strFact(k1,k2)
    enddo
enddo
close(10)

end subroutine writestg
!======================================!

!======================================!
subroutine initconf
!======================================!
use hyzer; implicit none

integer :: i
real(8) :: rndm

l=20; nh=0; nl=5

do i=1,nn
    st(i)=int(rndm()*2.0)
enddo

allocate(gstring(l))
gstring(:)=0

allocate(vert(0:l-1))
allocate(link(0:8*l-1))

end subroutine initconf
!=====================================!

!===================!
subroutine readconf
!===================!
use hyzer; implicit none

integer i

read(20,*)l,nh,nl
do i=1,nn
    read(20,*)st(i)
enddo
do i=1,l
    read(20,*)gstring(i)
enddo

end subroutine readconf
!=======================!

!====================!
subroutine writeconf
!====================!
use hyzer; implicit none

integer :: i,no
real :: mn

no=0
write(20,*)"l:",l,"nh:",nh,"nl:",nl
do i=1,nn
    if (st(i)==1) no=no+1
    write(20,"(i1,a1)",advance="no")st(i)," "
    if (mod(i,nx)==0) then
        write(20,*)
    endif
enddo
mn=no/dble(nn)
write(20,*)"Average density:",mn
write(20,*)"----------------------------------------"

end subroutine writeconf
!========================!
!=====================================!
subroutine read_params
!=====================================!
use hyzer;    implicit none

open (unit=10,file='read.in',status='old')
read(10,*)tt!1
read(10,*)tp!start with 0
read(10,*)vv!2,3
read(10,*)vv2
read(10,*)mu!2,3,4...
read(10,*)beta!8*nx
read(10,*)istep!10000
read(10,*)nruns!10
read(10,*)mstep!10000
read(10,*)equ!0
close(10)

end subroutine read_params
!===================================!

!======================================================!
real(8) function rndm()
!======================================================!
use hyzer, only: iir,jjr,kkr,nnr; implicit none

integer :: mzran

mzran=iir-kkr
if (mzran.lt.0) mzran=mzran+2147483579
iir=jjr
jjr=kkr
kkr=mzran
nnr=69069*nnr+1013904243
mzran=mzran+nnr
rndm=.5+.23283064e-9*mzran

return
end function rndm
!========================================================!
