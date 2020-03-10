!SSE QMC for effective Shastry-Sutherland lattice in dimer basis
!see Wessel et al. 10.1103/PhysRevB.98.174432 for details.

!Label local sites with 0 as singlet, 1 as tm (triplet with sz=-1), 2 as t0 (triplet with sz=0), 3 as tp (triplet with sz=+1)

!Bond number to bond mapping:
! 0: 0-0, 4: 0-1, 8: 0-2, 12: 0-3
! 1: 1-0, 5: 1-1, 9: 1-2, 13: 1-3
! 2: 2-0, 6: 2-1, 10: 2-2, 14: 2-3
! 3: 3-0, 7: 3-1, 11: 3-2, 15: 3-3

!==============================!
module hyzer
save

integer, parameter :: nx=8, ny=8, nn=nx*ny, nb=2*nn, nn2=nn*nn
real(8), parameter :: z=4.d0 ! 2*nb/nn

integer,parameter :: max_bond_num = 4**2 - 1, max_vx_num=4**4 - 1, op_num = 12 - 1

integer,parameter :: ntau=100

real(8) :: j,j2,beta
integer :: istep,mstep,nruns,equ

integer :: xy(2,nn),xy1(0:nx-1,0:ny-1)

integer :: l,mloop,nh
integer :: st(nn)
integer :: ns2iq(0:3,0:3),iq2ns(0:1,0:max_bond_num)
integer :: bond(0:1,nb),btyp(nb),phase(nn)
real(8) :: amax,wgt(0:max_bond_num),awgt(0:max_bond_num),dwgt(0:max_bond_num)

integer :: vxoper(max_vx_num),vxcode(0:op_num,0:max_bond_num),vxleg(0:3,max_vx_num)
integer :: vxnew(0:3,0:3,max_vx_num),op(0:op_num,0:max_bond_num) !6 allowed vertices for 2 site bonds, 16 total number of verties
integer :: ivx(0:max_vx_num),vxi(max_vx_num)
real(8) :: vxprb(0:3,0:3,max_vx_num)

integer, allocatable :: gstring(:)
integer, allocatable :: vert(:),link(:)

!wgt arrays tells you the eigenvalues of a local state after the operators acts on it. helps in calculating weights. for eg. tdp == T_d^+ operator
!act arrays give you the value of the resulting state after action of an operator. for eg. tdp(2) = 3 because T_d^+ !0> gives you !+>.
real(8) :: tdp_wgt(0:3), tdm_wgt(0:3), tdz_wgt(0:3), ddp_wgt(0:3), ddm_wgt(0:3), ddz_wgt(0:3)
integer :: act_tdp(0:3), act_tdm(0:3), act_tdz(0:3), act_ddp(0:3), act_ddm(0:3), act_ddz(0:3)


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

!==========================!
subroutine writeres (nmsr)
!==========================!

use bmsr; use hyzer;

implicit none

integer :: i,j,k,nmsr,x,y,dx,dy,k1,k2

umag=umag/dble(nmsr)
avu=-avu/(dble(beta)*dble(nmsr)*dble(nn))
avu=avu+amax*dble(nb/nn)
avk=-avk/(dble(beta)*dble(nmsr)*dble(nn))
avp=-avp/(dble(beta)*dble(nmsr)*dble(nn))
sxu=sxu/dble(nmsr)
ssa=ssa/dble(nmsr)
sxa=sxa/dble(nmsr)
! rhox=rhox/(dble(nmsr))
! rhoy=rhoy/(dble(nmsr))
! rhotx=rhotx/(dble(nmsr))
! rhoty=rhoty/(dble(nmsr))
! rhotpx=rhotpx/(dble(nmsr))
! rhotpy=rhotpy/(dble(nmsr))

open(UNIT=10,FILE='enr.dat',STATUS='unknown',ACCESS='append')
write(10,*)avu,avk,avp
close(10)

open(UNIT=10,FILE='uni.dat',STATUS='unknown',ACCESS='append')
write(10,*)sxu,umag
close(10)

! open(UNIT=10,FILE='stgpipi.dat',STATUS='unknown',ACCESS='append')
! write(10,*)ssa,sxa
! close(10)

! open(UNIT=10,FILE='rho.dat',STATUS='unknown',ACCESS='append')
! write(10,*)rhox,rhoy
! close(10)

! open(UNIT=10,FILE='rhot.dat',STATUS='unknown',ACCESS='append')
! write(10,*)rhotx,rhoty
! close(10)

! open(UNIT=10,FILE='rhotp.dat',STATUS='unknown',ACCESS='append')
! write(10,*)rhotpx,rhotpy
! close(10)

end subroutine writeres
!======================!

! !======================================!
! subroutine calcCorrStr(nmsr)
! !======================================!
! use hyzer; implicit none
! integer :: i,j,k1,k2,nmsr,ix1,iy1,ix2,iy2
! real(8) :: pi,cosk1,sink1,cosk2,sink2

! do i=1,nn
!     do j=1,nn
!         corr(i,j)=corr(i,j)/(dble(nmsr))
!     enddo
! enddo

! pi=3.141592654
! do k1=0,nx
!     do k2=0,nx
!         do i=1,nn
!             do j=1,nn
!                 ix1=xy(1,i)
!                 iy1=xy(2,i)
!                 ix2=xy(1,j)
!                 iy2=xy(2,j)
!                 cosk1=DCOS(dble(k1)*2.d0*pi*dble(ix1-ix2)/dble(nx))
!                 sink1=DSIN(dble(k1)*2.d0*pi*dble(ix1-ix2)/dble(nx))
!                 cosk2=DCOS(dble(k2)*2.d0*pi*dble(iy1-iy2)/dble(nx))
!                 sink2=DSIN(dble(k2)*2.d0*pi*dble(iy1-iy2)/dble(nx))
!                 strFactTemp(k1,k2)=strFactTemp(k1,k2)+(cosk1*cosk2-sink1*sink2)*corr(i,j)
!             enddo
!         enddo
!     enddo
! enddo

! end subroutine calcCorrStr
! !=======================================!

! !======================================!
! subroutine equatestg
! !======================================!
! use bmsr; use hyzer;

! implicit none

! integer :: k1,k2

! do k1=0,nx
!     do k2=0,nx
!         strFact(k1,k2) = strFact(k1,k2)+strFactTemp(k1,k2)/dble(nruns*nn*nn)
!     enddo
! enddo

! end subroutine equatestg
! !======================================!

! !======================================!
! subroutine writestg
! !======================================!
! use hyzer; implicit none
! integer :: k1,k2

! open(UNIT=10,FILE='stgfull.dat',STATUS='unknown',ACCESS='append')
! do k1=0,nx
!     do k2=0,nx
!         write(10,*)strFact(k1,k2)
!     enddo
! enddo
! close(10)

! end subroutine writestg
! !======================================!

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
subroutine initran
!======================================================!
use hyzer, only: iir,jjr,kkr,nnr;   implicit none

integer :: is,js,ks,ls
real(8) ::    rndm

open(10,file='rand.in',status='old')
read(10,*)is
read(10,*)js
read(10,*)ks
read(10,*)ls
close(10)
iir=1+abs(is)
jjr=1+abs(js)
kkr=1+abs(ks)
nnr=ls
open(10,file='rand.in',status='unknown')
write(10,*)abs(nint((rndm()-.5)/.23283064e-9))
write(10,*)abs(nint((rndm()-.5)/.23283064e-9))
write(10,*)abs(nint((rndm()-.5)/.23283064e-9))
write(10,*)abs(nint((rndm()-.5)/.23283064e-9))
close(10)

end subroutine initran
!=======================================================!

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
