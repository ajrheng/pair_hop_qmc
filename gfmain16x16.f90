!==============================!
module hyzer
save

integer, parameter :: nx=16, ny=16, nn=nx*ny, nq=nn, nb=2*nn, nn2=nn*nn
real(8), parameter :: z=4.d0 ! 2*nb/nn

integer,parameter :: nvx=52, ivmax=2**8 - 1

integer,parameter :: ntau=100

real(8) :: vv,tt,tp,mu,beta
integer :: istep,mstep,nruns,equ

integer :: xy(2,nn)

integer :: l,mloop,nh
integer :: st(nn),nst(0:15,0:7),ctra(1:6,0:1,0:1,0:1,0:1)
integer :: ns2iq(0:3,0:3,0:3,0:3),iq2ns(0:3,0:15)
integer :: plqt(0:3,nq),btyp(nb),phase(nn)
real(8) :: amax,wgt(0:15),awgt(0:15),dwgt(0:15)

integer :: vxoper(nvx),vxcode(0:6,0:15),vxleg(0:7,nvx)
integer :: vxnew(0:7,0:7,nvx),op(0:6,0:15)
integer :: vtyp(nvx),ivx(0:ivmax),vxi(nvx)
real(8) :: vxprb(0:7,0:7,nvx)

integer, allocatable :: gstring(:)
integer, allocatable :: vert(:),link(:)

integer :: nl
real(8) :: lopers,nloops

integer :: iir,jjr,kkr,nnr

end module hyzer
!==============================!

!==============================!
module bgfm
use hyzer, only: nn,ntau;

integer :: ngf
real(8) :: ccc(nn,nn),ggg(nn,nn,0:ntau)
 
end module bgfm
!=============================!

!=============================!
module bmsr

real(8),save :: avu,avk,avp,umag,sxu,ssa,sxa,rhox,rhoy
 
end module bmsr
!============================!

!============================!
program main
!============================!
  use hyzer; implicit none

  integer :: j,seed(4)

  call read_params
  write(*,*)'read params'
  write(*,*)"tt: ",tt,"tp: ",tp,"vv: ",vv,"mu: ",mu,"beta: ",beta

  call initran
  !if (equ.eq.0) then
     call initconf
  !else
     !open (UNIT=20,FILE='conf',STATUS='old')
     !call readconf
!     if (equ.eq.1) call readconf
!     if (equ.eq.2) call readconf2
     !close (20)
  !endif
  call simulation

end program main
!=============================!

!==========================!
subroutine writeres (nmsr)
!==========================!

 use bgfm; use bmsr; use hyzer;

 implicit none

 integer :: i,j,k,itau,nmsr,x,y
! real(8) :: fac(jmax),sum_amax
 real(8) :: gfm(0:nn/2,0:ntau),tau

! fac(1)=2.d0 ! ~nb/nn !
! sum_amax=0.d0
! do j=1,jmax
!    sum_amax=sum_amax+amax(j)*fac(j)
! enddo
 umag=umag/dble(nmsr)
 avu=-avu/(dble(beta)*dble(nmsr)*dble(nn))
! avu=avu+sum_amax
 avu=avu+amax*dble(nb/nn)
 avk=-avk/(dble(beta)*dble(nmsr)*dble(nn))
 avp=-avp/(dble(beta)*dble(nmsr)*dble(nn))
 sxu=sxu/dble(nmsr)
 ssa=ssa/dble(nmsr)
! ssx=ssx/dble(nmsr)
! ssy=ssy/dble(nmsr)
 sxa=sxa/dble(nmsr)
 rhox=rhox/(dble(nmsr))
 rhoy=rhoy/(dble(nmsr))

! open(UNIT=10,FILE='col.dat',STATUS='unknown',ACCESS='append')
! write(10,*)ssx,ssy
! close(10)

 open(UNIT=10,FILE='enr.dat',STATUS='unknown',ACCESS='append')
 write(10,*)avu,avk,avp
 close(10)

 open(UNIT=10,FILE='uni.dat',STATUS='unknown',ACCESS='append')
 write(10,*)sxu,umag
 close(10)

 open(UNIT=10,FILE='stg.dat',STATUS='unknown',ACCESS='append')
 write(10,*)ssa,sxa
 close(10)

 open(UNIT=10,FILE='rho.dat',STATUS='unknown',ACCESS='append')
 write(10,*)rhox,rhoy
 close(10)

 !gfm(:,:)=0.d0
 !do i=1,nn
 !   do k=0,nn/2
 !      j=i+k; if(j.gt. nn) j=j-nn
 !      gfm(k,:)=gfm(k,:)+ggg(i,j,:)
 !   enddo
 !enddo
 !gfm(:,:)=gfm(:,:)/(dble(nmsr)*dble(nl))

 !open(UNIT=10,FILE='gfm.dat',STATUS='unknown',ACCESS='append')
 !do itau=0,ntau/2
 !   do i=0,nn/2
 !      write(10,15)itau,i,0.5d0*(gfm(i,itau)+gfm(i,ntau-itau))
 !   enddo
 !enddo
 !close(10)
!15 format(I5,' ',I5,' ',F16.12)

end subroutine writeres
!======================!

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
 write(20,*)l,nh,nl
 do i=1,nn
    write(20,*)st(i)
    if (st(i)==1) no=no+1
 enddo
 mn=no/dble(nn)
 write(20,*)"Average n:",mn
 !do i=1,l
 !   write(20,*)gstring(i)
 !enddo

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
