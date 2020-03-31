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

integer, parameter :: nx=16, ny=16, nn=nx*ny, nb=2*nn, nn2=nn*nn
real(8), parameter :: z=4.d0 ! 2*nb/nn

integer,parameter :: max_bond_num = 4**2 - 1, max_vx_num=4**4 - 1, op_num = 6-1 !SIX types of opers (0:5), see vxweight

integer,parameter :: ntau=100

real(8) :: j1,j2,beta
integer :: istep,mstep,nruns,equ

integer :: xy(2,nn),xy1(0:nx-1,0:ny-1)

integer :: l,mloop,nh,nvx,num_op_tot !nh = num of operators in opstring, nvx = number of vertex, num_op_tot = total ops in opstring over all runs
integer :: st(nn)
integer :: ns2iq(0:3,0:3),iq2ns(0:1,0:max_bond_num)
integer :: bond(0:1,nb),btyp(nb),phase(nn)
real(8) :: amax,wgt(0:max_bond_num),awgt(0:max_bond_num),dwgt(0:max_bond_num)

integer :: vxoper(max_vx_num),vxcode(0:op_num,0:max_bond_num),vxleg(0:3,max_vx_num)
integer :: op(0:op_num,0:max_bond_num)
integer :: vxnew(0:3,0:3,0:3,max_vx_num) !vxnew(inleg, outleg, in_state_aft_flip, max_vx_num)
integer :: ivx(0:max_vx_num),vxi(max_vx_num) 
real(8) :: vxprb_t_worm(0:3,0:3,0:3,max_vx_num), vxprb_d_worm(0:3,0:3,0:3,max_vx_num) !vxprb(inleg, outleg, in_state_aft_flip ,max_vx_num)
real(8) :: vx_matrix_ele(max_vx_num) !given the vertex_num, gives the value of the matrix ele

integer, allocatable :: gstring(:)
integer, allocatable :: link(:)
integer, allocatable :: vert(:)

!wgt arrays tells you the eigenvalues of a local state after the operators acts on it. helps in calculating weights. for eg. tdp == T_d^+ operator
!act arrays give you the value of the resulting state after action of an operator. for eg. tdp(2) = 3 because T_d^+ !0> gives you !+>.
real(8) :: tdp_wgt(0:3), tdm_wgt(0:3), tdz_wgt(0:3), ddp_wgt(0:3), ddm_wgt(0:3), ddz_wgt(0:3)
integer :: act_tdp(0:3), act_tdm(0:3), act_tdz(0:3), act_ddp(0:3), act_ddm(0:3), act_ddz(0:3)

integer :: nl_t, nl_d
real(8) :: lopers_t,nloops_t,lopers_d, nloops_d

integer :: iir,jjr,kkr,nnr

real(8):: en !energy

end module hyzer
!==============================!

! !=============================!
! module bmsr

! real(8),save :: avu,avk,avp,umag,sxu,ssa,sxa,rhox,rhoy,rhotx,rhoty,rhotpx,rhotpy

! end module bmsr
! !============================!

!============================!
program main
!============================!
use hyzer; implicit none

call read_params
write(*,*)'read params'
write(*,*)"j1: ",j1,"j2: ",j2,"beta: ",beta

call initran
call initconf
call simulation

end program main
!=============================!

!==========================!
subroutine writeres (nmsr)
!==========================!
use hyzer;

implicit none

integer :: nmsr

en = - num_op_tot/(dble(nmsr)*dble(beta)) !energy = <n>/beta = n/(nmsr*beta)
en = en + (amax *dble(nb))
!en = en - ( 0.5*(j1+j2)* dble(nb) )!minus away 1/2(j1+j2)*nb, the diagonal shift
en = en/dble(nn) !energy PER dimer

! umag=umag/dble(nmsr)
! avu=-avu/(dble(beta)*dble(nmsr)*dble(nn))
! avu=avu+amax*dble(nb/nn)
! avk=-avk/(dble(beta)*dble(nmsr)*dble(nn))
! avp=-avp/(dble(beta)*dble(nmsr)*dble(nn))
! sxu=sxu/dble(nmsr)
! ssa=ssa/dble(nmsr)
! sxa=sxa/dble(nmsr)



open(UNIT=10,FILE='energy.txt',STATUS='unknown',ACCESS='append')
write(10,*)en
close(10)

! open(UNIT=10,FILE='uni.dat',STATUS='unknown',ACCESS='append')
! write(10,*)sxu,umag
! close(10)


end subroutine writeres

!======================================!
subroutine initconf
!======================================!
use hyzer; implicit none

integer :: i
real(8) :: rndm

l=20; nh=0; nl_t=5; nl_d=5; num_op_tot = 0

st(:) = -1

do i=1,nn
    st(i)=int(rndm()*4.0)
enddo

allocate(gstring(l))
gstring(:)=0

allocate(vert(0:l-1))
allocate(link(0:4*l-1))

end subroutine initconf
!=====================================!

!===================!
subroutine readconf
!===================!
use hyzer; implicit none

integer i

read(20,*)l,nh,nl_t
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

integer :: i

write(20,*)"l:",l,"nh:",nh,"nl_t:",nl_t, "nl_d", nl_d
do i=1,nn
    ! if (st(i)==1) no=no+1
    write(20,"(i1,a1)",advance="no")st(i)," "
    if (mod(i,nx)==0) then
        write(20,*)
    endif
enddo
! mn=no/dble(nn)
! write(20,*)"Average density:",mn
! write(20,*)"----------------------------------------"

end subroutine writeconf
!========================!
!=====================================!
subroutine read_params
!=====================================!
use hyzer;    implicit none

open (unit=10,file='read.in',status='old')
read(10,*)j1!1
read(10,*)j2!start with 0
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
use hyzer, only: iir,jjr,kkr,nnr;
implicit none

integer :: is,js,ks,ls
real(8) :: rndm

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
