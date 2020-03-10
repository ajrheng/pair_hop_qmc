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


program main
use hyzer; implicit none

j = 4
j2 = 0

call lattice
call matrix_ele
call pvect0
call vxweight

end program main

!===================!
subroutine lattice
!===================!
use hyzer; implicit none

integer :: i,q,ix,iy,ix1,iy1,ix2,iy2,ix3,iy3,ix4,iy4,k1,k2
integer :: iq,iiq,ns(0:1)

i=0
do iy=0,ny-1
    do ix=0,nx-1
        i=i+1
        xy(1,i)=ix !x coordinate of site number i
        xy(2,i)=iy !y coordinate of site number i
        xy1(ix,iy)=i !given coordinates ix, iy, what site number does it correspond to
    enddo
enddo

q=0
do iy=0,ny-1
    do ix=0,nx-1
        q = q+1
        ix1 = ix; iy1 = iy
        ix2 = mod(ix1+1,nx); iy2 = iy1
        bond(0,q)= xy1(ix1,iy1) !so site 0 of bond q is site index xy1(ix1,iy1)
        bond(1,q)= xy1(ix2,iy2) !etc etc
    enddo
enddo

do ix=0,ny-1
    do iy=0,nx-1
        q = q+1
        ix1 = ix; iy1 = iy
        ix2 = ix1; iy2 = mod(iy1+1,ny)
        bond(0,q)=xy1(ix1,iy1) 
        bond(1,q)=xy1(ix2,iy2)
    enddo
enddo

do iq=0,max_bond_num !there are a total of 4 possible plaquette configurations (2^2)
    iiq=iq
    ns(0)=mod(iiq,4); iiq=iiq/4 !we store these configurations as a 4-bit number
    ns(1)=mod(iiq,4); iiq=iiq/4 !eg. 7 corresponds to 3-1, so ns2iq(3,1) = 7
    ns2iq(ns(0),ns(1))=iq
    iq2ns(0,iq)=ns(0)
    iq2ns(1,iq)=ns(1)
enddo

end subroutine lattice
!====================!

!=====================!
subroutine matrix_ele
!=====================!
use hyzer; implicit none

tdp_wgt(:) = 0.d0; tdm_wgt(:) = 0.d0; tdz_wgt(:) = 0.d0; ddp_wgt(:) = 0.d0; ddm_wgt(:) = 0.d0; ddz_wgt(:) = 0.d0
act_tdp(:) = -1; act_tdm(:) = -1; act_tdz(:) = -1; act_ddp(:) = -1; act_ddm(:) = -1; act_ddz(:) = -1

!handle the weight of the action of each operator on a local state
ddz_wgt(0) = 1.d0
ddp_wgt(0) = SQRT(2.d0)
ddm_wgt(0) = SQRT(2.d0)

tdz_wgt(1) = -1.d0
tdp_wgt(1) = SQRT(2.d0)
ddp_wgt(1) = SQRT(2.d0)

tdp_wgt(2) = SQRT(2.d0)
tdm_wgt(2) = SQRT(2.d0)
ddz_wgt(2) = 1.d0

tdz_wgt(3) = 1.d0
tdm_wgt(3) = SQRT(2.d0)
ddm_wgt(3) = SQRT(2.d0)

!now handle the state resulting from the action of operator on state
act_ddz(0) = 2
act_ddp(0) = 3
act_ddm(0) = 1

act_tdz(1) = 1
act_tdp(1) = 2
act_ddp(1) = 0

act_tdp(2) = 3
act_tdm(2) = 1
act_ddz(2) = 0

act_tdz(3) = 1
act_tdm(3) = 2
act_ddm(3) = 0

end subroutine matrix_ele
!=====================!

!====================!
subroutine pvect0
!====================!
use hyzer; implicit none

integer :: iq,iiq,s1,s2
real(8) :: diag_weight

amax=0.d0; wgt(:)=0.d0; awgt(:)=0.d0; dwgt(:)=0.d0

do iq=0,max_bond_num
    iiq=iq
    s1=mod(iiq,4); iiq=iiq/4
    s2=mod(iiq,4); iiq=iiq/4
    diag_weight = 0.5*(j+j2)*tdz_wgt(s1)*tdz_wgt(s2) !only for the diagonal elemets which is Tdz Tdz
    wgt(iq)= diag_weight
    !if (wgt(iq).gt.amax) amax=wgt(iq) !set a maximum weight
enddo

!amax=amax+1.d0 !shift by 1/2(j+j2) to make certain diagonal vertices zero???
wgt(:) = wgt(:) + 0.5 * (j+j2)
do iq=0,max_bond_num
    !awgt(iq)=amax-wgt(iq)
    if (wgt(iq).gt.1.d-6) then
        dwgt(iq)=1.0/wgt(iq)
    else
        dwgt(iq)=1.d6
    endif
enddo

end subroutine pvect0
!======================!

!==========================!
subroutine vxweight
!==========================!
use hyzer; implicit none

integer :: i,k,m,iq,iiv,nv,opnum,iiq,jq
integer :: ns(0:3)

ivx(:)=-1; vxleg(:,:)=-1
nv=0

!=========================================!
! for opnum, 0 = TzTz (diagonal), 1 = TpTm, 2 = TmTp
! 3 = TzDz, 4 = TpDm, 5 = TmDp
!=========================================!

!=========================================!
! diagonal vertices:    n2  n3
!                       =======
!                       n0  n1
!represent Tdz*Tdz term
!=========================================!
opnum = 0
do iq = 0, max_bond_num
    if (wgt(iq) /= 0.d0) then
        ns(0:3) = 0
        iiq=iq
        do i=0,1
            ns(i)=mod(iiq,4); iiq=iiq/4
            ns(i+2) = ns(i)
        enddo
        iiv = 0
        do k=0,3
            iiv = iiv + ns(k)*(4**k)
        enddo
        nv=nv+1
        ivx(iiv) = nv; vxi(nv) = iiv
        op(0,iq) = iq; vxoper(nv) = opnum  
        vxcode(opnum,iq) = nv
        do k = 0,3
            vxleg(k,nv) = ns(k)
        enddo
    endif
enddo
write(*,*) nv

!========================================!
! Off diagonal vertices
! this represents Td+ Td-
!========================================!
do iq=0,max_bond_num
    ns(0:3) = 0
    iiq=iq
    do i=0,1
        ns(i)=mod(iiq,4); iiq=iiq/4
    enddo

    if (act_tdp(ns(0)) /= -1 .and. act_tdm(ns(1)) /= -1) then
        !meaning TpTm causes a valid vertex on this state
        opnum = 1
        ns(2) = act_tdp(ns(0))
        ns(3) = act_tdm(ns(1))
        iiv = 0
        do k=0,3
            iiv = iiv + ns(k)*(4**k)
        enddo
        nv=nv+1
        ivx(iiv) = nv; vxi(nv) = iiv
        jq=0
        do k=0,1
            jq=jq+ns(k+2)*(4**k) !encode the resulting vertex into an integer jq
        enddo
        op(opnum,iq) = jq; vxoper(nv) = opnum  
        vxcode(opnum,iq) = nv
        do k = 0,3
            vxleg(k,nv) = ns(k)
        enddo
    endif

    if (act_tdm(ns(0)) /= -1 .and. act_tdp(ns(1)) /= -1) then !if it is TmTp
    !two if statements in case one vertex number allows both TpTm and TmTp to act on it.
        opnum = 2
        ns(2) = act_tdm(ns(0))
        ns(3) = act_tdp(ns(1))
        iiv = 0
        do k=0,3
            iiv = iiv + ns(k)*(4**k)
        enddo
        nv=nv+1
        ivx(iiv) = nv; vxi(nv) = iiv
        jq=0
        do k=0,1
            jq=jq+ns(k+2)*(4**k)
        enddo
        op(opnum,iq) = jq; vxoper(nv) = opnum  
        vxcode(opnum,iq) = nv
        do k = 0,3
            vxleg(k,nv) = ns(k)
        enddo
    endif
enddo
write(*,*)nv
!======================================!
! Off diagonal vertices
! this represents TzDz, TpDm and TmDp
!======================================!

do iq=0,max_bond_num
    ns(0:3) = 0
    iiq=iq
    do i=0,1
        ns(i)=mod(iiq,4); iiq=iiq/4
    enddo

    if ( (act_tdz(ns(0)) /= -1 .and. act_ddz(ns(1)) /= -1)) then
        opnum = 3
        ns(2) = act_tdz(ns(0))
        ns(3) = act_ddz(ns(1))
        iiv = 0
        do k=0,3
            iiv = iiv + ns(k)*(4**k)
        enddo
        nv=nv+1
        ivx(iiv) = nv; vxi(nv) = iiv
        jq=0
        do k=0,1
            jq=jq+ns(k+2)*(4**k) !encode the resulting vertex into an integer jq
        enddo
        op(opnum,iq) = jq; vxoper(nv) = opnum  
        vxcode(opnum,iq) = nv
        do k = 0,3
            vxleg(k,nv) = ns(k)
        enddo
    endif

    if (act_tdp(ns(0)) /= -1 .and. act_ddm(ns(1)) /= -1) then
        opnum = 4
        ns(2) = act_tdp(ns(0))
        ns(3) = act_ddm(ns(1))
        iiv = 0
        do k=0,3
            iiv = iiv + ns(k)*(4**k)
        enddo
        nv=nv+1
        ivx(iiv) = nv; vxi(nv) = iiv
        jq=0
        do k=0,1
            jq=jq+ns(k+2)*(4**k) !encode the resulting vertex into an integer jq
        enddo
        op(opnum,iq) = jq; vxoper(nv) = opnum  
        vxcode(opnum,iq) = nv
        do k = 0,3
            vxleg(k,nv) = ns(k)
        enddo
    endif

    if (act_tdm(ns(0)) /= -1 .and. act_ddp(ns(1)) /= -1) then
        opnum = 5
        ns(2) = act_tdm(ns(0))
        ns(3) = act_ddp(ns(1))
        iiv = 0
        do k=0,3
            iiv = iiv + ns(k)*(4**k)
        enddo
        nv=nv+1
        ivx(iiv) = nv; vxi(nv) = iiv
        jq=0
        do k=0,1
            jq=jq+ns(k+2)*(4**k) !encode the resulting vertex into an integer jq
        enddo
        op(opnum,iq) = jq; vxoper(nv) = opnum  
        vxcode(opnum,iq) = nv
        do k = 0,3
            vxleg(k,nv) = ns(k)
        enddo
    endif
enddo
write(*,*) nv
end subroutine vxweight
!==================================!