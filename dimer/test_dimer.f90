!==============================!
module hyzer
save

integer, parameter :: nx=2, ny=2, nn=nx*ny, nb=2*nn, nn2=nn*nn
real(8), parameter :: z=4.d0 ! 2*nb/nn

integer,parameter :: max_bond_num = 4**2 - 1, max_vx_num=4**4 - 1, op_num = 6-1 !SIX types of opers (0:5), see vxweight

integer,parameter :: ntau=100

real(8) :: j1,j2,beta
integer :: istep,mstep,nruns,equ

integer :: xy(2,nn),xy1(0:nx-1,0:ny-1)

integer :: l,mloop,nh, nvx !nh = num of operators in opstring, nvx = number of vertex
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

integer :: nl
real(8) :: lopers,nloops

integer :: iir,jjr,kkr,nnr

real(8):: en !energy

end module hyzer
!==============================!
program main
use hyzer; implicit none

j1 = 4
j2 = 0

call lattice
call matrix_ele
call pvect0
call vxweight
call initvrtx



end program main
!===================!
subroutine lattice
!===================!
use hyzer; implicit none

integer :: i,q,ix,iy,ix1,iy1,ix2,iy2
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
        write(*,*)q,bond(0,q),bond(1,q)
    enddo
enddo

do ix=0,ny-1
    do iy=0,nx-1
        q = q+1
        ix1 = ix; iy1 = iy
        ix2 = ix1; iy2 = mod(iy1+1,ny)
        bond(0,q)=xy1(ix1,iy1) 
        bond(1,q)=xy1(ix2,iy2)
        write(*,*)q,bond(0,q),bond(1,q)
    enddo
enddo

do iq=0,max_bond_num 
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
!switch all off-diagonals to positive, but keep the negative signs of the diagonals
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

act_tdz(3) = 3
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
    diag_weight = 0.5*(j1+j2)*tdz_wgt(s1)*tdz_wgt(s2) !only for the diagonal elemets which is Tdz Tdz
    wgt(iq)= diag_weight
    !if (wgt(iq).gt.amax) amax=wgt(iq) !set a maximum weight
enddo

!shift by 1/2(j1+j2) to make certain diagonal vertices zero???
wgt(:) = wgt(:) + 0.5 * (j1+j2)
do iq=0,max_bond_num
    !awgt(iq)=amax-wgt(iq)
    if (wgt(iq).gt.1.d-6) then
        dwgt(iq)=1.0/wgt(iq)
    else
        dwgt(iq)=1.d6
    endif
enddo

awgt(:) = wgt(:)

end subroutine pvect0
!======================!

!==========================!
subroutine vxweight
!==========================!
use hyzer; implicit none

integer :: i,k,iq,iiv,opnum,iiq,jq
integer :: ns(0:3)

ivx(:)=-1; vxleg(:,:)=-1; vxcode(:,:) = -1
nvx=0

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
        nvx=nvx+1
        ivx(iiv) = nvx; vxi(nvx) = iiv
        op(opnum,iq) = iq; vxoper(nvx) = opnum  
        vxcode(opnum,iq) = nvx
        do k = 0,3
            vxleg(k,nvx) = ns(k)
        enddo
        vx_matrix_ele(iiv) = wgt(iq) !for diagonals
    endif
enddo


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
        !write(*,*)'iq',iq, 'ns0', ns(0), 'ns1', ns(1), 'opnum', opnum, act_tdp(ns(0)), act_tdm(ns(1))
        ns(2) = act_tdp(ns(0))
        ns(3) = act_tdm(ns(1))
        iiv = 0
        do k=0,3
            iiv = iiv + ns(k)*(4**k)
        enddo
        nvx=nvx+1
        ivx(iiv) = nvx; vxi(nvx) = iiv
        jq=0
        do k=0,1
            jq=jq+ns(k+2)*(4**k) !encode the resulting vertex into an integer jq
        enddo
        op(opnum,iq) = jq; vxoper(nvx) = opnum  
        vxcode(opnum,iq) = nvx
        do k = 0,3
            vxleg(k,nvx) = ns(k)
        enddo
        vx_matrix_ele(iiv) = 0.5*(j1+j2)*2.d0 !weight of T+T- matrix ele 
    endif

    if (act_tdm(ns(0)) /= -1 .and. act_tdp(ns(1)) /= -1) then !if it is TmTp
    !two if statements in case one vertex number allows both TpTm and TmTp to act on it.
        opnum = 2
        !write(*,*)'iq',iq, 'ns0', ns(0), 'ns1', ns(1), 'opnum', opnum, act_tdm(ns(0)), act_tdp(ns(1))
        ns(2) = act_tdm(ns(0))
        ns(3) = act_tdp(ns(1))
        iiv = 0
        do k=0,3
            iiv = iiv + ns(k)*(4**k)
        enddo
        nvx=nvx+1
        ivx(iiv) = nvx; vxi(nvx) = iiv
        jq=0
        do k=0,1
            jq=jq+ns(k+2)*(4**k)
        enddo
        op(opnum,iq) = jq; vxoper(nvx) = opnum  
        vxcode(opnum,iq) = nvx
        do k = 0,3
            vxleg(k,nvx) = ns(k)
        enddo
        vx_matrix_ele(iiv) = 0.5*(j1+j2)*2.d0 !weight of T+T- matrix ele 
    endif
enddo

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

    if ( act_tdz(ns(0)) /= -1 .and. act_ddz(ns(1)) /= -1) then !if Tz Dz
        opnum = 3
        !write(*,*)'iq',iq, 'ns0', ns(0), 'ns1', ns(1), 'opnum', opnum, act_tdz(ns(0)), act_ddz(ns(1))
        ns(2) = act_tdz(ns(0))
        ns(3) = act_ddz(ns(1))
        iiv = 0
        do k=0,3
            iiv = iiv + ns(k)*(4**k)
        enddo
        nvx=nvx+1
        ivx(iiv) = nvx; vxi(nvx) = iiv
        jq=0
        do k=0,1
            jq=jq+ns(k+2)*(4**k) !encode the resulting vertex into an integer jq
        enddo
        op(opnum,iq) = jq; vxoper(nvx) = opnum  
        vxcode(opnum,iq) = nvx
        do k = 0,3
            vxleg(k,nvx) = ns(k)
        enddo
        vx_matrix_ele(iiv) = 0.5*ABS(j1-j2)*1.d0 !weight of TzDz matrix ele (CHECK IF GOT MINUS SIGN)
    endif

    if (act_tdp(ns(0)) /= -1 .and. act_ddm(ns(1)) /= -1) then
        opnum = 4
        !write(*,*)'iq',iq, 'ns0', ns(0), 'ns1', ns(1), 'opnum', opnum, act_tdp(ns(0)), act_ddm(ns(1))
        ns(2) = act_tdp(ns(0))
        ns(3) = act_ddm(ns(1))
        iiv = 0
        do k=0,3
            iiv = iiv + ns(k)*(4**k)
        enddo
        nvx=nvx+1
        ivx(iiv) = nvx; vxi(nvx) = iiv
        jq=0
        do k=0,1
            jq=jq+ns(k+2)*(4**k) !encode the resulting vertex into an integer jq
        enddo
        op(opnum,iq) = jq; vxoper(nvx) = opnum  
        vxcode(opnum,iq) = nvx
        do k = 0,3
            vxleg(k,nvx) = ns(k)
        enddo
        vx_matrix_ele(iiv) = 0.5*ABS(j1-j2)*2.d0 !weight of T+D- matrix ele (CHECK IF GOT MINUS SIGN)
    endif

    if (act_tdm(ns(0)) /= -1 .and. act_ddp(ns(1)) /= -1) then
        opnum = 5
        !write(*,*)'iq',iq, 'ns0', ns(0), 'ns1', ns(1), 'opnum', opnum, act_tdm(ns(0)), act_ddp(ns(1))
        ns(2) = act_tdm(ns(0))
        ns(3) = act_ddp(ns(1))
        iiv = 0
        do k=0,3
            iiv = iiv + ns(k)*(4**k)
        enddo
        nvx=nvx+1
        ivx(iiv) = nvx; vxi(nvx) = iiv
        jq=0
        do k=0,1
            jq=jq+ns(k+2)*(4**k) !encode the resulting vertex into an integer jq
        enddo
        op(opnum,iq) = jq; vxoper(nvx) = opnum  
        vxcode(opnum,iq) = nvx
        do k = 0,3
            vxleg(k,nvx) = ns(k)
        enddo
        vx_matrix_ele(iiv) = 0.5*ABS(j1-j2)*2.d0 !weight of T-D+ matrix ele (CHECK IF GOT MINUS SIGN)
    endif
enddo

end subroutine vxweight
!==================================!

subroutine initvrtx
!==============================!
use hyzer; implicit none

integer :: ns(0:3),i,iiq,ic,oc,iiv,ns1(0:3),k,ns2(0:3),instate,outstate,vertex_num

vxprb_t_worm(:,:,:,:)=0.d0!vxprb stores the probability of accepting the change
vxprb_d_worm(:,:,:,:) =0.d0 
vxnew(:,:,:,:)=0
do i=1,nvx
    iiq=vxi(i) !retrieve the binary number representing the vertex
    do k = 0,3
        ns(k)=mod(iiq,4); iiq=iiq/4 !undo the binary number to retrieve the spins
    enddo

    do ic=0,3
        instate = ns(ic)
        !====================== START OF T+ T- WORM============================!
        !do TpTm/TmTp worm probabilities first, so check action of T+ or T-
        if (act_tdp(instate) /= -1) then !if can act T+ on inleg
            ns1(:)=ns(:) !copy a ns1 to flip in leg
            ns1(ic) = act_tdp(instate)
            do oc =0,3
                ns2(:) = ns1(:) !copy a ns2 to flip out leg
                outstate = ns2(oc)

                if (act_tdp(outstate) /= -1) then !if can act T+ on outleg
                    ns2(oc) = act_tdp(outstate)
                    iiv = 0
                    do k = 0,3
                        iiv = iiv + ns2(k)*(4**k)
                    enddo
                    if (ivx(iiv)/= -1) then ! if its not -1 then it is a valid vertex
                        vertex_num = ivx(iiv)
                        vxnew(ic,oc,ns1(ic),i) = vertex_num
                        vxprb_t_worm(ic,oc,ns1(ic),i) = vx_matrix_ele(iiv)
                    endif
                endif

                if (act_tdm(outstate) /= -1) then !if can act T- on outleg
                    ns2(oc) = act_tdm(outstate)
                    iiv = 0
                    do k = 0,3
                        iiv = iiv + ns2(k)*(4**k)
                    enddo
                    if (ivx(iiv)/= -1) then ! if its not -1 then it is a valid vertex
                        vertex_num = ivx(iiv)
                        vxnew(ic,oc,ns1(ic),i) = vertex_num
                        vxprb_t_worm(ic,oc,ns1(ic),i) = vx_matrix_ele(iiv)
                    endif
                endif
            enddo
        endif
        
        if (act_tdm(instate) /= -1) then !if can act T- on inleg
            ns1(:)=ns(:) !copy a ns1 to flip in leg
            ns1(ic) = act_tdm(instate)
            do oc =0,3
                ns2(:) = ns1(:) !copy a ns2 to flip out leg
                outstate = ns2(oc)

                if (act_tdp(outstate) /= -1) then !if can act T+ on outleg
                    ns2(oc) = act_tdp(outstate)
                    iiv = 0
                    do k = 0,3
                        iiv = iiv + ns2(k)*(4**k)
                    enddo
                    if (ivx(iiv)/= -1) then ! if its not -1 then it is a valid vertex
                        vertex_num = ivx(iiv)
                        vxnew(ic,oc,ns1(ic),i) = vertex_num
                        vxprb_t_worm(ic,oc,ns1(ic),i) = vx_matrix_ele(iiv)
                    endif
                endif

                if (act_tdm(outstate) /= -1) then !if can act T- on outleg
                    ns2(oc) = act_tdm(outstate)
                    iiv = 0
                    do k = 0,3
                        iiv = iiv + ns2(k)*(4**k)
                    enddo
                    if (ivx(iiv)/= -1) then ! if its not -1 then it is a valid vertex
                        vertex_num = ivx(iiv)
                        vxnew(ic,oc,ns1(ic),i) = vertex_num
                        vxprb_t_worm(ic,oc,ns1(ic),i) = vx_matrix_ele(iiv) 
                    endif
                endif
            enddo
        endif
        !===========END OF T+ T- WORM==============================!

        !============START OF Dz D+ D- WORM====================!
        if (act_ddp(instate) /= -1) then !if can act D+ on inleg
            ns1(:)=ns(:) !copy a ns1 to flip in leg
            ns1(ic) = act_ddp(instate)
            do oc =0,3
                ns2(:) = ns1(:) !copy a ns2 to flip out leg
                outstate = ns2(oc)

                if (act_ddp(outstate) /= -1) then !if can act D+ on outleg
                    ns2(oc) = act_ddp(outstate)
                    iiv = 0
                    do k = 0,3
                        iiv = iiv + ns2(k)*(4**k)
                    enddo
                    if (ivx(iiv)/= -1) then ! if its not -1 then it is a valid vertex
                        vertex_num = ivx(iiv)
                        vxnew(ic,oc,ns1(ic),i) = vertex_num
                        vxprb_d_worm(ic,oc,ns1(ic),i) = vx_matrix_ele(iiv)
                    endif
                endif

                if (act_ddm(outstate) /= -1) then !if can act D- on outleg
                    ns2(oc) = act_ddm(outstate)
                    iiv = 0
                    do k = 0,3
                        iiv = iiv + ns2(k)*(4**k)
                    enddo
                    if (ivx(iiv)/= -1) then ! if its not -1 then it is a valid vertex
                        vertex_num = ivx(iiv)
                        vxnew(ic,oc,ns1(ic),i) = vertex_num
                        vxprb_d_worm(ic,oc,ns1(ic),i) = vx_matrix_ele(iiv)
                    endif
                endif

                if (act_ddz(outstate) /= -1) then !if can act Dz on outleg
                    ns2(oc) = act_ddz(outstate)
                    iiv = 0
                    do k = 0,3
                        iiv = iiv + ns2(k)*(4**k)
                    enddo
                    if (ivx(iiv)/= -1) then ! if its not -1 then it is a valid vertex
                        vertex_num = ivx(iiv)
                        vxnew(ic,oc,ns1(ic),i) = vertex_num
                        vxprb_d_worm(ic,oc,ns1(ic),i) = vx_matrix_ele(iiv) 
                    endif
                endif

            enddo
        endif
        
        if (act_ddm(instate) /= -1) then !if can act D- on inleg
            ns1(:)=ns(:) !copy a ns1 to flip in leg
            ns1(ic) = act_ddm(instate)
            do oc =0,3
                ns2(:) = ns1(:) !copy a ns2 to flip out leg
                outstate = ns2(oc)

                if (act_ddp(outstate) /= -1) then !if can act D+ on outleg
                    ns2(oc) = act_ddp(outstate)
                    iiv = 0
                    do k = 0,3
                        iiv = iiv + ns2(k)*(4**k)
                    enddo
                    if (ivx(iiv)/= -1) then ! if its not -1 then it is a valid vertex
                        vertex_num = ivx(iiv)
                        vxnew(ic,oc,ns1(ic),i) = vertex_num
                        vxprb_d_worm(ic,oc,ns1(ic),i) = vx_matrix_ele(iiv) 
                    endif
                endif

                if (act_ddm(outstate) /= -1) then !if can act D- on outleg
                    ns2(oc) = act_ddm(outstate)
                    iiv = 0
                    do k = 0,3
                        iiv = iiv + ns2(k)*(4**k)
                    enddo
                    if (ivx(iiv)/= -1) then ! if its not -1 then it is a valid vertex
                        vertex_num = ivx(iiv)
                        vxnew(ic,oc,ns1(ic),i) = vertex_num
                        vxprb_d_worm(ic,oc,ns1(ic),i) = vx_matrix_ele(iiv) 
                    endif
                endif

                if (act_ddz(outstate) /= -1) then !if can act Dz on outleg
                    ns2(oc) = act_ddz(outstate)
                    iiv = 0
                    do k = 0,3
                        iiv = iiv + ns2(k)*(4**k)
                    enddo
                    if (ivx(iiv)/= -1) then ! if its not -1 then it is a valid vertex
                        vertex_num = ivx(iiv)
                        vxnew(ic,oc,ns1(ic),i) = vertex_num
                        vxprb_d_worm(ic,oc,ns1(ic),i) = vx_matrix_ele(iiv) 
                    endif
                endif

            enddo
        endif

        if (act_ddz(instate) /= -1) then !if can act Dz on inleg
            ns1(:)=ns(:) !copy a ns1 to flip in leg
            ns1(ic) = act_ddz(instate)
            do oc =0,3
                ns2(:) = ns1(:) !copy a ns2 to flip out leg
                outstate = ns2(oc)

                if (act_ddp(outstate) /= -1) then !if can act D+ on outleg
                    ns2(oc) = act_ddp(outstate)
                    iiv = 0
                    do k = 0,3
                        iiv = iiv + ns2(k)*(4**k)
                    enddo
                    if (ivx(iiv)/= -1) then ! if its not -1 then it is a valid vertex
                        vertex_num = ivx(iiv)
                        vxnew(ic,oc,ns1(ic),i) = vertex_num
                        vxprb_d_worm(ic,oc,ns1(ic),i) = vx_matrix_ele(iiv) 
                    endif
                endif

                if (act_ddm(outstate) /= -1) then !if can act D- on outleg
                    ns2(oc) = act_ddm(outstate)
                    iiv = 0
                    do k = 0,3
                        iiv = iiv + ns2(k)*(4**k)
                    enddo
                    if (ivx(iiv)/= -1) then ! if its not -1 then it is a valid vertex
                        vertex_num = ivx(iiv)
                        vxnew(ic,oc,ns1(ic),i) = vertex_num
                        vxprb_d_worm(ic,oc,ns1(ic),i) = vx_matrix_ele(iiv) 
                    endif
                endif

                if (act_ddz(outstate) /= -1) then !if can act Dz on outleg
                    ns2(oc) = act_ddz(outstate)
                    iiv = 0
                    do k = 0,3
                        iiv = iiv + ns2(k)*(4**k)
                    enddo
                    if (ivx(iiv)/= -1) then ! if its not -1 then it is a valid vertex
                        vertex_num = ivx(iiv)
                        vxnew(ic,oc,ns1(ic),i) = vertex_num
                        vxprb_d_worm(ic,oc,ns1(ic),i) = vx_matrix_ele(iiv) 
                    endif
                endif

            enddo
        endif
        !===========END OF Dz D+ D- WORM==============================!
    enddo
enddo

do i=1,nvx
    do ic=0,3
        do instate = 0,3
            do oc= 1,3
                vxprb_t_worm(ic,oc,instate,i)=vxprb_t_worm(ic,oc,instate,i)+vxprb_t_worm(ic,oc-1,instate,i)
                vxprb_d_worm(ic,oc,instate,i)=vxprb_d_worm(ic,oc,instate,i)+vxprb_d_worm(ic,oc-1,instate,i)
                !write(*,*) vxprb_t_worm(ic,oc,instate,i), vxprb_d_worm(ic,oc,instate,i)
                !write(*,*) 
            enddo
        enddo
    enddo
enddo

do i=1,nvx
    do ic=0,3
        do instate = 0,3
            do oc= 0,3
                if (vxprb_t_worm(ic,3,instate,i) /= 0.d0) then 
                !'if' statement necessary as sometimes it is 0, then you are dividing by 0
                    vxprb_t_worm(ic,oc,instate,i)=vxprb_t_worm(ic,oc,instate,i)/vxprb_t_worm(ic,3,instate,i)
                endif
                if (vxprb_d_worm(ic,3,instate,i) /= 0.d0) then
                    vxprb_d_worm(ic,oc,instate,i)=vxprb_d_worm(ic,oc,instate,i)/vxprb_d_worm(ic,3,instate,i)
                endif
                
                !set all probability = 0s to -1
                if (vxprb_t_worm(ic,oc,instate,i).lt.1.e-6) vxprb_t_worm(ic,oc,instate,i)=-1.
                if (vxprb_d_worm(ic,oc,instate,i).lt.1.e-6) vxprb_d_worm(ic,oc,instate,i)=-1. 
            enddo
        enddo
    enddo
enddo

do i=1,nvx
    do ic=0,3
        do instate = 0,3
            do oc= 0,3
                write(*,*)vxprb_t_worm(ic,oc,instate,i)
            enddo
            write(*,*)
        enddo
    enddo
enddo


end subroutine initvrtx
!===================================!
    