!===============!
module blink
!===============!
use hyzer, only: l,nn;

integer,save :: no
integer,save :: frst(nn),fspn(nn)
integer,allocatable :: bnd(:),lpos(:),lvtx(:)

end module blink
!===============!

!=====================!
subroutine mcstep (i)
!=====================!

implicit none

integer :: i
logical :: passed

if(i.eq.0)then
1   call dupdate
    call linkoper
    call updloop_t_worm(passed)
    call updloop_d_worm(passed)
    if (.not.passed) goto 1
else
2   call dupdate
    call linkoper
    call updloop_t_worm(passed)
    call updloop_d_worm(passed)
    if (.not.passed) goto 2
endif

end subroutine mcstep
!=====================!

!==================================!
subroutine dupdate
!==================================!
use hyzer; use bmsr; implicit none

integer :: i,ii,j,k,q,o,b,s1,s2,s3,s4,ss1,ss2,ss3,ss4,iiq
integer :: ns(0:1),nms
real(8) :: rndm,p,addp,delp
integer :: tt1,tt2,tt3,tt4,nu,nk,np,nh1,last,ph1,ph2,ph3,ph4,qi,qj
real(8) :: ssum,sstr,ssus,dl,scc(1:nn,1:nn)

nu=0;nk=0!nu=number of diagonal operators, nk=number of off-diagonal operators
last=0;nh1=0
do i=1,l
    ii=gstring(i)
    if (ii == 0) then !if identity operator
        q=min(int(rndm()*dble(nb))+1,nb) !pick a random bond
        ns(0)=st(bond(0,q)); ns(1) = st(bond(1,q))
        iiq=ns2iq(ns(0),ns(1)) !get the binary number
        addp=beta*dble(nb)/dble(l-nh) !l is like M, the expansion cutoff in the series expansion, nh is like n, the number of non-identity operators
        p=awgt(iiq)*addp
        if (p >= 1.d0) then
            gstring(i)=6*q !insert diagonal operator
            nh=nh+1 !increment num of operators by 1
        elseif (rndm() <= p) then !similarly insert diagonal operator
            gstring(i)=6*q
            nh=nh+1
        endif
    elseif (mod(ii,6) == 0) then !if there already exists a diagonal operator there
        q=ii/6 !get the bond
        do k=0,1
            ns(k)=st(bond(k,q))
        enddo
        iiq=ns2iq(ns(0),ns(1))
        delp=dble(l-nh+1)/(beta*dble(nb)) !probably to remove operator
        p=dwgt(iiq)*delp
        if (p >= 1.d0) then
            gstring(i)=0 !set to identity
            nh=nh-1 !reduce the power
        elseif (rndm() <= p) then
            gstring(i)=0
            nh=nh-1
        endif
        nu=nu+1
        nh1=nh1+1
    else !if off-diagonal operator
        q=ii/6
        o=mod(ii,6)
        ns(0)=st(bond(0,q)); ns(1) = st(bond(1,q))
        dl=dble(nh1-last)
        nh1=nh1+1
        nk = nk + 1 !increment off-diagonal operators
        iiq=ns2iq(ns(0),ns(1))
        jjq=op(o,iiq)
        do k=0,1
            st(bond(k,q))=iq2ns(k,jjq) !here we update the sites on the lattice due to the off diagonal operators
        enddo
    endif
enddo

dl=dble(nh1-last)
ssum=ssum+dl*dble(sa)
sstr=sstr+dl*dble(sa)**2
sstr=sstr/(dble(nh1+1)*dble(nn))
! if (nh1.ne.0) then
!     ssus=(ssum**2)/(dble(nh1)*dble(nh1+1)*dble(nn))
!     ssus=ssus+sstr/dble(nh1+1)
! else
!     ssus=sstr
! endif
! do qi=1,nn
!     do qj=1,nn
!         corr(qi,qj)=corr(qi,qj)+scc(qi,qj)/dble(nms) !boson-boson correlation function normalized.
!     enddo
! enddo

avu=avu+dble(nu)
umag=umag+dble(su)/dble(nn)
! sxu=sxu+beta*(dble(su)**2)/dble(nn)
! ssa=ssa+sstr
! sxa=sxa+beta*ssus
! rhox=rhox+(dble(jjx)**2)/(dble(nn)*beta)
! rhoy=rhoy+(dble(jjy)**2)/(dble(nn)*beta)
! rhotx=rhotx+(dble(jjtx)**2)/(dble(nn)*beta)
! rhoty=rhoty+(dble(jjty)**2)/(dble(nn)*beta)
! rhotpx=rhotpx+(dble(jjtpx)**2)/(dble(nn)*beta)
! rhotpy=rhotpy+(dble(jjtpy)**2)/(dble(nn)*beta)

end subroutine dupdate
!=======================================!

!====================!
subroutine linkoper
!====================!

use blink; use hyzer; implicit none

integer :: i,i0,i1,ii,b,k,o,iiq,jjq
integer :: p,q,p0,p1,p2,p3,s0,s1
integer :: ns(0:1)
integer :: last(nn)

if(allocated(bnd)) deallocate(bnd)
allocate(bnd(0:nh-1))

ii=0
i0=0
i1=1

frst(:)=-1;  last(:)=-1
vert(:)=-1;  link(:)=-1

do p=1,l
    if (gstring(p)/=0) then !if not identity operator
        o=mod(gstring(p),6)
        q=gstring(p)/6
        s0=bond(0,q); s1 = bond(1,q)
        ns(0)=st(bond(0,q)); ns(1) = st(bond(1,q))
        iiq=ns2iq(ns(0),ns(1))
        vert(ii)=vxcode(o,iiq) !store vertex number, we are labelling the operators sequentually from propagation 1, 2,3 ..
        jjq=op(o,iiq) !binary number for resulting plaqeutte state after action of operator
        do k=0,1
            st(bond(k,q))=iq2ns(k,jjq) !update. you need to update here so that the next vertex operator
            !that links to the PROPAGATED spins will 'see' the correct state of the spins.
            !this is necessary for the "vert(ii) = vxcode(o,iiq)" line for the subsequent operators
            !to store the correct vertex number.
        enddo
!       pos(ii)=q
        p0=last(s0)
        p1=last(s1)
        if (p0/=-1) then
            link(p0)=i0 !if the spin is last linked to something, then establish the link
            link(i0)=p0
        else
            frst(s0)=i0
        endif
        if (p1/=-1) then
            link(p1)=i1
            link(i1)=p1
        else
            frst(s1)=i1
        endif
        last(s0)=i0+2
        last(s1)=i1+2
        ii=ii+1
        i0=i0+4 !i0 i1 are labelling ONLY the sites in the linked list
        i1=i1+4 !in ascending order
    endif
enddo

do s0=1,nn !now link the first and last ones together
    i0=frst(s0) !before this, all in between vertices are linked, except for the first and last
    if (i0/=-1) then !check if this site even has an operator acting on it throughout the propagation, if not ignore
        p0=last(s0)
        link(p0)=i0
        link(i0)=p0
    endif
enddo
!no=nh

end subroutine linkoper
!-----------------------!

!=========================!
subroutine updloop_t_worm(passed)
!=========================!

!changes in the configuration is achieved with the loop update, when we traversed the linked list
!ACROSS periodic boundary conditions, i.e going from |\alpha(M)> to |\alpha(0)> state or vice-versa.
!this WILL change the intial and final configurations TOGETHER, but they will continue to remain the same.

use blink; use hyzer; implicit none

integer :: i
integer :: j,k,p0,p1,p2,n4,iiv
integer :: vx,vx0,nv,ml,instate_aft_flip, outstate_aft_flip,new_vx
integer :: ic,ic0,is,is0,vp,oc
real(8) :: r,rndm
logical :: passed

n4=4*nh
ml=100*l
nv=0
do i=1,nl
    ! since T+ or T- only acts on states 1 (-), 2 (0), 3 (+), find an instate that isnt 0
    is0 = 0
    do while (is0 == 0)
        p0=min(int(rndm()*n4),n4-1) !pick a random vertex leg
        p1=p0
        vx0=vert(p0/4) !get vertex number
        ic0=mod(p0,4) !get an in leg number (between 0 to 3)
        is0=vxleg(ic0,vx0) !get the site state of that leg, 0 to 3
    enddo

    do j=1,ml
        vp=p1/4
        vx=vert(vp)
        ic=mod(p1,4) !again calculate in leg number
        is=vxleg(ic,vx)

        if (j==1) then !if it is first time entering the loop
            if (is0 == 2) then
                r = rndm()
                if (r.le.0.5) then 
                    instate_aft_flip = act_tdp(is0)
                else 
                    instate_aft_flip = act_tdm(is0)
                endif
            else if (is0 == 1) then 
                instate_aft_flip = act_tdp(1)
            else
                instate_aft_flip = act_tdm(3)
            endif
        else !if this is not first time entering loop, the instate is just state of prev outstate in linked list
            instate_aft_flip = outstate_aft_flip
        endif

        r=rndm()
        do oc=0,3
            if (r.le.vxprb_t_worm(ic,oc,instate_aft_flip,vx)) then
                new_vx = vxnew(ic,oc,instate_aft_flip,vx)
                outstate_aft_flip = vxleg(oc,new_vx)
                vert(vp)=new_vx !if we accept a probability to exit from a certain leg, goto 10
                exit !if found a satisfactory out leg, exit this do loop
            endif
        enddo

        p1=4*vp+oc
        nv=nv+1
        if ((p1==p0) .and. (vxleg(oc,vert(vp)) == vxleg(ic0,vert(p0/4))) ) goto 20 !!!!!!need to check same state also
        p1=link(p1) !traverse the linked list to the next leg that is linked
        if ((p1==p0) .and. (vxleg(oc,vert(vp)) == vxleg(ic0,vert(p0/4))) ) goto 20
    enddo
    passed=.false.
    return
20   lopers=lopers+dble(nv)
!    if (lopers > 100*l) exit
enddo
j=0
do i=1,l !after loop update is done, update the gstring to reflect the change in operators
    if (gstring(i) /= 0) then
        gstring(i)=6*(gstring(i)/6)+vxoper(vert(j)) !here you can see that gstring(i) mod 6 gives the operator number (0-5)
        j=j+1 !this is going over all the operators in the propagation we labelled earlier, 1,2,3....
    endif
enddo

do i=1,nn
    if (frst(i) /= -1) then
        ic=mod(frst(i),2) !this step gets the in-leg index (0-2)
        vp=frst(i)/6 !gets the index of the operator (1,2,3....) if you look at how i0 and ii scales, this works
        st(i)=vxleg(ic,vert(vp)) !i think this doesnt change the value at st(i) at all..?
        if ((st(i) .ne. 0) .or. (st(i) .ne. 1) .or. (st(i) .ne. 2) .or. (st(i) .ne. 3)) then
            write(*,*)'wrong state',i,vp,ic,vert(vp), st(i)
        endif
    else
        r = rndm()
        if (r.le.0.25) then
            st(i) = 0
        else if (r.le.0.5) then
            st(i) = 1
        else if (r.le.0.75) then
            st(i) = 2
        else
            st(i) = 3
        endif
        !if (rndm().lt.0.5) st(i)= 1-st(i) !how does this ensure the initial and final config is same..?
        !ans: turns out even if this else statement is removed, the SSE algorithm still works.
        !ultimately it does not matter, because this else statement is only triggered if that spin state
        !has NO operators acting on it at all in the propagation, i.e even if you flip it,
        !the initial and final state will remain unchanged, hence |alpha(0)> = !alpha(M> periodicity is fulfilled
    endif
enddo
nloops=nloops+dble(nl)
passed=.true.

end subroutine updloop_t_worm
!====================================!

!=========================!
subroutine updloop_d_worm(passed)
!=========================!

!changes in the configuration is achieved with the loop update, when we traversed the linked list
!ACROSS periodic boundary conditions, i.e going from |\alpha(M)> to |\alpha(0)> state or vice-versa.
!this WILL change the intial and final configurations TOGETHER, but they will continue to remain the same.

use blink; use hyzer; implicit none

integer :: i
integer :: j,k,p0,p1,p2,n4
integer :: vx,vx0,nv,ml,instate_aft_flip
integer :: ic,ic0,is,is0,vp,oc,
real(8) :: r,rndm
logical :: passed

n4=4*nh
ml=100*l
nv=0
do i=1,nl
    ! D,D+ and D- acts on all states
    p0=min(int(rndm()*n4),n4-1) !pick a random vertex leg
    p1=p0
    vx0=vert(p0/4) !get vertex number
    ic0=mod(p0,4) !get an in leg number (between 0 to 3)
    is0=vxleg(ic0,vx0) !get the site state of that leg, 0 to 3

    do j=1,ml
        vp=p1/4
        vx=vert(vp)
        ic=mod(p1,4) !again calculate in leg number
        is=vxleg(ic,vx)

        if (is == 0) then
            r = rndm()
            if (r.le.0.5) then 
                instate_aft_flip = act_tdp(is)
            else 
                instate_aft_flip = act_tdm(is)
            endif
        endif

        r=rndm()
        do oc=0,3
            if (r.le.vxprb(ic,oc,instate_aft_flip,vx)) then
                vert(vp)=vxnew(ic,oc,instate_aft_flip,vx) !if we accept a probability to exit from a certain leg, goto 10
                exit !if found a satisfactory out leg, exit this do loop
            endif
        enddo

        p1=4*vp+oc
        nv=nv+1
        if (p1==p0) goto 20
        p1=link(p1) !traverse the linked list to the next leg that is linked
        if (p1==p0) goto 20 !if the link is closed, go to 20
        enddo
        passed=.false.
        return
20   lopers=lopers+dble(nv)
!    if (lopers > 100*l) exit
enddo
j=0
do i=1,l !after loop update is done, update the gstring to reflect the change in operators
    if (gstring(i) /= 0) then
        gstring(i)=6*(gstring(i)/6)+vxoper(vert(j)) !here you can see that gstring(i) mod 6 gives the operator number (0-5)
        j=j+1 !this is going over all the operators in the propagation we labelled earlier, 1,2,3....
    endif
enddo


do i=1,nn
    if (frst(i) /= -1) then
        ic=mod(frst(i),2) !this step gets the in-leg index (0-2)
        vp=frst(i)/6 !gets the index of the operator (1,2,3....) if you look at how i0 and ii scales, this works
        st(i)=vxleg(ic,vert(vp)) !i think this doesnt change the value at st(i) at all..?
        if ((st(i) .ne. 0) .or. (st(i) .ne. 1) .or. (st(i) .ne. 2) .or. (st(i) .ne. 3)) then
            write(*,*)'wrong state',i,vp,ic,vert(vp), st(i)
        endif
    else
        r = rndm()
        if (r.le.0.25) then
            st(i) = 0
        else if (r.le.0.5) then
            st(i) = 1
        else if (r.le.0.75) then
            st(i) = 2
        else
            st(i) = 3
        endif
        !if (rndm().lt.0.5) st(i)= 1-st(i) !how does this ensure the initial and final config is same..?
        !ans: turns out even if this else statement is removed, the SSE algorithm still works.
        !ultimately it does not matter, because this else statement is only triggered if that spin state
        !has NO operators acting on it at all in the propagation, i.e even if you flip it,
        !the initial and final state will remain unchanged, hence |alpha(0)> = !alpha(M> periodicity is fulfilled
        !endif
    endif
enddo
nloops=nloops+dble(nl)
passed=.true.

end subroutine updloop_d_worm
!====================================!

!================!
subroutine adjnl
!================!
use hyzer; implicit none

integer :: nl1

lopers=lopers/nloops

nl1=1+int(dble(2*l)/lopers)
nl=(nl+nl1)/2
lopers=0.d0
nloops=0.d0
if (nl<0) then
    write(*,*)'nl is negative'
    stop
endif

end subroutine adjnl
!===================!
