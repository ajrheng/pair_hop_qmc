!==============================!
module hyzer
save

integer, parameter :: nx=8, ny=8, nn=nx*ny, nq=nn, nb=2*nn, nn2=nn*nn
real(8), parameter :: z=4.d0 ! 2*nb/nn

integer,parameter :: nvx=52, ivmax=2**8 - 1

integer,parameter :: ntau=100

real(8) :: vv,tt,tp,mu,beta,vv2
integer :: istep,mstep,nruns,equ

integer :: xy(2,nn),xy1(0:nx-1,0:ny-1)

integer :: l,mloop,nh
integer :: st(nn),nst(0:15,0:7),ctra(1:6,0:1,0:1,0:1,0:1)
integer :: ns2iq(0:1,0:1,0:1,0:1),iq2ns(0:3,0:15)
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


!===================!
subroutine lattice
!===================!
use hyzer; implicit none

integer :: i,q,ix,iy,ix1,iy1,ix2,iy2,ix3,iy3,ix4,iy4,m,k1,k2
integer :: iq,iiq,ns(0:3)

i=0
m=0
do iy=0,ny-1
    do ix=0,nx-1
        i=i+1
        xy(1,i)=ix !x coordinate of site number i
        xy(2,i)=iy !y coordinate of site number i
        xy1(ix,iy)=i !given coordinates ix, iy, what site number does it correspond to
enddo
    enddo

do q=1,nn !iterate through plaquette number
    ix1=xy(1,q); iy1=xy(2,q) !get x and y coordinate of the 4 sites in a plaquette
    ix2=mod(ix1+1,nx); iy2=iy1 ! sites on plaquette go in anti-clockwise, if you imagine normal cartesian axes
    ix3=mod(ix1+1,nx); iy3=mod(iy1+1,ny) !the mod is to account for plaquettes wrapping around lattice, i.e PBC
    ix4=ix1; iy4=mod(iy1+1,ny)
    plqt(0,q)=xy1(ix1,iy1) !so site 0 of plaqeutte q is site index xy1(ix1,iy1)
    plqt(1,q)=xy1(ix2,iy2) !etc etc
    plqt(2,q)=xy1(ix3,iy3)
    plqt(3,q)=xy1(ix4,iy4)
    phase(q)=(-1)**(q+m) !this is for calculation of S(\pi,\pi) i think, not so important.
    if (mod(q,nx)==0) then
        m=m+1
    endif
enddo

do iq=0,15 !there are a total of 16 possible plaquette configurations (2^4)
    iiq=iq
    ns(0)=mod(iiq,2); iiq=iiq/2 !we store these configurations as a 4-bit number
    ns(1)=mod(iiq,2); iiq=iiq/2 !eg. 0100 corresponds to 2, so ns2iq(0,1,0,0) = 2
    ns(2)=mod(iiq,2); iiq=iiq/2 !iq2ns(1,2) = 1 for example.
    ns(3)=mod(iiq,2); iiq=iiq/2
    ns2iq(ns(0),ns(1),ns(2),ns(3))=iq
    iq2ns(0,iq)=ns(0)
    iq2ns(1,iq)=ns(1)
    iq2ns(2,iq)=ns(2)
    iq2ns(3,iq)=ns(3)
enddo

!intialize strFact to 0
strFact(:,:)=0.d0

end subroutine lattice
!====================!

!====================!
subroutine pvect0
!====================!
use hyzer; implicit none

integer :: iq,iiq,s1,s2,s3,s4

amax=0.d0; wgt(:)=0.d0; awgt(:)=0.d0; dwgt(:)=0.d0
vv=1; mu=1; vv2=1

do iq=0,15
    iiq=iq
    s1=mod(iiq,2); iiq=iiq/2
    s2=mod(iiq,2); iiq=iiq/2
    s3=mod(iiq,2); iiq=iiq/2
    s4=mod(iiq,2); iiq=iiq/2
    wgt(iq)= vv*dble(s1*s2 + s2*s3 + s3*s4 + s4*s1) !weight of a plaquette configuration, nn-repulsion term
    wgt(iq)= wgt(iq) + vv2*dble(s1*s3 + s2*s4) !diagonal repulsion term
    wgt(iq)= wgt(iq) - mu*dble(s1+s2+s3+s4)/z !chemical potential term
    if (wgt(iq).gt.amax) amax=wgt(iq) !set a maximum weight
enddo
!the larger the wgt(iq), the more unfavourable it is, or the higher energy the configuration is
amax=amax+1.d0
do iq=0,15
    awgt(iq)=amax-wgt(iq)
    if (awgt(iq).gt.1.d-6) then
        dwgt(iq)=1.0/awgt(iq)
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

integer :: i,j,k,m,iq,jq,iiv,nv
integer :: ns(0:7),ns1(0:7)
real(8) :: vvsum,musum

ivx(:)=-1; vxleg(:,:)=-1
nv=0

!=========================================!
! diagonal vertices:    n3  n2   n1  n0
!                       ===============
!                       n3  n2   n1  n0
!=========================================!
do iq=0,15
    !these are to keep track of the vertex with diagonal operators, i.e no change in the sites after action of operator
    ns(0:7)=0
    do i=0,3
        if(btest(iq,i)) ns(i)=1 !returns true if the bit is 1 at position i, for the integer iq
        ns(i+4)=ns(i) !eg. 3 = 2^0 + 2^1, so btest(3,0) and btest(3,1) are true, btest(3,2) etc. are false.
    enddo
    iiv=0
    do k=0,7
        iiv=iiv+ns(k)*(2**k) !converts the 8 site vertex to a unique integer iiv
    enddo
    nv=nv+1 !number the vertex from 1 onwards
    ivx(iiv)=nv; vxi(nv)=iiv !ivx(iiv) takes a vertex integer and returns the vertex number, vxi(nv) is the reverse
    jq=0
    op(0,iq)=iq; vxoper(nv)=0 !op tells you when operator 0 (off-diagonal) acts on plaquette iq, result is iq
    ! vxoper(nv) returns the operator type (0 = diagonal) for a vertex number nv
    vxcode(0,iq)=nv !given an diagonal operator and the state of the plaquette sites before action of operator, what is the vertex number
    do k=0,7
        vxleg(k,nv)=ns(k) !given a vertex number nv, what is the k'th site's state
    enddo
enddo

!========================================!
! single boson hopping vertices
!========================================!

do iq=0,15
    ns(0:7)=0
    do i=0,3
        if(btest(iq,i)) ns(i)=1
        ns(i+4)=ns(i)
    enddo

    do i=0,3
        ns1(:)=ns(:)
        j=mod(i+1,4) !find the neighboring site
        if(ns1(i)/= ns1(j)) then !if the next site is not the same as the previous site, single hop to that site
            ns1(i+4)=1-ns(i) !flip old site.
            ns1(j+4)=1-ns(j) !flip new site, so now boson hopped from old site to new site

            iiv=0
            do k=0,7
                iiv=iiv+ns1(k)*(2**k) !this encodes it as a binary 2 bit number. 2^0 + 2^1 + ...
            enddo
            nv=nv+1
            ivx(iiv)=nv; vxi(nv)=iiv
            jq=0
            do k=0,3
                jq=jq+ns1(k+4)*(2**k)
            enddo
            op(i+1,iq)=jq; vxoper(nv)=i+1 !act operator number i+1 on initial plaquette config iq gives jq
            vxcode(i+1,iq)=nv !given operator number (i+1) and the initial state iq, you know what vertex number it is referring to
            do k=0,7
                vxleg(k,nv)=ns1(k)
            enddo
        endif
    enddo
enddo

do iq=0,15
    ns(0:7)=0
    do i=0,3
        if(btest(iq,i)) ns(i)=1
        ns(i+4)=ns(i)
    enddo
    if((ns(0).eq.0).and.(ns(1).eq.0) &
        .and.(ns(2).eq.1).and.(ns(3).eq.1)) then
        ns(4)=1; ns(5)=1; ns(6)=0; ns(7)=0
        iiv=0
        do k=0,7
            iiv=iiv+ns(k)*(2**k)
        enddo
        nv=nv+1
        ivx(iiv)=nv; vxi(nv)=iiv
        jq=0
        do k=0,3
            jq=jq+ns(k+4)*(2**k)
        enddo
        op(5,iq)=jq; vxoper(nv)=5
        vxcode(5,iq)=nv
        do k=0,7
            vxleg(k,nv)=ns(k)
        enddo
        write(*,*) iiv, nv, jq
     elseif  ((ns(0).eq.1).and.(ns(1).eq.1) &
        .and.(ns(2).eq.0).and.(ns(3).eq.0))then
        ns(4)=0; ns(5)=0; ns(6)=1; ns(7)=1
        iiv=0
        do k=0,7
            iiv=iiv+ns(k)*(2**k)
        enddo
        nv=nv+1
        ivx(iiv)=nv; vxi(nv)=iiv
        jq=0
        do k=0,3
            jq=jq+ns(k+4)*(2**k)
        enddo
        op(5,iq)=jq; vxoper(nv)=5
        vxcode(5,iq)=nv
        do k=0,7
            vxleg(k,nv)=ns(k)
        enddo
        write(*,*) iiv, nv, jq
     elseif  ((ns(0).eq.1).and.(ns(1).eq.0) &
        .and.(ns(2).eq.0).and.(ns(3).eq.1))then
        ns(4)=0; ns(5)=1; ns(6)=1; ns(7)=0
        iiv=0
        do k=0,7
            iiv=iiv+ns(k)*(2**k)
        enddo
        nv=nv+1
        ivx(iiv)=nv; vxi(nv)=iiv
        jq=0
        do k=0,3
            jq=jq+ns(k+4)*(2**k)
        enddo
        op(6,iq)=jq; vxoper(nv)=6
        vxcode(6,iq)=nv
        do k=0,7
            vxleg(k,nv)=ns(k)
        enddo
        write(*,*) iiv, nv, jq
     elseif  ((ns(0).eq.0).and.(ns(1).eq.1) &
        .and.(ns(2).eq.1).and.(ns(3).eq.0))then
        ns(4)=1; ns(5)=0; ns(6)=0; ns(7)=1
        iiv=0
        do k=0,7
            iiv=iiv+ns(k)*(2**k)
        enddo
        nv=nv+1
        ivx(iiv)=nv; vxi(nv)=iiv
        jq=0
        do k=0,3
            jq=jq+ns(k+4)*(2**k)
        enddo
        op(6,iq)=jq; vxoper(nv)=6
        vxcode(6,iq)=nv
        do k=0,7
            vxleg(k,nv)=ns(k)
        enddo
        write(*,*) iiv, nv, jq
    endif
enddo


end subroutine vxweight

program main
use hyzer; implicit none
call lattice
call pvect0
call vxweight

end program main