!=========================================!
subroutine simulation
!=========================================!
use hyzer; implicit none

integer :: i,j

call lattice
call pvect0
call vxweight
call initvrtx

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

!===================!
subroutine zerodata
!===================!
use bmsr; use hyzer; implicit none
integer:: dx,dy,k1,k2

avu=0.d0
avk=0.d0
avp=0.d0
umag=0.d0
sxu=0.d0
ssa=0.d0
sxa=0.d0
rhox=0.d0
rhoy=0.d0
rhotx=0.d0
rhoty=0.d0
rhotpx=0.d0
rhotpy=0.d0
corr(:,:)=0.d0
strFactTemp(:,:)=0.d0

end subroutine zerodata
!=======================!

!=================!
subroutine adjstl
!=================!
use hyzer;  implicit none

integer :: i,j,p,p1,p2,dl,l1
integer,allocatable :: tstring(:)
real(8) :: r,rndm

dl=l/10+2
if (nh.lt.l-dl/2) return
l1=l+dl

allocate(tstring(l))
tstring(:)=gstring(:)

deallocate(gstring)
allocate(gstring(l1))
do i=1,l
    gstring(i)=tstring(i)
enddo
do i=l+1,l1
    gstring(i)=0
enddo
deallocate(tstring)
r=dble(nh)/dble(l1)
do i=l1,1,-1
    if(gstring(i) /= 0) then
        p1=i
        exit
    endif
enddo
outer: do p2=l1,1,-1
    if (p2 <= p1) then
        exit outer
    elseif (rndm() <= r) then
        gstring(p2)=gstring(p1)
        gstring(p1)=0
        innner: do i=p1-1,1,-1
            if (gstring(i) /= 0) then
                p1=i
                cycle outer
            endif
        enddo innner
    exit outer
    endif
enddo outer

l=l1
mloop=100*l

deallocate(vert)
deallocate(link)

allocate(vert(0:l-1))
allocate(link(0:8*l-1))

end subroutine adjstl
!=======================!


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
        xy(1,i)=ix
        xy(2,i)=iy
        xy1(ix,iy)=i
enddo
    enddo

do q=1,nn
    ix1=xy(1,q); iy1=xy(2,q)
    ix2=mod(ix1+1,nx); iy2=iy1
    ix3=mod(ix1+1,nx); iy3=mod(iy1+1,ny)
    ix4=ix1; iy4=mod(iy1+1,ny)
    plqt(0,q)=xy1(ix1,iy1)
    plqt(1,q)=xy1(ix2,iy2)
    plqt(2,q)=xy1(ix3,iy3)
    plqt(3,q)=xy1(ix4,iy4)
    phase(q)=(-1)**(q+m)
    if (mod(q,nx)==0) then
        m=m+1
    endif
enddo
do iq=0,15
    iiq=iq
    ns(0)=mod(iiq,2); iiq=iiq/2
    ns(1)=mod(iiq,2); iiq=iiq/2
    ns(2)=mod(iiq,2); iiq=iiq/2
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

do iq=0,15
    iiq=iq
    s1=mod(iiq,2); iiq=iiq/2
    s2=mod(iiq,2); iiq=iiq/2
    s3=mod(iiq,2); iiq=iiq/2
    s4=mod(iiq,2); iiq=iiq/2
    wgt(iq)= vv*dble(s1*s2 + s2*s3 + s3*s4 + s4*s1)
    wgt(iq)= wgt(iq) + vv2*dble(s1*s3 + s2*s4)
    wgt(iq)= wgt(iq) - mu*dble(s1+s2+s3+s4)/z
    if (wgt(iq).gt.amax) amax=wgt(iq)
enddo

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
    ns(0:7)=0
    do i=0,3
        if(btest(iq,i)) ns(i)=1
        ns(i+4)=ns(i)
    enddo
    iiv=0
    do k=0,7
        iiv=iiv+ns(k)*(2**k)
    enddo
    nv=nv+1
    ivx(iiv)=nv; vxi(nv)=iiv
    jq=0
    op(0,iq)=iq; vxoper(nv)=0
    vxcode(0,iq)=nv
    do k=0,7
        vxleg(k,nv)=ns(k)
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
        j=mod(i+1,4)
        if(ns1(i)/= ns1(j)) then
            ns1(i+4)=1-ns(i)
            ns1(j+4)=1-ns(j)
            iiv=0
            do k=0,7
                iiv=iiv+ns1(k)*(2**k)
            enddo
            nv=nv+1
            ivx(iiv)=nv; vxi(nv)=iiv
            jq=0
            do k=0,3
                jq=jq+ns1(k+4)*(2**k)
            enddo
            op(i+1,iq)=jq; vxoper(nv)=i+1
            vxcode(i+1,iq)=nv
            do k=0,7
                vxleg(k,nv)=ns1(k)
            enddo
        endif
    enddo
enddo
!======================================!
! pair hopping vertices
!======================================!
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
    endif
enddo

end subroutine vxweight
!==================================!

subroutine initvrtx
!==============================!
use hyzer; implicit none

integer::ns(0:7),i,iiq,ic,oc,iiv,ns1(0:7),k,ns2(0:7),iq,o

vxprb(:,:,:)=0
vxnew(:,:,:)=0
do i=1,nvx
    iiq=vxi(i)
    ns(0)=mod(iiq,2); iiq=iiq/2
    ns(1)=mod(iiq,2); iiq=iiq/2
    ns(2)=mod(iiq,2); iiq=iiq/2
    ns(3)=mod(iiq,2); iiq=iiq/2
    ns(4)=mod(iiq,2); iiq=iiq/2
    ns(5)=mod(iiq,2); iiq=iiq/2
    ns(6)=mod(iiq,2); iiq=iiq/2
    ns(7)=mod(iiq,2); iiq=iiq/2
    do ic=0,7
        ns1(:)=ns(:)
        ns1(ic)=1-ns1(ic)
        ns2(:)=ns1(:)
        do oc=0,7
            ns1(oc)=1-ns1(oc)
            iiv=0
            do k=0,7
                iiv=iiv+ns1(k)*(2**k)
            enddo
            if (ivx(iiv)/=-1) then
                vxnew(oc,ic,i)=ivx(iiv)
                o=vxoper(ivx(iiv))
                if (o==0) then
                    iq=0
                    do k=0,3
                        iq=iq+ns1(k)*(2**k)
                    enddo
                    vxprb(oc,ic,i)=awgt(iq)
                elseif (o==1 .or. o==2 .or. o==3 .or. o==4 ) then
                    vxprb(oc,ic,i)=tt
                    if (o==1 .and. (ns1(1)>ns1(0))) then
                        ctra(o,ns1(0),ns1(1),ns1(2),ns1(3))=-1
                    elseif (o==1 .and. (ns1(1)<ns1(0))) then
                        ctra(o,ns1(0),ns1(1),ns1(2),ns1(3))=1
                    elseif (o==2 .and. ns1(2)>ns1(1)) then
                        ctra(o,ns1(0),ns1(1),ns1(2),ns1(3))=-1
                    elseif (o==2 .and. ns1(2)<ns1(1)) then
                        ctra(o,ns1(0),ns1(1),ns1(2),ns1(3))=1
                    elseif (o==3 .and. ns1(3)>ns1(2)) then
                        ctra(o,ns1(0),ns1(1),ns1(2),ns1(3))=1
                    elseif (o==3 .and. ns1(3)<ns1(2)) then
                        ctra(o,ns1(0),ns1(1),ns1(2),ns1(3))=-1
                    elseif (o==4 .and. ns1(3)>ns1(0)) then
                        ctra(o,ns1(0),ns1(1),ns1(2),ns1(3))=-1
                    else
                        ctra(o,ns1(0),ns1(1),ns1(2),ns1(3))=1
                    endif
                else
                    vxprb(oc,ic,i)=tp
                    if (o==5 .and. ns1(1)<ns1(2)) then
                        ctra(o,ns1(0),ns1(1),ns1(2),ns1(3))=-2
                    elseif (o==5 .and. ns1(1)>ns1(2)) then
                        ctra(o,ns1(0),ns1(1),ns1(2),ns1(3))=2
                    elseif (o==6 .and. ns1(0)>ns1(1)) then
                        ctra(o,ns1(0),ns1(1),ns1(2),ns1(3))=2
                    else
                        ctra(o,ns1(0),ns1(1),ns1(2),ns1(3))=-2
                    endif
                endif
            endif
            ns1(:)=ns2(:)
        enddo
    enddo
enddo

do i=1,nvx
    do ic=0,7
        do oc=1,7
            vxprb(oc,ic,i)=vxprb(oc,ic,i)+vxprb(oc-1,ic,i)!sum the probabilities so it is easier to decide vertex change or not
            !eg, vxprb(0,0,1)=0.2,vxprb(1,0,1)=0.2,vxprb(2,0,1)=0.5,vxprb(3,0,1)=0.1. total probability = 1
            !vxprb(2,0,1)=0.2+0.2+0.5=0.9. Therefore to decide if you accept going in at leg 0, exiting at leg 2 of vertex number 1,
            !generate random number, if random number <0.9, accept it.
        enddo
    enddo
enddo
do i=1,nvx
    do ic=0,7
        do oc=0,7
            vxprb(oc,ic,i)=vxprb(oc,ic,i)/vxprb(7,ic,i)
            if (vxprb(oc,ic,i).lt.1.e-6) vxprb(oc,ic,i)=-1.
        enddo
    enddo
enddo


end subroutine initvrtx
!===================================!
