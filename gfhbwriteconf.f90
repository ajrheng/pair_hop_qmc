!=========================================!
subroutine simulation 
!=========================================!
 use hyzer; implicit none
      
 integer :: i,j

 call lattice
! call operators
 call pvect0
 call vxweight
 call initvrtx
 !call initvrtx_dirloop

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
    open(12,file='log.txt',status='unknown',access='append')
    write(12,*)'Completed run ',i
    close(12)
    open(UNIT=20,FILE='conf',STATUS='unknown',access='append')
    write(20,*)"Run",i,"conf: "
    call writeconf
    close(20)
 enddo

 deallocate(vert)
 deallocate(link)

end subroutine simulation
!=========================!

!===================!
subroutine zerodata
!===================!
 use bgfm; use bmsr;

 avu=0.d0
 avk=0.d0
 avp=0.d0
 umag=0.d0
 sxu=0.d0
 ssa=0.d0
 sxa=0.d0
 rhox=0.d0
 rhoy=0.d0

 ggg(:,:,:)=0.d0

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

 integer :: i,q,ix,iy,ix1,iy1,ix2,iy2,ix3,iy3,ix4,iy4
 integer :: iq,iiq,ns(0:3),xy1(0:nx-1,0:ny-1)

 i=0
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
    phase(q)=(-1)**q
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
    wgt(iq)=wgt(iq) - mu*dble(s1+s2+s3+s4)/z
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
!  xhop(:)=0; yhop(:)=0; xxhop(:)=0; yyhop(:)=0
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
!     do k=0,3
!        jq=jq+ns(k+4)*(2**k)
!     enddo
     op(0,iq)=iq; vxoper(nv)=0
     vxcode(0,iq)=nv
     do k=0,7
        vxleg(k,nv)=ns(k)
     enddo
     !wgt(nv)=awgt(iq)
!     wgt(nv)=0
!     do i=0,3
!        wgt(nv)=wgt(nv)+vv*ns(i)*ns(mod(i+1,4))
!        wgt(nv)=wgt(nv)-mu*dble(ns(i))/z
!     enddo
!     wgt(nv)=amax - wgt(nv)
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
          !wgt(nv)=tt
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
        !wgt(nv)=tp
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
        !wgt(nv)=tp
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
        !wgt(nv)=tp
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
        !wgt(nv)=tp
     endif
  enddo
  

     
  vtyp(1:16)=0; vtyp(17:48)=1; vtyp(49:52)=2

end subroutine vxweight
!==================================!

!============================!
subroutine initvrtx_dirloop
!============================!
 use hyzer; implicit none

 integer :: i,j,k,m,s1,s2,t1,t2,vx,ic,oc,is,os,vxn
 integer :: iiv,jjv,dleq,lst(0:7)
 real(8) :: vxwgt(nvx),legwgt(0:7,0:7,nvx),mwgt(0:7,0:7)

 integer :: jn(0:4),newv(4),dlistic(0:3),dlistis(0:3)
 integer :: dlistvx(0:3),outleg(4)
 real(8) :: jw(0:4)

 integer :: jdog(5),kdog(5),jcat(3),kcat(3)
 integer :: vxj,vxk,icn,ocn,isn,nk,nj

 data jdog/1,2,1,3,2/
 data kdog/4,3,2,4,3/
 data jcat/1,1,2/
 data kcat/3,2,3/

 legwgt(:,:,:)=0.d0
! vxcode(:,:,:)=0
 vxnew(:,:,:)=0
 vxprb(:,:,:)=0.d0


 !Create vxnew links
 do vx=1,nvx
    iiv=vxi(vx)
    do ic=0,7
       dleq=0
       do oc=0,7
          do k=0,7
             lst(k)=vxleg(vx,k)
          enddo
          lst(ic)=1-lst(ic)
          lst(oc)=1-lst(oc)
          jjv=0
          do k=0,7
             jjv=jjv+lst(k)*(2**k)
          enddo
          vxn=ivx(jjv)
          if (vxn.ne.0) then
             dleq=dleq+1
             newv(dleq)=vxn
             outleg(dleq)=oc
             vxnew(oc,ic,vx)=vxn
          endif
       enddo

 !Choose a vertex, in channel, and state change
       if (dleq.eq.4) then

          !Count number of allowed vertices in Directed Loop equations
          do m=1,4
             dlistic(m)=m
             dlistvx(m)=newv(m)
          enddo
          !Now we have list of vertices in Directed Loop equations
          !Initialize vertex order
          do j=1,4
             jn(j)=j
             jw(j)=0.d0
             vxn=dlistvx(j)
             if (vxn.ne.0) jw(j)=vxwgt(vxn)
          enddo
          !Five moves to absolute order
          do m=1,5
             j=jdog(m)
             k=kdog(m)
             if (jw(j).lt.jw(k)) then
                jn(0)=jn(j)
                jn(j)=jn(k)
                jn(k)=jn(0)
                jw(0)=jw(j)
                jw(j)=jw(k)
                jw(k)=jw(0)
             endif
          enddo
          !single-bounce or bounce-free (see Eq. 25 of Sylju\r{a}sen)
          mwgt(1,1)=max(0.d0,jw(1)-jw(2)-jw(3)-jw(4))
          mwgt(1,2)=min(jw(2),0.5d0*(jw(1)+jw(2)-jw(3)-jw(4)))
          mwgt(1,3)=min(jw(3),0.5d0*(jw(1)-jw(2)+jw(3)-jw(4)))
          mwgt(1,4)=jw(4)
          mwgt(2,3)=max(0.d0,0.5d0*(-jw(1)+jw(2)+jw(3)+jw(4)))
          mwgt(2,4)=0.d0
          mwgt(3,4)=0.d0
          do j=2,4
             mwgt(j,j)=0.d0
             do k=1,j-1
                mwgt(j,k)=mwgt(k,j)
             enddo
          enddo
          !Finally, assign weights to Directed Loop vertex elements
          do j=1,4
             nj=jn(j)
             vxn=dlistvx(nj)
             if (vxn.ne.0) then
                icn=dlistic(nj)
                do k=1,4
                   nk=jn(k)
                   vxk=dlistvx(nk)
                   if (vxk.ne.0) then
                      !Cycle through to find correct out channel
                      !(the one that changes vxj to vxk)
                      ocn=-1
                      do m=1,4
                         oc=outleg(m)
                         vxj=vxnew(oc,icn,vxn)
                         if (vxj.eq.vxk) ocn=oc
                      enddo
                      if (ocn.eq.-1) stop
                      legwgt(ocn,icn,vxn)=mwgt(j,k)
                   endif
                enddo
             endif
          enddo
       elseif (dleq.eq.3) then
          !Count number of allowed vertices in Directed Loop equations
          do m=1,3
             dlistic(m)=m
             dlistvx(m)=newv(m)
          enddo
          !Now we have list of vertices in Directed Loop equations
          !Initialize vertex order
          do j=1,3
             jn(j)=j
             jw(j)=0.d0
             vxn=dlistvx(j)
             if (vxn.ne.0) jw(j)=vxwgt(vxn)
          enddo
          !Three moves to absolute order
          
          do m=1,3
             j=jcat(m)
             k=kcat(m)
             if (jw(j).lt.jw(k)) then
                jn(0)=jn(j)
                jn(j)=jn(k)
                jn(k)=jn(0)
                jw(0)=jw(j)
                jw(j)=jw(k)
                jw(k)=jw(0)
             endif
          enddo
          !single-bounce or bounce-free (see Eq. 25 of Sylju\r{a}sen)
          mwgt(1,1)=max(0.d0,jw(1)-jw(2)-jw(3))
          mwgt(1,2)=min(jw(2),0.5d0*(jw(1)+jw(2)-jw(3)))
          mwgt(1,3)=min(jw(3),0.5d0*(jw(1)-jw(2)+jw(3)))
          mwgt(2,3)=max(0.d0,0.5d0*(-jw(1)+jw(2)+jw(3)))
          do j=2,3
             mwgt(j,j)=0.d0
             do k=1,j-1
                mwgt(j,k)=mwgt(k,j)
             enddo
          enddo
          !Finally, assign weights to Directed Loop vertex elements
          do j=1,3
             nj=jn(j)
             vxn=dlistvx(nj)
             if (vxn.ne.0) then
                icn=dlistic(nj)
                do k=1,3
                   nk=jn(k)
                   vxk=dlistvx(nk)
                   if (vxk.ne.0) then
                      !Cycle through to find correct out channel
                      !(the one that changes vxj to vxk)
                      ocn=-1
                      do m=1,3
                         oc=outleg(m)
                         vxj=vxnew(oc,icn,vxn)
                         if (vxj.eq.vxk) ocn=oc
                      enddo
                      if (ocn.eq.-1) stop
                      legwgt(ocn,icn,vxn)=mwgt(j,k)
                   endif
                enddo
             endif
          enddo
          
       endif
    enddo
 enddo


 !Convert matrix weights to probabilities
 do vx=1,nvx
    do ic=0,7
       do oc=0,7
          vxprb(oc,ic,vx)=legwgt(oc,ic,vx)/vxwgt(vx)
       enddo
    enddo
 enddo

 !Create cumulative probabilities
 do vx=1,nvx
    do ic=0,7
       do oc=1,7
          vxprb(oc,ic,vx)=vxprb(oc,ic,vx)+vxprb(oc-1,ic,vx)
       enddo
    enddo
 enddo

 !Normalize cumulative probabilities
 do vx=1,nvx
    do ic=0,7
       do oc=0,7
          vxprb(oc,ic,vx)=vxprb(oc,ic,vx)/vxprb(3,ic,vx)
          if (vxprb(oc,ic,vx).lt.1.d-6) vxprb(oc,ic,vx)=-1.d0
       enddo
    enddo
 enddo

end subroutine initvrtx_dirloop
!==============================!
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
            !write(*,*)i," ",ic," ",oc," ",vxprb(oc,ic,i)
            !write(*,*)ns(0)," ",ns(1)," ",ns(2)," ",ns(3)," ",ns(4)," ",ns(5)," ",ns(6)," ",ns(7)

        enddo
    enddo
enddo


end subroutine initvrtx







