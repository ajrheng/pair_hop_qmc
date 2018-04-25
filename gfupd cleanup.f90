!===============!
module blink
!===============!
use hyzer, only: l,nn;

integer,save :: no
integer,save :: frst(nn),fspn(nn)
integer,allocatable :: bnd(:),lpos(:),lvtx(:)

end module blink
!===============!

!===========!
!module btau
!===========!
!use hyzer, only: l,ntau;

!real(8),allocatable :: wwtaucdf(:,:)

!end module btau
!==============!

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
!2    call gdupdate
2   call dupdate
    call linkoper
    call updloop(passed)
    !call initw
    !call updloop_gfm(passed)
    if (.not.passed) goto 2
endif

end subroutine mcstep
!=====================!

!==================================!
subroutine dupdate
!==================================!
use hyzer; use bmsr; implicit none

integer :: i,ii,j,k,q,o,b,s1,s2,s3,s4,ss1,ss2,ss3,ss4,iiq,jjq,jjx,jjy,sa,su,jjtx,jjty,jjtpx,jjtpy
integer :: ns(0:3),nms
real(8) :: rndm,p,addp,delp
integer :: tt1,tt2,tt3,tt4,nu,nk,np,nh1,last,ph1,ph2,ph3,ph4,qi,qj
real(8) :: ssum,sstr,ssus,dl,scc(1:nn,1:nn)

su=0; sa=0; jjx=0; jjy=0; jjtx=0; jjty=0; jjtpx=0; jjtpy=0;nms=0
scc(:,:)=0.d0!intermediate array for corr function
do i=1,nn
    su=su+st(i)
    sa=sa+phase(i)*st(i)
enddo
ssum=0.d0; sstr=dble(sa)**2
nu=0;nk=0!nu=number of diagonal operators, nk=number of single hop off-diagonal operators
np=0!no of pair hop operators
last=0;nh1=0
do i=1,l
    ii=gstring(i)
    if (ii == 0) then
        q=min(int(rndm()*dble(nq))+1,nq)
        do k=0,3
            ns(k)=st(plqt(k,q))
        enddo
        iiq=ns2iq(ns(0),ns(1),ns(2),ns(3))
        addp=beta*dble(nq)/dble(l-nh)
        p=awgt(iiq)*addp
        if (p >= 1.d0) then
            gstring(i)=7*q
            nh=nh+1
        elseif (rndm() <= p) then
            gstring(i)=7*q
            nh=nh+1
        endif
    elseif (mod(ii,7) == 0) then
        q=ii/7
        do k=0,3
            ns(k)=st(plqt(k,q))
        enddo
        iiq=ns2iq(ns(0),ns(1),ns(2),ns(3))
        delp=dble(l-nh+1)/(beta*dble(nq))
        p=dwgt(iiq)*delp
        if (p >= 1.d0) then
            gstring(i)=0
            nh=nh-1
        elseif (rndm() <= p) then
            gstring(i)=0
            nh=nh-1
        endif
        nu=nu+1
        nh1=nh1+1
    else
        q=ii/7
        o=mod(ii,7)
        do k=0,3
            ns(k)=st(plqt(k,q))
        enddo
        s1=plqt(0,q); s2=plqt(1,q); s3=plqt(2,q); s4=plqt(3,q)
        ss1=st(s1); ss2=st(s2); ss3=st(s3); ss4=st(s4)
        ph1=phase(s1); ph2=phase(s2); ph3=phase(s3); ph4=phase(s4)
        dl=dble(nh1-last)
        nh1=nh1+1
        ssum=ssum+dl*dble(sa)
        sstr=sstr+dl*dble(sa)**2 !for every propagation level that has not changed (dl times)
        ! add sstr, since sstr should be calculated at every propagation level
        if (o==1) then
            tt1=1-ss1; tt2=1-ss2; tt3=ss3; tt4=ss4
            jjx=jjx+ctra(o,ss1,ss2,ss3,ss4)
            jjtx=jjtx+ctra(o,ss1,ss2,ss3,ss4)
            nk=nk+1
        elseif (o==2) then
            tt1=ss1; tt2=1-ss2; tt3=1-ss3; tt4=ss4
            jjy=jjy+ctra(o,ss1,ss2,ss3,ss4)
            jjty=jjty+ctra(o,ss1,ss2,ss3,ss4)
            nk=nk+1
        elseif (o==3) then
            tt1=ss1; tt2=ss2; tt3=1-ss3; tt4=1-ss4
            jjx=jjx+ctra(o,ss1,ss2,ss3,ss4)
            jjtx=jjtx+ctra(o,ss1,ss2,ss3,ss4)
            nk=nk+1
        elseif (o==4) then
            tt1=1-ss1; tt2=ss2; tt3=ss3; tt4=1-ss4
            jjy=jjy+ctra(o,ss1,ss2,ss3,ss4)
            jjty=jjty+ctra(o,ss1,ss2,ss3,ss4)
            nk=nk+1
        elseif (o==5) then
            tt1=1-ss1; tt2=1-ss2; tt3=1-ss3; tt4=1-ss4
            jjy=jjy+ctra(o,ss1,ss2,ss3,ss4)
            jjtpy=jjtpy+ctra(o,ss1,ss2,ss3,ss4)
            np=np+1
        else
            tt1=1-ss1; tt2=1-ss2; tt3=1-ss3; tt4=1-ss4
            jjx=jjx+ctra(o,ss1,ss2,ss3,ss4)
            jjtpx=jjtpx+ctra(o,ss1,ss2,ss3,ss4)
            np=np+1
        endif
        sa=sa+ph1*(tt1-ss1)+ph2*(tt2-ss2)+ph3*(tt3-ss3)+ph4*(tt4-ss4)!change sa only for spins that changed
        last=nh1
        iiq=ns2iq(ns(0),ns(1),ns(2),ns(3))
        jjq=op(o,iiq)
        do k=0,3
            st(plqt(k,q))=iq2ns(k,jjq)
        enddo
    endif
    if (mod(i,8*nn).eq.1) then!update structure factor
        nms=nms+1!separate counter
        do qi=1,nn
            do qj=1,nn
                ss1=st(qi)
                ss2=st(qj)
                scc(qi,qj)=scc(qi,qj)+dble(ss1*ss2)
            enddo
        enddo
    endif
enddo

dl=dble(nh1-last)
ssum=ssum+dl*dble(sa)
sstr=sstr+dl*dble(sa)**2
sstr=sstr/(dble(nh1+1)*dble(nn))
if (nh1.ne.0) then
    ssus=(ssum**2)/(dble(nh1)*dble(nh1+1)*dble(nn))
    ssus=ssus+sstr/dble(nh1+1)
else
    ssus=sstr
endif
do qi=1,nn
    do qj=1,nn
        corr(qi,qj)=corr(qi,qj)+scc(qi,qj)/dble(nms)
    enddo
enddo

avu=avu+dble(nu)
avk=avk+dble(nk)!kinetic of single hop
avp=avp+dble(np)!kinetic of double hop

umag=umag+dble(su)/dble(nn)
sxu=sxu+beta*(dble(su)**2)/dble(nn)
ssa=ssa+sstr
sxa=sxa+beta*ssus
rhox=rhox+3.d0/4*(dble(jjx)**2)/(dble(nn)*beta)
rhoy=rhoy+3.d0/4*(dble(jjy)**2)/(dble(nn)*beta)
rhotx=rhotx+3.d0/4*(dble(jjtx)**2)/(dble(nn)*beta)
rhoty=rhoty+3.d0/4*(dble(jjty)**2)/(dble(nn)*beta)
rhotpx=rhotpx+3.d0/4*(dble(jjtpx)**2)/(dble(nn)*beta)
rhotpy=rhotpy+3.d0/4*(dble(jjtpy)**2)/(dble(nn)*beta)

end subroutine dupdate
!=======================================!

!====================!
subroutine linkoper
!====================!

 use blink; use hyzer; implicit none

 integer :: i,i0,i1,i2,i3,ii,b,k,o,iiq,jjq
 integer :: p,q,p0,p1,p2,p3,s0,s1,s2,s3
 integer :: ns(0:3)
 integer :: last(nn)

 if(allocated(bnd)) deallocate(bnd)
 allocate(bnd(0:nh-1))

  ii=0
  i0=0
  i1=1
  i2=2
  i3=3

  frst(:)=-1;  last(:)=-1
  vert(:)=-1;  link(:)=-1

  do p=1,l
     if (gstring(p)/=0) then
        o=mod(gstring(p),7)
        q=gstring(p)/7
        s0=plqt(0,q)
        s1=plqt(1,q)
        s2=plqt(2,q)
        s3=plqt(3,q)
        do k=0,3
           ns(k)=st(plqt(k,q))
        enddo
        iiq=ns2iq(ns(0),ns(1),ns(2),ns(3))
        vert(ii)=vxcode(o,iiq)
        jjq=op(o,iiq)
        do k=0,3
           st(plqt(k,q))=iq2ns(k,jjq)
        enddo
!        pos(ii)=q
        p0=last(s0)
        p1=last(s1)
        p2=last(s2)
        p3=last(s3)
        if (p0/=-1) then
           link(p0)=i0
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
        if (p2/=-1) then
           link(p2)=i2
           link(i2)=p2
        else
           frst(s2)=i2
        endif
        if (p3/=-1) then
           link(p3)=i3
           link(i3)=p3
        else
           frst(s3)=i3
        endif
        last(s0)=i0+4
        last(s1)=i1+4
        last(s2)=i2+4
        last(s3)=i3+4
        ii=ii+1
        i0=i0+8
        i1=i1+8
        i2=i2+8
        i3=i3+8
     endif
  enddo
  do s0=1,nn
     i0=frst(s0)
     if (i0/=-1) then
        p0=last(s0)
        link(p0)=i0
        link(i0)=p0
     endif
  enddo
  no=nh

end subroutine linkoper
!-----------------------!

!=========================!
subroutine updloop(passed)
!=========================!
 use blink; use hyzer; implicit none

 integer :: i
 integer :: j,k,p0,p1,p2,n8
 integer :: vx,vx0,nv,nv1,ml
 integer :: ic,ic0,is,is0,vp,oc,nop,nop1
 real(8) :: r,rndm
 logical :: passed

  n8=8*nh
  ml=100*l
  nv=0
  do i=1,nl
     nv1=0
     p0=min(int(rndm()*n8),n8-1)
     p1=p0
     vx0=vert(p0/8)
     ic0=mod(p0,8)
     is0=vxleg(ic0,vx0)
     do j=1,ml
        vp=p1/8
        vx=vert(vp)
        ic=mod(p1,8)
        r=rndm()
        do oc=0,6
           if (r.le.vxprb(oc,ic,vx)) goto 10
        enddo
        oc=7
10      vert(vp)=vxnew(oc,ic,vx)
        p1=8*vp+oc
        nv=nv+1
        if (p1==p0) goto 20
        p1=link(p1)
        if (p1==p0) goto 20
     enddo
     passed=.false.
     return
20   lopers=lopers+dble(nv)
!    if (lopers > 100*l) exit
  enddo
  j=0
  do i=1,l
     if (gstring(i) /= 0) then
        gstring(i)=7*(gstring(i)/7)+vxoper(vert(j))
        j=j+1
     endif
  enddo


  do i=1,nn
     if (frst(i) /= -1) then
        ic=mod(frst(i),4)
        vp=frst(i)/8
        st(i)=vxleg(ic,vert(vp))
        if ((st(i) .lt. 0) .or. (st(i) .gt. 1)) then
           write(*,*)'wrong state',i,vp,ic,vert(vp)
        endif
     else
        if (rndm().lt.0.5) st(i)= 1-st(i)
     endif
  enddo
  nloops=nloops+dble(nl)
  passed=.true.

end subroutine updloop
!====================================!


!================!
subroutine adjnl
!================!
 use hyzer; implicit none

 integer :: nl1

 lopers=lopers/nloops
 !write(*,*)'lopers: ',lopers,'nloops: ',nloops
 nl1=1+int(dble(2*l)/lopers)
 nl=(nl+nl1)/2
 lopers=0.d0
 nloops=0.d0
 !write(*,*)'nl:',nl
  if (nl<0) then
    write(*,*)'nl is negative'
    stop
  endif

end subroutine adjnl
!===================!

!==============================!
subroutine updloop_gfm(passed)
!==============================!
 use bgfm; use blink; use btau; use hyzer; implicit none

 integer :: j
 integer :: i,k,p0,p1,p2,vx,vx0,vp,ic,ic0,oc,is,is0,nop,nop1
 real(8) :: r,rndm

 integer :: el0,el1,el2,i0,ij,dn0,nv,ml
 integer :: dirn,idir,vxn
 real(8) :: wsum(0:ntau)
 logical :: cmsr,passed

 ml=100*l
 nv=0
 do j=1,nl
    el0=min(int(rndm()*no),no)
    i0=min(int(nn*rndm())+1,nn)
    idir=1; if(rndm().le.0.5d0) idir=-idir
    dn0=1; if(rndm().le.0.5d0) dn0=-1
    p1=frst(i0)
    if (p1.eq.-1) then
       is=st(i0)+dn0
       if((is.ge.-1).and.(is.le.1)) then
          ccc(i0,i0)=ccc(i0,i0)+1.d0
          ggg(i0,i0,:)=ggg(i0,i0,:)+1.d0
          st(i0)=is
       endif
    else
       el1=p1/8
       ic=mod(p1,8)
       do ! find vertex leg above i0,p0 !
          el1=p1/8
          if (el1.ge.el0) exit
          p1=p1+4   !oc=ic+4
          p2=link(p1)
          el2=p2/8
!          ic=mod(p2,4)
          p1=p2
          if (el2.le.el1) exit
       enddo
       ! set starting vertex and entry leg !
       if (idir.eq.1) then
          p0=p1
       else
          p0=link(p1)
       endif
       vx0=vert(p0/8)
       ic0=mod(p0,8)
       is=vxleg(ic0,vx0)+dn0
       if ((is.ge.-1).and.(is.le.1)) then
          p1=p0
          el1=p1/0
          call gfm_begin(el1,el0,idir,wsum)
          ccc(i0,i0)=ccc(i0,i0)+1.d0 ! equal-time self-correlation
          ggg(i0,i0,:)=ggg(i0,i0,:)+wsum(:) ! self-correlation
          nop1=1
          vp=p1/8; ic=mod(p1,8)
          vx=vert(vp)
          do k=1,ml
!             vp=p1/4
!             vx=vert(p1/4)
!             ic=mod(p1,4)
             r=rndm()
             do oc=0,6
                if (r.le.vxprb(oc,ic,vx)) goto 10
             enddo
             oc=7
10           vert(vp)=vxnew(oc,ic,vx)
             p1=8*vp+oc
             nv=nv+1
             if (oc.gt.1) then
                el1=mod(vp+1,no)
                idir=1
             else
                el1=vp
                idir=-1
             endif
             if (p1==p0) then
                call gfm_end(el1,el0,idir,wsum)
                ggg(i0,i0,:)=ggg(i0,i0,:)+wsum(:)
                go to 20
             endif
             p1=link(p1)
             if (p1==p0) then
                call gfm_end(el1,el0,idir,wsum)
                ggg(i0,i0,:)=ggg(i0,i0,:)+wsum(:)
                go to 20
             endif
             vp=p1/8; ic=mod(p1,8)
             el2=vp; if (ic.gt.1) el2=mod(vp+1,no)
             call gfm_segment(el1,el2,el0,idir,wsum,cmsr)
             ij = vxleg(ic,vert(vp))
!             ij=bst(mod(ic,2)+1,bnd(vp))
             ggg(i0,ij,:)=ggg(i0,ij,:)+wsum(:)
             if (cmsr) ccc(i0,ij)=ccc(i0,ij)+1.d0
             vx=vert(vp)
             nv=nv+1
          enddo
          passed=.false.
          return
20        lopers=lopers+dble(nv)
       endif
    endif
 enddo
 j=0
 do i=1,l
    if (gstring(i) /= 0) then
       gstring(i)=7*(gstring(i)/7)+vxoper(vert(j))
       j=j+1
    endif
 enddo



 do i=1,nn
    if (frst(i) /= -1) then
       ic=mod(frst(i),4)
       vp=frst(i)/8
       st(i)=vxleg(ic,vert(vp))
       if ((st(i) .lt. -1) .or. (st(i) .gt. 1)) then
          write(*,*)'wrong state',i,vp,ic
       endif
    else
       st(i)= int(rndm()*3.0)-1
    endif
 enddo

 nloops=nloops+dble(nl)
 passed=.true.

 deallocate(wwtaucdf)
 deallocate(bnd)

end subroutine updloop_gfm
!=========================!

!=====================!
subroutine initw
!=====================!
 use blink, only: no; use btau; use hyzer, only: l,ntau; implicit none

 integer :: itau,il,i,j
 real(8) :: asum,lnw,sig,tau
 real(8),allocatable :: wwtau(:,:)

 if(allocated(wwtaucdf)) deallocate(wwtaucdf)
 allocate(wwtaucdf(0:ntau,0:l))
 allocate(wwtau(0:ntau,0:l))
 ! assume tau in units of beta !

 do itau=1,ntau-1
    tau=dble(itau)/dble(ntau)
    sig=1.d0-tau
    lnw=no*log(sig)
    wwtau(itau,0)=exp(lnw)
    asum=wwtau(itau,0)
    do il=1,no
       lnw=lnw+log(dble(no-il+1)*tau/(dble(il)*sig))
       wwtau(itau,il)=exp(lnw)
       asum=asum+wwtau(itau,il)
    enddo
    if(abs(asum-1.d0).gt.1.d-6)stop 'un-normalized binomial distribution'
 enddo

 wwtaucdf(:,0)=wwtau(:,0)
 do il=1,no
    wwtaucdf(:,il)=wwtau(:,il)+wwtaucdf(:,il-1)
 enddo
 do il=0,no
    wwtaucdf(0,il)=1.d0
    wwtaucdf(ntau,il)=0.d0
 enddo
 wwtaucdf(ntau,no)=1.d0

 deallocate(wwtau)

end subroutine initw
!===================!

!=====================================================!
subroutine gfm_segment(lone,ltwo,lzero,dirn,tsum,ccm)
!=====================================================!
 use blink, only: no; use btau; use hyzer, only: ntau; implicit none

 integer :: lone,ltwo,lzero,dirn
 real(8) :: tsum(0:ntau)
 logical :: ccm

 ccm=.false.
 tsum(:)=0.d0

 if (dirn.eq.1) then
    if (lone.le.lzero.and.lzero.le.ltwo) then
       tsum(:)=tsum(:)+wwtaucdf(:,no)
       tsum(:)=tsum(:)-wwtaucdf(:,no-lzero+lone-1)
       tsum(:)=tsum(:)+wwtaucdf(:,ltwo-lzero)
       ccm=.true.
    elseif (lzero.lt.lone.and.lone.le.ltwo) then
       tsum(:)=tsum(:)+wwtaucdf(:,ltwo-lzero)
       tsum(:)=tsum(:)-wwtaucdf(:,lone-lzero-1)
    elseif (lone.le.ltwo.and.ltwo.lt.lzero) then
       tsum(:)=tsum(:)+wwtaucdf(:,no-lzero+ltwo)
       tsum(:)=tsum(:)-wwtaucdf(:,no-lzero+lone-1)
    elseif (ltwo.lt.lzero.and.lzero.lt.lone) then
       tsum(:)=tsum(:)+wwtaucdf(:,no-lzero+ltwo)
       tsum(:)=tsum(:)-wwtaucdf(:,lone-lzero-1)
    elseif (lzero.le.ltwo.and.ltwo.lt.lone) then
       tsum(:)=tsum(:)+wwtaucdf(:,no)
       tsum(:)=tsum(:)-wwtaucdf(:,lone-lzero-1)
       tsum(:)=tsum(:)+wwtaucdf(:,ltwo-lzero)
       ccm=.true.
    elseif (ltwo.lt.lone.and.lone.le.lzero) then
       tsum(:)=tsum(:)+wwtaucdf(:,no)
       tsum(:)=tsum(:)-wwtaucdf(:,no-lzero+lone-1)
       tsum(:)=tsum(:)+wwtaucdf(:,no-lzero+ltwo)
       ccm=.true.
    else
       stop 'gfm_segment 1'
    endif
 else if (dirn.eq.-1) then
    if (lone.lt.lzero.and.lzero.lt.ltwo) then
       tsum(:)=tsum(:)+wwtaucdf(:,no-lzero+lone)
       tsum(:)=tsum(:)-wwtaucdf(:,ltwo-lzero-1)
    elseif (lzero.le.lone.and.lone.lt.ltwo) then
       tsum(:)=tsum(:)+wwtaucdf(:,no)
       tsum(:)=tsum(:)-wwtaucdf(:,ltwo-lzero-1)
       tsum(:)=tsum(:)+wwtaucdf(:,lone-lzero)
       ccm=.true.
    elseif (lone.lt.ltwo.and.ltwo.le.lzero) then
       tsum(:)=tsum(:)+wwtaucdf(:,no)
       tsum(:)=tsum(:)-wwtaucdf(:,no-lzero+ltwo-1)
       tsum(:)=tsum(:)+wwtaucdf(:,no-lzero+lone)
       ccm=.true.
    elseif (ltwo.le.lzero.and.lzero.le.lone) then
       tsum(:)=tsum(:)+wwtaucdf(:,no)
       tsum(:)=tsum(:)-wwtaucdf(:,no-lzero+ltwo-1)
       tsum(:)=tsum(:)+wwtaucdf(:,lone-lzero)
       ccm=.true.
    elseif (lzero.lt.ltwo.and.ltwo.le.lone) then
       tsum(:)=tsum(:)+wwtaucdf(:,lone-lzero)
       tsum(:)=tsum(:)-wwtaucdf(:,ltwo-lzero-1)
    elseif (ltwo.le.lone.and.lone.lt.lzero) then
       tsum(:)=tsum(:)+wwtaucdf(:,no-lzero+lone)
       tsum(:)=tsum(:)-wwtaucdf(:,no-lzero+ltwo-1)
    else
       stop 'gfm_segment 2'
    endif
 else
    stop 'gfm_segment 3'
 endif

end subroutine gfm_segment
!=========================!

!==========================================!
subroutine gfm_begin(lone,lzero,dirn,tsum)
!==========================================!
 use blink, only: no; use btau; use hyzer, only: ntau; implicit none

 integer :: lzero,lone,dirn
 real(8) :: tsum(0:ntau)

 tsum(:)=0.d0

 if (dirn.eq.1) then
    if (lone.ge.lzero) then
       tsum(:)=tsum(:)+wwtaucdf(:,lone-lzero)
    else if (lone.lt.lzero) then
       tsum(:)=tsum(:)+wwtaucdf(:,lone-lzero+no)
    else
       stop 'gfm_begin 1'
    endif
 else if (dirn.eq.-1) then
    if (lone.le.lzero) then
       tsum(:)=tsum(:)+wwtaucdf(:,no)
       tsum(:)=tsum(:)-wwtaucdf(:,lone-lzero-1+no)
    else if (lone.gt.lzero) then
       tsum(:)=tsum(:)+wwtaucdf(:,no)
       tsum(:)=tsum(:)-wwtaucdf(:,lone-lzero-1)
    else
       stop 'gfm_begin 2'
    endif
 else
    stop 'gfm_begin 3'
 endif

end subroutine gfm_begin
!========================!

!========================================!
subroutine gfm_end(lone,lzero,dirn,tsum)
!========================================!
 use blink, only: no; use btau; use hyzer, only: ntau; implicit none

 integer :: lone,lzero,dirn
 real(8) :: tsum(0:ntau)

 tsum(:)=0.d0

 if (dirn.eq.1) then
    if (lone.le.lzero) then
       tsum(:)=tsum(:)+wwtaucdf(:,no)
       tsum(:)=tsum(:)-wwtaucdf(:,lone-lzero-1+no)
    else if (lone.gt.lzero) then
       tsum(:)=tsum(:)+wwtaucdf(:,no)
       tsum(:)=tsum(:)-wwtaucdf(:,lone-lzero-1)
    else
       stop 'gfm_end 1'
    endif
 else if (dirn.eq.-1) then
    if (lone.ge.lzero) then
       tsum(:)=tsum(:)+wwtaucdf(:,lone-lzero)
    else if (lone.lt.lzero) then
       tsum(:)=tsum(:)+wwtaucdf(:,lone-lzero+no)
    else
       stop 'gfm_end 2'
    endif
 else
    stop 'gfm_end 3'
 endif

 end subroutine gfm_end
!======================!
