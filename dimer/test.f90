!==============================!
module hyzer
    save
    
    integer, parameter :: nx=4, ny=4, nn=nx*ny, nb=2*nn, nn2=nn*nn
    real(8), parameter :: z=4.d0 ! 2*nb/nn
    
    integer,parameter :: nvx=6, ivmax=2**4 - 1
    
    integer,parameter :: ntau=100
    
    real(8) :: j1,j2,beta
    integer :: istep,mstep,nruns,equ
    
    integer :: xy(2,nn),xy1(0:nx-1,0:ny-1)
    
    integer :: l,mloop,nh
    integer :: st(nn)!,ctra(1:6,0:1,0:1,0:1,0:1)
    integer :: ns2iq(0:1,0:1),iq2ns(0:1,0:3)
    integer :: bond(0:1,nb),btyp(nb),phase(nn)
    real(8) :: amax,wgt(0:3),awgt(0:3),dwgt(0:3)
    
    integer :: vxoper(nvx),vxcode(0:5,0:15),vxleg(0:7,nvx)
    integer :: vxnew(0:3,0:3,nvx),op(0:5,0:3) !6 allowed vertices for 2 site bonds, 16 total number of verties
    integer :: ivx(0:ivmax),vxi(nvx)
    real(8) :: vxprb(0:3,0:3,nvx)!,corr(1:nn,1:nn),strFactTemp(0:nx,0:nx),strFact(0:nx,0:nx)
    !corr and strFactTemp survives between msteps, strFact writes averaged
    !strFact for whole simulation.
    
    integer, allocatable :: gstring(:)
    integer, allocatable :: vert(:),link(:)
    
    integer :: nl
    real(8) :: lopers,nloops
    
    integer :: iir,jjr,kkr,nnr
    
    end module hyzer
    !==============================!
    


program main
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

do iq=0,3 !there are a total of 4 possible plaquette configurations (2^2)
    iiq=iq
    ns(0)=mod(iiq,2); iiq=iiq/2 !we store these configurations as a 2-bit number
    ns(1)=mod(iiq,2); iiq=iiq/2 !eg. 10 corresponds to 2, so ns2iq(1,0) = 2
    ns2iq(ns(0),ns(1))=iq
    iq2ns(0,iq)=ns(0)
    iq2ns(1,iq)=ns(1)
enddo

do i = 1,nb
    write(*,*) bond(0,i), bond(1,i)
enddo

end program main