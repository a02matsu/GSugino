!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ver.01: correct the normalization of the PCSC relation 
!! ver.04: bug fix for mass part of PCSC
!! ver.05: include WT id. in naive quench
!! ver.06: added compensator for SU(2) (triple cover version)
!module observables
!use global_parameters
!use global_subroutines
!implicit none
!
!contains

#include "WT_identities.f90"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_bosonic_action(Sb,UMAT,PhiMat)
!use hamiltonian
implicit none

double precision, intent(out) :: Sb
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

double precision :: SB_S,SB_L,SB_F,SB_local

SB_S=0d0
SB_L=0d0
SB_F=0d0
SB=0d0
call bosonic_action_site(SB_S,PhiMat)
call bosonic_action_link(SB_L,UMAT,PhiMat)
call bosonic_action_face(SB_F,UMAT)

Sb_local=SB_S+SB_L+SB_F

#ifdef PARALLEL
call MPI_REDUCE(SB_local,SB,1,MPI_DOUBLE_PRECISION, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)
#else
SB=SB_local
#endif

end subroutine calc_bosonic_action

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_TrX2(TrX2, PhiMat)
implicit none

double precision, intent(out) :: TrX2
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
double precision :: tmp

integer :: s,i,j

tmp=0d0
do s=1,num_sites
  do i=1,NMAT
    do j=1,NMAT
      tmp=tmp+ dble( PhiMat(i,j,s)*conjg( PhiMat(i,j,s) ) )
    enddo
  enddo
enddo
#ifdef PARALLEL
call MPI_REDUCE(tmp,TrX2,1,MPI_DOUBLE_PRECISION, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)
#else
TrX2=tmp
#endif

TrX2 = TrX2 / (2d0 * dble(global_num_sites) )

end subroutine calc_TrX2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_trace_compensator(Acomp,PhiMat)
use parallel
implicit none

complex(kind(0d0)), intent(out) :: Acomp
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) tmp,A_tmp
double precision :: radius, phase, ratio
integer :: s,i,j,eular


eular=global_num_sites-global_num_links+global_num_faces 
ratio=dble(-(NMAT*NMAT-1)*eular)/4d0 
Acomp=(0d0,0d0)
A_tmp=(0d0,0d0)
do s=1,num_sites
  tmp=(0d0,0d0)
  do i=1,NMAT
    do j=1,NMAT
      tmp=tmp+PhiMat(i,j,s)*dconjg(PhiMat(i,j,s))
    enddo
  enddo
  tmp=(tmp/dcmplx(dble(NMAT)))
  radius=cdabs(tmp)
  phase=atan2(dble(tmp),dble(tmp*(0d0,-1d0)))

  A_tmp=A_tmp + dcmplx(radius**ratio) * cdexp( (0d0,1d0)*dcmplx(phase*ratio) )
enddo

call MPI_REDUCE(A_tmp,Acomp,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)
  
Acomp=Acomp/dcmplx(dble(global_num_sites))

end subroutine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_VM_compensator(Acomp,PhiMat)
use parallel
implicit none

complex(kind(0d0)), intent(out) :: Acomp
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: Atmp,cartan(1:NMAT-1),VMdet
integer :: s,i,j,eular

eular=global_num_sites-global_num_links+global_num_faces
Acomp=(0d0,0d0)
Atmp=(0d0,0d0)

do s=1, num_sites
  call calc_cartan(cartan,PhiMat(:,:,s))
  VMdet=(1d0,0d0)
  do I=1,NMAT-2
    do j=I+1,NMAT-1
      VMdet=VMdet*(cartan(i)-cartan(j))
    enddo
  enddo
  do i=1,NMAT-1
    VMdet=VMdet * cartan(i)**0.5d0
  enddo
  VMdet=VMdet**(-eular)

  Atmp=Atmp+VMdet
enddo

call MPI_REDUCE(Atmp,Acomp,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)
  
Acomp=Acomp/dcmplx(dble(global_num_sites))

end subroutine calc_VM_compensator


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_eigenvalues_Dirac(eigenvalues,UMAT,PhiMat)
use Dirac_operator, only : make_Dirac
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
!complex(kind(0d0)), intent(inout) :: eigenvalues(1:sizeD)
complex(kind(0d0)), intent(inout) :: eigenvalues(:)

!complex(kind(0d0)) :: Dirac(1:sizeD,1:sizeD)
!complex(kind(0d0)) VL(1,1:sizeD), VR(1,1:sizeD)
!double precision RWORK(1:2*sizeD)
!complex(kind(0d0)) WORK(2*(sizeD)),tako1
complex(kind(0d0)), allocatable :: Dirac(:,:)
complex(kind(0d0)), allocatable :: VL(:,:), VR(:,:)
double precision, allocatable :: RWORK(:)
complex(kind(0d0)), allocatable :: WORK(:)
complex(kind(0d0)) tako1
character JOBVL,JOBVR
integer info,lwork
integer i,j
integer :: sizeD

sizeD=size(eigenvalues,1)

allocate( Dirac(1:sizeD,1:sizeD) )
allocate( VL(1,1:sizeD), VR(1,1:sizeD) )
allocate( RWORK(1:2*sizeD) )
allocate( WORK(2*(sizeD)) )

lwork=2*sizeD
JOBVL='N'
JOBVR='N'

call make_Dirac(Dirac,UMAT,PhiMat)
!if( MYRANK==0 ) then 
!do i=1,sizeD
!  do j=1,sizeD
!    write(*,*) i,j,Dirac(i,j)
!  enddo
!enddo
!endif

#ifdef PARALLEL
if( MYRANK==0 ) then 
#endif
call ZGEEV(JOBVL,JOBVR,sizeD,&
     DIRAC,sizeD,eigenvalues,VL,1,VR,1,WORK,lwork,RWORK,info)
#ifdef PARALLEL
endif
#endif

! sort the eigenvalues
do i=1,sizeD
 do j=i+1,sizeD
  tako1 = eigenvalues(i)
  if(abs(eigenvalues(j)).LT.abs(eigenvalues(i))) then 
    eigenvalues(i) = eigenvalues(j)
    eigenvalues(j) = tako1
  endif
 enddo
enddo

end subroutine calc_eigenvalues_Dirac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_eigenvalues_DdagD(eigenvalues,UMAT,PhiMat)
use Dirac_operator, only : make_DdagD
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
!complex(kind(0d0)), intent(inout) :: eigenvalues(1:sizeD)
complex(kind(0d0)), intent(inout) :: eigenvalues(:)

!complex(kind(0d0)) :: Dirac(1:sizeD,1:sizeD)
!complex(kind(0d0)) VL(1,1:sizeD), VR(1,1:sizeD)
!double precision RWORK(1:2*sizeD)
!complex(kind(0d0)) WORK(2*(sizeD)),tako1
complex(kind(0d0)), allocatable :: DdagD(:,:)
complex(kind(0d0)), allocatable :: VL(:,:), VR(:,:)
double precision, allocatable :: RWORK(:)
complex(kind(0d0)), allocatable :: WORK(:)
complex(kind(0d0)) tako1
character JOBVL,JOBVR
integer info,lwork
integer i,j
integer :: sizeD

sizeD=size(eigenvalues,1)

allocate( DdagD(1:sizeD,1:sizeD) )
allocate( VL(1,1:sizeD), VR(1,1:sizeD) )
allocate( RWORK(1:2*sizeD) )
allocate( WORK(2*(sizeD)) )

lwork=2*sizeD
JOBVL='N'
JOBVR='N'

call make_DdagD(DdagD,UMAT,PhiMat)
!if( MYRANK==0 ) then 
!do i=1,sizeD
!  do j=1,sizeD
!    write(*,*) i,j,DdagD(i,j)
!  enddo
!enddo
!endif

#ifdef PARALLEL
if( MYRANK==0 ) then 
#endif
call ZGEEV(JOBVL,JOBVR,sizeD,&
     DdagD,sizeD,eigenvalues,VL,1,VR,1,WORK,lwork,RWORK,info)
#ifdef PARALLEL
endif
#endif

! sort the eigenvalues
do i=1,sizeD
 do j=i+1,sizeD
  tako1 = eigenvalues(i)
  if(abs(eigenvalues(j)).LT.abs(eigenvalues(i))) then 
    eigenvalues(i) = eigenvalues(j)
    eigenvalues(j) = tako1
  endif
 enddo
enddo

end subroutine calc_eigenvalues_DdagD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_min_and_max_of_eigenvalues_Dirac(minimal,maximal,UMAT,PhiMat)
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(inout) :: minimal, maximal

!complex(kind(0d0)) :: eigenvalues(1:sizeD)
complex(kind(0d0)), allocatable :: eigenvalues(:)
integer :: sizeD

#ifdef PARALLEL
sizeD=dimG*(global_num_sites+global_num_links+global_num_faces)
#else
sizeD=dimG*(num_sites+num_links+num_faces)
#endif
allocate(eigenvalues(1:sizeD))

call calc_eigenvalues_Dirac(eigenvalues,UMAT,PhiMat)

minimal=eigenvalues(1)
maximal=eigenvalues(sizeD)

end subroutine calc_min_and_max_of_eigenvalues_Dirac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Sf part 
subroutine partial_Sf(obs_site, obs_link, obs_face,Dirac,Dirac_inv)
implicit none

complex(kind(0d0)), intent(out) :: obs_site, obs_link, obs_face
complex(kind(0d0)), intent(in) :: Dirac(:,:), Dirac_inv(:,:)
integer :: s, t, l, l2, f, f2, a, b

!!!
obs_site=(0d0,0d0)
do s=1,num_sites
  do t=1,num_sites
    do a=1,dimG
      do b=1,dimG
        obs_site = obs_site + (0.5d0,0d0) &
          * Dirac( site_index(a,s,NMAT), site_index(b,t,NMAT) ) &
          * Dirac_inv( site_index(b,t,NMAT), site_index(a,s,NMAT) )
      enddo
    enddo
  enddo
enddo
!!!
obs_link=(0d0,0d0)
do l=1,num_links
  do t=1,num_sites
    do a=1,dimG
      do b=1,dimG
        obs_link = obs_link + (0.5d0,0d0) &
          * Dirac( link_index(a,l,NMAT,num_sites), site_index(b,t,NMAT) ) &
          * Dirac_inv( site_index(b,t,NMAT), link_index(a,l,NMAT,num_sites) )
        !!
        obs_link = obs_link + (0.5d0,0d0) &
          * Dirac( site_index(a,t,NMAT), link_index(b,l,NMAT,num_sites) ) &
          * Dirac_inv( link_index(b,l,NMAT,num_sites), site_index(a,t,NMAT) )
      enddo
    enddo
  enddo
enddo
do l=1,num_links
  do l2=1,num_links
    do a=1,dimG
      do b=1,dimG
        obs_link = obs_link + (0.5d0,0d0) &
          * Dirac( link_index(a,l,NMAT,num_sites), link_index(b,l2,NMAT,num_sites) ) &
          * Dirac_inv( link_index(b,l2,NMAT,num_sites), link_index(a,l,NMAT,num_sites) )
      enddo
    enddo
  enddo
enddo
!!!
obs_face=(0d0,0d0)
do f=1,num_faces
  do l=1,num_links
    do a=1,dimG
      do b=1,dimG
        obs_face = obs_face + (0.5d0,0d0) &
          * Dirac( face_index(a,f,NMAT,num_sites,num_links), link_index(b,l,NMAT,num_sites) ) &
          * Dirac_inv( link_index(b,l,NMAT,num_sites), face_index(a,f,NMAT,num_sites,num_links) )
        !!
        obs_face = obs_face + (0.5d0,0d0) &
          * Dirac( link_index(a,l,NMAT,num_sites), face_index(b,f,NMAT,num_sites,num_links) ) &
          * Dirac_inv( face_index(b,f,NMAT,num_sites,num_links), link_index(a,l,NMAT,num_sites) )
      enddo
    enddo
  enddo
enddo
do f=1,num_faces
  do f2=1,num_faces
    do a=1,dimG
      do b=1,dimG
        obs_face = obs_face + (0.5d0,0d0) &
          * Dirac( face_index(a,f,NMAT,num_sites,num_links), face_index(b,f2,NMAT,num_sites,num_links) ) &
          * Dirac_inv( face_index(b,f2,NMAT,num_sites,num_links), face_index(a,f,NMAT,num_sites,num_links) )
      enddo
    enddo
  enddo
enddo
end subroutine partial_Sf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Maximal eigenvalue of (D^dag D)
subroutine max_eigen_DdagD(eigen,Umat,PhiMat)
use Dirac_operator, only :prod_DdagD
implicit none

complex(kind(0d0)), intent(out) :: eigen
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) :: eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: chi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)) :: Xeta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: Xlambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: Xchi(1:NMAT,1:NMAT,1:num_faces)

complex(kind(0d0)) :: gauss(1:(NMAT*NMAT-1)*(num_sites+num_links+num_faces) ) 
double precision :: gauss2(1:2*(NMAT*NMAT-1)*(num_sites+num_links+num_faces) ) 

complex(kind(0d0)) :: val,tmp_val, num,denom,eigen_pre
integer :: i,j,s,l,f,info,ite

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! random eta/lambda/chi
call BoxMuller2(gauss2, (NMAT*NMAT-1)*(num_sites+num_links+num_faces) )

val=(0d0,0d0)
do i=1,(NMAT*NMAT-1)*(num_sites+num_links+num_faces)
  gauss(i)=(gauss2(2*i-1) + gauss2(2*i)*(0d0,1d0)) * dcmplx(dsqrt(0.5d0))
  val=val+gauss(i)*dconjg(gauss(i))
enddo
val=dcmplx(dsqrt( dble(val) ))
do i=1,(NMAT*NMAT-1)*(num_sites+num_links+num_faces)
  gauss(i)=gauss(i)/val
enddo
call vec_to_mat(eta(:,:,1:num_sites),lambda(:,:,1:num_links),chi(:,:,1:num_faces),gauss)
#ifdef PARALLEL
  call syncronize_sites(eta)
  call syncronize_links(lambda)
  call syncronize_faces(chi)
#endif

eigen=(0d0,0d0)
do ite=1,10000
  !!
  eigen_pre=eigen
  !!
  call prod_DdagD(Xeta,Xlambda,Xchi,eta,lambda,chi,Umat,PhiMat)
  call InnerProd(num, &
    Xeta, Xlambda, Xchi, &
    eta, lambda, chi)
  call InnerProd(denom, &
    eta, lambda, chi, &
    eta, lambda, chi)
  eigen=num/denom
  !!
  if( cdabs(eigen-eigen_pre) < epsilon) exit
  eta(:,:,1:num_sites)=Xeta
  lambda(:,:,1:num_links)=Xlambda
  chi(:,:,1:num_faces)=Xchi
  !! normalization
  tmp_val=(0d0,0d0)
  do s=1,num_sites
    do j=1,NMAT
      do i=1,NMAT
        tmp_val=tmp_val+eta(i,j,s)*dconjg( eta(i,j,s) )
      enddo
    enddo
  enddo
  do l=1,num_links
    do j=1,NMAT
      do i=1,NMAT
        tmp_val=tmp_val+lambda(i,j,l)*dconjg( lambda(i,j,l) )
      enddo
    enddo
  enddo
  do f=1,num_faces
    do j=1,NMAT
      do i=1,NMAT
        tmp_val=tmp_val+chi(i,j,f)*dconjg( chi(i,j,f) )
      enddo
    enddo
  enddo
  val=(0d0,0d0)
#ifdef PARALLEL
  call MPI_ALLREDUCE(tmp_val,val,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
#else
  val=tmp_val
#endif
  val=dcmplx(dsqrt(dble(val)))
  do s=1,num_sites
    do j=1,NMAT
      do i=1,NMAT
        eta(i,j,s)=eta(i,j,s)/val
      enddo
    enddo
  enddo
  do l=1,num_links
    do j=1,NMAT
      do i=1,NMAT
        lambda(i,j,l)=lambda(i,j,l)/val
      enddo
    enddo
  enddo
  do f=1,num_faces
    do j=1,NMAT
      do i=1,NMAT
        chi(i,j,f)=chi(i,j,f)/val
      enddo
    enddo
  enddo
#ifdef PARALLEL
  call syncronize_sites(eta)
  call syncronize_links(lambda)
  call syncronize_faces(chi)
#endif
enddo

end subroutine max_eigen_DdagD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Minimal eigenvalue of (D^dag D)
subroutine min_eigen_DdagD(eigen,Umat,PhiMat)
use Dirac_operator, only :prod_DdagD
implicit none

complex(kind(0d0)), intent(out) :: eigen
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) :: eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: chi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)) :: Xeta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: Xlambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: Xchi(1:NMAT,1:NMAT,1:num_faces)

complex(kind(0d0)) :: gauss(1:(NMAT*NMAT-1)*(num_sites+num_links+num_faces) ) 
double precision :: gauss2(1:2*(NMAT*NMAT-1)*(num_sites+num_links+num_faces) ) 

complex(kind(0d0)) :: val,tmp_val, num,denom,eigen_pre
integer :: i,j,s,l,f,info,ite,counter

counter=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! random eta/lambda/chi
1111 call BoxMuller2(gauss2, (NMAT*NMAT-1)*(num_sites+num_links+num_faces) )
val=(0d0,0d0)
do i=1,(NMAT*NMAT-1)*(num_sites+num_links+num_faces)
  gauss(i)=(gauss2(2*i-1) + gauss2(2*i)*(0d0,1d0)) * dcmplx(dsqrt(0.5d0))
  val=val+gauss(i)*dconjg(gauss(i))
enddo
val=dcmplx(dsqrt( dble(val) ))
do i=1,(NMAT*NMAT-1)*(num_sites+num_links+num_faces)
  gauss(i)=gauss(i)/val
enddo
call vec_to_mat(eta(:,:,1:num_sites),lambda(:,:,1:num_links),chi(:,:,1:num_faces),gauss)
#ifdef PARALLEL
  call syncronize_sites(eta)
  call syncronize_links(lambda)
  call syncronize_faces(chi)
#endif

eigen=(0d0,0d0)
do ite=1,10000
  !!
  eigen_pre=eigen
  !!
  call calc_DdagDinvF(Xeta,Xlambda,Xchi,eta,lambda,chi,Umat,Phimat,info)
  if( info == 1 .and. counter<5 ) then
    if( counter < 5 ) then
      counter=counter+1
      goto 1111
    else
      exit
    endif
  endif
  call InnerProd(num, &
    Xeta, Xlambda, Xchi, &
    eta, lambda, chi)
  call InnerProd(denom, &
    eta, lambda, chi, &
    eta, lambda, chi)
  eigen=num/denom
  !!
  if( cdabs(eigen-eigen_pre) < epsilon) then
    !if(MYRANK==0) write(*,*) cdabs(eigen-eigen_pre)
    exit
  endif

  eta(:,:,1:num_sites)=Xeta
  lambda(:,:,1:num_links)=Xlambda
  chi(:,:,1:num_faces)=Xchi
  !! normalization
  tmp_val=(0d0,0d0)
  do s=1,num_sites
    do j=1,NMAT
      do i=1,NMAT
        tmp_val=tmp_val+eta(i,j,s)*dconjg( eta(i,j,s) )
      enddo
    enddo
  enddo
  do l=1,num_links
    do j=1,NMAT
      do i=1,NMAT
        tmp_val=tmp_val+lambda(i,j,l)*dconjg( lambda(i,j,l) )
      enddo
    enddo
  enddo
  do f=1,num_faces
    do j=1,NMAT
      do i=1,NMAT
        tmp_val=tmp_val+chi(i,j,f)*dconjg( chi(i,j,f) )
      enddo
    enddo
  enddo
  val=(0d0,0d0)
#ifdef PARALLEL
  call MPI_ALLREDUCE(tmp_val,val,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
#else
  val=tmp_val
#endif
  val=dcmplx(dsqrt(dble(val)))
  do s=1,num_sites
    do j=1,NMAT
      do i=1,NMAT
        eta(i,j,s)=eta(i,j,s)/val
      enddo
    enddo
  enddo
  do l=1,num_links
    do j=1,NMAT
      do i=1,NMAT
        lambda(i,j,l)=lambda(i,j,l)/val
      enddo
    enddo
  enddo
  do f=1,num_faces
    do j=1,NMAT
      do i=1,NMAT
        chi(i,j,f)=chi(i,j,f)/val
      enddo
    enddo
  enddo
#ifdef PARALLEL
  call syncronize_sites(eta)
  call syncronize_links(lambda)
  call syncronize_faces(chi)
#endif
enddo
if( info==1 ) write(*,*) "Dirac is singular"
eigen=(1d0,0d0)/eigen


end subroutine min_eigen_DdagD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to take inner product for matrix value
!!  (V1,V2) = V1^\dagger_i V2_i
subroutine InnerProd(v1v2, &
   vec1_eta, vec1_lambda, vec1_chi, &
   vec2_eta, vec2_lambda, vec2_chi)
#ifdef PARALLEL
use parallel
#endif
implicit none

complex(kind(0d0)), intent(in) :: vec1_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(in) :: vec1_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: vec1_chi(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)), intent(in) :: vec2_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(in) :: vec2_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: vec2_chi(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)), intent(out) :: v1v2

integer :: s,l,f,i,j
complex(kind(0d0)) :: v1v2_tmp

v1v2=(0d0,0d0)
v1v2_tmp=(0d0,0d0)
do s=1,num_sites
  do j=1,NMAT
    do i=1,NMAT
      v1v2_tmp=v1v2_tmp+dconjg( vec1_eta(i,j,s) ) *vec2_eta(i,j,s) 
    enddo
  enddo
enddo
do l=1,num_links
  do j=1,NMAT
    do i=1,NMAT
      v1v2_tmp=v1v2_tmp+dconjg( vec1_lambda(i,j,l) ) *vec2_lambda(i,j,l) 
    enddo
  enddo
enddo
do f=1,num_faces
  do j=1,NMAT
    do i=1,NMAT
      v1v2_tmp=v1v2_tmp+dconjg( vec1_chi(i,j,f) ) *vec2_chi(i,j,f)
    enddo
  enddo
enddo
#ifdef PARALLEL
call MPI_REDUCE(v1v2_tmp,v1v2,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(v1v2,1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
#else
v1v2=v1v2_tmp
#endif 

end subroutine InnerProd



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to test Q-transformation of S
subroutine check_QS(Umat,PhiMat)
use matrix_functions, only : matrix_commutator, matrix_3_product
use Dirac_operator, only : Prod_Dirac
#ifdef PARALLEL
use parallel
#endif
implicit none

complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) :: Qeta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: Qlambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: Qchi(1:NMAT,1:NMAT,1:num_necessary_faces)

complex(kind(0d0)) :: tmp_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: tmp_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: tmp_chi(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: Omega(1:NMAT,1:NMAT,1:num_faces)

complex(kind(0d0)) :: Bforce_s(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: Bforce_l(1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)) :: Fforce_s(1:NMAT,1:NMAT,1:num_sites)
!complex(kind(0d0)) :: Fforce_l(1:NMAT,1:NMAT,1:num_links)

complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
double precision :: tmp,QS
integer :: info,s,l,f,i,j,triger

call make_Qfermion(Qeta,Qlambda,Qchi,Omega,Umat,PhiMat)
call make_bosonic_force_nomass(Bforce_s,Bforce_l,Umat,PhiMat)

call Prod_Dirac(tmp_eta,tmp_lambda,tmp_chi,Qeta,Qlambda,Qchi,UMAT,Phimat)
! Q^2 \Omega を care する
!do f=1,num_faces
!  call matrix_commutator(tmpmat,PhiMat(:,:,sites_in_f(f)%label_(1)),Omega(:,:,f))
!  tmp_chi(:,:,f)=tmp_chi(:,:,f)+(0d0,1d0)*dcmplx(beta_f(f))*tmpmat
!enddo

if(MYRANK==0) write(*,*) "# QS = 0 ?"
!write(*,*) tmp_chi
QS=0d0
tmp=0d0
do s=1,num_sites
  !tmp=0d0
  do i=1,NMAT
    do j=1,NMAT
      tmp_eta(i,j,s)=-tmp_eta(i,j,s)+dconjg(Bforce_s(j,i,s))
      tmp=tmp+dble( tmp_eta(i,j,s)*dconjg(tmp_eta(i,j,s)) )
    enddo
  enddo
  !write(*,*) "# site",global_site_of_local(s),tmp
enddo
call MPI_REDUCE(tmp,QS,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
if(MYRANK==0) write(*,*) "#   site:",QS

QS=0d0
tmp=0d0
do l=1,num_links
  !tmp=0d0
  do i=1,NMAT
    do j=1,NMAT
      tmp_lambda(i,j,l) = -tmp_lambda(i,j,l)+Bforce_l(i,j,l)
      tmp=tmp+dble( tmp_lambda(i,j,l)*dconjg(tmp_lambda(i,j,l)) )
    enddo
  enddo
  !write(*,*) "# link",global_link_of_local(l),tmp
enddo
call MPI_REDUCE(tmp,QS,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
if(MYRANK==0) write(*,*) "#   link:",QS

QS=0d0
tmp=0d0
do f=1,num_faces
  !tmp=0d0
  do i=1,NMAT
    do j=1,NMAT
      tmp_chi(i,j,f)=-tmp_chi(j,i,f)
      tmp=tmp+dble( tmp_chi(i,j,f)*dconjg(tmp_chi(i,j,f)) )
    enddo
  enddo
  !write(*,*) "# face",global_face_of_local(f),tmp
enddo
call MPI_REDUCE(tmp,QS,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
if(MYRANK==0) write(*,*) "#   face:",QS

end subroutine check_QS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to make Q\Psi
subroutine make_Qfermion(Qeta,Qlambda,Qchi,Omega,Umat,PhiMat)
use matrix_functions, only : matrix_commutator, matrix_3_product
#ifdef PARALLEL
use parallel
#endif
implicit none

complex(kind(0d0)), intent(out) :: Qeta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(out) :: Qlambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(out) :: Qchi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), intent(out) :: Omega(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
!complex(kind(0d0)), intent(in) :: PF_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)) :: comm(1:NMAT,1:NMAT)

integer :: s,l,f,i,j
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT)

do s=1,num_sites
  call matrix_commutator(Qeta(:,:,s),PhiMat(:,:,s),PhiMat(:,:,s),'N','C')
enddo

do l=1,num_links
  Qlambda(:,:,l)=(0d0,-1d0)*PhiMat(:,:,link_org(l))
  call matrix_3_product(Qlambda(:,:,l),&
    Umat(:,:,l),PhiMat(:,:,link_tip(l)),Umat(:,:,l),&
    'N','N','C',(0d0,1d0),'ADD')
enddo

do f=1,num_faces
  call Make_face_variable(Uf,f,UMAT)
  if(m_omega == 0) then 
    call Make_moment_map0(Omega(:,:,f),Uf)
  elseif(m_omega == -1) then
    call Make_moment_map_adm(Omega(:,:,f),Uf)
  endif
  Qchi(:,:,f)=(0d0,-0.5d0)*dcmplx(beta_f(f))*Omega(:,:,f)
enddo

#ifdef PARALLEL
call syncronize_sites(Qeta)
call syncronize_links(Qlambda)
call syncronize_faces(Qchi)
#endif

end subroutine make_Qfermion


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Cartan elements
!!   cartan(1:NMAT-1)
!!   MAT(1:NMAT,1:NMAT)
subroutine calc_cartan(cartan,MAT)
use SUN_generators, only : trace_MTa
implicit none

complex(kind(0d0)), intent(out) :: cartan(:)
complex(kind(0d0)), intent(in) :: MAT(:,:)
integer :: n,i

n=size(MAT,1)

do i=1,n-1
  call trace_MTa(cartan(i),MAT,i,n)
enddo

end subroutine calc_cartan




!end module observables
