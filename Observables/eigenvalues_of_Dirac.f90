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


