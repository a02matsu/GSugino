!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! module to compute (MM^\dagger)^r by the rational approximation
module rational_algorithm
implicit none
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compute (MM^\dagger)^r . Bvec
subroutine calc_matrix_rational_power(&
    r_eta, r_lambda, r_chi, & ! 1..num_sites,links,faces
    B_eta, B_lambda, B_chi, & ! 1..num_necessary_sites,links,faces
    epsilon, MaxIte, info, CGite, alpha_r, beta_r,UMAT,PhiMat, ProdMat)
!use global_subroutines, only : mat_to_vec, vec_to_mat
use global_subroutines, only :  site_abs,link_abs,face_abs
#ifdef PARALLEL
use parallel
!use global_subroutines, only : site_abs, face_abs, link_abs, stop_for_test
#endif
implicit none
interface 
!  subroutine ProdMat(R_VEC,VEC,MSIZE,UMAT,PhiMat)
!    integer, intent(in) :: MSIZE
!    complex(kind(0d0)), intent(in) :: VEC(1:MSIZE)
!    complex(kind(0d0)), intent(inout) :: R_VEC(1:MSIZE)
!    complex(kind(0d0)), intent(in) :: UMAT(:,:,:),PhiMat(:,:,:)
!  end subroutine ProdMat
  subroutine ProdMat(D_eta,D_lambda,D_chi,eta,lambda,chi,UMAT,PhiMat)
    complex(kind(0d0)), intent(out) :: D_eta(:,:,:),D_lambda(:,:,:),D_chi(:,:,:)
    complex(kind(0d0)), intent(in) :: eta(:,:,:),lambda(:,:,:),chi(:,:,:)
    complex(kind(0d0)), intent(in) :: UMAT(:,:,:),PhiMat(:,:,:)
  end subroutine ProdMat
end interface

complex(kind(0d0)), intent(out) :: r_eta(:,:,:) ! 1:num_sites
complex(kind(0d0)), intent(out) :: r_lambda(:,:,:) ! 1:num_links
complex(kind(0d0)), intent(out) :: r_chi(:,:,:) ! 1:num_faces
complex(kind(0d0)), intent(in) :: B_eta(:,:,:) ! 1:num_necessary_sites
complex(kind(0d0)), intent(in) :: B_lambda(:,:,:) ! 1:num_necessary_links
complex(kind(0d0)), intent(in) :: B_chi(:,:,:) ! 1:num_necessary_faces
double precision, intent(in) :: alpha_r(0:)
double precision, intent(in) :: beta_r(:) 
double precision, intent(in) :: epsilon
complex(kind(0d0)), intent(in) :: UMAT(:,:,:)
complex(kind(0d0)), intent(in) :: PhiMat(:,:,:) 
integer, intent(in) :: MaxIte
integer, intent(inout) :: CGite
integer, intent(inout) :: info 
integer :: num_sites, num_links, num_faces, NumVec, NMAT
integer :: num_necessary_sites, num_necessary_links, num_necessary_faces

complex(kind(0d0)), allocatable :: X_eta(:,:,:,:)
complex(kind(0d0)), allocatable :: X_lambda(:,:,:,:)
complex(kind(0d0)), allocatable :: X_chi(:,:,:,:)
integer :: r
double precision :: rtmp

NMAT=size(Umat,1)
NumVec=size(beta_r,1)
num_sites=size(r_eta,3)
num_links=size(r_lambda,3)
num_faces=size(r_chi,3)
num_necessary_sites=size(B_eta,3)
num_necessary_links=size(B_lambda,3)
num_necessary_faces=size(B_chi,3)

allocate( X_eta(1:NMAT,1:NMAT,1:num_sites,1:NumVec) )
allocate( X_lambda(1:NMAT,1:NMAT,1:num_links,1:NumVec) )
allocate( X_chi(1:NMAT,1:NMAT,1:num_faces,1:NumVec) )

!rtmp=site_abs(B_eta(:,:,1:num_sites))
!if( MYRANK == 0 ) then
!  write(*,'(a,E25.18)') "|B_eta|^2=",rtmp
!endif
!rtmp=link_abs(B_lambda(:,:,1:num_links))
!if( MYRANK == 0 ) then
!  write(*,'(a,E25.18)') "|B_lambda|^2=",rtmp
!endif
!rtmp=face_abs(B_chi(:,:,1:num_faces))
!if( MYRANK == 0 ) then
!  write(*,'(a,E25.18)') "|B_chi|^2=",rtmp
!endif




call mmBiCG( &
  X_eta, X_lambda, X_chi, &
  B_eta, B_lambda, B_chi, &
  beta_r, epsilon, MaxIte, info, CGite, UMAT, PhiMat, ProdMat)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! if CG is failed, we should reject
if( info==1 ) return
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

r_eta = dcmplx(alpha_r(0)) * B_eta(:,:,1:num_sites)
r_lambda = dcmplx(alpha_r(0)) * B_lambda(:,:,1:num_links)
r_chi = dcmplx(alpha_r(0)) * B_chi(:,:,1:num_faces)
do r=1,NumVec
  r_eta = r_eta + dcmplx(alpha_r(r))*X_eta(:,:,:,r)
  r_lambda = r_lambda + dcmplx(alpha_r(r))*X_lambda(:,:,:,r)
  r_chi = r_chi + dcmplx(alpha_r(r))*X_chi(:,:,:,r)
enddo
  
end subroutine calc_matrix_rational_power

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mmBiCG_vec( Xvec, Bvec, sigma, VSize, NumVec, &
    epsilon, MaxIte, info, CGite, UMAT,PhiMat, ProdMat, &
    num_sites, num_links, num_faces )
use global_subroutines, only : mat_to_vec, vec_to_mat
#ifdef PARALLEL
use parallel
#endif
implicit none
interface
!! R_VEC = MATRIX . VEC を計算するサブルーチン
!! ただし、MATRIX は明示的に指定しない。
!! 呼び出す側で定義する必要がある。
!  subroutine ProdMat(R_VEC,VEC,MSIZE,UMAT,PhiMat)
!    integer, intent(in) :: MSIZE
!    complex(kind(0d0)), intent(in) :: VEC(1:MSIZE)
!    complex(kind(0d0)), intent(inout) :: R_VEC(1:MSIZE)
!    complex(kind(0d0)), intent(in) :: UMAT(:,:,:),PhiMat(:,:,:)
!  end subroutine ProdMat
  subroutine ProdMat(D_eta,D_lambda,D_chi,eta,lambda,chi,UMAT,PhiMat)
    complex(kind(0d0)), intent(out) :: D_eta(:,:,:),D_lambda(:,:,:),D_chi(:,:,:)
    complex(kind(0d0)), intent(in) :: eta(:,:,:),lambda(:,:,:),chi(:,:,:)
    complex(kind(0d0)), intent(in) :: UMAT(:,:,:),PhiMat(:,:,:)
  end subroutine ProdMat
end interface

integer, intent(in) :: VSize, NumVec
complex(kind(0d0)), intent(inout) :: Xvec(1:VSize,1:NumVec)
complex(kind(0d0)), intent(in) :: Bvec(1:VSize)
double precision, intent(in) :: sigma(1:NumVec)
double precision, intent(in) :: epsilon
integer, intent(in) :: MaxIte
integer, intent(inout) :: info 
integer, intent(inout) :: CGite
integer, intent(in) :: num_sites, num_links, num_faces
complex(kind(0d0)), intent(in) :: UMAT(:,:,:), PhiMat(:,:,:)

complex(kind(0d0)), allocatable :: X_eta(:,:,:,:),X_lambda(:,:,:,:),X_chi(:,:,:,:)
complex(kind(0d0)), allocatable :: B_eta(:,:,:),B_lambda(:,:,:),B_chi(:,:,:)

integer :: NMAT,s 

NMAT=size(UMAT,1)

allocate( X_eta(1:NMAT,1:NMAT,1:num_sites,1:NumVec) )
allocate( X_lambda(1:NMAT,1:NMAT,1:num_links,1:NumVec) )
allocate( X_chi(1:NMAT,1:NMAT,1:num_faces,1:NumVec) )

allocate( B_eta(1:NMAT,1:NMAT,1:num_sites) )
allocate( B_lambda(1:NMAT,1:NMAT,1:num_links) )
allocate( B_chi(1:NMAT,1:NMAT,1:num_faces) )

call vec_to_mat(B_eta,B_lambda,B_chi,Bvec)

call mmBiCG( &
  X_eta,X_lambda,X_chi, &
  B_eta,B_lambda,B_chi, &
  sigma, epsilon, MaxIte, info, CGite, UMAT,PhiMat, ProdMat)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! if CG is failed, we should reject
if( info==1 ) return
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do s=1,NumVec
  call mat_to_vec(Xvec(:,s),X_eta(:,:,:,s),X_lambda(:,:,:,s),X_chi(:,:,:,s))
enddo

end subroutine mmBiCG_vec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Multi-shift CG Solver
!!
!!  MAT: NxN positive hermitian matrix
!!  b: N-vector
!! に対して、
!!  (Mat + \sigma) x = b
!! を解くサブルーチン
!!
!! ＜設計＞
!! 将来的にMPIでの並列化も視野に入れるので、
!! ベクトル演算は出来るだけサブルーチンにしたい。
!! 並列化の際は、サブルーチンをいじればいいようにする。
!! ver.00

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Multi-mass BiCG
!! solve (MAT + sigma) Xvec = Bvec
!! 
!! Usage:
!!  mmBiCG( Xvec, MAT, Bvec, sigma, VSize, NumVec, epsilon, MaxIte, info )
!!   Xvec(1:VSize,1:NumVec) : double complex   :: solution 
!!   Bvec(1:VSize)          : double complex   :: given vector 
!!   sigma(1:NumVec)       : double presision :: mass shift 
!!   VSize                  : integer          :: size of vectors and MAT 
!!   NumVec                : integer          :: number of mass shifts 
!!   epsilon               : double presision :: error range 
!!   MaxIte                : integer          :: maximal iteration
!!   info                  : integer :: take 1 if iteration reaches MaxIte
!!   ProdMat : subroutine to compute MAT.VEC :: ProdMat(R_VEC,VEC,VSize)
!!               R_VEC = "MATRIX" . VEC
!!             "MATRIX" must be defined globally
!!
!!  (1) "MAT" is only used in the external subroutine ProdMat, 
!!      so the structure of the data may be different depending 
!!      the subroutine. 
!!      For instance, "MAT" may be a set of array with different sizes. 
!!      One should rewrite it depending on "ProdMat".  
!! 
!! This subroutine calls external subroutine 
!!   ProdMat( result_vec, MAT, given_vec, VSize )
!!    result_vec: MAT * v
!!    MAT: a given matrix
!!    given_vec: a vector v 
!!    VSize: size of the vector
subroutine mmBiCG( &
    X_eta, X_lambda, X_chi, &
    B_eta, B_lambda, B_chi, &
    sigma, epsilon, MaxIte, info, CGite, UMAT,PhiMat, ProdMat)
!use global_subroutines, only : mat_to_vec, vec_to_mat
#ifdef PARALLEL
use parallel
use global_subroutines, only : syncronize_sites, syncronize_links, syncronize_faces, site_abs, link_abs, face_abs
#endif
implicit none
interface
!! R_VEC = MATRIX . VEC を計算するサブルーチン
!! ただし、MATRIX は明示的に指定しない。
!! 呼び出す側で定義する必要がある。
!  subroutine ProdMat(R_VEC,VEC,MSIZE,UMAT,PhiMat)
!    integer, intent(in) :: MSIZE
!    complex(kind(0d0)), intent(in) :: VEC(1:MSIZE)
!    complex(kind(0d0)), intent(inout) :: R_VEC(1:MSIZE)
!    complex(kind(0d0)), intent(in) :: UMAT(:,:,:),PhiMat(:,:,:)
!  end subroutine ProdMat
  subroutine ProdMat(D_eta,D_lambda,D_chi,eta,lambda,chi,UMAT,PhiMat)
    complex(kind(0d0)), intent(out) :: D_eta(:,:,:),D_lambda(:,:,:),D_chi(:,:,:)
    complex(kind(0d0)), intent(in) :: eta(:,:,:),lambda(:,:,:),chi(:,:,:)
    complex(kind(0d0)), intent(in) :: UMAT(:,:,:),PhiMat(:,:,:)
  end subroutine ProdMat
end interface

complex(kind(0d0)), intent(out) :: X_eta(:,:,:,:),X_lambda(:,:,:,:),X_chi(:,:,:,:)
complex(kind(0d0)), intent(in) :: B_eta(:,:,:),B_lambda(:,:,:),B_chi(:,:,:)
complex(kind(0d0)), intent(in) :: UMAT(:,:,:), PhiMat(:,:,:)
!complex(kind(0d0)), intent(inout) :: Xvec(1:VSize,1:NumVec)
!complex(kind(0d0)), intent(in) :: Bvec(1:VSize)
double precision, intent(in) :: sigma(:)
double precision, intent(in) :: epsilon
integer, intent(in) :: MaxIte
integer, intent(inout) :: info 
integer, intent(inout) :: CGite
integer :: num_sites, num_links, num_faces
integer :: num_necessary_sites, num_necessary_links, num_necessary_faces
integer :: NumVec
!! 行列を構成する変数
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!== DEPEND ON ProdMat YOUR MADE ==
!complex(kind(0d0)), intent(in) :: MAT(1:VSize,1:VSize)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!! private variables
!double precision :: zeta_k(1:NumVec), zeta_km1(1:NumVec)
double precision, allocatable :: zeta_k(:), zeta_km1(:)
double precision :: alpha_k, alpha_km1
double precision :: beta_k, beta_km1
double precision :: rtmp

complex(kind(0d0)), allocatable :: r_eta(:,:,:)
complex(kind(0d0)), allocatable :: r_lambda(:,:,:)
complex(kind(0d0)), allocatable :: r_chi(:,:,:)
!!
complex(kind(0d0)), allocatable :: p_eta(:,:,:)
complex(kind(0d0)), allocatable :: p_lambda(:,:,:)
complex(kind(0d0)), allocatable :: p_chi(:,:,:)
!!
complex(kind(0d0)), allocatable :: psigma_eta(:,:,:,:)
complex(kind(0d0)), allocatable :: psigma_lambda(:,:,:,:)
complex(kind(0d0)), allocatable :: psigma_chi(:,:,:,:)
!!
complex(kind(0d0)), allocatable :: tmp_eta(:,:,:)
complex(kind(0d0)), allocatable :: tmp_lambda(:,:,:)
complex(kind(0d0)), allocatable :: tmp_chi(:,:,:)
!!

!integer :: flag(NumVec), ff
integer, allocatable :: flag(:)
integer :: ff


!! for iterations
integer :: i,ite, s,j
double precision :: tmp_r1
complex(kind(0d0)) :: tmp_c1, tmp_c2, rkrk
integer :: NMAT


NMAT=size(UMAT,1)
num_sites=size(X_eta,3)
num_links=size(X_lambda,3)
num_faces=size(X_chi,3)
num_necessary_sites=size(B_eta,3)
num_necessary_links=size(B_lambda,3)
num_necessary_faces=size(B_chi,3)
NumVec=size(sigma,1)

allocate( zeta_k(1:NumVec), zeta_km1(1:NumVec) )
allocate( flag(1:NumVec) )
!!
allocate( r_eta(1:NMAT,1:NMAT,1:num_sites) )
allocate( p_eta(1:NMAT,1:NMAT,1:num_necessary_sites) )
allocate( tmp_eta(1:NMAT,1:NMAT,1:num_sites) )
allocate( psigma_eta(1:NMAT,1:NMAT,1:num_sites,1:NumVec) )
!!
allocate( r_lambda(1:NMAT,1:NMAT,1:num_links) )
allocate( p_lambda(1:NMAT,1:NMAT,1:num_necessary_links) )
allocate( tmp_lambda(1:NMAT,1:NMAT,1:num_links) )
allocate( psigma_lambda(1:NMAT,1:NMAT,1:num_links,1:NumVec) )
!!
allocate( r_chi(1:NMAT,1:NMAT,1:num_faces) )
allocate( p_chi(1:NMAT,1:NMAT,1:num_necessary_faces) )
allocate( tmp_chi(1:NMAT,1:NMAT,1:num_faces) )
allocate( psigma_chi(1:NMAT,1:NMAT,1:num_faces,1:NumVec) )



!write(*,*) Bvec
!! initialization
info = 0
flag = 0

!Xvec=(0d0,0d0)
X_eta=(0d0,0d0)
X_lambda=(0d0,0d0)
X_chi=(0d0,0d0)

!r_vec = Bvec
r_eta = B_eta(:,:,1:num_sites)
r_lambda = B_lambda(:,:,1:num_links)
r_chi = B_chi(:,:,1:num_faces)

!p_vec = Bvec
p_eta = B_eta
p_lambda = B_lambda
p_chi = B_chi

do s=1,NumVec
  !psigma_vec(:,s)=Bvec
  psigma_eta(:,:,:,s)=B_eta
  psigma_lambda(:,:,:,s)=B_lambda
  psigma_chi(:,:,:,s)=B_chi
enddo

zeta_k = 1d0
zeta_km1 = 1d0
alpha_km1 = 1d0
beta_km1 = 0d0


!! iteration start
do ite=1,MaxIte

! (1) construct \alpha_k
 ! rkrk = (r_k, r_k)
 call InnerProd2(rkrk, &
   r_eta(:,:,1:num_sites), r_lambda(:,:,1:num_links), r_chi(:,:,1:num_faces),&
   r_eta(:,:,1:num_sites), r_lambda(:,:,1:num_links), r_chi(:,:,1:num_faces))
!if( MYRANK == 0 ) then
!  write(*,*) "#########",ite,"##########"
!  write(*,'(a,E25.18)') "rkrk=",dble(rkrk)
!endif
!rtmp=site_abs(B_eta(:,:,1:num_sites))
!if( MYRANK == 0 ) then
!  !write(*,*) "#########",ite,"##########"
!  write(*,'(a,E25.18)') "|B_eta|^2=",rtmp
!endif
!rtmp=link_abs(B_lambda(:,:,1:num_links))
!if( MYRANK == 0 ) then
!  write(*,'(a,E25.18)') "|B_lambda|^2=",rtmp
!endif
!rtmp=face_abs(B_chi(:,:,1:num_faces))
!if( MYRANK == 0 ) then
!  write(*,'(a,E25.18)') "|B_chi|^2=",rtmp
!endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!== DEPEND ON ProdMat YOU MADE ==
! call ProdMat(tmp_vec, p_vec, VSize,UMAT,PhiMat)
!if(MYRANK==0) write(*,*) "## prod_ddagd in bicg ##"
 call ProdMat(&
   tmp_eta,tmp_lambda,tmp_chi,&
   p_eta,p_lambda,p_chi,&
   UMAT,PhiMat)

!rtmp=site_abs(PhiMat(:,:,1:num_sites))
!if( MYRANK == 0 ) then
!  write(*,*) "|PhiMat|^2=",rtmp
!endif
!rtmp=link_abs(UMat(:,:,1:num_links))
!if( MYRANK == 0 ) then
!  write(*,*) "|UMat|^2=",rtmp
!endif
!rtmp=site_abs(p_eta(:,:,1:num_sites))
!if( MYRANK == 0 ) then
!  write(*,*) "|p_eta|^2=",rtmp
!endif
!rtmp=link_abs(p_lambda(:,:,1:num_links))
!if( MYRANK == 0 ) then
!  write(*,*) "|p_lambda|^2=",rtmp
!endif
!rtmp=face_abs(p_chi(:,:,1:num_faces))
!if( MYRANK == 0 ) then
!  write(*,*) "|p_chi|^2=",rtmp
!endif
!rtmp=site_abs(tmp_eta(:,:,1:num_sites))
!if( MYRANK == 0 ) then
!  write(*,*) "|tmp_eta|^2=",rtmp
!endif
!rtmp=link_abs(tmp_lambda(:,:,1:num_links))
!if( MYRANK == 0 ) then
!  write(*,*) "|tmp_lambda|^2=",rtmp
!endif
!rtmp=face_abs(tmp_chi(:,:,1:num_faces))
!if( MYRANK == 0 ) then
!  write(*,*) "|tmp_chi|^2=",rtmp
!endif
!call stop_for_test
 !call syncronize_sites(tmp_eta)
 !call syncronize_links(tmp_lambda)
 !call syncronize_faces(tmp_chi)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !call InnerProd(tmp_c1, p_vec, tmp_vec, VSize)
 call InnerProd2(tmp_c1, &
   p_eta(:,:,1:num_sites), p_lambda(:,:,1:num_links), p_chi(:,:,1:num_faces),&
   tmp_eta(:,:,1:num_sites), tmp_lambda(:,:,1:num_links), tmp_chi(:,:,1:num_faces))
 !alpha_k = dble(rkrk) / dble(tmp_c1)
 alpha_k = dexp(dlog(dble(rkrk))-dlog(dble(tmp_c1)))

! (2) update r_k --> r_{k+1}
  r_eta = r_eta - dcmplx(alpha_k) * tmp_eta
  r_lambda = r_lambda - dcmplx(alpha_k) * tmp_lambda
  r_chi = r_chi - dcmplx(alpha_k) * tmp_chi

!rtmp=site_abs(r_eta(:,:,1:num_sites))
!if( MYRANK == 0 ) then
!  !write(*,*) "#########",ite,"##########"
!  write(*,*) "alpha_k=",alpha_k
!  write(*,*) "|r_eta|^2=",rtmp
!endif
!rtmp=link_abs(r_lambda(:,:,1:num_links))
!if( MYRANK == 0 ) then
!  write(*,*) "|r_lambda|^2=",rtmp
!endif
!rtmp=face_abs(r_chi(:,:,1:num_faces))
!if( MYRANK == 0 ) then
!  write(*,*) "|r_chi|^2=",rtmp
!endif

! (3) construct beta_k
 call InnerProd2(tmp_c1, &
   r_eta, r_lambda, r_chi,&
   r_eta, r_lambda, r_chi)
 !beta_k = dble(tmp_c1) / dble(rkrk)
 beta_k = dexp(dlog(dble(tmp_c1)) - dlog(dble(rkrk)))

! (4) update p_k --> p_{k+1}
  p_eta(:,:,1:num_sites) = r_eta(:,:,1:num_sites) + dcmplx(beta_k) * p_eta(:,:,1:num_sites)
  p_lambda(:,:,1:num_links) = r_lambda(:,:,1:num_links) + dcmplx(beta_k) * p_lambda(:,:,1:num_links)
  p_chi(:,:,1:num_faces) = r_chi(:,:,1:num_faces) + dcmplx(beta_k) * p_chi(:,:,1:num_faces)
#ifdef PARALLEL
  call syncronize_sites(p_eta)
  call syncronize_faces(p_chi)
  call syncronize_links(p_lambda)
#endif


! (5) update zeta^{sigma}_{k} --> zeta^{sigma}_{k+1}
 do s=1,NumVec
   if ( flag(s) .ne. 1 ) then
     !tmp_r1 = zeta_k(s)*zeta_km1(s)*alpha_km1 &
     !    / ( zeta_km1(s)*alpha_km1*(1d0+sigma(s)*alpha_k) &
     !        + ( zeta_km1(s) - zeta_k(s) )*alpha_k*beta_km1 )
     tmp_r1 = dexp(dlog(zeta_k(s)*zeta_km1(s)*alpha_km1) &
         - dlog( zeta_km1(s)*alpha_km1*(1d0+sigma(s)*alpha_k) &
             + ( zeta_km1(s) - zeta_k(s) )*alpha_k*beta_km1 ) )
     zeta_km1(s)=zeta_k(s)
     zeta_k(s)=tmp_r1
   endif
 enddo

! (6) update Xvec and psigma_vec
 do s=1, NumVec
   if ( flag(s) .ne. 1 ) then
     tmp_c1 = dcmplx(zeta_k(s)/zeta_km1(s)*alpha_k)
     tmp_c2 = dcmplx((zeta_k(s)*zeta_k(s))/(zeta_km1(s)*zeta_km1(s))*beta_k)
! do i=1, VSize
!     Xvec(i,s) = Xvec(i,s) + tmp_c1 * psigma_vec(i,s)
!     psigma_vec(i,s) = dcmplx(zeta_k(s))*r_vec(i) + tmp_c2*psigma_vec(i,s) 
! enddo
     X_eta(:,:,:,s) =  X_eta(:,:,:,s) + tmp_c1 * psigma_eta(:,:,:,s)
     X_lambda(:,:,:,s) = X_lambda(:,:,:,s) + tmp_c1 * psigma_lambda(:,:,:,s)
     X_chi(:,:,:,s) =  X_chi(:,:,:,s) + tmp_c1 * psigma_chi(:,:,:,s)
     !!
     psigma_eta(:,:,:,s) = dcmplx(zeta_k(s))*r_eta(:,:,:) &
       + tmp_c2*psigma_eta(:,:,:,s) 
     psigma_lambda(:,:,:,s) = dcmplx(zeta_k(s))*r_lambda(:,:,:) &
       + tmp_c2*psigma_lambda(:,:,:,s) 
     psigma_chi(:,:,:,s) = dcmplx(zeta_k(s))*r_chi(:,:,:) &
       + tmp_c2*psigma_chi(:,:,:,s) 
     !!
   endif
 enddo
   
! conservation check
 ff=1
 call InnerProd2(tmp_c1, &
   r_eta, r_lambda, r_chi,&
   r_eta, r_lambda, r_chi)
 tmp_r1 = dsqrt( dble(tmp_c1) )
 !write(*,*) "##########",ite,"###########"
 do s=1,NumVec
 !write(*,*) s,dabs(zeta_k(s)*tmp_r1)
     if ( dabs(zeta_k(s)*tmp_r1) < epsilon ) then
         flag(s) = 1
     endif
     ff=ff*flag(s)
 enddo

 
 if ( ff==1 ) then 
     !write(*,*) "iteration=",ite
     info = 0
     CGite=ite
     !do s=1,NumVec
       !call mat_to_vec(Xvec(:,s),X_eta(:,:,:,s),X_lambda(:,:,:,s),X_chi(:,:,:,s))
     !enddo
     return
 else
     alpha_km1 = alpha_k
     beta_km1 = beta_k
 endif
enddo

!write(*,*) "iteration=",ite
 
 CGite=Maxite
 info = 1
 write(0,*) "# [mmBiCG] The iteration reached to the maximum."
! do s=1,NumVec
!   call mat_to_vec(Xvec(:,s),X_eta(:,:,:,s),X_lambda(:,:,:,s),X_chi(:,:,:,s))
! enddo

return
end subroutine mmBiCG

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Multi-shift CG Solver
!!
!!  MAT: NxN positive hermitian matrix
!!  b: N-vectorcall syncronize_sites(tmp2_eta)
!! に対して、
!!  (Mat + \sigma) x = b
!! を解くサブルーチン
!!
!! ＜設計＞
!! 将来的にMPIでの並列化も視野に入れるので、
!! ベクトル演算は出来るだけサブルーチンにしたい。
!! 並列化の際は、サブルーチンをいじればいいようにする。
!! ver.00

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Multi-mass BiCG
!! solve (MAT + sigma) Xvec = Bvec
!! 
!! Usage:
!!  mmBiCG( Xvec, MAT, Bvec, sigma, VSize, NumVec, epsilon, MaxIte, info )
!!   Xvec(1:VSize,1:NumVec) : double complex   :: solution 
!!   Bvec(1:VSize)          : double complex   :: given vector 
!!   sigma(1:NumVec)       : double presision :: mass shift 
!!   VSize                  : integer          :: size of vectors and MAT 
!!   NumVec                : integer          :: number of mass shifts 
!!   epsilon               : double presision :: error range 
!!   MaxIte                : integer          :: maximal iteration
!!   info                  : integer :: take 1 if iteration reaches MaxIte
!!   ProdMat : subroutine to compute MAT.VEC :: ProdMat(R_VEC,VEC,VSize)
!!               R_VEC = "MATRIX" . VEC
!!             "MATRIX" must be defined globally
!!
!!  (1) "MAT" is only used in the external subroutine ProdMat, 
!!      so the structure of the data may be different depending 
!!      the subroutine. 
!!      For instance, "MAT" may be a set of array with different sizes. 
!!      One should rewrite it depending on "ProdMat".  
!! 
!! This subroutine calls external subroutine 
!!   ProdMat( result_vec, MAT, given_vec, VSize )
!!    result_vec: MAT * v
!!    MAT: a given matrix
!!    given_vec: a vector v 
!!    VSize: size of the vector
subroutine mmBiCG_org( Xvec, Bvec, sigma, VSize, NumVec, epsilon, MaxIte, info, CGite, UMAT,PhiMat, ProdMat )
implicit none
interface
!! R_VEC = MATRIX . VEC を計算するサブルーチン
!! ただし、MATRIX は明示的に指定しない。
!! 呼び出す側で定義する必要がある。
  subroutine ProdMat(R_VEC,VEC,MSIZE,UMAT,PhiMat)
    integer, intent(in) :: MSIZE
    complex(kind(0d0)), intent(in) :: VEC(1:MSIZE)
    complex(kind(0d0)), intent(inout) :: R_VEC(1:MSIZE)
    complex(kind(0d0)), intent(in) :: UMAT(:,:,:),PhiMat(:,:,:)
  end subroutine ProdMat
end interface

integer, intent(in) :: VSize, NumVec
complex(kind(0d0)), intent(inout) :: Xvec(1:VSize,1:NumVec)
complex(kind(0d0)), intent(in) :: Bvec(1:VSize)
double precision, intent(in) :: sigma(1:NumVec)
double precision, intent(in) :: epsilon
integer, intent(in) :: MaxIte
integer, intent(inout) :: info 
integer, intent(inout) :: CGite
!! 行列を構成する変数
complex(kind(0d0)), intent(in) :: UMAT(:,:,:), PhiMat(:,:,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!== DEPEND ON ProdMat YOUR MADE ==
!complex(kind(0d0)), intent(in) :: MAT(1:VSize,1:VSize)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!! private variables
complex(kind(0d0)) :: r_vec(1:VSize)
complex(kind(0d0)) :: p_vec(1:VSize)
double precision :: zeta_k(1:NumVec), zeta_km1(1:NumVec)
double precision :: alpha_k, alpha_km1
double precision :: beta_k, beta_km1
complex(kind(0d0)) :: psigma_vec(1:VSize,1:NumVec)

integer :: flag(NumVec), ff


!! for iterations
integer :: i,ite, s
double precision :: tmp_r1
complex(kind(0d0)) :: tmp_c1, tmp_c2, rkrk
complex(kind(0d0)) :: tmp_vec(1:VSize)




!write(*,*) Bvec
!! initialization
info = 0
flag = 0

Xvec=(0d0,0d0)
r_vec = Bvec
p_vec = Bvec
zeta_k = 1d0
do s=1,NumVec
    psigma_vec(:,s)=Bvec
enddo

zeta_km1 = 1d0
alpha_km1 = 1d0
beta_km1 = 0d0

!! iteration start
do ite=1,MaxIte

! (1) construct \alpha_k
 ! rkrk = (r_k, r_k)
 call InnerProd(rkrk, r_vec, r_vec, VSize) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!== DEPEND ON ProdMat YOU MADE ==
 call ProdMat(tmp_vec, p_vec, VSize,UMAT,PhiMat)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 call InnerProd(tmp_c1, p_vec, tmp_vec, VSize)
 alpha_k = rkrk / dble(tmp_c1)

! (2) update r_k --> r_{k+1}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!== DEPEND ON ProdMat YOU MADE ==
 call ProdMat(tmp_vec, p_vec, VSize,UMAT,PhiMat)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 do i=1,VSize
   r_vec(i) = r_vec(i) - dcmplx(alpha_k) * tmp_vec(i)
 enddo

! (3) construct beta_k
 call InnerProd(tmp_c1, r_vec, r_vec, VSize) 
 beta_k = dble(tmp_c1) / rkrk

! (4) update p_k --> p_{k+1}
 do i=1,VSize
   p_vec(i) = r_vec(i) + dcmplx(beta_k) * p_vec(i)
 enddo
! (5) update zeta^{sigma}_{k} --> zeta^{sigma}_{k+1}
 do s=1,NumVec
   if ( flag(s) .ne. 1 ) then
   tmp_r1 = zeta_k(s)*zeta_km1(s)*alpha_km1 &
       / ( zeta_km1(s)*alpha_km1*(1d0+sigma(s)*alpha_k) &
           + ( zeta_km1(s) - zeta_k(s) )*alpha_k*beta_km1 )
   zeta_km1(s)=zeta_k(s)
   zeta_k(s)=tmp_r1
 endif
 enddo

! (6) update Xvec and psigma_vec
 do s=1, NumVec
   if ( flag(s) .ne. 1 ) then
     tmp_c1 = dcmplx(zeta_k(s)/zeta_km1(s)*alpha_k)
     tmp_c2 = dcmplx((zeta_k(s)*zeta_k(s))/(zeta_km1(s)*zeta_km1(s))*beta_k)
 do i=1, VSize
     Xvec(i,s) = Xvec(i,s) + tmp_c1 * psigma_vec(i,s)
     psigma_vec(i,s) = dcmplx(zeta_k(s))*r_vec(i) + tmp_c2*psigma_vec(i,s) 
 enddo
 endif
 !write(*,*) ite,s,zeta_k(s),zeta_km1(s)
 enddo
   
! conservation check
 ff=1
 call InnerProd( tmp_c1, r_vec, r_vec, VSize)
 tmp_r1 = dsqrt( dble(tmp_c1) )
 do s=1,NumVec
     if ( dabs(zeta_k(s)*tmp_r1) < epsilon ) then
         flag(s) = 1
     endif
     ff=ff*flag(s)
 enddo

 
 if ( ff==1 ) then 
     !write(*,*) "iteration=",ite
     info = 0
     CGite=ite
     return
 else
     alpha_km1 = alpha_k
     beta_km1 = beta_k
 endif
enddo

!write(*,*) "iteration=",ite
 
 CGite=Maxite
 info = 1

return
end subroutine mmBiCG_org

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to take inner product
!!  (V1,V2) = V1^\dagger_i V2_i
subroutine  InnerProd(v1v2, vec1, vec2, VSize)
implicit none

integer, intent(in) :: VSize
complex(kind(0d0)), intent(in) :: vec1(1:VSize), vec2(1:VSize)
complex(kind(0d0)), intent(inout) :: v1v2

integer i

v1v2=(0d0,0d0)
do i=1,VSize
    v1v2 = v1v2 + dconjg(vec1(i)) * vec2(i)
enddo
end subroutine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   ProdMat( result_vec, MAT, given_vec, VSize )
!!    result_ve
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to take inner product for matrix value
!!  (V1,V2) = V1^\dagger_i V2_i
subroutine InnerProd2(v1v2, &
   vec1_eta, vec1_lambda, vec1_chi, &
   vec2_eta, vec2_lambda, vec2_chi)
#ifdef PARALLEL
use parallel
use global_subroutines, only : syncronize_sites, syncronize_links, syncronize_faces
#endif
implicit none

integer :: NMAT,num_sites,num_links,num_faces
complex(kind(0d0)), intent(in) :: vec1_eta(:,:,:)
complex(kind(0d0)), intent(in) :: vec1_lambda(:,:,:)
complex(kind(0d0)), intent(in) :: vec1_chi(:,:,:)
complex(kind(0d0)), intent(in) :: vec2_eta(:,:,:)
complex(kind(0d0)), intent(in) :: vec2_lambda(:,:,:)
complex(kind(0d0)), intent(in) :: vec2_chi(:,:,:)
complex(kind(0d0)), intent(out) :: v1v2

integer :: s,l,f,i,j
complex(kind(0d0)) :: v1v2_tmp

NMAT=size(vec1_eta,1)
num_sites=size(vec1_eta,3)
num_links=size(vec1_lambda,3)
num_faces=size(vec1_chi,3)

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
#endif 

end subroutine InnerProd2

!!    MAT: a given matrix
!!    given_vec: a vector v 
!!    VSize: size of the vector
!subroutine ProdMat( result_vec, MAT, given_vec, VSize )
!implicit none
!
!integer, intent(in) :: VSize
!complex(kind(0d0)), intent(in) :: MAT(1:VSize,1:VSize), given_vec(1:VSize)
!complex(kind(0d0)), intent(inout) :: result_vec(1:VSize)
!
!integer i,j
!
!result_vec=(0d0,0d0)
!do i=1,VSize
!do j=1,VSize
    !result_vec(i)=result_vec(i)+MAT(i,j)*given_vec(j)
!enddo
!enddo
!
!end subroutine

end module rational_algorithm
