!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! module to compute (MM^\dagger)^r by the rational approximation
module rational_algorithm
implicit none
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compute (MM^\dagger)^r . Bvec
subroutine calc_matrix_rational_power(rvec, Bvec, VSize, &
        NumVec, epsilon, MaxIte, info, CGite, alpha_r, beta_r,UMAT,PhiMat, ProdMat,num_sites, num_links, num_faces)
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

integer, intent(in) :: VSize, NumVec
complex(kind(0d0)), intent(inout) :: rvec(1:VSize)
complex(kind(0d0)), intent(in) :: Bvec(1:VSize)
double precision, intent(in) :: alpha_r(0:NumVec), beta_r(1:NumVec)
double precision, intent(in) :: epsilon
integer, intent(in) :: MaxIte
integer, intent(inout) :: CGite
integer, intent(inout) :: info 
integer, intent(in) :: num_sites, num_links, num_faces

complex(kind(0d0)) :: Xvec(1:VSize,1:NumVec)
integer :: r

!! 行列を構成する変数
complex(kind(0d0)), intent(in) :: UMAT(:,:,:), PhiMat(:,:,:) 


call mmBiCG( Xvec, Bvec, beta_r, VSize, NumVec, epsilon, MaxIte, info, CGite, UMAT, PhiMat, ProdMat,num_sites, num_links, num_faces )


rvec = dcmplx(alpha_r(0)) * Bvec
do r=1,NumVec
  rvec = rvec + dcmplx(alpha_r(r))*Xvec(:,r)
enddo
  
end subroutine calc_matrix_rational_power

subroutine mmBiCG_dev( Xvec, Bvec, sigma, VSize, NumVec, &
    epsilon, MaxIte, info, CGite, UMAT,PhiMat, ProdMat, &
    num_sites, num_links, num_faces )
use global_subroutines, only : mat_to_vec, vec_to_mat
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

integer :: NMAT

NMAT=size(UMAT,1)

allocate( X_eta(1:NMAT,1:NMAT,1:num_sites,1:NumVec) )
allocate( X_lambda(1:NMAT,1:NMAT,1:num_links,1:NumVec) )
allocate( X_chi(1:NMAT,1:NMAT,1:num_faces,1:NumVec) )

allocate( B_eta(1:NMAT,1:NMAT,1:num_sites) )
allocate( B_lambda(1:NMAT,1:NMAT,1:num_links) )
allocate( B_chi(1:NMAT,1:NMAT,1:num_faces) )

call vec_to_mat(B_eta,B_lambda,B_chi,Bvec)

!call mmBiCG_mat( &
!  X_eta,X_lambda,X_chi, &
!  B_eta,B_lambda,B_chi, &
!  sigma, epsilon, MaxIte, info, CGite, UMAT,PhiMat, ProdMat)

end subroutine mmBiCG_dev

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
subroutine mmBiCG( Xvec, Bvec, sigma, VSize, NumVec, &
    epsilon, MaxIte, info, CGite, UMAT,PhiMat, ProdMat, &
    num_sites, num_links, num_faces )
use global_subroutines, only : mat_to_vec, vec_to_mat
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
!! 行列を構成する変数
complex(kind(0d0)), intent(in) :: UMAT(:,:,:), PhiMat(:,:,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!== DEPEND ON ProdMat YOUR MADE ==
!complex(kind(0d0)), intent(in) :: MAT(1:VSize,1:VSize)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!! private variables
!complex(kind(0d0)) :: r_vec(1:VSize)
!complex(kind(0d0)) :: p_vec(1:VSize)
double precision :: zeta_k(1:NumVec), zeta_km1(1:NumVec)
double precision :: alpha_k, alpha_km1
double precision :: beta_k, beta_km1
!complex(kind(0d0)) :: psigma_vec(1:VSize,1:NumVec)

complex(kind(0d0)) :: r_vec(1:Vsize)
complex(kind(0d0)), allocatable :: r_eta(:,:,:)
complex(kind(0d0)), allocatable :: r_lambda(:,:,:)
complex(kind(0d0)), allocatable :: r_chi(:,:,:)
!!
complex(kind(0d0)) :: p_vec(1:Vsize)
complex(kind(0d0)), allocatable :: p_eta(:,:,:)
complex(kind(0d0)), allocatable :: p_lambda(:,:,:)
complex(kind(0d0)), allocatable :: p_chi(:,:,:)
!!
complex(kind(0d0)) :: psigma_vec(1:Vsize,1:NumVec)
complex(kind(0d0)), allocatable :: psigma_eta(:,:,:,:)
complex(kind(0d0)), allocatable :: psigma_lambda(:,:,:,:)
complex(kind(0d0)), allocatable :: psigma_chi(:,:,:,:)
!!
complex(kind(0d0)) :: tmp_vec(1:Vsize)
complex(kind(0d0)), allocatable :: tmp_eta(:,:,:)
complex(kind(0d0)), allocatable :: tmp_lambda(:,:,:)
complex(kind(0d0)), allocatable :: tmp_chi(:,:,:)
!!
complex(kind(0d0)), allocatable :: X_eta(:,:,:,:),X_lambda(:,:,:,:),X_chi(:,:,:,:)
!complex(kind(0d0)), allocatable :: B_eta(:,:,:),B_lambda(:,:,:),B_chi(:,:,:)

integer :: flag(NumVec), ff


!! for iterations
integer :: i,ite, s
double precision :: tmp_r1
complex(kind(0d0)) :: tmp_c1, tmp_c2, rkrk
integer :: NMAT

NMAT=size(UMAT,1)

!!
allocate( X_eta(1:NMAT,1:NMAT,1:num_sites,1:NumVec) )
allocate( r_eta(1:NMAT,1:NMAT,1:num_sites) )
allocate( p_eta(1:NMAT,1:NMAT,1:num_sites) )
allocate( tmp_eta(1:NMAT,1:NMAT,1:num_sites) )
allocate( psigma_eta(1:NMAT,1:NMAT,1:num_sites,1:NumVec) )
!!
allocate( X_lambda(1:NMAT,1:NMAT,1:num_links,1:NumVec) )
allocate( r_lambda(1:NMAT,1:NMAT,1:num_links) )
allocate( p_lambda(1:NMAT,1:NMAT,1:num_links) )
allocate( tmp_lambda(1:NMAT,1:NMAT,1:num_links) )
allocate( psigma_lambda(1:NMAT,1:NMAT,1:num_links,1:NumVec) )
!!
allocate( X_chi(1:NMAT,1:NMAT,1:num_faces,1:NumVec) )
allocate( r_chi(1:NMAT,1:NMAT,1:num_faces) )
allocate( p_chi(1:NMAT,1:NMAT,1:num_faces) )
allocate( tmp_chi(1:NMAT,1:NMAT,1:num_faces) )
allocate( psigma_chi(1:NMAT,1:NMAT,1:num_faces,1:NumVec) )



!write(*,*) Bvec
!! initialization
info = 0
flag = 0

Xvec=(0d0,0d0)
r_vec = Bvec
p_vec = Bvec
do s=1,NumVec
    psigma_vec(:,s)=Bvec
enddo

call vec_to_mat(r_eta(:,:,:),r_lambda(:,:,:),r_chi(:,:,:),r_vec)
do s=1,NumVec
  call vec_to_mat(psigma_eta(:,:,:,s),psigma_lambda(:,:,:,s),psigma_chi(:,:,:,s),psigma_vec(:,s))
call vec_to_mat(X_eta(:,:,:,s),X_lambda(:,:,:,s),X_chi(:,:,:,s),Xvec(:,s))
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
   r_eta, r_lambda, r_chi,&
   r_eta, r_lambda, r_chi)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!== DEPEND ON ProdMat YOU MADE ==
! call ProdMat(tmp_vec, p_vec, VSize,UMAT,PhiMat)
 call vec_to_mat(p_eta,p_lambda,p_chi,p_vec)
 call ProdMat(&
   tmp_eta,tmp_lambda,tmp_chi,&
   p_eta,p_lambda,p_chi,&
   UMAT,PhiMat)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !call InnerProd(tmp_c1, p_vec, tmp_vec, VSize)
 call InnerProd2(tmp_c1, &
   p_eta, p_lambda, p_chi,&
   tmp_eta, tmp_lambda, tmp_chi)
 alpha_k = rkrk / dble(tmp_c1)

! (2) update r_k --> r_{k+1}
  r_eta = r_eta - dcmplx(alpha_k) * tmp_eta
  r_lambda = r_lambda - dcmplx(alpha_k) * tmp_lambda
  r_chi = r_chi - dcmplx(alpha_k) * tmp_chi

! (3) construct beta_k
 call InnerProd2(tmp_c1, &
   r_eta, r_lambda, r_chi,&
   r_eta, r_lambda, r_chi)
 beta_k = dble(tmp_c1) / rkrk

! (4) update p_k --> p_{k+1}
  p_eta = r_eta + cmplx(beta_k) * p_eta
  p_lambda = r_lambda + cmplx(beta_k) * p_lambda
  p_chi = r_chi + cmplx(beta_k) * p_chi


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
! ここは要チェック
! このブロックを下のブロックに交換すると結果が変わる。
! do i=1, VSize
!     Xvec(i,s) = Xvec(i,s) + tmp_c1 * psigma_vec(i,s)
!     psigma_vec(i,s) = dcmplx(zeta_k(s))*r_vec(i) + tmp_c2*psigma_vec(i,s) 
! enddo
     X_eta(:,:,:,s) =  X_eta(:,:,:,s) + tmp_c1 * psigma_eta(:,:,:,s)
     X_lambda(:,:,:,s) = X_lambda(:,:,:,s) + tmp_c1 * psigma_lambda(:,:,:,s)
     X_chi(:,:,:,s) =  X_chi(:,:,:,s) + tmp_c1 * psigma_chi(:,:,:,s)
     !!
     psigma_eta(:,:,:,s) = cmplx(zeta_k(s))*r_eta(:,:,:) &
       + tmp_c2*psigma_eta(:,:,:,s) 
     psigma_lambda(:,:,:,s) = cmplx(zeta_k(s))*r_lambda(:,:,:) &
       + tmp_c2*psigma_lambda(:,:,:,s) 
     psigma_chi(:,:,:,s) = cmplx(zeta_k(s))*r_chi(:,:,:) &
       + tmp_c2*psigma_chi(:,:,:,s) 
   endif
 enddo
   
! conservation check
 ff=1
 call InnerProd2(tmp_c1, &
   r_eta, r_lambda, r_chi,&
   r_eta, r_lambda, r_chi)
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
end subroutine mmBiCG

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

NMAT=size(vec1_eta,1)
num_sites=size(vec1_eta,3)
num_links=size(vec1_lambda,3)
num_faces=size(vec1_chi,3)

v1v2=(0d0,0d0)
do s=1,num_sites
  do j=1,NMAT
    do i=1,NMAT
      v1v2=v1v2+conjg( vec1_eta(i,j,s) ) *vec2_eta(i,j,s) 
    enddo
  enddo
enddo
do l=1,num_links
  do j=1,NMAT
    do i=1,NMAT
      v1v2=v1v2+conjg( vec1_lambda(i,j,l) ) *vec2_lambda(i,j,l) 
    enddo
  enddo
enddo
do f=1,num_faces
  do j=1,NMAT
    do i=1,NMAT
      v1v2=v1v2+conjg( vec1_chi(i,j,f) ) *vec2_chi(i,j,f) 
    enddo
  enddo
enddo

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
