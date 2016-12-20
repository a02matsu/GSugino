!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! module to compute (MM^\dagger)^r by the rational approximation
module rational_algorithm
implicit none
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compute (MM^\dagger)^r . Bvec
subroutine calc_matrix_rational_power(rvec, Bvec, NMAT, &
        NumVec, epsilon, MaxIte, info, CGite, alpha_r, beta_r,UMAT,PhiMat, ProdMat)
implicit none
interface 
  subroutine ProdMat(R_VEC,VEC,MSIZE,UMAT,PhiMat)
    integer, intent(in) :: MSIZE
    complex(kind(0d0)), intent(in) :: VEC(1:MSIZE)
    complex(kind(0d0)), intent(inout) :: R_VEC(1:MSIZE)
    complex(kind(0d0)), intent(in) :: UMAT(:,:,:),PhiMat(:,:,:)
  end subroutine ProdMat
end interface

complex(kind(0d0)), intent(inout) :: rvec(1:NMAT)
integer, intent(in) :: NMAT, NumVec
complex(kind(0d0)), intent(in) :: Bvec(1:NMAT)
double precision, intent(in) :: alpha_r(0:NumVec), beta_r(1:NumVec)
double precision, intent(in) :: epsilon
integer, intent(in) :: MaxIte
integer, intent(inout) :: CGite
integer, intent(inout) :: info 

complex(kind(0d0)) :: Xvec(1:NMAT,1:NumVec)
integer :: r

!! 行列を構成する変数
complex(kind(0d0)), intent(in) :: UMAT(:,:,:), PhiMat(:,:,:) 


call mmBiCG( Xvec, Bvec, beta_r, NMAT, NumVec, epsilon, MaxIte, info, CGite, UMAT, PhiMat, ProdMat )

rvec = dcmplx(alpha_r(0)) * Bvec
do r=1,NumVec
  rvec = rvec + dcmplx(alpha_r(r))*Xvec(:,r)
enddo
  
end subroutine calc_matrix_rational_power


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
!!  mmBiCG( Xvec, MAT, Bvec, sigma, NMAT, NumVec, epsilon, MaxIte, info )
!!   Xvec(1:NMAT,1:NumVec) : double complex   :: solution 
!!   Bvec(1:NMAT)          : double complex   :: given vector 
!!   sigma(1:NumVec)       : double presision :: mass shift 
!!   NMAT                  : integer          :: size of vectors and MAT 
!!   NumVec                : integer          :: number of mass shifts 
!!   epsilon               : double presision :: error range 
!!   MaxIte                : integer          :: maximal iteration
!!   info                  : integer :: take 1 if iteration reaches MaxIte
!!   ProdMat : subroutine to compute MAT.VEC :: ProdMat(R_VEC,VEC,NMAT)
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
!!   ProdMat( result_vec, MAT, given_vec, NMAT )
!!    result_vec: MAT * v
!!    MAT: a given matrix
!!    given_vec: a vector v 
!!    NMAT: size of the vector
subroutine mmBiCG( Xvec, Bvec, sigma, NMAT, NumVec, epsilon, MaxIte, info, CGite, UMAT,PhiMat, ProdMat )
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

integer, intent(in) :: NMAT, NumVec
complex(kind(0d0)), intent(inout) :: Xvec(1:NMAT,1:NumVec)
complex(kind(0d0)), intent(in) :: Bvec(1:NMAT)
double precision, intent(in) :: sigma(1:NumVec)
double precision, intent(in) :: epsilon
integer, intent(in) :: MaxIte
integer, intent(inout) :: info 
integer, intent(inout) :: CGite
!! 行列を構成する変数
complex(kind(0d0)), intent(in) :: UMAT(:,:,:), PhiMat(:,:,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!== DEPEND ON ProdMat YOUR MADE ==
!complex(kind(0d0)), intent(in) :: MAT(1:NMAT,1:NMAT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!! private variables
complex(kind(0d0)) :: r_vec(1:NMAT)
complex(kind(0d0)) :: p_vec(1:NMAT)
double precision :: zeta_k(1:NumVec), zeta_km1(1:NumVec)
double precision :: alpha_k, alpha_km1
double precision :: beta_k, beta_km1
complex(kind(0d0)) :: psigma_vec(1:NMAT,1:NumVec)

integer :: flag(NumVec), ff


!! for iterations
integer :: i,ite, s
double precision :: tmp_r1
complex(kind(0d0)) :: tmp_c1, tmp_c2, rkrk
complex(kind(0d0)) :: tmp_vec(1:NMAT)




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
 call InnerProd(rkrk, r_vec, r_vec, NMAT) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!== DEPEND ON ProdMat YOU MADE ==
 call ProdMat(tmp_vec, p_vec, NMAT,UMAT,PhiMat)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 call InnerProd(tmp_c1, p_vec, tmp_vec, NMAT)
 alpha_k = rkrk / dble(tmp_c1)

! (2) update r_k --> r_{k+1}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!== DEPEND ON ProdMat YOU MADE ==
 call ProdMat(tmp_vec, p_vec, NMAT,UMAT,PhiMat)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 do i=1,NMAT
   r_vec(i) = r_vec(i) - dcmplx(alpha_k) * tmp_vec(i)
 enddo

! (3) construct beta_k
 call InnerProd(tmp_c1, r_vec, r_vec, NMAT) 
 beta_k = dble(tmp_c1) / rkrk

! (4) update p_k --> p_{k+1}
 do i=1,NMAT
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
 do i=1, NMAT
     Xvec(i,s) = Xvec(i,s) + tmp_c1 * psigma_vec(i,s)
     psigma_vec(i,s) = dcmplx(zeta_k(s))*r_vec(i) + tmp_c2*psigma_vec(i,s) 
 enddo
 endif
 !write(*,*) ite,s,zeta_k(s),zeta_km1(s)
 enddo
   
! conservation check
 ff=1
 call InnerProd( tmp_c1, r_vec, r_vec, NMAT)
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to take inner product
!!  (V1,V2) = V1^\dagger_i V2_i
subroutine  InnerProd(v1v2, vec1, vec2, NMAT)
implicit none

integer, intent(in) :: NMAT
complex(kind(0d0)), intent(in) :: vec1(1:NMAT), vec2(1:NMAT)
complex(kind(0d0)), intent(inout) :: v1v2

integer i

v1v2=(0d0,0d0)
do i=1,NMAT
    v1v2 = v1v2 + dconjg(vec1(i)) * vec2(i)
enddo
end subroutine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   ProdMat( result_vec, MAT, given_vec, NMAT )
!!    result_vec: MAT * v
!!    MAT: a given matrix
!!    given_vec: a vector v 
!!    NMAT: size of the vector
!subroutine ProdMat( result_vec, MAT, given_vec, NMAT )
!implicit none
!
!integer, intent(in) :: NMAT
!complex(kind(0d0)), intent(in) :: MAT(1:NMAT,1:NMAT), given_vec(1:NMAT)
!complex(kind(0d0)), intent(inout) :: result_vec(1:NMAT)
!
!integer i,j
!
!result_vec=(0d0,0d0)
!do i=1,NMAT
!do j=1,NMAT
    !result_vec(i)=result_vec(i)+MAT(i,j)*given_vec(j)
!enddo
!enddo
!
!end subroutine

end module rational_algorithm
