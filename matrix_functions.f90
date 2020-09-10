!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! module including matrix funcitons
module matrix_functions
implicit none

contains

!***********************************************************
!***********************************************************
!  make hermitian transpose
SUBROUTINE Hermitian_Conjugate(MATDAG,MAT)
implicit none

complex(kind(0d0)), intent(in) :: MAT(:,:)
complex(kind(0d0)), intent(out) :: MATDAG(:,:)

integer :: NMAT,i,j

NMAT=size(MAT,1)
do i=1,NMAT
  do j=1,NMAT
    MATDAG(i,j)=dconjg(MAT(j,i))
  enddo
enddo


end SUBROUTINE Hermitian_Conjugate


!***********************************************************
!***********************************************************
!  Trace of matrix  
SUBROUTINE MATRIX_TRACE(trace,MAT)
implicit none

complex(kind(0d0)), intent(in) :: MAT(:,:)
complex(kind(0d0)), intent(out) :: trace

integer :: NMAT,i,j

NMAT=size(MAT,1)
trace=(0d0,0d0)
do i=1,NMAT
  trace=trace+MAT(i,i)
enddo

end SUBROUTINE MATRIX_TRACE
  
!***********************************************************
!***********************************************************
!  Trace of matrix product
SUBROUTINE TRACE_MM(trace,MAT1,MAT2)
implicit none

complex(kind(0d0)), intent(in) :: MAT1(:,:)
complex(kind(0d0)), intent(in) :: MAT2(:,:)
complex(kind(0d0)), intent(out) :: trace

integer :: NMAT,i,j

NMAT=size(MAT1,1)
trace=(0d0,0d0)
do i=1,NMAT
  do j=1,NMAT
    trace=trace+MAT1(i,j)*MAT2(j,i)
  enddo
enddo

end SUBROUTINE TRACE_MM

!***********************************************************
!***********************************************************
!   Log of a matrix : P -> Log(P)
SUBROUTINE MATRIX_LOG(MATLOG,MAT)

  implicit none 

  complex(kind(0d0)), intent(in) :: MAT(:,:)
  complex(kind(0d0)), intent(out) :: MATLOG(:,:)

  integer MatSize

  integer i,j,k
  
  complex(kind(0d0)), allocatable :: eig_log(:)
  !double complex eig_log(1:MatSize)
  !lapack
  character jobvl,jobvr
  integer lda,ldvl,ldvr,info,lwork
  complex(kind(0d0)), allocatable :: w(:),VL(:),VR(:,:), work(:),rwork(:)
  !double complex w(1:MatSize),VL(1:MatSize),VR(1:MatSize,1:MATSIZE),&
       !work(1:2*MatSize),rwork(1:2*MatSize)

  MatSize=size(MAT,1)
  allocate( &
    eig_log(1:MatSize), &
    w(1:MatSize), &
    VL(1:MatSize), &
    VR(1:MatSize,1:MATSIZE), &
    work(1:2*MatSize), &
    rwork(1:2*MatSize) )


  jobvl='N'
  jobvr='V'
  lda=MatSize
  ldvl=1
  ldvr=MatSize
  lwork=2*MatSize

  !MAT2=MAT

  call ZGEEV(jobvl,jobvr,MatSize,MAT,lda,w,vl,ldvl,vr,ldvr,work,lwork,rwork,info)
  ! in this case w(i) is real. 
  do i=1,MatSize
     !write(*,*)w(i)
     eig_log(i)=dcmplx(dlog(dble(w(i))))
     !write(*,*)eig_log(i) 
 end do
  MATLOG=(0d0,0d0)
  do i=1,MatSize
     do j=1,MatSize
        do k=1,MatSize
           MATLOG(i,j)=MATLOG(i,j)+VR(i,k)*eig_log(k)*dconjg(VR(j,k))
       end do
     end do
  end do

  return

END SUBROUTINE MATRIX_LOG

!***********************************************************
!***********************************************************
!   exp of a matrix : P -> Exp(P)
SUBROUTINE MATRIX_EXP(MATEXP,MAT)

  implicit none 

  complex(kind(0d0)), intent(in) :: MAT(:,:)
  complex(kind(0d0)), intent(out) :: MATEXP(:,:)
  integer MatSize

  integer i,j,k
  complex(kind(0d0)), allocatable :: eig_exp(:)
  !lapack
  character jobvl,jobvr
  integer lda,ldvl,ldvr,info,lwork
  complex(kind(0d0)), allocatable :: w(:),VL(:),VR(:,:),work(:),rwork(:)
  !double complex eig_exp(1:MatSize)
  !lapack
  !character jobvl,jobvr
  !integer lda,ldvl,ldvr,info,lwork
  !double complex w(1:MatSize),VL(1:MatSize),VR(1:MatSize,1:MATSIZE),&
       !work(1:2*MatSize),rwork(1:2*MatSize)

  MatSize=size(MAT,1)
  allocate( &
    eig_exp(1:MatSize), &
    w(1:MatSize), &
    VL(1:MatSize), &
    VR(1:MatSize,1:MATSIZE), &
    work(1:2*MatSize), &
    rwork(1:2*MatSize) )

  jobvl='N'
  jobvr='V'
  lda=MatSize
  ldvl=1
  ldvr=MatSize
  lwork=2*MatSize

  !MAT2=MAT

  call ZGEEV(jobvl,jobvr,MatSize,MAT,lda,w,vl,ldvl,vr,ldvr,work,lwork,rwork,info)
  do i=1,MatSize
     !write(*,*)w(i)
     eig_exp(i)=cdexp(w(i)) !dcmplx(dlog(dble(w(i))))
     !write(*,*)eig_exp(i) 
 end do
  MATEXP=(0d0,0d0)
  do i=1,MatSize
     do j=1,MatSize
        do k=1,MatSize
           MATEXP(i,j)=MATEXP(i,j)+VR(i,k)*eig_exp(k)*dconjg(VR(j,k))
       end do
     end do
  end do

  return

END SUBROUTINE MATRIX_EXP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! norm of matrix
!! || A || = \sqrt{ A.A^\dag }
!subroutine matrix_norm(MAT,NMAT,NORM)
subroutine matrix_norm(NORM,MAT)
implicit none

double precision, intent(out) :: NORM
complex(kind(0d0)), intent(in) :: MAT(:,:)
integer :: NMAT
integer i,j,k

NMAT=size(MAT,1)
NORM=0d0
do j=1,NMAT
  do i=1,NMAT
    norm=norm+dble(MAT(i,j)*dconjg(MAT(i,j)))
  enddo
enddo
!norm=dsqrt(norm/dble(NMAT))
norm=dsqrt(norm)
end subroutine matrix_norm



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! make unit matrix with a given size
subroutine make_unit_matrix(UNITMAT)
implicit none

complex(kind(0d0)), intent(inout) :: UNITMAT(:,:)
integer :: NMAT
integer i

NMAT=size(UNITMAT,1)
UNITMAT=(0d0,0d0)
do i=1,NMAT
    UNITMAT(i,i)=(1d0,0d0)
enddo

end subroutine make_unit_matrix


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MAT^m 
!subroutine matrix_power(NMAT,MAT,m,pMAT)
subroutine matrix_power(pMAT,MAT,m)
implicit none

complex(kind(0d0)), intent(in) :: MAT(:,:)
integer, intent(in) :: m
complex(kind(0d0)), intent(out) :: pMAT(:,:)
complex(kind(0d0)), allocatable :: tmpMAT(:,:)
integer :: NMAT
integer i

NMAT=size(MAT,1)
allocate( tmpMAT(1:NMAT,1:NMAT) )
if( m < 0 ) then 
  write(*,*) "input positive power"
  stop
elseif( m==0 ) then
  call make_unit_matrix(pMAT)
  return
elseif( m==1 ) then
  pMAT=MAT
  return
else
  pMAT=MAT
  do i=2,m
    tmpMAT=pMAT
    call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
      tmpMAT,NMAT,&
      MAT,NMAT,&
      (0d0,0d0),pMAT,NMAT)
  enddo
endif


end subroutine matrix_power

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! (MAT^\dagger)^m 
subroutine matrix_dagger_power(pMAT,MAT,m)
implicit none

complex(kind(0d0)), intent(in) :: MAT(:,:)
integer, intent(in) :: m
complex(kind(0d0)), intent(out) :: pMAT(:,:)
complex(kind(0d0)), allocatable :: tmpMAT(:,:)
integer :: NMAT
integer i,j

NMAT=size(MAT,1)
allocate( tmpMAT(1:NMAT,1:NMAT) )
if( m < 0 ) then 
  write(*,*) "input positive power"
  stop
elseif( m==0 ) then
  call make_unit_matrix(pMAT)
  return
elseif( m==1 ) then
  do i=1,NMAT
    do j=1,NMAT
      pMAT(i,j)=dconjg(MAT(j,i))
    enddo
  enddo
  return
else
  do i=1,NMAT
    do j=1,NMAT
      pMAT(i,j)=dconjg(MAT(j,i))
    enddo
  enddo
  do i=2,m
    tmpMAT=pMAT
    call ZGEMM('N','C',NMAT,NMAT,NMAT,(1d0,0d0), &
      tmpMAT,NMAT,&
      MAT,NMAT,&
      (0d0,0d0),pMAT,NMAT)
  enddo
endif
end subroutine matrix_dagger_power

!***********************************************************
!***********************************************************
!  (MAT)**r (r: real)
SUBROUTINE MATRIX_RATIONAL_POWER(M_POWER, MAT, r)
implicit none 

double precision, intent(in) :: r
complex(kind(0d0)), intent(in) :: MAT(:,:)
complex(kind(0d0)), intent(inout) :: M_Power(:,:)
integer :: MatSize

  integer i,j,k
  double complex, allocatable :: eig_power(:)
  !lapack
  character jobvl,jobvr
  integer lda,ldvl,ldvr,info,lwork
  !double complex w(1:MatSize),VL(1:MatSize),VR(1:MatSize,1:MATSIZE),&
       !work(1:2*MatSize),rwork(1:2*MatSize)
  complex(kind(0d0)), allocatable :: w(:),VL(:),VR(:,:),&
       work(:),rwork(:)

  MatSize=size(MAT,1)
  allocate(eig_power(1:MatSize),&
       w(1:MatSize),VL(1:MatSize),VR(1:MatSize,1:MATSIZE),&
       work(1:2*MatSize),rwork(1:2*MatSize))

  jobvl='N'
  jobvr='V'
  lda=MatSize
  ldvl=1
  ldvr=MatSize
  lwork=2*MatSize

  !MAT2=MAT

  call ZGEEV(jobvl,jobvr,MatSize,MAT,lda,w,vl,ldvl,vr,ldvr,work,lwork,rwork,info)
  do i=1,MatSize
     eig_power(i)=w(i)**r
 end do
  M_Power=(0d0,0d0)
  do i=1,MatSize
     do j=1,MatSize
        do k=1,MatSize
           M_Power(i,j)=M_Power(i,j)+VR(i,k)*eig_power(k)*dconjg(VR(j,k))
       end do
     end do
  end do

  return

END SUBROUTINE MATRIX_RATIONAL_POWER

!***********************************************************
!***********************************************************
! output -> (MAT)^{-1}
!SUBROUTINE Matrix_Inverse(NMAT,MAT)
SUBROUTINE Matrix_Inverse(MAT)

  implicit none

  complex(kind(0d0)), intent(in) ::  MAT(:,:)
  integer NMAT
  integer lwork,info 
  integer, allocatable :: ipiv(:)
  !work(lwork);lwork must be equal to or larger than 2*(NMAT).    
  !doublecomplex WORK(16*(NMAT)) 
  complex(kind(0d0)), allocatable :: WORK(:)

  NMAT=size(MAT,1)
  lwork=16*(NMAT)
  allocate( ipiv(1:NMAT) )
  allocate( WORK( 16*NMAT ) )

  call zgetrf(NMAT,NMAT,MAT,NMAT,ipiv,info)
!  write(*,*)"info=",info
  call zgetri(NMAT,MAT,NMAT,ipiv,work,lwork,info)
!  write(*,*)"info=",info
  return

END SUBROUTINE Matrix_Inverse


!***********************************************************
!***********************************************************
! output -> (MAT)^{-1}
! MAT must be hermitian matrix
SUBROUTINE HermitianMatrix_Inverse(MAT)

  implicit none

  complex(kind(0d0)), intent(in) :: MAT(:,:)
  integer, allocatable :: ipiv(:)
  integer :: info
  !work(lwork);lwork must be equal to or larger than 2*(NMAT).    
  !complex(kind(0d0)) :: WORK(16*(NMAT)),work2(NMAT)
  complex(kind(0d0)), allocatable :: WORK(:),work2(:)
  integer NMAT
  integer :: lwork

  NMAT=size(MAT,1)
  lwork=16*(NMAT)
  allocate( ipiv(1:NMAT) )
  allocate( WORK(16*(NMAT)),work2(NMAT) )

  call zhetrf('U',NMAT,MAT,NMAT,ipiv,work,lwork,info)
!  write(*,*)"info=",info
  call zhetri('L',NMAT,MAT,NMAT,ipiv,work2,info)
!  write(*,*)"info=",info
  return

END SUBROUTINE HermitianMatrix_Inverse



!***********************************************************
!***********************************************************
! output -> (MAT)^{-1}
SUBROUTINE Matrix_Determinant(logDet,argDet,MAT)
implicit none

double precision, parameter :: PI=dacos(-1d0)
double precision, intent(out) :: logDet
complex(kind(0d0)), intent(out) :: argDet
complex(kind(0d0)), intent(in) ::  MAT(:,:)
integer ::  MatSize
double precision :: phase


integer i,j,k
!lapack
character jobvl,jobvr
integer lda,ldvl,ldvr,info,lwork
complex(kind(0d0)), allocatable :: w(:),VL(:),VR(:,:),&
     work(:),rwork(:)

MatSize=size(MAT,1)
allocate( w(1:MatSize),VL(1:MatSize),VR(1:MatSize,1:MATSIZE),&
     work(1:2*MatSize),rwork(1:2*MatSize) )
jobvl='N'
jobvr='V'
lda=MatSize
ldvl=1
ldvr=MatSize
lwork=2*MatSize


call ZGEEV(jobvl,jobvr,MatSize,MAT,lda,w,vl,ldvl,vr,ldvr,work,lwork,rwork,info)
logDet=0d0
phase=0d0
do i=1,MatSize
    !write(*,*) i,w(i)
  !logr = logr + 0.5d0 * dlog( dble( w(i) * dconjg(w(i)) ) )
  !phase = phase + datan( dble(w(i))/dble( (0d0,-1d0)*w(i) ) ) 
  logdet = logDet + dlog( cdabs(w(i)) ) 
  phase = phase + datan2( dble((0d0,-1d0)*w(i)), dble(w(i)) )  
  !if ( phase >= 2*PI ) then 
    !phase = phase - 2*PI
  !elseif (phase < 0) then 
    !phase = phase + 2*PI
  !endif 
enddo

argDet = dcmplx( dcos(phase) ) + (0d0,1d0)*dcmplx( dsin(phase) ) 

!Det= dcmplx(dexp(logr) * dcos(phase)) &
    !+ (0d0,1d0) * dcmplx(dexp(logr) * dsin(phase))

end subroutine Matrix_Determinant

!***********************************************************
!***********************************************************
! prod = alpha*MAT1.MAT2 
!  or
! prod = alpha*MAT1.MAT2 + prod  with ADD=A
!  * the size of comm, MAT1 and MAT2 must be the same.
!  * beta is potional and the default is (0d0,0d0) 
!  * conj1 and conj2 are optional and 
!    conj1 = N or C or T (default: N)
!    conj2 = N or C or T (default: N)
SUBROUTINE Matrix_Product(prod,MAT1,MAT2,char1,char2,alpha_input,ADD)
implicit none

complex(kind(0d0)), intent(in) :: MAT1(:,:), MAT2(:,:)
complex(kind(0d0)), intent(inout) :: prod(:,:)
character, optional :: char1,char2
character(3), optional :: ADD
complex(kind(0d0)), optional :: alpha_input
character :: C1,C2
complex(kind(0d0)) :: beta,alpha

integer :: NMAT

NMAT=size(MAT1,1)

if (present(char1)) then
  C1=char1
else
  C1='N'
endif

if (present(char2)) then
  C2=char2
else
  C2='N'
endif

alpha=(1d0,0d0)
beta=(0d0,0d0)
if (present(alpha_input)) then
  alpha=alpha_input
endif

if (present(ADD)) then
  if( ADD == 'ADD' .or. ADD == 'A' ) then 
    beta=(1d0,0d0)
  endif
endif

call ZGEMM(C1,C2,NMAT,NMAT,NMAT,alpha, &
  MAT1, NMAT, &
  MAT2, NMAT, &
  beta, prod, NMAT)

end subroutine Matrix_Product

!***********************************************************
!***********************************************************
! prod = alpha*MAT1.MAT2 
!  or
! prod = alpha*MAT1.MAT2 + prod  with ADD=A
!  * the size of comm, MAT1 and MAT2 must be the same.
!  * beta is potional and the default is (0d0,0d0) 
!  * conj1 and conj2 are optional and 
!    conj1 = N or C or T (default: N)
!    conj2 = N or C or T (default: N)
SUBROUTINE Matrix_Product2(prod,MAT1,MAT2,NMAT,char1,char2,alpha_input,ADD)
implicit none


integer, intent(in) :: NMAT
complex(kind(0d0)), intent(in) :: MAT1(NMAT,NMAT), MAT2(NMAT,NMAT)
complex(kind(0d0)), intent(inout) :: prod(NMAT,NMAT)
character, optional :: char1,char2
character(3), optional :: ADD
complex(kind(0d0)), optional :: alpha_input
character :: C1,C2
complex(kind(0d0)) :: beta,alpha


if (present(char1)) then
  C1=char1
else
  C1='N'
endif

if (present(char2)) then
  C2=char2
else
  C2='N'
endif

alpha=(1d0,0d0)
beta=(0d0,0d0)
if (present(alpha_input)) then
  alpha=alpha_input
endif

if (present(ADD)) then
  if( ADD == 'ADD' .or. ADD == 'A' ) then 
    beta=(1d0,0d0)
  endif
endif

call ZGEMM(C1,C2,NMAT,NMAT,NMAT,alpha, &
  MAT1, NMAT, &
  MAT2, NMAT, &
  beta, prod, NMAT)

end subroutine Matrix_Product2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! product of 3 matrices
! prod = alpha * MAT1.MAT2.Mat3 
!   or
! prod = alpha * MAT1.MAT2.Mat3 + prod for ADD="A"
!  * the size of comm, MAT1, MAT2 and MAT3 must be the same.
!  * beta is potional and the default is (0d0,0d0) 
!  * conj1 and conj2 and conj3 are optional and 
!    conj1 = N or C or T (default: N)
!    conj2 = N or C or T (default: N)
!    conj3 = N or C or T (default: N)
subroutine matrix_3_product(prod,mat1,mat2,mat3,char1,char2,char3,alpha_input,ADD)
implicit none

complex(kind(0d0)), intent(in) :: MAT1(:,:), MAT2(:,:), MAT3(:,:)
complex(kind(0d0)), intent(inout) :: prod(:,:)
complex(kind(0d0)), allocatable :: tmp(:,:)
character, optional :: char1,char2,char3
character(3), optional :: ADD
complex(kind(0d0)), optional :: alpha_input
character :: C1,C2,C3
complex(kind(0d0)) :: beta,alpha

integer :: NMAT

NMAT=size(MAT1,1)
allocate( tmp(1:NMAT,1:NMAT) )

if (present(char1)) then
  C1=char1
else
  C1='N'
endif

if (present(char2)) then
  C2=char2
else
  C2='N'
endif

if (present(char3)) then
  C3=char3
else
  C3='N'
endif

alpha=(1d0,0d0)
if (present(alpha_input)) then
  alpha=alpha_input
endif

beta=(0d0,0d0)
if (present(ADD)) then
  if (ADD == "ADD" .or. ADD == 'A') beta=(1d0,0d0)
endif

call ZGEMM(C1,C2,NMAT,NMAT,NMAT,(1d0,0d0), &
  MAT1, NMAT, &
  MAT2, NMAT, &
  (0d0,0d0), tmp, NMAT)

call ZGEMM('N',C3,NMAT,NMAT,NMAT,alpha, &
  tmp, NMAT, &
  MAT3, NMAT, &
  beta, prod, NMAT)

end subroutine Matrix_3_Product

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! product of 4 matrices
! prod = alpha * MAT1.MAT2.Mat3.Mat4
!   or
! prod = alpha * MAT1.MAT2.Mat3.Mat4 + prod for ADD="A"
!  * the size of comm, MAT1, MAT2, Mat3 and MAT4 must be the same.
!  * beta is potional and the default is (0d0,0d0) 
!  * conj1 and conj2 and conj3 and conj4 are optional and 
!    conj1 = N or C or T (default: N)
!    conj2 = N or C or T (default: N)
!    conj3 = N or C or T (default: N)
!    conj4 = N or C or T (default: N)
subroutine matrix_4_product(prod,mat1,mat2,mat3,mat4,char1,char2,char3,char4,alpha_input,ADD)
implicit none

complex(kind(0d0)), intent(in) :: MAT1(:,:), MAT2(:,:), MAT3(:,:),MAT4(:,:)
complex(kind(0d0)), intent(inout) :: prod(:,:)
complex(kind(0d0)), allocatable :: tmp(:,:),tmp2(:,:)
character, optional :: char1,char2,char3,char4
character(3), optional :: ADD
complex(kind(0d0)), optional :: alpha_input
character :: C1,C2,C3,C4
complex(kind(0d0)) :: beta,alpha

integer :: NMAT

NMAT=size(MAT1,1)
allocate( tmp(1:NMAT,1:NMAT) )
allocate( tmp2(1:NMAT,1:NMAT) )

if (present(char1)) then
  C1=char1
else
  C1='N'
endif

if (present(char2)) then
  C2=char2
else
  C2='N'
endif

if (present(char3)) then
  C3=char3
else
  C3='N'
endif

if (present(char4)) then
  C4=char4
else
  C4='N'
endif

alpha=(1d0,0d0)
if (present(alpha_input)) then
  alpha=alpha_input
endif

beta=(0d0,0d0)
if (present(ADD)) then
  if (ADD == "ADD" .or. ADD == 'A') beta=(1d0,0d0)
endif

call ZGEMM(C1,C2,NMAT,NMAT,NMAT,(1d0,0d0), &
  MAT1, NMAT, &
  MAT2, NMAT, &
  (0d0,0d0), tmp, NMAT)

call ZGEMM('N',C3,NMAT,NMAT,NMAT,(1d0,0d0), &
  tmp, NMAT, &
  MAT3, NMAT, &
  (0d0,0d0), tmp2, NMAT)

call ZGEMM('N',C4,NMAT,NMAT,NMAT,alpha, &
  tmp2, NMAT, &
  MAT4, NMAT, &
  beta, prod, NMAT)

end subroutine Matrix_4_Product


!***********************************************************
!***********************************************************
! comm = [MAT1, MAT2]
!  the size of comm, MAT1 and MAT2 must be the same
!  conj1 and conj2 are optional and 
!    conj1 = N or C or T (default: N)
!    conj2 = N or C or T (default: N)
SUBROUTINE Matrix_Commutator(comm,MAT1,MAT2,char1,char2)
implicit none

complex(kind(0d0)), intent(in) :: MAT1(:,:), MAT2(:,:)
complex(kind(0d0)), intent(inout) :: comm(:,:)
character, optional :: char1,char2
character :: C1,C2 

integer :: NMAT

NMAT=size(MAT1,1)

if (present(char1)) then
  C1=char1
else
  C1='N'
endif

if (present(char2)) then
  C2=char2
else
  C2='N'
endif

call ZGEMM(C1,C2,NMAT,NMAT,NMAT,(1d0,0d0), &
  MAT1, NMAT, &
  MAT2, NMAT, &
  (0d0,0d0), comm, NMAT)

call ZGEMM(C2,C1,NMAT,NMAT,NMAT,(-1d0,0d0), &
  MAT2, NMAT, &
  MAT1, NMAT, &
  (1d0,0d0), comm, NMAT)

end subroutine Matrix_Commutator

!***********************************************************
!***********************************************************
! anticomm = {MAT1, MAT2}
!  the size of comm, MAT1 and MAT2 must be the same
!  conj1 and conj2 are optional and 
!    conj1 = N or C or T (default: N)
!    conj2 = N or C or T (default: N)
SUBROUTINE Matrix_antiCommutator(acomm,MAT1,MAT2,char1,char2)
implicit none

complex(kind(0d0)), intent(in) :: MAT1(:,:), MAT2(:,:)
complex(kind(0d0)), intent(inout) :: acomm(:,:)
character, optional :: char1,char2
character :: C1,C2 

integer :: NMAT

NMAT=size(MAT1,1)

if (present(char1)) then
  C1=char1
else
  C1='N'
endif

if (present(char2)) then
  C2=char2
else
  C2='N'
endif

call ZGEMM(C1,C2,NMAT,NMAT,NMAT,(1d0,0d0), &
  MAT1, NMAT, &
  MAT2, NMAT, &
  (0d0,0d0), acomm, NMAT)

call ZGEMM(C2,C1,NMAT,NMAT,NMAT,(1d0,0d0), &
  MAT2, NMAT, &
  MAT1, NMAT, &
  (1d0,0d0), acomm, NMAT)

end subroutine Matrix_AntiCommutator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Matrix_Eigenvalues
subroutine Matrix_Eigenvalues(eigen,MAT)
implicit none

complex(kind(0d0)), intent(in) :: MAT(:,:)
complex(kind(0d0)), intent(out) :: eigen(:)

integer :: NMAT, i, j
complex(kind(0d0)) :: tako1

!lapack
character jobvl,jobvr
integer lda,ldvl,ldvr,info,lwork
complex(kind(0d0)), allocatable :: VL(:),VR(:,:),work(:),rwork(:) 
!double complex w(1:MatSize),VL(1:MatSize),VR(1:MatSize,1:MATSIZE),&
!     work(1:2*MatSize),rwork(1:2*MatSize)

NMAT=size(MAT,1)
if ( size(eigen,1) .ne. NMAT ) then
  write(*,*) "The size of the first array must be the same with that of the imput matrix."
  stop
endif
allocate( VL(1:NMAT) )
allocate( VR(1:NMAT,1:NMAT) )
allocate( work(1:2*NMAT) ) 
allocate( rwork(1:2*NMAT) )

jobvl='N'
jobvr='V'
lda=NMAT
ldvl=1
ldvr=NMAT
lwork=2*NMAT

call ZGEEV(jobvl,jobvr,NMAT,MAT,lda,eigen,vl,ldvl,vr,ldvr,work,lwork,rwork,info)

! sort the eigenvalues
do i=1,NMAT
 do j=i+1,NMAT
 tako1 = eigen(i)
  if(cdabs(eigen(j)).LT.cdabs(eigen(i))) then 
    eigen(i) = eigen(j)
    eigen(j) = tako1
  endif
 enddo
enddo


end subroutine Matrix_Eigenvalues

!CC    CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC    After C code by Simon Catterall
!CC    Phys. Rev. D 68, 014503 (2003) 
!CC
!CC    The meaning of rows has been interchanged with columns
!CC    in the comments because of the translation from C->Fortran
!CC    for more efficient programming
!CC
!CC    K. Anagnostopoulos, NTUA Feb 2007
!CC
!CC    changed by S.Matsuura Sep 2015
!CC    CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!C     If Phase pfaffian is 0.0 then the Pfaffian is zero:
subroutine PfaffianLog(MAT,logPfaffian,phasePfaffian)
  implicit none
  !complex(kind(0d0)) M(N,N)
  complex(kind(0d0)), intent(in) :: MAT(:,:)
  double precision, intent(out)  ::   logPfaffian
  complex(kind(0d0)), intent(out) :: phasePfaffian

  integer N
  complex(kind(0d0)), allocatable :: M(:,:)
  
  integer i,j,k,jpiv,interchange
  double precision pivot
  complex(kind(0d0)) scale,f,phasePf,mphase,tmppf
  complex(kind(0d0)), allocatable :: dum(:)
  double precision logPf,mmod
  
  interchange=1
  N=size( MAT,1 )
  allocate( M(N,N) )
  allocate( dum(N) )
  M=MAT

  !Loop ovel all rows in steps of 2
  do i=1,N-1,2
     !first row i
     !find column whose ith component is biggest to use as pivot
     pivot=CDABS(M(i+1,i))
     jpiv=i+1
     do j=i+2,N
        if(CDABS(M(j,i)).GT.pivot)then
           pivot=CDABS(M(j,i))
           jpiv=j
        endif
     enddo

     !interchange col(i+1) with col(jpiv)
     do j=1,N
        dum(j)=M(i+1,j)
     enddo
     do j=1,N
        M(i+1 ,j)=M(jpiv,j)
        M(jpiv,j)=dum(j)
     enddo

     !interchange row(i+1) with row(jpiv)
     do j=1,N
        dum(j)=M(j,i+1)
     enddo
     do j=1,N
        M(j,i+1 )=M(j,jpiv)
        M(j,jpiv)=dum(j)
     enddo

     if(jpiv .NE. i+1) interchange = -interchange

     !using this,zero progressively elements of row M(j,i),j=i+2,N
     do j=i+2,N
        scale=M(j,i)/M(i+1,i)
        do k=1,N
           M(j,k)=M(j,k)-scale*M(i+1,k)
        enddo
        !zero out elements along corresponding column M(i,j) too
        do k=1,N
           M(k,j)=M(k,j)-scale*M(k,i+1)
        enddo
     enddo!do j=i+2,N

     !next row i+1
     !using this,zero progressively elements of row M(j,i),j=i+2,N
     do j=i+2,N
        scale=M(j,i+1)/M(i,i+1)

        do k=1,N
           M(j,k)=M(j,k)-scale*M(i,k)
        enddo
        !zero out elements along corresponding column too
        do k=1,N
           M(k,j)=M(k,j)-scale*M(k,i)
        enddo
     enddo!do j=i+2,N

  enddo !do i=1,N-1,2,Loop ovel all rows in steps of 2

  phasePf= DCMPLX(1.0D0,0.0D0)
  logPf  = 0.0D0
  !c      tmppf  = (0.0D0,0.0D0)
  do i=1,N,2
     mmod=CDABS(M(i+1,i))
     if(mmod .le. 0.0D0)then
        phasePfaffian=DCMPLX(0.0D0    ,0.0D0)
        logPfaffian  =0.0D0
        RETURN
     endif
     mphase  = M(i+1,i) / mmod
     phasePf = phasePf  * mphase
     logPf   = logPf    + DLOG(mmod)
     !c       tmppf   = tmppf + DCMPLX(DLOG(CDABS(M(i+1,i))), 
     !c     $      DATAN2(DIMAG(M(i+1,i)),DREAL(M(i+1,i))))
  enddo

  !since we interchanged rows with columns and Pf(A^T)=(-1)^(N/2)
  !we have to change the sign appropriately:
  if(MOD(N/2,2).EQ.1)  interchange=-interchange
  phasePf = DBLE(interchange)*phasePf

  logPfaffian   = logPf
  phasePfaffian = phasePf
  !c      logPfaffian   = DREAL(tmppf)
  !c      phasePfaffian = DCMPLX(DCOS(DIMAG(tmppf)),DSIN(DIMAG(tmppf)))
  !c      phasePfaffian = DBLE(interchange)*phasePfaffian
  !c      print *,'++++++++++++++++++++++++++++++'
  !c      print *,'logPfaffian   = ',logPfaffian
  !c      print *,'phasePfaffian = ',phasePfaffian
  !c      print *,'++++++++++++++++++++++++++++++'
end subroutine PfaffianLog

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! abs_pf = |Pfaffian|
!! Pfaffian/|Pfaffian| = exp( i arg_pf ) -Pi< arg_pf < Pi
!subroutine CalcPfaffian(abs_pf,arg_pf,MAT)
subroutine CalcPfaffian(abs_pf,phase_pf,MAT)
implicit none

real(8), intent(out) :: abs_pf!, arg_pf
complex(kind(0d0)), intent(out) :: phase_pf
complex(kind(0d0)), intent(in) :: MAT(:,:)

integer :: NMAT
integer, allocatable :: IWORK(:)
integer  LWORK, LDA,NB, INFO
complex(kind(0d0)), allocatable :: WORK(:)
real(8), allocatable :: RWORK(:)
real(8) LogPf
!complex(kind(0d0)) phase

NMAT=size(MAT,1)

NB = NMAT*NMAT
LWORK = NMAT*NB
allocate ( IWORK(1:NMAT) )
allocate ( RWORK(1:NMAT-1) )
allocate ( WORK(LWORK) )
LDA = NMAT

call ZSKPFA('U','P',NMAT,MAT,LDA,LogPf,phase_pf,&
IWORK,WORK,LWORK,RWORK,INFO)

abs_pf=dexp(LogPf)
!arg_pf=datan2( dble((0d0,-1d0)*phase), dble(phase) )
end subroutine CalcPfaffian


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine make_matrix_traceless(MAT)
implicit none

complex(kind(0d0)), intent(inout) :: MAT(:,:)
integer :: NMAT
integer :: i
complex(kind(0d0)) :: trace

NMAT=size(MAT,1)

trace=(0d0,0d0)
do i=1,NMAT
  trace=trace+MAT(i,i)
enddo
do i=1,NMAT
  MAT(i,i)=MAT(i,i)-trace/dcmplx(dble(NMAT))
enddo

end subroutine make_matrix_traceless


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! CGで逆行列。16x16くらいまではこちらの方が速い
subroutine HermitianMatrix_Inverse_CG( invMAT, MAT, epsilon, info )
implicit none

complex(kind(0d0)), intent(out) :: invMat(:,:)
complex(kind(0d0)), intent(in) :: Mat(:,:)
double precision, intent(in) :: epsilon
integer, intent(out) :: info

integer :: NMAT,i,j
complex(kind(0d0)), allocatable :: unit_vec(:)
complex(kind(0d0)), allocatable :: tmpvec(:)

NMAT=size(MAT,1)
allocate( unit_vec(1:NMAT))
allocate( tmpvec(1:NMAT))

do i=1,NMAT
  unit_vec=(0d0,0d0)
  unit_vec(i)=(1d0,0d0)
  call CG( invmat(:,i), MAT, unit_vec, epsilon, 100*NMAT, info )
enddo

end subroutine HermitianMatrix_Inverse_CG


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Conjugate Graduate Solver 
! evaluate Xvec of 
!   MAT.Xvec = Bvec
subroutine CG( Xvec, MAT, Bvec, epsilon, MaxIte, info)
implicit none

complex(kind(0d0)), intent(inout) :: Xvec(:)
complex(kind(0d0)), intent(in) :: Bvec(:)
real(8), intent(in) :: epsilon
integer, intent(in) :: MaxIte
integer, intent(out) :: info 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!== DEPEND ON ProdMat YOUR MADE ==
complex(kind(0d0)), intent(in) :: MAT(:,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! private variables
integer :: NMAT
complex(kind(0d0)), allocatable :: r_vec(:)
complex(kind(0d0)), allocatable :: p_vec(:)
real(8) :: alpha_k!, alpha_km1
real(8) :: beta_k!, beta_km1


!! for iterations
integer :: i,ite
real(8) :: tmp_r1
complex(kind(0d0)) :: tmp_c1, tmp_c2, rkrk
complex(kind(0d0)), allocatable :: tmp_vec(:)

!! initialization
info = 0
NMAT=size(Mat,1)

allocate(r_vec(1:NMAT))
allocate(p_vec(1:NMAT))
allocate(tmp_vec(1:NMAT))

!Xvec=(0d0,0d0)
do i=1,NMAT
  if( cdabs(MAT(i,i)) > epsilon ) then
    Xvec(i)=Bvec(i)/MAT(i,i)
  else
    Xvec(i)=(0d0,0d0)
  endif
enddo
call Prod_MatVec(tmp_vec, MAT, Xvec)
r_vec = Bvec-tmp_vec
!r_vec = Bvec
p_vec = r_vec

!alpha_km1 = 1d0
!beta_km1 = 0d0

!! iteration start
do ite=1,MaxIte

  ! (1) construct \alpha_k
  ! rkrk = (r_k, r_k)
  call InnerProd(rkrk, r_vec, r_vec)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!== DEPEND ON ProdMat YOUR MADE ==
  call Prod_MatVec(tmp_vec, MAT, p_vec)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call InnerProd(tmp_c1, p_vec, tmp_vec)
  alpha_k = rkrk / dble(tmp_c1)
 
 ! (2) update Xvec and r_vec
  Xvec = Xvec + dcmplx(alpha_k)*p_vec
  r_vec = r_vec - dcmplx(alpha_k)*tmp_vec
    
 ! (3) conservation check
  call InnerProd( tmp_c1, r_vec, r_vec)
  if ( dsqrt( dble(tmp_c1) ) < epsilon ) then
    return
  endif

! (4) update p_k --> p_{k+1}
!construct beta_k
  !call InnerProd(tmp_c1, r_vec, r_vec)
  beta_k = dble(tmp_c1) / rkrk
  p_vec = r_vec + dcmplx(beta_k) * p_vec
enddo

info = 1
return
end subroutine CG

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to take inner product
!!  MatVec = Mat . Vec
subroutine Prod_MatVec(MatVec, MAT, vec)
implicit none

complex(kind(0d0)), intent(out) :: MatVec(:)
complex(kind(0d0)), intent(in) :: Mat(:,:)
complex(kind(0d0)), intent(in) :: Vec(:)

integer :: NMAT,i,j

NMAT=size(Mat,2)

MatVec=(0d0,0d0)
do i=1,NMAT
  do j=1,NMAT
    Matvec(i) = MatVec(i) + Mat(i,j)*Vec(j)
  enddo
enddo

end subroutine Prod_MatVec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to take inner product
!!  (V1,V2) = V1^\dagger_i V2_i
subroutine InnerProd(v1v2,v1,v2)
implicit none

complex(kind(0d0)), intent(out) ::  v1v2
complex(kind(0d0)), intent(in) :: v1(:),v2(:)

integer :: Nvec,i

Nvec=size(v1,1)
v1v2=(0d0,0d0)
do i=1,Nvec
  v1v2=v1v2+dconjg(v1(i))*v2(i)
enddo


end subroutine InnerProd



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
complex(kind(0d0)) function tr_MdagM(MAT1,MAT2)
implicit none


complex(kind(0d0)), intent(in) :: MAT1(:,:)
complex(kind(0d0)), intent(in) :: MAT2(:,:)

integer :: NMAT
integer :: i,j

NMAT=size(MAT1,1)
tr_MdagM=(0d0,0d0)
if( size(MAT1,2)/=NMAT .or. size(MAT2,1)/=NMAT .or. size(MAT2,2)/NMAT) then
  write(*,*) "### Check matrix size in tr_MdagM ###"
  return
endif

do i=1,NMAT
  do j=1,NMAT
    tr_MdagM=tr_MdagM+dconjg(MAT1(j,i))*MAT2(i,j)
  enddo
enddo

end function tr_MdagM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine make_random_matrix(MAT)
implicit none

complex(kind(0d0)), intent(out) :: MAT(:,:)
integer :: NMAT
double precision, allocatable :: tmpmat(:,:)

NMAT=size(MAT,1)
allocate(tmpmat(1:NMAT,1:NMAT))

call random_number(tmpmat)

MAT=dcmplx(tmpmat)

call random_number(tmpmat)

MAT=MAT + (0d0,1d0)*dcmplx(tmpmat)

end subroutine make_random_matrix


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine CalcPfaffian2(abs_pf,arg_pf,MAT)
!implicit none
!
!complex(kind(0d0)) :: PF
!real(8), intent(out) :: abs_pf, arg_pf
!complex(kind(0d0)), intent(in) :: MAT(:,:)
!
!integer :: NMAT
!integer, allocatable :: IWORK(:)
!integer  LWORK, LDA,NB, INFO
!complex(kind(0d0)), allocatable :: WORK(:)
!real(8), allocatable :: RWORK(:)
!real(8) LogPf
!complex(kind(0d0)) phase
!
!NMAT=size(MAT,1)
!
!NB = NMAT*NMAT
!LWORK = NMAT*NB
!allocate ( IWORK(1:NMAT) )
!allocate ( RWORK(1:NMAT-1) )
!allocate ( WORK(LWORK) )
!LDA = NMAT
!
!call ZSKPFA('U','P',NMAT,MAT,LDA,PF, &
!IWORK,WORK,LWORK,RWORK,INFO)
!
!abs_pf=abs(PF)
!arg_pf=atan2( dble((0d0,-1d0)*PF), dble(PF) )
!end subroutine CalcPfaffian2

end module matrix_functions
