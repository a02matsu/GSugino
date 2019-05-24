!!!!!!!!!!!!!!!!!!!!
!! 
!! 1) construct_fabc(f,num_nonzerof, NMAT)
!!  subroutine to create the structure constant of SU(N)
!! 
!!  f: out : complex(kind(0d0)) : f(1:NMAT**2-1,1:NMAT**2-1,1:NMAT**2-1)
!!  num_nonzerof : out : integer 
!!  NMAT : in : integer
!!
!! 2) get_nonzerof( num_nonzerof, nonzerof_index, nonzerof_value, NMAT )
!!  subroutine to obtain data of nonzero f(a,b,c)
!!
!!  num_nonzerof : out : integer 
!!  nonzerof_index : out : integer(1:3,1:num_nonzerof)
!!  nonzerof_value : out : complex(kind(0d0))(1:num_nonzerof)
!!  NMAT : in : integer
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module SUN_generators
implicit none
!private
!public set_NZF
contains

!subroutine set_NZF(NZF,NZF_index,NZF_value,NMAT)
!implicit none
!integer :: NZF      ! number of nonzero elements of f(a,b,c)
!integer, allocatable :: NZF_index(:,:)  ! index of nonzero f(a,b,c)
!double precision, allocatable :: NZF_value(:) ! value of nonzero f(a,b,c)
!integer NMAT
! 
!call get_num_nonzerof(NZF,NMAT)
!allocate( NZF_index(1:3,1:NZF) )
!allocate( NZF_value(1:NZF) )
!call get_nonzerof(NZF_index,NZF_value,NZF,NMAT)
!
!end subroutine set_NZF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to create the SU(N) generators
subroutine make_SUN_generators(T,NMAT)
implicit none

integer, intent(in) :: NMAT
complex(kind(0d0)), intent(out) :: T(1:NMAT,1:NMAT,1:NMAT**2-1) ! generators
integer a,i,j

T=(0d0,0d0)
a=0
do i=1,NMAT-1
  do j=i+1,NMAT
    a=a+1
    T(i,j,a)=(1d0,0d0)/dcmplx(dsqrt(2d0))
    T(j,i,a)=(1d0,0d0)/dcmplx(dsqrt(2d0))
    a=a+1
    T(i,j,a)=(0d0,-1d0)/dcmplx(dsqrt(2d0))
    T(j,i,a)=(0d0,1d0)/dcmplx(dsqrt(2d0))
  enddo
enddo
do i=1,NMAT-1
  a=a+1
  do j=1,i
    T(j,j,a)=dcmplx(1d0/dsqrt(dble(i**2+i)))
  enddo
  T(i+1,i+1,a)=dcmplx(-dble(i)/dsqrt(dble(i**2+i)))
enddo
end subroutine make_SUN_generators


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to create the structure constant of SU(N)
subroutine construct_fabc(f, NMAT)
implicit none

integer, intent(in) :: NMAT
double precision, intent(out) :: f(1:NMAT**2-1,1:NMAT**2-1,1:NMAT**2-1)

double precision, parameter :: epsilon=1d-10 ! zero indeed
complex(kind(0d0)) T(1:NMAT,1:NMAT,1:NMAT**2-1) ! generators
complex(kind(0d0)) tmp
complex(kind(0d0)) tmpmat(1:NMAT,1:NMAT)


integer i,j,k
integer a,b,c

call Make_SUN_generators(T,NMAT)

f=0d0
do a=1,NMAT**2-1
do b=1,NMAT**2-1
  tmpmat=(0d0,0d0)
  do i=1,NMAT
  do j=1,NMAT
  do k=1,NMAT
    tmpmat(i,j)=tmpmat(i,j)+T(i,k,a)*T(k,j,b)
    tmpmat(i,j)=tmpmat(i,j)-T(i,k,b)*T(k,j,a)
  enddo
  enddo
  enddo
  
  ! f_{a,b,c}=-i Tr( [T_a,T_b] T_c )
  do c=1,NMAT**2-1
    tmp=(0d0,0d0)
    do i=1,NMAT
    do j=1,NMAT
      tmp=tmp+tmpmat(i,j)*T(j,i,c)*(0d0,-1d0)
    enddo
    enddo
    f(a,b,c)=dble(tmp)
  !! check
    !if ( abs(f(a,b,c)) > epsilon ) then
      !write(*,*) a,b,c,f(a,b,c)
    !endif
  enddo
enddo
enddo

!! check
!do a=1,NMAT**2-1
!do b=1,NMAT**2-1
!do c=1,NMAT**2-1
!  tmp=f(a,b,c)-f(c,a,b) 
!  if ( abs(tmp) > epsilon ) then
!  write(*,*) a,b,c,tmp
!endif
!enddo
!enddo
!enddo
end subroutine construct_fabc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to obtain data of nonzero f(a,b,c)
subroutine get_num_nonzerof( num_nonzerof, NMAT )
implicit none

double precision, parameter :: epsilon=1d-10 ! zero indeed
integer, intent(in) :: NMAT
integer, intent(out) :: num_nonzerof ! number of nonzero components of f(a,b,c)
double precision f(1:NMAT**2-1,1:NMAT**2-1,1:NMAT**2-1)
integer a,b,c

call construct_fabc(f, NMAT)

num_nonzerof=0
do a=1,NMAT**2-1
do b=1,NMAT**2-1
do c=1,NMAT**2-1
    if( abs(f(a,b,c)) > epsilon ) then
        num_nonzerof=num_nonzerof+1
    endif
enddo
enddo
enddo


end subroutine get_num_nonzerof


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to obtain data of nonzero f(a,b,c)
subroutine get_nonzerof( nonzerof_index, nonzerof_value, num_nonzerof, NMAT )
implicit none

double precision, parameter :: epsilon=1d-10 ! zero indeed
integer, intent(in) :: NMAT
integer, intent(in) :: num_nonzerof ! number of nonzero components of f(a,b,c)
integer, intent(out) :: nonzerof_index(1:3,1:num_nonzerof) ! index of nonzero f(a,b,c)
double precision, intent(out) :: nonzerof_value(1:num_nonzerof) ! value of nonzero f(a,b,c)

double precision f(1:NMAT**2-1,1:NMAT**2-1,1:NMAT**2-1)
integer a,b,c
integer i


call construct_fabc(f,NMAT)

! find non-zero components
i=0
do a=1,NMAT**2-1
do b=1,NMAT**2-1
do c=1,NMAT**2-1
  if ( abs(f(a,b,c)) > epsilon ) then
    i=i+1
    nonzerof_index(1,i)=a
    nonzerof_index(2,i)=b
    nonzerof_index(3,i)=c
    nonzerof_value(i)=f(a,b,c)
  endif
enddo
enddo
enddo

end subroutine get_nonzerof

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to obtain data of nonzero f(a,b,e)*f(c,d,e)
subroutine get_num_nonzero_ff( num_nonzero_ff, NMAT )
implicit none

double precision, parameter :: epsilon=1d-10 ! zero indeed
integer, intent(in) :: NMAT
integer, intent(out) :: num_nonzero_ff ! number of nonzero components of f*f
double precision f(1:NMAT**2-1,1:NMAT**2-1,1:NMAT**2-1),tmp
integer a,b,c,d,e

call construct_fabc(f, NMAT)

num_nonzero_ff=0
do a=1,NMAT**2-1
do b=1,NMAT**2-1
do c=1,NMAT**2-1
do d=1,NMAT**2-1
  tmp=0d0
  do e=1,NMAT**2-1
    tmp=tmp+f(a,b,e)*f(c,d,e)
  enddo
    if( dabs(tmp) > epsilon ) then
        num_nonzero_ff=num_nonzero_ff+1
    endif
enddo
enddo
enddo
enddo

end subroutine get_num_nonzero_ff


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to obtain data of nonzero f(a,b,e)*f(c,d,e)
subroutine get_nonzero_ff( nonzeroff_index, nonzeroff_value, num_nonzero_ff, NMAT )
implicit none

double precision, parameter :: epsilon=1d-10 ! zero indeed
integer, intent(in) :: NMAT
integer, intent(in) :: num_nonzero_ff ! number of nonzero components of f(a,b,c)
integer, intent(out) :: nonzeroff_index(1:4,1:num_nonzero_ff) ! index of nonzero f(a,b,e)*f(c,d,e)
double precision, intent(out) :: nonzeroff_value(1:num_nonzero_ff) ! value of nonzero f(a,b,e)*f(c,d,e)

double precision f(1:NMAT**2-1,1:NMAT**2-1,1:NMAT**2-1), tmp
integer a,b,c,d,e
integer i

call construct_fabc(f,NMAT)

i=0
do a=1,NMAT**2-1
do b=1,NMAT**2-1
do c=1,NMAT**2-1
do d=1,NMAT**2-1
  tmp=0d0
  do e=1,NMAT**2-1
    tmp=tmp+f(a,b,e)*f(c,d,e)
  enddo
    if( dabs(tmp) > epsilon ) then
      i=i+1
      nonzeroff_index(1,i)=a
      nonzeroff_index(2,i)=b
      nonzeroff_index(3,i)=c
      nonzeroff_index(4,i)=d
      nonzeroff_value(i)=tmp
    endif
enddo
enddo
enddo
enddo

end subroutine get_nonzero_ff




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! obtain T^p_{(ij)} m\in[ 1, NMAT**2-1 ] 
!!
!! i<j, p=1,2
!!  T^1_{(ij)}_{kl} = 1/sqrt(2) (\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk}
!!  T^2_{(ij}}_{kl} = i/sqrt(2) (-\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk}
!! i=j, p=3
!!  T^3_{(ii}}_{kl} = 1/sqrt(i^2+i) ( \delta_{kl} \theta(i-k) -i \delta_{ki}\delta_{li}
subroutine get_Taij(NMAT,a,i,j,p)
implicit none

integer, intent(out) :: i,j,p
integer, intent(in) :: a,NMAT
integer :: n,k


if(a <= NMAT*(NMAT-1)) then
  if( mod(a,2) == 1 ) then 
    p=1
    n=(a+1)/2
  else
    p=2
    n=a/2
  endif
    i=1
  do while ( n > (i-1)*(2*NMAT-i)/2 )
    i=i+1
  enddo
  i=i-1
  j=n-(i-1)*(2*NMAT-i)/2+i
else
  p=3
  i = a - NMAT*(NMAT-1)
  j=i
endif

end subroutine get_Taij


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MAT1.T_a
subroutine tmpMtimesT(MT,MAT,a,NMAT)
implicit none

complex(kind(0d0)), intent(out) :: MT(1:NMAT,1:NMAT)
integer, intent(in) :: a, NMAT
complex(kind(0d0)), intent(in) :: MAT(1:NMAT,1:NMAT)
complex(kind(0d0)) :: T(1:NMAT,1:NMAT,1:NMAT**2-1)

call make_SUN_generators(T,NMAT)

call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
  MAT, NMAT, &
  T(:,:,a), NMAT, &
  (0d0,0d0), MT, NMAT)
end subroutine tmpMtimesT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MAT1.T_a
subroutine tmpTtimesM(MT,MAT,a,NMAT)
implicit none

complex(kind(0d0)), intent(out) :: MT(1:NMAT,1:NMAT)
integer, intent(in) :: a, NMAT
complex(kind(0d0)), intent(in) :: MAT(1:NMAT,1:NMAT)
complex(kind(0d0)) :: T(1:NMAT,1:NMAT,1:NMAT**2-1)

call make_SUN_generators(T,NMAT)

call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
  T(:,:,a), NMAT, &
  MAT, NMAT, &
  (0d0,0d0), MT, NMAT)
end subroutine tmpTtimesM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MAT1.T_a.MAT2
subroutine tmpMTN(MTNMAT,MAT1,MAT2,a,NMAT)
implicit none

complex(kind(0d0)), intent(out) :: MTNMAT(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: MAT1(1:NMAT,1:NMAT), MAT2(1:NMAT,1:NMAT)
integer, intent(in) :: a, NMAT
complex(kind(0d0)) :: T(1:NMAT,1:NMAT,1:NMAT**2-1),tmpmat(1:NMAT,1:NMAT)

call make_SUN_generators(T,NMAT)

call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
  MAT1, NMAT, &
  T(:,:,a), NMAT, &
  (0d0,0d0), tmpmat, NMAT)

call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
  tmpmat, NMAT, &
  MAT2, NMAT, &
  (0d0,0d0), MTNMAT, NMAT)

end subroutine tmpMTN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MAT1.T_a
subroutine MtimesT(MT,MAT,a,NMAT)
implicit none

complex(kind(0d0)), intent(out) :: MT(1:NMAT,1:NMAT)
integer, intent(in) :: a, NMAT
complex(kind(0d0)), intent(in) :: MAT(1:NMAT,1:NMAT)
integer :: i,j,p
integer :: k,l

MT=(0d0,0d0)
call get_Taij(NMAT,a,i,j,p)


if( p==1 ) then
  do k=1,NMAT
    MT(k,j)=MAT(k,i)/dcmplx(dsqrt(2d0))
    MT(k,i)=MAT(k,j)/dcmplx(dsqrt(2d0))
  enddo
elseif( p==2 ) then
  do k=1,NMAT
    MT(k,j)=(0d0,-1d0)*MAT(k,i)/dcmplx(dsqrt(2d0))
    MT(k,i)=(0d0,1d0)*MAT(k,j)/dcmplx(dsqrt(2d0))
  enddo
else
  do l=1,i
    do k=1,NMAT
      MT(k,l)=MAT(k,l)/dcmplx(dsqrt(dble(i*i+i)))
    enddo
  enddo
  do k=1,NMAT
    MT(k,i+1)=MAT(k,i+1)*dcmplx(dble(-i))/dcmplx(dsqrt(dble(i*i+i)))
  enddo
endif
end subroutine MtimesT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MAT1.T_a
subroutine TtimesM(TM,MAT,a,NMAT)
implicit none

complex(kind(0d0)), intent(out) :: TM(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: MAT(1:NMAT,1:NMAT)
integer, intent(in) :: a, NMAT
integer :: i,j,p
integer :: k,l


TM=(0d0,0d0)
call get_Taij(NMAT,a,i,j,p)

if( p==1 ) then
  do l=1,NMAT
    TM(i,l)=MAT(j,l)/dcmplx(dsqrt(2d0))
    TM(j,l)=MAT(i,l)/dcmplx(dsqrt(2d0))
  enddo
elseif( p==2 ) then
  do l=1,NMAT
    TM(i,l)=(0d0,-1d0)*MAT(j,l)/dcmplx(dsqrt(2d0))
    TM(j,l)=(0d0,1d0)*MAT(i,l)/dcmplx(dsqrt(2d0))
  enddo
else
  do l=1,NMAT
    do k=1,i
      TM(k,l)=MAT(k,l)/dcmplx(dsqrt(dble(i*i+i)))
    enddo
  enddo
  do l=1,NMAT
    TM(i+1,l)=MAT(i+1,l)*dcmplx(dble(-i))/dcmplx(dsqrt(dble(i*i+i)))
  enddo
endif
end subroutine TtimesM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! [T_a, MAT]
subroutine commutator_TaM(comm,MAT,a,NMAT)
implicit none

complex(kind(0d0)), intent(out) :: comm(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: MAT(1:NMAT,1:NMAT)
integer, intent(in) :: NMAT,a

complex(kind(0d0)) :: TM(1:NMAT,1:NMAT),MT(1:NMAT,1:NMAT)

call TtimesM(TM,MAT,a,NMAT)
call MtimesT(MT,MAT,a,NMAT)
comm=TM-MT

end subroutine commutator_TaM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! [T_a, MAT]
subroutine commutator_MTa(comm,MAT,a,NMAT)
implicit none

complex(kind(0d0)), intent(out) :: comm(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: MAT(1:NMAT,1:NMAT)
integer, intent(in) :: NMAT,a

complex(kind(0d0)) :: TM(1:NMAT,1:NMAT),MT(1:NMAT,1:NMAT)

call TtimesM(TM,MAT,a,NMAT)
call MtimesT(MT,MAT,a,NMAT)
comm=-TM+MT

end subroutine commutator_MTa


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MAT1.T_a.MAT2
subroutine MTN(MTNMAT,MAT1,MAT2,a,NMAT)
implicit none

complex(kind(0d0)), intent(out) :: MTNMAT(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: MAT1(1:NMAT,1:NMAT), MAT2(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: tmp(1:NMAT,1:NMAT)
integer, intent(in) :: a, NMAT
integer :: i,j,p
integer :: k,l,m

MTNMAT=(0d0,0d0)
call get_Taij(NMAT,a,i,j,p)
!write(*,*) a,i,j,p
!call tmpMTN(tmp,MAT1,MAT2,a,NMAT)

if( p==1 ) then
  do l=1,NMAT
    do k=1,NMAT
      MTNMAT(k,l)=MAT1(k,j)*MAT2(i,l)/dcmplx(dsqrt(2d0))&
                 +MAT1(k,i)*MAT2(j,l)/dcmplx(dsqrt(2d0))
    enddo
  enddo
elseif( p==2 ) then
  do l=1,NMAT
    do k=1,NMAT
      MTNMAT(k,l)=MAT1(k,j)*MAT2(i,l)*(0d0,1d0)/dcmplx(dsqrt(2d0))&
                 +MAT1(k,i)*MAT2(j,l)*(0d0,-1d0)/dcmplx(dsqrt(2d0))
    enddo
  enddo
else
  do l=1,NMAT
    do k=1,NMAT
      do m=1,i
        MTNMAT(k,l)=MTNMAT(k,l)&
          +MAT1(k,m)*MAT2(m,l)/dcmplx(dsqrt(dble(i*i+i)))
      enddo
      MTNMAT(k,l)=MTNMAT(k,l)&
          -dcmplx(dble(i))*MAT1(k,i+1)*MAT2(i+1,l)/dcmplx(dsqrt(dble(i*i+i)))
    enddo
  enddo
endif

!write(*,*) "=---================="
!write(*,*) p
!write(*,*) MTNMAT-tmp

end subroutine MTN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MAT1.T_a.MAT2
subroutine MTMd(MTNMAT,MAT1,a,NMAT)
implicit none

complex(kind(0d0)), intent(out) :: MTNMAT(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: MAT1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: MAT2(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: tmp(1:NMAT,1:NMAT)
integer, intent(in) :: a, NMAT
integer :: i,j,p
integer :: k,l,m

do i=1,NMAT
  do j=1,NMAT
    MAT2(i,j)=dconjg(MAT1(j,i))
  enddo
enddo

call MTN(MTNMAT,MAT1,MAT2,a,NMAT)

end subroutine MTMd


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MAT^\dagger .T_a.MAT
subroutine MdTM(MTNMAT,MAT1,a,NMAT)
implicit none

complex(kind(0d0)), intent(out) :: MTNMAT(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: MAT1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: MAT2(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: tmp(1:NMAT,1:NMAT)
integer, intent(in) :: a, NMAT
integer :: i,j,p
integer :: k,l,m

do i=1,NMAT
  do j=1,NMAT
    MAT2(i,j)=dconjg(MAT1(j,i))
  enddo
enddo

call MTN(MTNMAT,MAT2,MAT1,a,NMAT)

end subroutine MdTM



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Tr( MAT.T_a )
subroutine trace_MTa(trace,MAT,a,NMAT)
implicit none

complex(kind(0d0)), intent(out) :: trace
complex(kind(0d0)), intent(in) :: MAT(1:NMAT,1:NMAT)
integer, intent(in) :: a, NMAT
integer :: i,j,p
integer :: k,l

call get_Taij(NMAT,a,i,j,p)
!write(*,*) a,p,i,j

if( p==1 ) then
  trace = ( MAT(j,i) + MAT(i,j) ) / dcmplx(dsqrt(2d0))
  return
elseif( p==2 ) then
  trace = ( MAT(j,i) - MAT(i,j) ) * (0d0,-1d0) / dcmplx(dsqrt(2d0))
  return
elseif( p==3 ) then
  trace = (0d0,0d0)
  do k=1,i
    trace = trace + MAT(k,k)
  enddo
  trace=trace-dcmplx(dble(i))*MAT(i+1,i+1)

  trace=trace/dcmplx(dsqrt(dble(i*i+i)))
  return
endif

end subroutine trace_MTa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  tr( T_a.MAT1.T_b.MAT2 )
subroutine trace_TMTN(trace,MAT1,MAT2,a,b,NMAT)
implicit none

complex(kind(0d0)), intent(out) :: trace
integer, intent(in) :: NMAT,a,b
complex(kind(0d0)), intent(in) :: MAT1(1:NMAT,1:NMAT),MAT2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
integer i

call MTN(tmpmat1,MAT1,MAT2,b,NMAT)
call trace_MTa(trace,tmpmat1,a,NMAT)
!call TtimesM(tmpmat2,tmpmat1,a,NMAT)

!trace=(0d0,0d0)
!do i=1,NMAT
!  trace=trace + tmpmat2(i,i)
!enddo
end subroutine trace_TMTN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  tr( T_a.MAT.T_b.MAT^\dagger )
subroutine trace_TMTMd(trace,MAT,a,b,NMAT)
implicit none

complex(kind(0d0)), intent(out) :: trace
integer, intent(in) :: NMAT,a,b
complex(kind(0d0)), intent(in) :: MAT(1:NMAT,1:NMAT)
complex(kind(0d0)) :: MATd(1:NMAT,1:NMAT)
integer i,j

do j=1,NMAT
  do i=1,NMAT
    MATd(i,j)=dconjg(MAT(j,i))
  enddo
enddo

call trace_TMTN(trace,MAT,MATd,a,b,NMAT)

end subroutine trace_TMTMd


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Create traceless matrix from MODES(1:dimG)
!!  MODES(1:dimG): complex modes M_a
!!  dimG = NMAT*NMAT-1
!!  
!!  MAT=T^a M_a (T^a: SU(N) generator)
subroutine Make_traceless_matrix_from_modes(MAT,NMAT,MODES)
implicit none

integer, intent(in) :: NMAT
complex(kind(0d0)), intent(in) :: MODES(1:NMAT*NMAT-1)
complex(kind(0d0)), intent(out) :: MAT(1:NMAT,1:NMAT)
integer a,i,j

a=0
MAT=(0d0,0d0)
do i=1,NMAT-1
  do j=i+1,NMAT
    a=a+1
    MAT(i,j)=MAT(i,j)+MODES(a)/dcmplx(dsqrt(2d0))
    MAT(j,i)=MAT(j,i)+MODES(a)/dcmplx(dsqrt(2d0))
    a=a+1
    MAT(i,j)=MAT(i,j)+(0d0,-1d0)*MODES(a)/dcmplx(dsqrt(2d0))
    MAT(j,i)=MAT(j,i)+(0d0,1d0)*MODES(a)/dcmplx(dsqrt(2d0))
  enddo
enddo
do i=1,NMAT-1
  a=a+1
  do j=1,i
    MAT(j,j)=MAT(j,j)+MODES(a)*dcmplx(1d0/dsqrt(dble(i**2+i)))
  enddo
  MAT(i+1,i+1)=MAT(i+1,i+1)+MODES(a)*dcmplx(-dble(i)/dsqrt(dble(i**2+i)))
enddo

end subroutine Make_traceless_matrix_from_modes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compute
!!  MODES(a,b)*Ta_{ij}*Tb_{kl}
subroutine make_Mijkl_from_modes(DP,MODES,NMAT)
integer, intent(in) :: NMAT
complex(kind(0d0)), intent(in) :: MODES(1:NMAT*NMAT-1,1:NMAT*NMAT-1)
complex(kind(0d0)), intent(out) :: DP(1:NMAT,1:NMAT,1:NMAT,1:NMAT)

complex(kind(0d0)) :: MED(1:NMAT*NMAT-1,1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
integer :: b,i,j,k,l

do b=1,NMAT*NMAT-1
  call Make_traceless_matrix_from_modes(tmpmat,NMAT,MODES(:,b))
  do j=1,NMAT
    do i=1,NMAT
      MED(b,i,j)=tmpmat(i,j)
    enddo
  enddo
enddo

do j=1,NMAT
  do i=1,NMAT
    call Make_traceless_matrix_from_modes(tmpmat,NMAT,MED(:,i,j))
    do l=1,NMAT
      do k=1,NMAT
        DP(i,j,k,l)=tmpmat(k,l)
      enddo
    enddo
  enddo
enddo

end subroutine make_Mijkl_from_modes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compute
!!  MODES(a,b)*Ta_{ij}*Tb_{kl}
subroutine make_Mijkl_from_modes_org(DP,MODES,NMAT)
implicit none

integer, intent(in) :: NMAT
complex(kind(0d0)), intent(in) :: MODES(1:NMAT*NMAT-1,1:NMAT*NMAT-1)
complex(kind(0d0)), intent(out) :: DP(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
integer :: a,b,ai,aj,bi,bj

DP=(0d0,0d0)
b=0
do bi=1,NMAT-1
  do bj=bi+1,NMAT
    b=b+1
    !!!!!!!!!!
    a=0
    do ai=1,NMAT-1
      do aj=ai+1,NMAT
        a=a+1
        DP(ai,aj,bi,bj)= DP(ai,aj,bi,bj) &
          + MODES(a,b)/dcmplx(2d0)
        DP(aj,ai,bi,bj)= DP(aj,ai,bi,bj) &
          + MODES(a,b)/dcmplx(2d0)
        DP(ai,aj,bj,bi)= DP(ai,aj,bj,bi) &
          + MODES(a,b)/dcmplx(2d0)
        DP(aj,ai,bj,bi)= DP(aj,ai,bj,bi) &
          + MODES(a,b)/dcmplx(2d0)
        a=a+1
        DP(ai,aj,bi,bj)= DP(ai,aj,bi,bj) &
          + MODES(a,b)*(0d0,-1d0)/dcmplx(2d0)
        DP(aj,ai,bi,bj)= DP(aj,ai,bi,bj) &
          + MODES(a,b)*(0d0,1d0)/dcmplx(2d0)
        DP(ai,aj,bj,bi)= DP(ai,aj,bj,bi) &
          + MODES(a,b)*(0d0,-1d0)/dcmplx(2d0)
        DP(aj,ai,bj,bi)= DP(aj,ai,bj,bi) &
          + MODES(a,b)*(0d0,1d0)/dcmplx(2d0)
      enddo
    enddo
    do ai=1,NMAT-1
      a=a+1
      do aj=1,ai
        DP(aj,aj,bi,bj)= DP(aj,aj,bi,bj) &
          + MODES(a,b)*dcmplx(1d0/dsqrt(2d0*dble(ai**2+ai)))
        DP(aj,aj,bj,bi)= DP(aj,aj,bj,bi) &
          + MODES(a,b)*dcmplx(1d0/dsqrt(2d0*dble(ai**2+ai)))
      enddo
      DP(ai+1,ai+1,bi,bj)= DP(aj,aj,bi,bj) &
        + MODES(a,b)*dcmplx(-dble(ai)/dsqrt(2d0*dble(ai**2+ai)))
      DP(ai+1,ai+1,bj,bi)= DP(ai+1,ai+1,bj,bi) &
        + MODES(a,b)*dcmplx(-dble(ai)/dsqrt(2d0*dble(ai**2+ai)))
    enddo
    !!!!!!!!!
    b=b+1
    a=0
    do ai=1,NMAT-1
      do aj=ai+1,NMAT
        a=a+1
        DP(ai,aj,bi,bj)= DP(ai,aj,bi,bj) &
          + MODES(a,b)*(0d0,-1d0)/dcmplx(2d0)
        DP(aj,ai,bi,bj)= DP(aj,ai,bi,bj) &
          + MODES(a,b)*(0d0,-1d0)/dcmplx(2d0)
        DP(ai,aj,bj,bi)= DP(ai,aj,bj,bi) &
          + MODES(a,b)*(0d0,1d0)/dcmplx(2d0)
        DP(aj,ai,bj,bi)= DP(aj,ai,bj,bi) &
          + MODES(a,b)*(0d0,1d0)/dcmplx(2d0)
        a=a+1
        DP(ai,aj,bi,bj)= DP(ai,aj,bi,bj) &
          + MODES(a,b)*(-1d0,0d0)/dcmplx(2d0)
        DP(aj,ai,bi,bj)= DP(aj,ai,bi,bj) &
          + MODES(a,b)*(1d0,0d0)/dcmplx(2d0)
        DP(ai,aj,bj,bi)= DP(ai,aj,bj,bi) &
          + MODES(a,b)*(1d0,0d0)/dcmplx(2d0)
        DP(aj,ai,bj,bi)= DP(aj,ai,bj,bi) &
          + MODES(a,b)*(-1d0,0d0)/dcmplx(2d0)
      enddo
    enddo
    do ai=1,NMAT-1
      a=a+1
      do aj=1,ai
        DP(aj,aj,bi,bj)= DP(aj,aj,bi,bj) &
          + MODES(a,b)*(0d0,-1d0)*dcmplx(1d0/dsqrt(2d0*dble(ai**2+ai)))
        DP(aj,aj,bj,bi)= DP(aj,aj,bj,bi) &
          + MODES(a,b)*(0d0,1d0)*dcmplx(1d0/dsqrt(2d0*dble(ai**2+ai)))
      enddo
      DP(ai+1,ai+1,bi,bj)= DP(ai+1,ai+1,bi,bj) &
        + MODES(a,b)*(0d0,-1d0)*dcmplx(-dble(ai)/dsqrt(2d0*dble(ai**2+ai)))
      DP(ai+1,ai+1,bj,bi)= DP(ai+1,ai+1,bj,bi) &
        + MODES(a,b)*(0d0,1d0)*dcmplx(-dble(ai)/dsqrt(2d0*dble(ai**2+ai)))
    enddo
  enddo
enddo
do bi=1,NMAT-1
  b=b+1
  do bj=1,bi
    a=0
    do ai=1,NMAT-1
      do aj=ai+1,NMAT
        a=a+1
        DP(ai,aj,bj,bj)= DP(ai,aj,bj,bj) &
          + MODES(a,b)*dcmplx(1d0/dsqrt(2d0*dble(bi**2+bi)))
        DP(aj,ai,bj,bj)= DP(aj,ai,bj,bj) &
          + MODES(a,b)*dcmplx(1d0/dsqrt(2d0*dble(bi**2+bi)))
        a=a+1
        DP(ai,aj,bj,bj)= DP(ai,aj,bj,bj) &
          + MODES(a,b)*(0d0,-1d0)*dcmplx(1d0/dsqrt(2d0*dble(bi**2+bi)))
        DP(aj,ai,bj,bj)= DP(aj,ai,bj,bj) &
          + MODES(a,b)*(0d0,1d0)*dcmplx(1d0/dsqrt(2d0*dble(bi**2+bi)))
      enddo
    enddo
    do ai=1,NMAT-1
      a=a+1
      do aj=1,ai
        DP(aj,aj,bj,bj)= DP(aj,aj,bj,bj) &
          + MODES(a,b)&
             *dcmplx(1d0/dsqrt(dble(ai**2+ai))) &
             *dcmplx(1d0/dsqrt(dble(bi**2+bi)))
      enddo
      DP(ai+1,ai+1,bj,bj)= DP(ai+1,ai+1,bj,bj) &
        + MODES(a,b)&
             *dcmplx(-dble(ai)/dsqrt(dble(ai**2+ai))) &
             *dcmplx(1d0/dsqrt(dble(bi**2+bi)))
    enddo
  enddo
  !!!!!!!!!!!!!!!!!!!!
  a=0
  do ai=1,NMAT-1
    do aj=ai+1,NMAT
      a=a+1
      DP(ai,aj,bi+1,bi+1)= DP(ai,aj,bi+1,bi+1) &
        + MODES(a,b) &
           *dcmplx(1d0/dsqrt(2d0)) &
           *dcmplx(-dble(bi)/dsqrt(dble(bi**2+bi)))
      DP(aj,ai,bi+1,bi+1)= DP(aj,ai,bi+1,bi+1) &
        + MODES(a,b) &
           *dcmplx(1d0/dsqrt(2d0)) &
           *dcmplx(-dble(bi)/dsqrt(dble(bi**2+bi)))
      a=a+1
      DP(ai,aj,bi+1,bi+1)= DP(ai,aj,bi+1,bi+1) &
        + MODES(a,b) &
           *(0d0,-1d0)*dcmplx(1d0/dsqrt(2d0)) &
           *dcmplx(-dble(bi)/dsqrt(dble(bi**2+bi)))
      DP(aj,ai,bi+1,bi+1)= DP(aj,ai,bi+1,bi+1) &
        + MODES(a,b) &
           *(0d0,1d0)*dcmplx(1d0/dsqrt(2d0)) &
           *dcmplx(-dble(bi)/dsqrt(dble(bi**2+bi)))
    enddo
  enddo
  do ai=1,NMAT-1
    a=a+1
    do aj=1,ai
      DP(aj,aj,bi+1,bi+1)= DP(aj,aj,bi+1,bi+1) &
        + MODES(a,b)&
           *dcmplx(1d0/dsqrt(dble(ai**2+ai))) &
           *dcmplx(-dble(bi)/dsqrt(dble(bi**2+bi)))
    enddo
    DP(ai+1,ai+1,bi+1,bi+1)= DP(ai+1,ai+1,bi+1,bi+1) &
      + MODES(a,b)&
           *dcmplx(1d0/dsqrt(dble(ai**2+ai))) &
           *dcmplx(1d0/dsqrt(dble(bi**2+bi)))
  enddo
enddo


end subroutine make_Mijkl_from_modes_org





    
end module SUN_generators
