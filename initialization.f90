module initialization
use global_parameters
use mt95
#ifdef PARALLEL
use parallel
#endif
implicit none

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! read the preevious configuration
subroutine read_config(total_ite,UMAT,PhiMat,state_mt95)
implicit none

integer, intent(inout) :: total_ite
type(genrand_state), intent(inout) :: state_mt95
complex(kind(0d0)), intent(inout) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(inout) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)

type(genrand_srepr) :: char_mt95

open(IN_CONF_FILE, file=Fconfigin, status='OLD',action='READ',form='unformatted')
read(IN_CONF_FILE) total_ite
read(IN_CONF_FILE) UMAT
read(IN_CONF_FILE) PHIMat
read(IN_CONF_FILE) char_mt95
close(IN_CONF_FILE)
state_mt95=char_mt95

end subroutine read_config

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! set random configuration
!! 
!! For SU(N) case, the initial UMAT must be nearer than
!! all the Z_N centers, Omega_n=diag( exp( 2\pi i n / N ) ) from 1_N. 
!! We see
!!   min_n( || 1 - Omega_n || ) = 2 sin( \pi/N ). 
!! On the other hand, 
!!   || 1 - U || = 4/N Tr( sin^2( \theat T / 2 ) ) \sim < \theta^2 > 
!! Thus the random number must satisfy 
!!   <\theta^2> < 2\pi/N
!! 
subroutine set_random_config(UMAT,PhiMat)
use SUN_generators, only : Make_SUN_generators
use matrix_functions, only : matrix_exp
use global_subroutines, only : BoxMuller2
implicit none

complex(kind(0d0)), intent(inout) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(inout) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)

complex(kind(0d0)) :: TMAT(1:NMAT,1:NMAT,1:NMAT**2-1)
double precision :: rsite(1:2*NMAT*NMAT*num_sites) ! for PHI
double precision :: rlink(1:dimG,1:num_links) ! for UMAT
complex(kind(0d0)) :: AMAT(1:NMAT,1:NMAT,1:num_links)
integer :: s,l,a,f,i,j,num

call make_SUN_generators(TMAT,NMAT)

!call genrand_real3(rsite)
call BoxMuller2(rsite,2*num_sites*NMAT*NMAT)
call genrand_real3(rlink)

rsite=rsite * 0.01d0 !/ mass_square_phi
num=0
do s=1,num_sites
  do i=1,NMAT
    do j=1,NMAT
      num=num+1
      if ( i.ne.NMAT .or. j.ne.NMAT ) then
        PHIMAT(i,j,s)=dcmplx(rsite(2*num-1))+(0d0,1d0)*dcmplx(rsite(2*num))
      endif
    enddo
  enddo
  PhiMat(NMAT,NMAT,s)=(0d0,0d0)
  do i=1,NMAT-1
    PhiMat(NMAT,NMAT,s)=PhiMat(NMAT,NMAT,s)-PhiMat(i,i,s)
  enddo
enddo


! random number must be sufficiently small
if( m_omega == 0 ) then 
  rlink=rlink * ( 1d0/dble(NMAT*NMAT) )
else
  rlink=rlink * ( 1d0/dble(NMAT*NMAT*m_omega) )
endif
AMAT=(0d0,0d0)
do l=1,num_links
  do a=1,dimG
    AMAT(:,:,l)=AMAT(:,:,l)+rlink(a,l)*TMAT(:,:,a)
  enddo
enddo


do l=1,num_links
!call matrix_exp(NMAT,(0d0,1d0)*AMAT(:,:,l),UMAT(:,:,l))
call matrix_exp(UMAT(:,:,l),(0d0,1d0)*AMAT(:,:,l))
enddo

end subroutine set_random_config


end module initialization
