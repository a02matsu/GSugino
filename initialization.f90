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

#ifdef PARATEST
complex(kind(0d0)) :: g_UMAT(1:NMAT,1:NMAT,1:global_num_links)
complex(kind(0d0)) :: g_PhiMat(1:NMAT,1:NMAT,1:global_num_sites)
complex(kind(0d0)) :: tmp_UMAT(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmp_PhiMat(1:NMAT,1:NMAT)
double precision :: g_rsite(1:2*NMAT*NMAT*global_num_sites) ! for PHI
double precision :: g_rlink(1:dimG,1:global_num_links) ! for UMAT
complex(kind(0d0)) :: g_AMAT(1:NMAT,1:NMAT,1:global_num_links)
integer :: ll,ls,rank
#endif

call make_SUN_generators(TMAT,NMAT)

!! テスト用に、シングルコアと同じconfigurationを用意する。
#ifdef PARATEST
!call genrand_real3(rsite)
if( MYRANK == 0 ) then 
  call BoxMuller2(g_rsite,2*global_num_sites*NMAT*NMAT)
  call genrand_real3(g_rlink)

  g_rsite=g_rsite * 0.01d0 !/ mass_square_phi
  num=0
  do s=1,global_num_sites
    do i=1,NMAT
      do j=1,NMAT
        num=num+1
        if ( i.ne.NMAT .or. j.ne.NMAT ) then
          G_PHIMAT(i,j,s)=dcmplx(g_rsite(2*num-1))+(0d0,1d0)*dcmplx(g_rsite(2*num))
        endif
      enddo
    enddo
    G_PhiMat(NMAT,NMAT,s)=(0d0,0d0)
    do i=1,NMAT-1
      G_PhiMat(NMAT,NMAT,s)=G_PhiMat(NMAT,NMAT,s)-G_PhiMat(i,i,s)
    enddo
  enddo

  ! random number must be sufficiently small
  if( m_omega == 0 ) then 
    g_rlink=g_rlink * ( 1d0/dble(NMAT*NMAT) )
  else
    g_rlink=g_rlink * ( 1d0/dble(NMAT*NMAT*m_omega) )
  endif
  G_AMAT=(0d0,0d0)
  do l=1,global_num_links
    do a=1,dimG
      G_AMAT(:,:,l)=G_AMAT(:,:,l)+g_rlink(a,l)*TMAT(:,:,a)
    enddo
  enddo
  
  do l=1,global_num_links
    call matrix_exp(G_UMAT(:,:,l),(0d0,1d0)*G_AMAT(:,:,l))
  enddo
endif

do s=1,global_num_sites
  if( MYRANK == 0 ) then
    tmp_PhiMat = g_PhiMat(:,:,s)
  endif
  rank=local_site_of_global(s)%rank_ 
  ls=local_site_of_global(s)%label_
  if( MYRANK == 0 ) then
    if( rank == 0 ) then 
      PhiMat(:,:,ls) = tmp_PhiMat
    else
      call MPI_SEND(tmp_PhiMat,NMAT*NMAT,MPI_DOUBLE_COMPLEX,rank,s,MPI_COMM_WORLD,IERR)
    endif
  else
    if( rank /= 0 ) then 
      call MPI_RECV(PhiMat(:,:,ls),NMAT*NMAT,MPI_DOUBLE_COMPLEX,0,s,MPI_COMM_WORLD,ISTATUS,IERR)
    endif
  endif
enddo


do l=1,global_num_links
  if( MYRANK == 0 ) then
    tmp_UMat = g_UMat(:,:,s)
  endif
  rank=local_link_of_global(l)%rank_ 
  ll=local_link_of_global(l)%label_
  if( MYRANK == 0 ) then
    if( rank == 0 ) then 
      UMat(:,:,ll) = tmp_UMat
    else
      call MPI_SEND(tmp_UMat,NMAT*NMAT,MPI_DOUBLE_COMPLEX,rank,global_num_sites+l,MPI_COMM_WORLD,IERR)
    endif
  else
    if( rank /= 0 ) then 
      call MPI_RECV(UMat(:,:,ll),NMAT*NMAT,MPI_DOUBLE_COMPLEX,0,global_num_sites+l,MPI_COMM_WORLD,ISTATUS,IERR)
    endif
  endif
enddo
      
  

#else
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
#endif

end subroutine set_random_config


end module initialization
