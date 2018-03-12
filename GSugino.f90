program main
use mt95
use global_parameters
!use initialization
use check_routines
use SUN_generators, only : make_traceless_matrix_from_modes
use simulation
#ifdef PARALLEL
use parallel
#endif

implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! variables
complex(kind(0d0)), allocatable :: UMAT(:,:,:) ! unitary link variables
complex(kind(0d0)), allocatable :: PHIMAT(:,:,:) ! complex scalar at sites

integer :: total_ite ! total iteration
integer :: seed, time
type(genrand_state) :: state_mt95 ! state of mt95
double precision :: tmp

integer :: s,l,i,j

integer num_para, iargc

#ifdef PARALLEL
! for type genrand_srepr
!integer :: MPI_GENRAND_SREPR, MPI_LOCAL_LABEL, MPI_GENRAND_STATE
integer IBLOCK1(1),IDISP1(1),ITYPE1(1)
integer IBLOCK2(1),IDISP2(1),ITYPE2(1)

call MPI_INIT(IERR)
call MPI_COMM_SIZE(MPI_COMM_WORLD,NPROCS,IERR)
call MPI_COMM_RANK(MPI_COMM_WORLD,MYRANK,IERR)

! setup MPI_GENRAND_SREPR
IBLOCK1(1)=1
IDISP1(1)=0
ITYPE1(1)=MPI_CHARACTER
call MPI_TYPE_STRUCT(1,IBLOCK1,IDISP1,ITYPE1,MPI_GENRAND_SREPR,IERR)
call MPI_TYPE_COMMIT(MPI_GENRAND_SREPR,IERR)

! setup MPI_LOCAL_LABEL
IBLOCK2(1)=2
IDISP2(1)=0
ITYPE2(1)=MPI_INTEGER
call MPI_TYPE_STRUCT(1,IBLOCK2,IDISP2,ITYPE2,MPI_LOCAL_LABEL,IERR)
call MPI_TYPE_COMMIT(MPI_LOCAL_LABEL,IERR)

#endif

num_para = iargc()
if(num_para > 0) then
 call getarg(1,PAR_FILE_NAME)
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! set the parameters in the module "global_parameters"
! read input parameter from the parameter file
  call set_parameters(seed)
!  write(*,*) NMAT
!  write(*,*) LatticeSpacing
!  write(*,*) SC_FILE_NAME
!  write(*,*) FsubSC
!  write(*,*) ALPHA_BETA
!  write(*,*) test_mode
!  write(*,*) force_measurement
!  write(*,*) new_config
!  write(*,*) fix_seed
!  write(*,*) reset_ite
!  write(*,*) read_alpha
!  write(*,*) save_med_step
!  write(*,*) save_config_step
!  write(*,*) obs_step
!  write(*,*) seed
!  write(*,*) m_omega
!  write(*,*) phys_mass_square_phi
!  write(*,*) mass_f
!  write(*,*) Remez_factor4
!  write(*,*) Remez_factor8
!  write(*,*) epsilon
!  write(*,*) CG_max
!  write(*,*) num_ite
!  write(*,*) Ntau_base
!  write(*,*) FBPhi_ratio
!  write(*,*) FBA_ratio
!  write(*,*) Dtau
!  write(*,*) Fconfigin
!  write(*,*) Foutput
!  write(*,*) Fconfigout
!  write(*,*) Fmedconf

  !call stop_for_test
! set the structure constant of SU(NMAT)
  call set_NZF
! set the data of the simplicial complex
!  "sizeD" is set here
  !call set_sc 
#ifdef PARALLEL
  call set_sc_alpha_beta_parallel
#else
  call set_sc_alpha_beta
#endif
! remove comment if test of the module is needed
!    call test_module_simplicial_complex(sc)

! set necesarry data in each node
#ifdef PARALLEL
  call set_sub_SC
  call switch_globalnum_to_localnum
  call set_local_alpha_beta
#endif

if( check_sub_sc == 1 ) then
#ifdef PARALLEL
  call check_local_sc
  call stop_for_test
#endif
endif
 

!write(*,*) MYRANK, num_sites,num_links,num_faces,num_necessary_sites,num_necessary_links, num_necessary_faces
! initialize the size of the variables
  allocate( UMAT(1:NMAT,1:NMAT,1:num_necessary_links) )
  allocate( PHIMAT(1:NMAT,1:NMAT, 1:num_necessary_sites) )

! set the seed of the random number generators
! The default value is set in the parameter file.
! When fix_seed=2, seed is set by the system time.
  if (fix_seed==2 .or. fix_seed==0) then
    call system_clock(count=time)
    seed=time
#ifdef PARALLEL
    seed=seed/(MYRANK+1)
    call genrand_init( put=seed )
#endif
  else
#ifdef PARALLEL
    seed=seed/(MYRANK+1)
    call genrand_init( put=seed )
#endif
  endif
! at this stage, seed is fixed to some value

! set the variables depending on simulation_mode and test_mode
  !if (test_mode==1 .or. new_config==1) then
  if (new_config==1) then
    call set_random_config(UMAT,PHIMAT) 
    total_ite=0
  else
    call read_config(total_ite,UMAT,PHIMat,state_mt95,seed)
  endif
  if( fix_seed == 0 ) then
    call genrand_init( put=state_mt95 )
  else
    call genrand_init( put=seed )
  endif
  if( reset_ite == 1 ) then
    total_ite=0
  endif
! set Remez data
    call set_Remez_data
!! check unitaryty
    !call check_unitary
!! check the distance from 1_N 
    !call check_distance
!! check mode expansion
    !call check_Ta_expansion
    !stop
!! check anti-symmetricity of Dirac and Hermiticity of D\dag D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!do s=1,num_sites
!call make_traceless_matrix_from_modes(PhiMat(:,:,s),NMAT,Phi(:,s))
!enddo

  !if (writedown_mode==1) then
    !call writedown_config_action_and_fores(UMAT,PhiMat,seed)
    !stop
  !endif

  if (test_mode==1) then
    !call check_Dirac(UMAT,Phi)
    call test_hamiltonian(UMAT,PhiMat,seed)
  else
    call HybridMonteCarlo(UMAT,PhiMat,seed,total_ite)
  endif

#ifdef PARALLEL
  call MPI_FINALIZE(IERR)
#endif

end program main


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! read the preevious configuration
subroutine read_config(total_ite,UMAT,PhiMat,state_mt95,seed)
use global_parameters
use mt95
#ifdef PARALLEL
use parallel
use global_subroutines, only : syncronize_bosons,site_abs,link_abs,face_abs
#endif
implicit none

integer, intent(inout) :: total_ite
type(genrand_state), intent(inout) :: state_mt95
integer :: num_core
complex(kind(0d0)), intent(inout) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(inout) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
integer, intent(in) :: seed
#ifdef PARALLEL
complex(kind(0d0)) ::tmpMat(1:NMAT,1:NMAT)
integer :: s,l,tag,counter,rank
type(genrand_srepr) :: tmp_char_mt95
integer :: ISEND(1:NPROCS-1),err(1:NPROCS-1),err0
!double precision :: rtmp
#endif

type(genrand_srepr) :: char_mt95

#ifdef PARALLEL
if( MYRANK == 0 ) then
  open(IN_CONF_FILE, file=Fconfigin, status='OLD',action='READ',form='unformatted')
  read(IN_CONF_FILE) total_ite
endif
call MPI_BCAST(total_ite,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
!!!!!!!!!!!!!!!!!!!!
!! Umat
do l=1,global_num_links
  rank=local_link_of_global(l)%rank_
  tag=l
  if( MYRANK == 0 ) then
    read(IN_CONF_FILE) tmpmat
    if( rank == 0 ) then
      UMAT(:,:,local_link_of_global(l)%label_) = tmpmat
    else
      call MPI_SEND(tmpmat,NMAT*NMAT,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,IERR)
    endif
  elseif( MYRANK == rank ) then
    call MPI_RECV(UMAT(:,:,local_link_of_global(l)%label_),&
      NMAT*NMAT,MPI_DOUBLE_COMPLEX, 0,tag,MPI_COMM_WORLD,ISTATUS,IERR)
  endif
enddo
!!!!!!!!!!!!!!!!!!!!
!! PhiMat
do s=1,global_num_sites
  tag=global_num_links+s
  rank=local_site_of_global(s)%rank_ 
  if( MYRANK == 0 ) then
    read(IN_CONF_FILE) tmpmat
    if( rank == 0 ) then
      PHIMAT(:,:,local_site_of_global(s)%label_) = tmpmat
    else
      call MPI_SEND(tmpmat,NMAT*NMAT,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,IERR)
    endif
  elseif( MYRANK == rank ) then
    call MPI_RECV(PHIMAT(:,:,local_site_of_global(s)%label_),&
      NMAT*NMAT,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,ISTATUS,IERR)
  endif
enddo

!rtmp=site_abs(PhiMat(:,:,1:num_sites))
!if( MYRANK == 0 ) then
!  write(*,'(a,E25.18)') "|PhiMat|^2=",rtmp
!endif
!rtmp=link_abs(UMat(:,:,1:num_links))
!if( MYRANK == 0 ) then
!  write(*,'(a,E25.18)') "|UMat|^2=",rtmp
!endif
!call stop_for_test
!!!!!!!!!!!!!!!!!!!!!!!!1
!! syncronize 
call syncronize_bosons(UMAT,Phimat)

!!!!!!!!!!!!!!!!!!!!!!!!!
!! random seed
if( fix_seed==0 ) then 
  if( MYRANK == 0 ) then
    read(IN_CONF_FILE,iostat=err0) tmp_char_mt95
    do counter=1,NPROCS-1
      read(IN_CONF_FILE,iostat=err(counter)) tmp_char_mt95
      call MPI_SEND(err(counter),1,MPI_INTEGER,counter,counter,MPI_COMM_WORLD,IERR)
      if( err(counter) == 0 ) then 
        if( counter == 0 ) then
          char_mt95=tmp_char_mt95
        else
          call MPI_ISEND(tmp_char_mt95,1,MPI_GENRAND_SREPR,counter,counter,MPI_COMM_WORLD,ISEND(counter),IERR)
        endif
      endif
    enddo
    do counter=1,NPROCS-1
      if( err(counter)==0 ) call MPI_WAIT(ISEND(counter),ISTATUS,IERR)
    enddo
  else
    call MPI_RECV(err0,1,MPI_INTEGER,0,MYRANK,MPI_COMM_WORLD,ISTATUS,IERR)
    if( err0 == 0 ) then
      call MPI_RECV(char_mt95,1,MPI_GENRAND_SREPR,0,MYRANK,MPI_COMM_WORLD,ISTATUS,IERR)
    endif
  endif
!  if( err0==0 ) then
!    state_mt95=char_mt95
!    call genrand_init( put=state_mt95 )
!  else
!    call genrand_init( put=seed )
!  endif
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if( MYRANK == 0 ) then
  close(IN_CONF_FILE)
endif

#else
open(IN_CONF_FILE, file=Fconfigin, status='OLD',action='READ',form='unformatted')
read(IN_CONF_FILE) total_ite
read(IN_CONF_FILE) UMAT
read(IN_CONF_FILE) PHIMat
read(IN_CONF_FILE) char_mt95
close(IN_CONF_FILE)
if( fix_seed==0 ) then
  state_mt95=char_mt95
  call genrand_init( put=state_mt95 )
else
  call genrand_init( put=seed )
endif
#endif

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
use global_parameters
use mt95
#ifdef PARALLEL
use parallel
use global_subroutines, only : syncronize_bosons
#endif
implicit none

complex(kind(0d0)), intent(inout) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(inout) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) :: TMAT(1:NMAT,1:NMAT,1:NMAT**2-1)
double precision :: rsite(1:2*NMAT*NMAT*num_sites) ! for PHI
double precision :: rlink(1:dimG,1:num_links) ! for UMAT
complex(kind(0d0)) :: AMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: trace
integer :: s,l,a,f,i,j,num

#ifdef PARALLEL
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

#ifdef PARALLEL
!! テスト用に、シングルコアと同じconfigurationを用意する。
if( PARATEST== 1 ) then 
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
          !if ( i.ne.NMAT .or. j.ne.NMAT ) then
            G_PHIMAT(i,j,s)=dcmplx(g_rsite(2*num-1))+(0d0,1d0)*dcmplx(g_rsite(2*num))
          !endif
        enddo
      enddo
      trace=(0d0,0d0)
      do i=1,NMAT
        trace=trace+G_PHIMAT(i,i,s)
      enddo
      do i=1,NMAT
        G_PHIMAT(i,i,s)=G_PHIMAT(i,i,s)-trace/dcmplx(dble(NMAT))
      enddo
    enddo
  
    ! random number must be sufficiently small
    if( m_omega == 0 .or. m_omega == -1 ) then 
      g_rlink=g_rlink * ( 1d0/dble(NMAT*NMAT) )
    else
      g_rlink=g_rlink * ( 1d0/dble(NMAT*NMAT*m_omega) )
    endif
    G_AMAT=(0d0,0d0)
    do l=1,global_num_links
      !trace=(0d0,0d0)
      !do i=1,NMAT
      !  trace=trace+g_rlink(i,l)
      !enddo
      !do i=1,NMAT
      !  G_AMat(i,i,l)=g_rlink(i,l)-trace/dble(NMAT)
      !enddo
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
      if( MYRANK == rank ) then 
        call MPI_RECV(PhiMat(:,:,ls),NMAT*NMAT,MPI_DOUBLE_COMPLEX,0,s,MPI_COMM_WORLD,ISTATUS,IERR)
      endif
    endif
  enddo
  
  
  do l=1,global_num_links
    if( MYRANK == 0 ) then
      tmp_UMat = g_UMat(:,:,l)
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
      if( MYRANK == rank  ) then 
        call MPI_RECV(UMat(:,:,ll),NMAT*NMAT,MPI_DOUBLE_COMPLEX,0,global_num_sites+l,MPI_COMM_WORLD,ISTATUS,IERR)
      endif
    endif
  enddo
        
  call syncronize_bosons(UMAT,Phimat)
  return

else ! PARATEST==0 case
  call BoxMuller2(g_rsite,2*num_sites*NMAT*NMAT)
  call genrand_real3(g_rlink)
  
  g_rsite=g_rsite * 0.01d0 !/ mass_square_phi
  num=0
  do s=1,num_sites
    do i=1,NMAT
      do j=1,NMAT
        num=num+1
        !if ( i.ne.NMAT .or. j.ne.NMAT ) then
          PHIMAT(i,j,s)=dcmplx(g_rsite(2*num-1))+(0d0,1d0)*dcmplx(g_rsite(2*num))
        !endif
      enddo
    enddo
    trace=(0d0,0d0)
    do i=1,NMAT
      trace=trace+PhiMat(i,i,s)
    enddo
    do i=1,NMAT
      PhiMat(i,i,s)=PhiMat(i,i,s)-trace/dcmplx(dble(NMAT))
    enddo
  enddo
  
  ! random number must be sufficiently small
  if( m_omega == 0 .or. m_omega == -1 ) then 
    g_rlink=g_rlink * ( 1d0/dble(NMAT*NMAT) )
  else
    g_rlink=g_rlink * ( 1d0/dble(NMAT*NMAT*m_omega) )
  endif
  AMAT=(0d0,0d0)
  do l=1,num_links
    do a=1,dimG
      AMAT(:,:,l)=AMAT(:,:,l)+g_rlink(a,l)*TMAT(:,:,a)
    enddo
  enddo
    
  do l=1,num_links
    call matrix_exp(UMAT(:,:,l),(0d0,1d0)*AMAT(:,:,l))
  enddo
  
  !write(*,*) "test"
  call syncronize_bosons(UMAT,Phimat)
endif

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
  if( m_omega == 0 .or. m_omega == -1 ) then 
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




