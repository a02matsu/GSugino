program main
use mt95
use global_parameters
use initialization
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
!integer :: MPI_GENRAND_SREPR, MPI_LOCAL_LABEL
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
! set the structure constant of SU(NMAT)
  call set_NZF
! set the data of the simplicial complex
!  "sizeD" is set here
  call set_sc 
  call set_alpha_beta
! remove comment if test of the module is needed
!    call test_module_simplicial_complex(sc)

! set necesarry data in each node
  call set_sub_SC

#ifdef PARALLEL
! set local_alpha_s/l/f and local_beta_f
  call set_local_alpha_beta

! 以下の量をlocalな量で再定義
!   num_sites/links/faces
!   linktip_from_s, linkorg_to_s
!   links_in_f
!   face_in_l
!   sites_in_f
!   alpha_s/l/f
!   beta_f
  call switch_globalnum_to_localnum
#endif
 

! initialize the size of the variables
  allocate( UMAT(1:NMAT,1:NMAT,1:num_necessary_links) )
  allocate( PHIMAT(1:NMAT,1:NMAT, 1:num_necessary_sites) )

! set the seed of the random number generators
! The default value is set in the parameter file.
! When fix_seed=2, seed is set by the system time.
  if (fix_seed==2 .or. fix_seed==0) then
    call system_clock(count=time)
    seed=time
  endif
  if( fix_seed == 0 ) then
    call genrand_init( put=state_mt95 )
  else
#ifdef PARALLEL
    seed=seed/(MYRANK+1)
    !call MPI_BCAST(seed,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
#endif
    call genrand_init( put=seed )
  endif
! at this stage, seed is fixed to some value

! set the variables depending on simulation_mode and test_mode
  !if (test_mode==1 .or. new_config==1) then
  if (new_config==1) then
    call set_random_config(UMAT,PHIMAT) 
    !do s=1,num_sites
      !write(*,*) global_site_of_local(s),PhiMat(:,:,s)
    !enddo
    !call stop_for_test

    !do s=1,num_sites
    !tmp=(0d0,0d0)
      !do i=1,NMAT
        !do j=1,NMAT
          !tmp=tmp+abs(PhiMat(i,j,s)*conjg(PhiMat(i,j,s)))
        !enddo
      !enddo
      !write(*,*) MYRANK,s,tmp
    !enddo
!
    !do l=1,1
      !do i=1,NMAT
        !do j=1,NMAT
          !write(*,*) MYRANK,l,i,j,UMAT(i,j,l)
        !enddo
      !enddo
    !enddo
    total_ite=0
  else
    call read_config(total_ite,UMAT,PHIMat,state_mt95) 
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

