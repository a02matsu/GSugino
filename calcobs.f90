module global_calcobs
implicit none

character(128), parameter :: PARAFILE="parameters_calcobs.dat"
character(128) :: MEDFILE
integer, parameter :: num_calcobs=2 ! 考えているobservableの数
integer :: trig_obs(1:num_calcobs)

double precision :: Sb, TrX2

integer, parameter :: N_MEDFILE=100
integer, parameter :: N_PARAFILE=101

integer :: trig_trx2
integer :: trig_Sb

end module global_calcobs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
use global_parameters
use global_calcobs
use simulation
use parallel
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! variables
complex(kind(0d0)), allocatable :: UMAT(:,:,:) ! unitary link variables
complex(kind(0d0)), allocatable :: PHIMAT(:,:,:) ! complex scalar at sites
complex(kind(0d0)), allocatable :: tmpMAT(:,:) 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! FOR PARALLEL
type(SITE_DIST), allocatable,save :: local_site_list(:) !(0:NPROCS-1)

! for type genrand_srepr
!integer :: MPI_GENRAND_SREPR, MPI_LOCAL_LABEL, MPI_GENRAND_STATE
integer IBLOCK1(1),IDISP1(1),ITYPE1(1)
integer IBLOCK2(1),IDISP2(1),ITYPE2(1)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: seed
integer :: l,ll,s,ls,tag,rank,i, ite

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

allocate(local_site_list(0:NPROCS-1) )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call set_theory_parameters(seed)
! set the structure constant of SU(NMAT)
call set_NZF
! set global data of the simplicial complex
call set_global_simplicial_complex 
! set global parameters
overall_factor=dble(NMAT)/(2d0*LatticeSpacing*LatticeSpacing)
mass_square_phi = phys_mass_square_phi * (LatticeSpacing*latticeSpacing)
! set simulation data
!   ここで、inputfileを読み込み、各コアに必要な情報を設定して、
!   Rank=0にたまっているalphaを振り分ける。
!   この段階で、
!      num_sites,num_links,num_faces
!      local_{site,link,face}_of_globalの(1:num_{sites,links,faces})が確定する。
call set_simulation_parameters(local_site_list)

! set local data
call set_local_data(local_site_list)

allocate( UMAT(1:NMAT,1:NMAT,1:num_necessary_links) )
allocate( PHIMAT(1:NMAT,1:NMAT, 1:num_necessary_sites) )
allocate( tmpMAT(1:NMAT,1:NMAT) )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! read parameters_calcobs.dat
if( MYRANK == 0 ) then
  open(N_PARAFILE, file=PARAFILE, status='OLD',action='READ')
  read(N_PARAFILE,*) MEDFILE
  do i=1,num_calcobs
    read(N_PARAFILE,*) trig_obs(i)
  enddo
endif
call MPI_BCAST(trig_obs, num_calcobs, MPI_INTEGER,0,MPI_COMM_WORLD,IERR)


if( MYRANK == 0 ) then
  open(N_MEDFILE, file=MEDFILE, status='OLD',action='READ',form='unformatted')
  !call read_basic_info_from_medfile(N_MEDFILE) ! in output.f90
endif

do
  call read_config_from_medfile(Umat,PhiMat,ite,N_MEDFILE)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! calculate observables
  if( MYRANK == 0 ) then
    write(*,*) ite
  endif
  if( trig_obs(1) == 0 ) then 
    call calc_TrX2(TrX2,PhiMat)
    if( MYRANK == 0 ) write(*,*) TrX2
  endif
  if( trig_obs(2) == 0 ) then 
    call calc_bosonic_action(Sb,Umat,PhiMat)
    if( MYRANK == 0 ) write(*,*) Sb
  endif


enddo

end program main




