!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
module initialization_calcobs
use global_parameters
!use global_calcobs
implicit none

!type A_in_B
!  integer :: num_ ! number of elements
!  integer, allocatable :: label_(:) ! global link label
!  complex(kind(0d0)), allocatable :: val_(:) ! 1~num_
!end type A_in_B


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! variables
complex(kind(0d0)), allocatable :: UMAT(:,:,:) ! unitary link variables
complex(kind(0d0)), allocatable :: PHIMAT(:,:,:) ! complex scalar at sites
!complex(kind(0d0)), allocatable :: tmpMAT(:,:) 
!!!
complex(kind(0d0)), allocatable :: Geta_eta(:,:,:,:,:,:)
!(1:dimG,1:num_global_sites,1:dimG,1:num_sites)
complex(kind(0d0)), allocatable :: Glambda_eta(:,:,:,:,:,:)
!(1:dimG,1:num_global_links,1:dimG,1:num_sites)
complex(kind(0d0)), allocatable :: Gchi_eta(:,:,:,:,:,:)
!(1:dimG,1:num_global_faces,1:dimG,1:num_sites)
!!!
complex(kind(0d0)), allocatable :: Geta_lambda(:,:,:,:,:,:)
!(1:dimG,1:num_global_sites,1:dimG,1:num_links)
complex(kind(0d0)), allocatable :: Glambda_lambda(:,:,:,:,:,:)
!(1:dimG,1:num_global_links,1:dimG,1:num_links)
complex(kind(0d0)), allocatable :: Gchi_lambda(:,:,:,:,:,:)
!(1:dimG,1:num_global_faces,1:dimG,1:num_links)
!!!
complex(kind(0d0)), allocatable :: Geta_chi(:,:,:,:,:,:)
!(1:dimG,1:num_global_sites,1:dimG,1:num_faces)
complex(kind(0d0)), allocatable :: Glambda_chi(:,:,:,:,:,:)
!(1:dimG,1:num_global_links,1:dimG,1:num_faces)
complex(kind(0d0)), allocatable :: Gchi_chi(:,:,:,:,:,:)
!(1:dimG,1:num_global_faces,1:dimG,1:num_faces)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! FOR PARALLEL
type(SITE_DIST), allocatable,save :: local_site_list(:) !(0:NPROCS-1)

! for type genrand_srepr
!integer :: MPI_GENRAND_SREPR, MPI_LOCAL_LABEL, MPI_GENRAND_STATE
integer IBLOCK1(1),IDISP1(1),ITYPE1(1)
integer IBLOCK2(1),IDISP2(1),ITYPE2(1)

contains
!!!
subroutine initialization
use global_parameters
!use global_calcobs
implicit none

integer seed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!call set_simulation_parameters(local_site_list)

! (2020/06/17変更)
! 分割情報は別ファイルにまとめて、siteの割り振りを決定
call set_mpi_distribution(local_site_list)
!   Rank=0にたまっているsiteのデータを使って、
!      num_sites,num_links,num_faces
!      local_{site,link,face}_of_globalの(1:num_{sites,links,faces})
!   を確定させる
call set_local_data(local_site_list)

allocate( UMAT(1:NMAT,1:NMAT,1:num_necessary_links) )
allocate( PHIMAT(1:NMAT,1:NMAT, 1:num_necessary_sites) )
!allocate( tmpMAT(1:NMAT,1:NMAT) )

!!!
allocate( Geta_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_necessary_sites) )
allocate( Glambda_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_necessary_sites) )
allocate( Gchi_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_necessary_sites) )
!!!
allocate( Geta_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_necessary_links) ) 
allocate( Glambda_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_necessary_links) )
allocate( Gchi_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_necessary_links) )
!!!
allocate( Geta_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_necessary_faces) ) 
allocate( Glambda_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_necessary_faces) )
allocate( Gchi_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_necessary_faces) )

!num_fermion=(global_num_sites+global_num_links+global_num_faces)*(NMAT*NMAT-1)
!num_sitelink=(global_num_sites+global_num_links)*(NMAT*NMAT-1)

end subroutine initialization


end module initialization_calcobs

