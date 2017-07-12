!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! module to define variables in main routine
module global_parameters
use simplicial_complex
#ifdef PARALLEL
use parallel
#endif
implicit none
type SITE_LINKDAT
  integer :: num_
  integer, allocatable :: labels_(:)
  integer, allocatable :: sites_(:)
end type SITE_LINKDAT

type INDEXDAT
  integer :: num_
  integer, allocatable :: b_(:)
  integer, allocatable :: c_(:)
  double precision, allocatable :: value_(:)
end type INDEXDAT

type FACEDAT
  integer :: num_ ! 要素の数
  integer, allocatable :: label_(:) ! faceのラベル label_(1:num_)
end type FACEDAT

type LINKDAT
  integer :: num_ ! linkの数
  integer, allocatable :: link_labels_(:) ! linkのラベル
  integer, allocatable :: link_dirs_(:) ! linkの方向
end type LINKDAT

type SITEDAT
  integer :: num_ ! 要素の数
  integer, allocatable :: label_(:) ! siteのラベル label_(1:num_)
end type SITEDAT

type VECT
  integer :: num_ ! 要素の数
  integer, allocatable :: labels_(:) ! ラベル label_(1:num_)
end type VECT

!type diffdiff_by_linkvals_in_face
!  integer :: num_links
!  complex(kind(0d0)), allocatable :: diffdiff_val(:,:,:,:,:,:)
!end type diffdiff_by_linkvals_in_face

! parameter for test hamiltonian
integer, parameter :: p1=0
integer, parameter :: p2=0
integer, parameter :: p3=0
integer, parameter :: p4=0
integer, parameter :: p5=0




! writedown_mode=1 : write down config, hamiltonian and forces
integer, parameter :: writedown_mode=0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! parameter file as an argument
!integer num_para, iargc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! constants
double precision, parameter :: PI=dacos(-1d0)
complex, parameter :: im_unit=(0d0,1d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! structure constant
integer :: NZF      ! number of nonzero elements of f(a,b,c)
! index of nonzero f(a,b,c)
integer, allocatable :: NZF_index(:,:)  ! NZF_index(1:3,1:NZF)
! value of nonzero f(a,b,c)
double precision, allocatable :: NZF_value(:) ! NZF_value(1:NZF)
! number of indices (b,c) with nonzero f(a,b,c) for given a
type(INDEXDAT), allocatable :: NZF_a(:) 
! number of nonzero elementes of f(a,b,e)*f(c,d,e)
integer :: NZFF

! index of nonzero f(a,b,e)*f(c,d,e)
integer, allocatable :: NZFF_index(:,:)  ! NZF_index(1:4,1:NZFF)
! value of nonzero f(a,b,e)*f(c,d,e)
double precision, allocatable :: NZFF_value(:) ! NZF_value(1:NZFF)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! data of the simplicial complex
integer :: global_num_sites ! number of sites
integer :: global_num_links ! number of links
integer :: global_num_faces ! number of faces
integer :: num_sites ! number of sites
integer :: num_links ! number of links
integer :: num_faces ! number of faces
type(SmpCom) :: SC
integer, allocatable :: link_org(:) ! l番目のリンクのorigin
integer, allocatable :: link_tip(:) ! l番目のリンクのtop
!! <s,・>
type(SITE_LINKDAT), allocatable :: linktip_from_s(:) ! sから出るリンクのtip
!! <・,t>
type(SITE_LINKDAT), allocatable :: linkorg_to_s(:) ! sに入るリンクのorigin
!!
type(LINKDAT), allocatable :: links_in_f(:) ! fに含まれるリンクのデータ
!! 
type(FACEDAT), allocatable :: face_in_l(:) ! lに含まれるface
!!
type(SITEDAT), allocatable :: sites_in_f(:) ! face f に含まれるサイト

double precision, allocatable :: alpha_s(:)
double precision, allocatable :: alpha_l(:)
double precision, allocatable :: alpha_f(:)
double precision, allocatable :: beta_f(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! data of the LOCAL simplicial complex
integer :: num_sub_SC
integer :: num_local_sites ! number of sites for each sub_SC
integer :: num_local_links ! number of links for each sub_SC
integer :: num_local_faces ! number of faces for each sub_SC


integer :: num_necessary_sites ! number of sites needed for each rank
integer :: num_necessary_links ! number of links needed for each rank
integer :: num_necessary_faces ! number of faces needed for each rank


integer, allocatable :: local_link_org(:) ! l番目のリンクのorigin
integer, allocatable :: local_link_tip(:) ! l番目のリンクのtop
!! <s,・>
type(VECT), allocatable :: local_linktip_from_s(:) ! sから出るリンクのtip(実はSITE_LINKDATである必要はない。)
!! <・,t>
type(VECT), allocatable :: local_linkorg_to_s(:) ! sに入るリンクのorigin
!!
type(LINKDAT), allocatable :: local_links_in_f(:) ! fに含まれるリンクのデータ
!! 
type(FACEDAT), allocatable :: local_face_in_l(:) ! lに含まれるface
!!
type(SITEDAT), allocatable :: local_sites_in_f(:) ! face f に含まれるサイト

double precision, allocatable :: local_alpha_s(:)
double precision, allocatable :: local_alpha_l(:)
double precision, allocatable :: local_alpha_f(:)
double precision, allocatable :: local_beta_f(:)

integer, allocatable :: global_site_of_local(:) ! map from local site to global site
integer, allocatable :: global_link_of_local(:) ! map from local link to global link
integer, allocatable :: global_face_of_local(:) ! map from local face to global face

! (rank, local label) corresponding to a global label
#ifdef PARALLEL
type LOCAL_LABEL
  integer :: rank_ ! nodes
  integer :: label_ ! label
end type LOCAL_LABEL

type(LOCAL_LABEL), allocatable :: local_site_of_global(:) ! map from global site to local site
type(LOCAL_LABEL), allocatable :: local_link_of_global(:) ! map from global link to local link
type(LOCAL_LABEL), allocatable :: local_face_of_global(:) ! map from global face to local face


integer :: num_send_sites, num_recv_sites
integer :: num_send_links, num_recv_links
integer :: num_send_faces, num_recv_faces

type(LOCAL_LABEL), allocatable :: send_sites(:) 
type(LOCAL_LABEL), allocatable :: send_links(:) 
type(LOCAL_LABEL), allocatable :: send_faces(:) 

type(LOCAL_LABEL), allocatable :: recv_sites(:) 
type(LOCAL_LABEL), allocatable :: recv_links(:) 
type(LOCAL_LABEL), allocatable :: recv_faces(:) 
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! controle paramaeters of simulation 
integer :: test_mode !! 0: Simulation mode, 1: Test mode, 
integer :: fix_seed !! 0: use the value at the end of the previous simulation
                    !! 1: fix seed to the value in the parameter file
                    !! 2: determine by the system time
integer :: new_config !! 0: Old config, 1: New config
integer :: read_alpha !! 0: alpha,beta=1 1: read ALPHA_BETA_FILE
integer :: save_med_step !! steps to write out configurations
integer :: obs_step !! steps to compute observables
integer :: reset_ite !! 0:iteration number continue, 1:iteration number reset 
integer :: save_config_step !! save config by this step. 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! paramaeters in simulation 
integer :: NMAT     ! rank of the gauge group
integer :: dimG     ! dimension of su(NMAT)=NMAT**2-1
integer :: sizeD    ! size of the Dirac matrix
integer :: m_omega ! the integer m to avoid degenerated vacua
!integer :: seed !! seed for random number generator
integer :: num_ite ! number of iteration
integer :: Ntau ! Ntau
double precision :: Dtau     ! Dtau base 
double precision :: Dtau_phi ! Dtau for Phi
double precision :: Dtau_A   ! Dtau for A
double precision :: R_phi    ! Dtau_phi = R_phi * Dtau
double precision :: R_A      ! Dtau_A = R_A * Dtau
double precision :: phys_mass_square_phi ! mass^2 of Phi in mass-dim 2
double precision :: mass_square_phi ! dimensionless mass^2 of Phi
double precision :: mass_f   ! mass of fermions
double precision :: LatticeSpacing ! lattice spacing


double precision :: one_ov_2g2N ! 1/2g^2N=1/2a^2
double precision :: overall_factor ! N/2g^2N=N/2a^2
double precision :: maximal_dist ! the largest allowed value of |U_f|



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! parameters for Remez algorithm
character(128), parameter :: Remez_1ovminus4="remez_Q24_1ov-4.dat"
character(128), parameter :: Remez_1ov8="remez_Q24_1ov8.dat" 
integer :: N_Remez4 ! number of the rational function in rational approximation 
integer :: N_Remez8 ! number of the rational function in rational approximation 
double precision :: Remez_min4
double precision :: Remez_max4
double precision :: Remez_min8
double precision :: Remez_max8
double precision, allocatable :: Remez_alpha4(:)
double precision, allocatable :: Remez_beta4(:)
double precision, allocatable :: Remez_alpha8(:)
double precision, allocatable :: Remez_beta8(:)
double precision :: Remez_factor4
double precision :: Remez_factor8
double precision :: epsilon ! range to stop CG solver
integer :: CG_max ! maximum number of CG iteration

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Files
character(128) :: PAR_FILE_NAME="parameters.dat" ! parameter file
character(128) :: SC_FILE_NAME ! structure of simplicial complex
character(128) :: ALPHA_BETA ! parameterfile of alpha_s, alpha_l, alpha_f, beta_f
character(128) :: Fconfigin ! configuration input 
character(128) :: Fconfigout ! configuration output in the end
character(128) :: Foutput ! output file
character(128) :: Fmedconf ! configurations during the simulation 
character(128) :: FsubSC ! file to describe sub simplicial complex


integer, parameter :: PAR_FILE = 10 ! parameter file
integer, parameter :: SC_FILE = 11 ! data of simplicial complex
integer, parameter :: ALPHA_BETA_FILE = 12 ! factor of the action
integer, parameter :: Remez4_FILE = 13 ! remez parameter for -1/4
integer, parameter :: Remez8_FILE = 14 ! remez parameter for 1/8
integer, parameter :: IN_CONF_FILE = 20 ! input configuration file
integer, parameter :: OUTPUT_FILE = 21 ! output file
integer, parameter :: OUT_CONF_FILE = 22 ! output configuration file
integer, parameter :: MED_CONF_FILE = 23 ! configurations during the simulation 
integer, parameter :: SUBSC_FILE = 24 ! file to describe sub simplicial complex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Control Algorighms
!integer, parameter :: USE_CG = 0 ! 0: do not use CG, 1: use CG

contains 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11!!!!
!! SUBROUTINES for initialization
subroutine set_parameters(seed)
implicit none

integer, intent(inout) :: seed
! Open parameter file
#ifdef PARALLEL
if (MYRANK==0) then
#endif
open(PAR_FILE, file=PAR_FILE_NAME, status='old',action='READ')
  read(PAR_FILE,*) NMAT
  read(PAR_FILE,*) LatticeSpacing
  read(PAR_FILE,*) SC_FILE_NAME
  read(PAR_FILE,*) FsubSC
  read(PAR_FILE,*) ALPHA_BETA
  read(PAR_FILE,*) test_mode
  read(PAR_FILE,*) new_config
  read(PAR_FILE,*) fix_seed
  read(PAR_FILE,*) reset_ite
  read(PAR_FILE,*) read_alpha
  read(PAR_FILE,*) save_med_step
  read(PAR_FILE,*) save_config_step
  read(PAR_FILE,*) obs_step
  read(PAR_FILE,*) seed
  read(PAR_FILE,*) m_omega
  read(PAR_FILE,*) phys_mass_square_phi
  read(PAR_FILE,*) mass_f
  read(PAR_FILE,*) Remez_factor4
  read(PAR_FILE,*) Remez_factor8
  read(PAR_FILE,*) epsilon
  read(PAR_FILE,*) CG_max
  read(PAR_FILE,*) num_ite
  read(PAR_FILE,*) Ntau
  read(PAR_FILE,*) Dtau
  read(PAR_FILE,*) R_phi
  read(PAR_FILE,*) R_A
  read(PAR_FILE,*) Fconfigin
  read(PAR_FILE,*) Foutput
  read(PAR_FILE,*) Fconfigout
  read(PAR_FILE,*) Fmedconf
close(PAR_FILE)

dimG=NMAT*NMAT-1
Dtau_phi = R_phi * Dtau
Dtau_A = R_A * Dtau
mass_square_phi = phys_mass_square_phi * (LatticeSpacing*latticeSpacing)

!one_ov_2g2N=1d0/(2d0*LatticeSpacing*LatticeSpacing)
overall_factor=dble(NMAT)/(2d0*LatticeSpacing*LatticeSpacing)

if (NMAT<=4) then
  maximal_dist = 2d0*dsqrt(2d0)
else
  maximal_dist = 2d0*sqrt(dble(NMAT))*sin(3.1415926535898d0/dble(NMAT))
endif

#ifdef PARALLEL
endif

! send all parameters to all the other rank
!  read(PAR_FILE,*) NMAT
call MPI_BCAST(NMAT,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) LatticeSpacing
call MPI_BCAST(LatticeSpacing,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) SC_FILE_NAME
call MPI_BCAST(SC_FILE_NAME,128,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) FsubSC
call MPI_BCAST(FsubSC,128,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) ALPHA_BETA
call MPI_BCAST(ALPHA_BETA,128,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) test_mode
call MPI_BCAST(test_mode,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) new_config
call MPI_BCAST(new_config,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) fix_seed
call MPI_BCAST(fix_seed,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) reset_ite
call MPI_BCAST(reset_ite,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) read_alpha
call MPI_BCAST(read_alpha,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) save_med_step
call MPI_BCAST(save_med_step,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) save_config_step
call MPI_BCAST(save_config_step,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) obs_step
call MPI_BCAST(obs_step,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) seed
call MPI_BCAST(seed,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) m_omega
call MPI_BCAST(m_omega,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) phys_mass_square_phi
call MPI_BCAST(phys_mass_square_phi,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) mass_f
call MPI_BCAST(mass_f,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) Remez_factor4
call MPI_BCAST(Remez_factor4,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) Remez_factor8
call MPI_BCAST(Remez_factor8,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) epsilon
call MPI_BCAST(epsilon,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) CG_max
call MPI_BCAST(CG_max,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) num_ite
call MPI_BCAST(num_ite,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) Ntau
call MPI_BCAST(Ntau,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) Dtau
call MPI_BCAST(Dtau,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) R_phi
call MPI_BCAST(R_phi,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) R_A
call MPI_BCAST(R_A,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) Fconfigin
call MPI_BCAST(Fconfigin,128,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) Foutput
call MPI_BCAST(Foutput,128,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) Fconfigout
call MPI_BCAST(Fconfigout,128,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) Fmedconf
call MPI_BCAST(Fmedconf,128,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
  !dimG=NMAT*NMAT-1
call MPI_BCAST(dimG,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  !Dtau_phi = R_phi * Dtau
call MPI_BCAST(Dtau_phi,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !Dtau_A = R_A * Dtau
call MPI_BCAST(Dtau_A,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !mass_square_phi = phys_mass_square_phi * (LatticeSpacing*latticeSpacing)
call MPI_BCAST(mass_square_phi,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !overall_factor=dble(NMAT)/(2d0*LatticeSpacing*LatticeSpacing)
call MPI_BCAST(overall_factor,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !maximal_dist = 2d0*dsqrt(2d0)
call MPI_BCAST(maximal_dist,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

#endif

end subroutine set_parameters


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! set nonzero values of the structure constant
subroutine set_NZF
use SUN_generators
implicit none
integer :: a,b,c,i
integer :: num(1:dimG)
 
call get_num_nonzerof(NZF,NMAT)
allocate( NZF_index(1:3,1:NZF) )
allocate( NZF_value(1:NZF) )
call get_nonzerof(NZF_index,NZF_value,NZF,NMAT)

allocate( NZF_a(1:dimG) )

!! for f(a,b,e)*f(c,d,e)
call get_num_nonzero_ff(NZFF,NMAT)
allocate( NZFF_index(1:4,1:NZFF) )
allocate( NZFF_value(1:NZFF) )
call get_nonzero_ff( NZFF_index, NZFF_value, NZFF, NMAT)

!! count the number of nonzero components of f(a,b,c) with a fixed a
num=0
do i=1,NZF
  a=NZF_index(1,i)
  num(a)=num(a)+1
enddo

do a=1,dimG
!write(*,*) a,num(a)
NZF_a(a)%num_=num(a)
allocate( NZF_a(a)%b_(1:num(a)) )
allocate( NZF_a(a)%c_(1:num(a)) )
allocate( NZF_a(a)%value_(1:num(a)) )
enddo

num=0
do i=1,NZF
  a=NZF_index(1,i)
  b=NZF_index(2,i)
  c=NZF_index(3,i)
  num(a)=num(a)+1
  NZF_a(a)%b_(num(a))=b
  NZF_a(a)%c_(num(a))=c
  NZF_a(a)%value_(num(a))=NZF_value(i)
enddo
end subroutine set_NZF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_sc
use simplicial_complex
implicit none
integer l,origin,tip
integer s,f
integer :: fsize
integer, allocatable :: fsites(:), faces_l(:), sites_f(:)
integer :: FaceSize
integer, allocatable :: sites(:)
integer k
character(128) tmp
! open SC_FILE 
#ifdef PARALLEL
if(MYRANK==0) then
  open(SC_FILE, file=SC_FILE_NAME, status='old',action='READ')
  read(SC_FILE,*) num_sites
  read(SC_FILE,*) num_links
  read(SC_FILE,*) num_faces
endif
  call MPI_BCAST(num_sites,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(num_links,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(num_faces,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
#else
open(SC_FILE, file=SC_FILE_NAME, status='old',action='READ')
  read(SC_FILE,*) num_sites
  read(SC_FILE,*) num_links
  read(SC_FILE,*) num_faces
#endif

  allocate ( alpha_s(1:num_sites) )
  allocate ( alpha_l(1:num_links) )
  allocate ( alpha_f(1:num_faces) )
  allocate ( beta_f(1:num_faces) )

! initialize the simplicial complex
call init_smpcom(SC,num_sites,num_links,num_faces)

do l=1,num_links
#ifdef PARALLEL
if(MYRANK==0) then
  read(SC_FILE,*) origin,tip
endif
  call MPI_BCAST(origin,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(tip,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  !write(*,*) MYRANK,l,origin,tip
#else
  read(SC_FILE,*) origin,tip
#endif
  call put_link_sc(SC,l,origin,tip)
enddo

do f=1,num_faces
#ifdef PARALLEL
if(MYRANK==0) then
  read(SC_FILE,"(I1)",advance="no") fsize
endif
call MPI_BCAST(fsize,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
#else
  read(SC_FILE,"(I1)",advance="no") fsize
#endif

  allocate( fsites(1:fsize ) )

  ! set sites constructing i'th face
#ifdef PARALLEL
  if (MYRANK==0) then
    read(SC_FILE,*) (fsites(k),k=1,fsize)
  endif
  call MPI_BCAST(fsites,fsize,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  !write(*,*) MYRANK, f,fsites
#else
  read(SC_FILE,*) (fsites(k),k=1,fsize)
#endif
  call put_face_sc(SC,f,fsites)

  deallocate( fsites )
enddo
#ifdef PARALLEL
if (MYRANK==0) then
  close(SC_FILE)
endif
#else 
  close(SC_FILE)
#endif



!write(*,*) "==="
!write(*,*) sc%FaceSet_(1)%FaceSize_
!write(*,*) "==="

! set links
allocate(link_org(1:num_links))
allocate(link_tip(1:num_links))
do l=1,num_links
  call get_link_sc(sc,l,link_org(l),link_tip(l))
enddo

! tips of links from s
allocate(linktip_from_s(1:num_sites))
do s=1,num_sites
call get_links_from_s_sc(sc,s,&
  linktip_from_s(s)%labels_,&
  linktip_from_s(s)%sites_)
linktip_from_s(s)%num_=size(linktip_from_s(s)%labels_)
enddo

! origins of links to s
allocate(linkorg_to_s(1:num_sites))
do s=1,num_sites
call get_links_to_s_sc(sc,s,&
  linkorg_to_s(s)%labels_,&
  linkorg_to_s(s)%sites_)
linkorg_to_s(s)%num_=size(linkorg_to_s(s)%labels_)
enddo

! links included in a face f
allocate(links_in_f(1:num_faces))
do f=1,num_faces
  call get_links_in_face_sc(sc,f,&
    links_in_f(f)%num_,&
    sites,&
    links_in_f(f)%link_labels_,&
    links_in_f(f)%link_dirs_)
enddo

! faces included in a link l
allocate(face_in_l(1:num_links))
do l=1,num_links
  call get_faces_in_link_sc(sc,l,faces_l)
  face_in_l(l)%num_=size(faces_l)
  allocate(face_in_l(l)%label_(1:size(faces_l)))
  face_in_l(l)%label_=faces_l
enddo
  
! sites included in a face f
allocate(sites_in_f(1:num_faces))
do f=1,num_faces
  call get_face_components_sc(sc,f,sites_f)
  sites_in_f(f)%num_=size(sites_f)
  allocate(sites_in_f(f)%label_(1:size(sites_f)))
  sites_in_f(f)%label_=sites_f
enddo

  
! set size of Dirac matrix
sizeD=(dimG)*(num_sites+num_links+num_faces)

end subroutine set_sc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! set alpha and beta
subroutine set_alpha_beta
implicit none
integer s,l,f
#ifdef PARALLEL
if(MYRANK==0) then
#endif
if (read_alpha==1) then
open(ALPHA_BETA_FILE, file=ALPHA_BETA, status='old',action='READ')
  read(ALPHA_BETA_FILE,'()') ! skip 1 line
  do s=1,num_sites
    read(ALPHA_BETA_FILE,*) alpha_s(s)
  enddo
  read(ALPHA_BETA_FILE,'()') !! skip 1 line 
  do l=1,num_links
    read(ALPHA_BETA_FILE,*) alpha_l(l)
  enddo
  read(ALPHA_BETA_FILE,'()')  !! skip 1 line 
  do f=1,num_faces
    read(ALPHA_BETA_FILE,*) alpha_f(f)
  enddo
  read(ALPHA_BETA_FILE,'()') !! skip 1 line 
  do f=1,num_faces
    read(ALPHA_BETA_FILE,*) beta_f(f)
  enddo
close(ALPHA_BETA_FILE)
else
    alpha_s=1d0
    alpha_l=1d0
    alpha_f=1d0
    beta_f=1d0
endif
#ifdef PARALLEL
endif
call MPI_BCAST(alpha_s,num_sites,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(alpha_l,num_links,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(alpha_f,num_faces,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(beta_f,num_faces,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
#endif

end subroutine set_alpha_beta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! set Remez data
subroutine set_Remez_data
implicit none

integer i

#ifdef PARALLEL
if (MYRANK==0) then
#endif

open(Remez4_FILE, file=Remez_1ovminus4, status='OLD',action='READ')
read(Remez4_FILE,*) N_Remez4
read(Remez4_FILE,*) Remez_min4
read(Remez4_FILE,*) Remez_max4
allocate( Remez_alpha4(0:N_Remez4) )
allocate( Remez_beta4(1:N_Remez4) )
do i=0,N_Remez4 
  read(Remez4_FILE,*) Remez_alpha4(i)
enddo
do i=1,N_Remez4 
  read(Remez4_FILE,*) Remez_beta4(i)
enddo
close(Remez4_FILE)

#ifdef PARALLEL
endif
call MPI_BCAST(N_Remez4,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(Remez_min4,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(Remez_max4,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
if( MYRANK .ne. 0) then
  allocate( Remez_alpha4(0:N_Remez4) )
  allocate( Remez_beta4(1:N_Remez4) )
endif
call MPI_BCAST(Remez_alpha4,N_Remez4+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(Remez_beta4,N_Remez4,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
#endif

!! rescaling
Remez_alpha4(0)=Remez_alpha4(0)*Remez_factor4**(-0.25d0)
do i=1,N_Remez4
  Remez_alpha4(i)=Remez_alpha4(i)*Remez_factor4**(0.75d0)
  Remez_beta4(i)=Remez_beta4(i)*Remez_factor4
enddo



!!
#ifdef PARALLEL
if (MYRANK==0) then
#endif
open(Remez8_FILE, file=Remez_1ov8, status='OLD',action='READ')
read(Remez8_FILE,*) N_Remez8
read(Remez8_FILE,*) Remez_min8
read(Remez8_FILE,*) Remez_max8
allocate( Remez_alpha8(0:N_Remez8) )
allocate( Remez_beta8(1:N_Remez8) )
do i=0,N_Remez8 
  read(Remez8_FILE,*) Remez_alpha8(i)
enddo
do i=1,N_Remez8 
  read(Remez8_FILE,*) Remez_beta8(i)
enddo
close(Remez8_FILE)

#ifdef PARALLEL
endif
call MPI_BCAST(N_Remez8,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(Remez_min8,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(Remez_max8,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
if( MYRANK .ne. 0) then
  allocate( Remez_alpha8(0:N_Remez8) )
  allocate( Remez_beta8(1:N_Remez8) )
endif
call MPI_BCAST(Remez_alpha8,N_Remez8+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(Remez_beta8,N_Remez8,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
#endif

!! rescaling
Remez_alpha8(0)=Remez_alpha8(0)*Remez_factor8**(0.125d0)
do i=1,N_Remez8
  Remez_alpha8(i)=Remez_alpha8(i)*Remez_factor8**(1.125d0)
  Remez_beta8(i)=Remez_beta8(i)*Remez_factor8
enddo
end subroutine set_Remez_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! set sub simplicial complex
subroutine set_sub_SC
#ifdef PARALLEL
use parallel
#endif
implicit none
integer :: s,l,f,i


#ifdef PARALLEL
  ! sub_SC file is read only by MYRANK==0. 
  num_sub_SC=0
  !! set (rank, label) for a global site/link/face label
  allocate( local_site_of_global(1:num_sites) )
  allocate( local_link_of_global(1:num_links) )
  allocate( local_face_of_global(1:num_faces) )

  !! SUBSC_FILEを読み取って、local_{site/link/face}_of_global(:)を設定
  !! すなわち、各ノードが担当するsite/link/faceを設定
  !!   num_local_sites/links/faces
  !!   local_{site/link/face}_of_global(:)
  call set_map_from_global_to_local_data

  !call stop_for_test
  !! 各ノードが計算に必要なsiteの情報と、site情報の送り先、受信元を設定する
  !!   num_necessary_sites
  !!   global_site_of_local(s)
  !!   send_sites(s)
  !!   recv_sites(s)
  call set_local_sites
  !! 各ノードが計算に必要なlinkの情報と、link情報の送り先、受信元を設定する
  !!   num_necessary_links
  !!   global_link_of_local(l)
  !!   send_links(l)
  !!   recv_links(l)
  call set_local_links
  !! 各ノードが計算に必要なfaceの情報を設定する
  !!   num_necessary_faces
  !!   global_face_of_local(l)
  call set_local_faces

  !! localなlinkのorgとtipを設定
  call set_local_link_org_and_tip
  !! localなsiteから出入りするlocal linkのラベルを設定
  call set_local_link_fromto_site
  !! localなfaceに含まれるlinkのlocalラベルを設定
  call set_local_links_in_f
  !! localなlinkを共有するfaceのlocalラベルを設定
  call set_local_face_in_l
  !! localなfaceに含まれるsiteのlocalラベルを設定
  !! ただし、使うのは sites_in_f(f)%label_(1)だけなので、
  !! はじめから size=1 にしておく。
  call set_local_sites_in_f

#else
  num_local_sites=num_sites
  num_local_links=num_links
  num_local_faces=num_faces
  num_necessary_sites=num_sites
  num_necessary_links=num_links
  num_necessary_faces=num_faces
  allocate( global_site_of_local(1:num_sites) )
  allocate( global_link_of_local(1:num_links) )
  allocate( global_face_of_local(1:num_faces) )

  do s=1,num_sites
    global_site_of_local(s)=s
  enddo
  do l=1,num_links
    global_link_of_local(l)=l
  enddo
  do f=1,num_faces
    global_face_of_local(f)=f
  enddo
#endif

end subroutine set_sub_SC

#ifdef PARALLEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_map_from_global_to_local_data
use parallel
implicit none

integer :: tmp_num_local_sites,tmp_num_local_links,tmp_num_local_faces
integer,allocatable :: tmp_global_site_of_local(:)
integer :: s,l,f,rank,i,part

if (MYRANK == 0) then
  open(SUBSC_FILE, file=FsubSC, status='OLD',action='READ')
  read(SUBSC_FILE,*) num_sub_SC
  ! send num_sub_SC to all the other nodes
  call MPI_BCAST(num_sub_SC,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)

  if (num_sub_SC .ne. NPROCS) then
    write(*,*) "number of core is mismatch."
    close(SUBSC_FILE)
    call MPI_FINALIZE(IERR)
  stop
  endif
  do part=1,num_sub_SC
    ! sub_SCの構造は、
    !   ・Num_Sub_SC
    ! の後に、
    !   ・RANK
    !   ・number of sites in the node
    !   ・そのnodeが担当するsiteのリスト
    ! が順番に並んでいる。
    
    !! RANKとそれが担当するノードを読み込む
    read(SUBSC_FILE,*) rank
    read(SUBSC_FILE,*) tmp_num_local_sites
    allocate( tmp_global_site_of_local(1:tmp_num_local_sites) )
    read(SUBSC_FILE,*) (tmp_global_site_of_local(i),i=1,tmp_num_local_sites)
    !  
    !! set local_site_of_global(:) to a temporary variable
    do i=1,tmp_num_local_sites
      local_site_of_global(tmp_global_site_of_local(i))%rank_=rank
      local_site_of_global(tmp_global_site_of_local(i))%label_=i
    enddo
    tmp_num_local_links=0
    do l=1,num_links
      do s=1,tmp_num_local_sites
        if( tmp_global_site_of_local(s)==link_org(l) ) then
          tmp_num_local_links=tmp_num_local_links+1
          local_link_of_global(l)%rank_=rank
          local_link_of_global(l)%label_=tmp_num_local_links
        endif
      enddo
    enddo
    tmp_num_local_faces=0
    do f=1,num_faces
      do s=1,tmp_num_local_sites
        if( tmp_global_site_of_local(s) == sites_in_f(f)%label_(1)) then
          tmp_num_local_faces=tmp_num_local_faces+1
          i=i+1
          local_face_of_global(f)%rank_=rank
          local_face_of_global(f)%label_=tmp_num_local_faces
        endif
      enddo
    enddo
    !!
    if( rank == 0 ) then
      num_local_sites=tmp_num_local_sites
      num_local_links=tmp_num_local_links
      num_local_faces=tmp_num_local_faces
    else
      call MPI_SEND(tmp_num_local_sites,1,MPI_INTEGER,rank,1,MPI_COMM_WORLD,IERR)
      call MPI_SEND(tmp_num_local_links,1,MPI_INTEGER,rank,2,MPI_COMM_WORLD,IERR)
      call MPI_SEND(tmp_num_local_faces,1,MPI_INTEGER,rank,3,MPI_COMM_WORLD,IERR)
    endif
    deallocate( tmp_global_site_of_local ) 
  enddo
  close(SUBSC_FILE)
  else
    call MPI_BCAST(num_sub_SC,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
    call MPI_RECV(num_local_sites,1,MPI_INTEGER,0,1,MPI_COMM_WORLD,ISTATUS,IERR)
    call MPI_RECV(num_local_links,1,MPI_INTEGER,0,2,MPI_COMM_WORLD,ISTATUS,IERR)
    call MPI_RECV(num_local_faces,1,MPI_INTEGER,0,3,MPI_COMM_WORLD,ISTATUS,IERR)
endif
call MPI_BCAST(local_site_of_global,num_sites,MPI_LOCAL_LABEL,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(local_link_of_global,num_links,MPI_LOCAL_LABEL,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(local_face_of_global,num_faces,MPI_LOCAL_LABEL,0,MPI_COMM_WORLD,IERR)

end subroutine set_map_from_global_to_local_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! 計算に必要なsiteの情報を設定
subroutine set_local_sites
use parallel
implicit none

integer s,l,i,nsend,nrecv,info
integer :: tmp_global_site_of_local(1:num_sites)
!integer :: tmp_send_sites(1:num_sites)
!integer :: tmp_recv_sites(1:num_sites)
type(LOCAL_LABEL) :: tmp_send_sites(1:num_sites)
type(LOCAL_LABEL) :: tmp_recv_sites(1:num_sites)

!! 計算を担当するsiteを tmp_global_site_of_local の最初の方に設定
num_necessary_sites=0
tmp_global_site_of_local=0
do s=1,num_sites
  if( local_site_of_global(s)%rank_ == MYRANK ) then
    num_necessary_sites = num_necessary_sites + 1
    tmp_global_site_of_local(num_necessary_sites)=s
  endif
enddo
!! この段階で、num_necesarry_sites=num_local_sites
!! linkの始点と終点が別のnodeに属していたら、
!! tipの位置を、orgを担当するnodeに送る必要がある
!! 重複があり得るので、逐一チェックしながら進める
nsend=0
nrecv=0
do l=1,num_links
  if( local_site_of_global(link_org(l))%rank_ &
        /= local_site_of_global(link_tip(l))%rank_ ) then
    !!!!!!!!!!!
    if( local_site_of_global(link_tip(l))%rank_ == MYRANK ) then
      !! 重複チェック
      info=0
      do i=1,nsend
        if( tmp_send_sites(i)%rank_ == local_site_of_global(link_org(l))%rank_ &
          .and. &
            tmp_send_sites(i)%label_ == link_tip(l) ) then
          info=1
          exit
        endif
      enddo
      if (info==0) then 
        nsend=nsend+1
        tmp_send_sites(nsend)%rank_ = local_site_of_global(link_org(l))%rank_ 
        tmp_send_sites(nsend)%label_ = link_tip(l) 
      endif
    !!!!!!!!!!!
    elseif( local_site_of_global(link_org(l))%rank_ == MYRANK ) then
      !! 重複チェック
      info=0
      do i=1,nrecv
        if( tmp_recv_sites(i)%rank_ == local_site_of_global(link_tip(l))%rank_ &
          .and. &
            tmp_recv_sites(i)%label_ == link_tip(l) ) then
          info=1
          exit
        endif
      enddo
      if (info==0) then 
        nrecv=nrecv+1
        tmp_recv_sites(nrecv)%rank_ = local_site_of_global(link_tip(l))%rank_ 
        tmp_recv_sites(nrecv)%label_ = link_tip(l) 
        tmp_global_site_of_local(num_local_sites+nrecv)=link_tip(l)
      endif
    endif
  endif
enddo


!! 重複度を除いて配列を定義
num_send_sites=nsend
num_recv_sites=nrecv
allocate( send_sites(num_send_sites) )
allocate( recv_sites(num_recv_sites) )
num_necessary_sites = num_local_sites + num_recv_sites
allocate( global_site_of_local(1:num_necessary_sites) )

!! 代入
do s=1,num_send_sites
  send_sites(s)%rank_ = tmp_send_sites(s)%rank_
  send_sites(s)%label_ = local_site_of_global(tmp_send_sites(s)%label_)%label_
enddo
do s=1,num_recv_sites
  recv_sites(s)%rank_ = tmp_recv_sites(s)%rank_
  recv_sites(s)%label_ = num_local_sites+s
enddo
do s=1,num_necessary_sites
  global_site_of_local(s) = tmp_global_site_of_local(s)
enddo
end subroutine set_local_sites

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! 計算に必要なlinkの情報を設定
subroutine set_local_links
use parallel
implicit none

integer s,l,f,i,j,k,ll,ll_label
integer :: tmp_global_link_of_local(1:num_links)
type(LOCAL_LABEL) :: tmp_send_links(1:num_links)
type(LOCAL_LABEL) :: tmp_recv_links(1:num_links)
integer :: nsend,nrecv,info


do l=1,num_links
  tmp_global_link_of_local(l)=0
  tmp_send_links(l)%rank_=-1
  tmp_send_links(l)%label_=-1
  tmp_recv_links(l)%rank_=-1
  tmp_recv_links(l)%label_=-1
enddo

!! 計算を担当するlinkを tmp_global_link_of_local の最初の方に設定
num_necessary_links=0
tmp_global_link_of_local=0
do l=1,num_links
  if( local_link_of_global(l)%rank_ == MYRANK ) then
    num_necessary_links = num_necessary_links + 1
    tmp_global_link_of_local(num_necessary_links)=l
  endif
enddo
!! この段階で、num_necesarry_links=num_local_links
!! 計算に必要なのは、
!!   (1) 担当するsiteをtipとするlink
!!   (2) 担当するfaceを構成するlink
!!   (3) 担当するlinkを共有するfaceに含まれる全link
!! 重複があり得るので、逐一チェックしながら進める

!!   (1) 担当するsiteをtipとするlink
nsend=0
nrecv=0
do l=1,num_links
  if( local_site_of_global(link_tip(l))%rank_ /= &
      local_site_of_global(link_org(l))%rank_ ) then
    if( local_site_of_global(link_org(l))%rank_ == MYRANK ) then
      !! 重複チェック
      info=0
      do i=1,nsend
        if( tmp_send_links(i)%rank_ == local_site_of_global(link_tip(l))%rank_ &
          .and. &
            tmp_send_links(i)%label_ == l ) then
          info=1
          exit
        endif
      enddo
      if( info == 0 ) then
        nsend=nsend+1
        tmp_send_links(nsend)%rank_ = local_site_of_global(link_tip(l))%rank_
        tmp_send_links(nsend)%label_ = l
      endif
    !!!!!!!!!!!
    elseif( local_site_of_global(link_tip(l))%rank_ == MYRANK ) then
      !! 重複チェック
      info=0
      do i=1,nrecv
        if( tmp_recv_links(i)%rank_ == local_site_of_global(link_org(l))%rank_ &
          .and. &
            tmp_recv_links(i)%label_ == l ) then
          info=1
          exit
        endif
      enddo
      if (info==0) then 
        nrecv=nrecv+1
        tmp_recv_links(nrecv)%rank_ = local_site_of_global(link_org(l))%rank_ 
        tmp_recv_links(nrecv)%label_ = l
        tmp_global_link_of_local(num_local_links+nrecv)=l
      endif
    endif
  endif
enddo

!!   (2) 担当するfaceを構成するlink
do f=1,num_faces
  do j=1,links_in_f(f)%num_
    l=links_in_f(f)%link_labels_(j)
    if( local_face_of_global(f)%rank_ /= local_link_of_global(l)%rank_ ) then
      if( local_link_of_global(l)%rank_ == MYRANK ) then
        !! 重複チェック
        info=0
        do i=1,nsend
          if( tmp_send_links(i)%rank_ == local_face_of_global(f)%rank_ &
            .and. &
              tmp_send_links(i)%label_ == l ) then
            info=1
            exit
          endif
        enddo
        if( info == 0 ) then
          nsend=nsend+1
          tmp_send_links(nsend)%rank_ = local_face_of_global(f)%rank_ 
          tmp_send_links(nsend)%label_ = l
        endif
      !!!!!!!!!!!
      elseif( local_face_of_global(f)%rank_ == MYRANK ) then
        !! 重複チェック
        info=0
        do i=1,nrecv
          if( tmp_recv_links(i)%rank_ == local_link_of_global(l)%rank_ &
            .and. &
              tmp_recv_links(i)%label_ == l ) then
            info=1
            exit
          endif
        enddo
        if (info==0) then 
          nrecv=nrecv+1
          tmp_recv_links(nrecv)%rank_ = local_link_of_global(l)%rank_
          tmp_recv_links(nrecv)%label_ = l
          tmp_global_link_of_local(num_local_links+nrecv)=l
        endif
      endif
    endif
  enddo
enddo

!!   (3) 担当するlinkを共有するfaceに含まれる全link
do l=1,num_links
  do k=1,face_in_l(l)%num_
    f=face_in_l(l)%label_(k)
    do ll_label=1,links_in_f(f)%num_
      ll=links_in_f(f)%link_labels_(ll_label)
      if( local_link_of_global(l)%rank_ /= local_link_of_global(ll)%rank_ ) then 
        if( local_link_of_global(ll)%rank_ == MYRANK ) then
          !! 重複チェック
          info=0
          do i=1,nsend
            if( tmp_send_links(i)%rank_ == local_link_of_global(l)%rank_ &
                .and. &
                tmp_send_links(i)%label_ == ll ) then 
              info=1
              exit
            endif
          enddo
          if (info == 0 ) then
            nsend=nsend+1
            tmp_send_links(nsend)%rank_ = local_link_of_global(l)%rank_ 
            tmp_send_links(nsend)%label_ = ll
          endif
        !!!!!!!!!!!
        elseif( local_link_of_global(l)%rank_ == MYRANK ) then
          !! 重複チェック
          info=0
          do i=1,nrecv
            if( tmp_recv_links(i)%rank_ == local_link_of_global(ll)%rank_ &
                .and. &
                tmp_recv_links(i)%label_ == ll ) then 
              info=1
              exit
            endif
          enddo
          if (info == 0 ) then
            nrecv=nrecv+1
            tmp_recv_links(nrecv)%rank_ = local_link_of_global(ll)%rank_ 
            tmp_recv_links(nrecv)%label_ = ll
            tmp_global_link_of_local(num_local_links+nrecv)=ll
          endif
        endif
      endif
    enddo
  enddo
enddo

!! 重複度を除いて配列を定義
num_send_links=nsend
num_recv_links=nrecv
allocate( send_links(num_send_links) )
allocate( recv_links(num_recv_links) )
num_necessary_links = num_local_links + num_recv_links
allocate( global_link_of_local(1:num_necessary_links) )

!! 代入
do l=1,num_send_links
  send_links(l)%rank_ = tmp_send_links(l)%rank_
  send_links(l)%label_ = local_link_of_global(tmp_send_links(l)%label_)%label_
enddo
do l=1,num_recv_links
  recv_links(l)%rank_ = tmp_recv_links(l)%rank_
  recv_links(l)%label_ = num_local_links+l
enddo
do l=1,num_necessary_links
  global_link_of_local(l) = tmp_global_link_of_local(l)
enddo

end subroutine set_local_links


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! 計算に必要なfaceの情報を設定
subroutine set_local_faces
use parallel
implicit none

integer f,l,i,k
type(LOCAL_LABEL) :: tmp_send_faces(1:num_faces)
type(LOCAL_LABEL) :: tmp_recv_faces(1:num_faces)
integer :: nsend,nrecv,info
integer :: tmp_global_face_of_local(1:num_faces)

!! 計算を担当するfaceを tmp_global_face_of_local の最初の方に設定
num_necessary_faces=0
do f=1,num_faces
  if( local_face_of_global(f)%rank_ == MYRANK ) then
    num_necessary_faces = num_necessary_faces + 1
    tmp_global_face_of_local(num_necessary_faces)=f
  endif
enddo

nsend=0
nrecv=0
!! 担当するlinkを共有するfaceが全て必要
do l=1,num_links
  do k=1,face_in_l(l)%num_
    f=face_in_l(l)%label_(k)
    if( local_link_of_global(l)%rank_ &
        /= local_face_of_global(f)%rank_ ) then
      if( local_face_of_global(f)%rank_ == MYRANK ) then
        !重複チェック
        info=0
        do i=1,nsend
          if( tmp_send_faces(i)%rank_ == local_link_of_global(l)%rank_ &
              .and. &
              tmp_send_faces(i)%label_ == f ) then 
            info=1
            exit
          endif
        enddo
        if (info == 0 ) then
          nsend=nsend+1
          tmp_send_faces(nsend)%rank_ = local_link_of_global(l)%rank_ 
          tmp_send_faces(nsend)%label_ = f
        endif
      !!!!!!!!!!!
      elseif( local_link_of_global(l)%rank_ == MYRANK ) then
        !! 重複チェック
        info=0
        do i=1,nrecv
          if( tmp_recv_faces(i)%rank_ == local_face_of_global(f)%rank_ &
              .and. &
              tmp_recv_faces(i)%label_ == f ) then 
            info=1
            exit
          endif
        enddo
        if (info == 0 ) then
          nrecv=nrecv+1
          tmp_recv_faces(nrecv)%rank_ = local_face_of_global(f)%rank_ 
          tmp_recv_faces(nrecv)%label_ = f
          tmp_global_face_of_local(num_local_faces+nrecv)=f
        endif
      endif
    endif
  enddo
enddo


!! 重複度を除いて配列を定義
num_send_faces=nsend
num_recv_faces=nrecv
allocate( send_faces(num_send_faces) )
allocate( recv_faces(num_recv_faces) )
num_necessary_faces = num_local_faces + num_recv_faces
allocate( global_face_of_local(1:num_necessary_links) )

!! 代入
do i=1,num_send_faces
  send_faces(i)%rank_ = tmp_send_faces(i)%rank_
  send_faces(i)%label_ = local_face_of_global(tmp_send_faces(i)%label_)%label_
enddo
do i=1,num_recv_faces
  recv_faces(i)%rank_ = tmp_recv_faces(i)%rank_
  recv_faces(i)%label_ = num_local_faces+i
enddo
do f=1,num_necessary_faces
  global_face_of_local(f) = tmp_global_face_of_local(f)
enddo

end subroutine set_local_faces


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! それぞれのnode上でのlink_tipとlink_orgを設定
subroutine set_local_link_org_and_tip
use parallel
implicit none

integer info_o, info_t
integer l,s

allocate( local_link_org(1:num_local_links) )
allocate( local_link_tip(1:num_local_links) )

do l=1,num_local_links
  info_o=0
  info_t=0
  do s=1,num_necessary_sites
    if( global_site_of_local(s) == link_org(global_link_of_local(l)) ) then
      local_link_org(l) = s
      info_o=1
    elseif( global_site_of_local(s) == link_tip(global_link_of_local(l)) ) then
      local_link_tip(l) = s
      info_t=1
    endif
    if( info_o * info_t == 1 ) exit
  enddo
  if( info_o * info_t == 0 ) then 
    write(*,*) "something is wrong in local link", l, "in RANK",MYRANK
    call stop_for_test
  endif
enddo


end subroutine set_local_link_org_and_tip



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! それぞれのnode上で、siteに出入りするlinkのラベルを設定
subroutine set_local_link_fromto_site
use parallel
implicit none

integer s,ss,l,num,info,i

allocate(local_linktip_from_s(1:num_local_sites))
allocate(local_linkorg_to_s(1:num_local_sites))
do s=1,num_local_sites
  !! local_linktip
  num=linktip_from_s(global_site_of_local(s))%num_
  local_linktip_from_s(s)%num_=num
  allocate(local_linktip_from_s(s)%labels_(1:num) )
  do i=1,num
    info=0
    do l=1,num_necessary_links
      if( linktip_from_s(global_site_of_local(s))%labels_(i) &
          == global_link_of_local(l) ) then
        local_linktip_from_s(s)%labels_(i)=l
        info=1
        exit
      endif
    enddo
    if( info==0 ) then 
      write(*,*) MYRANK,"Something happpend in set_local_linktip_from_s"
      call stop_for_test
    endif
  enddo
  !! local_linkorg
  num=linkorg_to_s(global_site_of_local(s))%num_
  local_linkorg_to_s(s)%num_= num
  allocate(local_linkorg_to_s(s)%labels_(1:num) )
  do i=1,num
    do l=1,num_necessary_links
      info=0
      if( linkorg_to_s(global_site_of_local(s))%labels_(i) &
          == global_link_of_local(l) ) then
        local_linkorg_to_s(s)%labels_(i)=l
        info=1
        exit
      endif
    enddo
    if( info==0 ) then 
      write(*,*) MYRANK,"Something happpend in set_local_linkorg_to_s 1"
      call stop_for_test
    endif
  enddo
enddo

end subroutine set_local_link_fromto_site

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! それぞれのnode上で、faceを構成するlinkのラベルを設定
subroutine  set_local_links_in_f
use parallel
implicit none

integer f,l,ll, num, g_link, dir

allocate( local_links_in_f(1:num_local_faces) )
do f=1,num_local_faces
  num = links_in_f(global_face_of_local(f))%num_
  local_links_in_f(f)%num_ = num
  allocate( local_links_in_f(f)%link_labels_(1:num) )
  allocate( local_links_in_f(f)%link_dirs_(1:num))

  do l=1,num
    g_link = links_in_f(global_face_of_local(f))%link_labels_(l)
    dir = links_in_f(global_face_of_local(f))%link_dirs_(l)
    do ll=1,num_necessary_links
      if( global_link_of_local(ll) == g_link ) then
        local_links_in_f(f)%link_labels_(l) = ll
        local_links_in_f(f)%link_dirs_(l) = dir
      endif
    enddo
  enddo
enddo

end subroutine  set_local_links_in_f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! それぞれのnode上で、linkに共有されるfaceのラベルを設定
subroutine  set_local_face_in_l
use parallel
implicit none

integer l,num ,i,f,g_face,info

allocate( local_face_in_l(1:num_local_links) )
do l=1,num_local_links
  num = face_in_l(global_link_of_local(l))%num_
  local_face_in_l(l)%num_ = num
  allocate( local_face_in_l(l)%label_(1:num) )

  do i=1,num
    info=0
    g_face = face_in_l(global_link_of_local(l))%label_(i)
    do f=1,num_necessary_faces
      if( global_face_of_local(f) == g_face ) then
        local_face_in_l(l)%label_(i) = f
        info=1
      endif
    enddo
    if( info == 0 ) then 
      write(*,*) "something happened in setting local_face_in_l."
      write(*,*) MYRANK, f,i
      call stop_for_test
    endif
  enddo
enddo


end subroutine  set_local_face_in_l


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! それぞれのnode上で担当するfaceのlocalな代表点を設定
subroutine set_local_sites_in_f
use parallel
implicit none

integer f,s,global_s,info

allocate( local_sites_in_f(1:num_local_faces) )
do f=1,num_local_faces
  local_sites_in_f(f)%num_=1 
  allocate( local_sites_in_f(f)%label_(1:1) )

  global_s=sites_in_f(global_face_of_local(f))%label_(1)
  do s=1,num_necessary_sites
    
    if( global_s == global_site_of_local(s) ) then
      local_sites_in_f(f)%label_(1) = s
      info=1
      exit
    endif
  enddo
  if( info == 0 ) then
    write(*,*) "something happened in set_local_sites_in_f"
    call stop_for_test
  endif
enddo

end subroutine set_local_sites_in_f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! それぞれのnode上で必要な alpha, beta を設定する
subroutine set_local_alpha_beta
use parallel
implicit none

integer s,l,f

allocate(local_alpha_s(1:num_necessary_sites) )
allocate(local_alpha_l(1:num_necessary_links) )
allocate(local_alpha_f(1:num_necessary_faces) )
allocate(local_beta_f(1:num_necessary_faces) )
do s=1,num_necessary_sites
  local_alpha_s(s) = alpha_s(global_site_of_local(s) )
enddo

do l=1,num_necessary_links
  local_alpha_l(l) = alpha_l(global_link_of_local(l) )
enddo

do f=1,num_necessary_faces
  local_alpha_f(f) = alpha_f(global_face_of_local(f) )
  local_beta_f(f) = beta_f(global_face_of_local(f) )
enddo


end subroutine set_local_alpha_beta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! globalに定義していた諸々の数をlocalに置き換える
subroutine switch_globalnum_to_localnum
use parallel
implicit none

integer s,l,f,i,num

global_num_sites = num_sites
global_num_links = num_links
global_num_faces = num_faces

num_sites = num_local_sites
num_links = num_local_links
num_faces = num_local_faces

deallocate( link_org, link_tip )
allocate(link_org(1:num_local_links))
allocate(link_tip(1:num_local_links))
link_org = local_link_org
link_tip = local_link_tip

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
deallocate( linktip_from_s, linkorg_to_s )
allocate(linktip_from_s(1:num_local_sites)) 
allocate(linkorg_to_s(1:num_local_sites))
do s=1,num_sites
  num = local_linktip_from_s(s)%num_
  linktip_from_s(s)%num_ = num
  allocate(linktip_from_s(s)%labels_(1:num) )
  allocate(linktip_from_s(s)%sites_(1:num) )
  do i=1,num
    linktip_from_s(s)%labels_(i) = local_linktip_from_s(s)%labels_(i)
  enddo
  !!!!!
  num = local_linkorg_to_s(s)%num_
  allocate(linkorg_to_s(s)%labels_(1:num) )
  allocate(linkorg_to_s(s)%sites_(1:num) )
  linkorg_to_s(s)%num_ = num
  do i=1,num
    linkorg_to_s(s)%labels_(i) = local_linkorg_to_s(s)%labels_(i)
  enddo
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
deallocate( links_in_f )
allocate( links_in_f(1:num_faces) )
do f=1,num_faces
  num=local_links_in_f(f)%num_
  links_in_f(f)%num_ = num
  allocate( links_in_f(f)%link_labels_(1:num) )
  allocate( links_in_f(f)%link_dirs_(1:num) )
  links_in_f(f)%link_labels_=local_links_in_f(f)%link_labels_
  links_in_f(f)%link_dirs_=local_links_in_f(f)%link_dirs_
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
deallocate( face_in_l )
allocate( face_in_l(1:num_links) )
do l=1,num_links
  num=local_face_in_l(l)%num_
  face_in_l(l)%num_ = num
  allocate( face_in_l(l)%label_(1:num) )
  face_in_l(l)%label_=local_face_in_l(l)%label_
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
deallocate( sites_in_f )
allocate( sites_in_f(1:num_faces) )
do f=1,num_faces
  sites_in_f(f)%num_=1
  allocate( sites_in_f(f)%label_(1:1) )
  sites_in_f(f)%label_ = local_sites_in_f(f)%label_
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
deallocate( alpha_s )
allocate( alpha_s(1:num_necessary_sites) )
alpha_s = local_alpha_s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
deallocate( alpha_l )
allocate( alpha_l(1:num_necessary_links) )
alpha_l = local_alpha_l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
deallocate( alpha_f )
allocate( alpha_f(1:num_necessary_faces) )
alpha_f = local_alpha_f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
deallocate( beta_f )
allocate( beta_f(1:num_necessary_faces) )
beta_f = local_beta_f

end subroutine switch_globalnum_to_localnum

#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! テスト用に止めるルーチン
subroutine stop_for_test
#ifdef PARALLEL
use parallel
#endif
implicit none

#ifdef PARALLEL
  call MPI_FINALIZE(IERR)
#endif
  stop
end subroutine stop_for_test



end module global_parameters
