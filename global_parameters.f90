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

type SITE_DIST
  integer :: rank_
  integer :: num_sites_
  integer, allocatable :: site_list_(:)
end type SITE_DIST

! parameter for test hamiltonian
integer, parameter :: pb_mass=0
integer, parameter :: pb_site=0
integer, parameter :: pb_link=0
integer, parameter :: pb_face=0 

integer, parameter :: pf=0 ! set 1 if you do NOT wanna include fermion
integer, parameter :: p1=0
integer, parameter :: p2=0
integer, parameter :: p3=0
integer, parameter :: p4=0
integer, parameter :: p5=0
integer, parameter :: p_mass=1

integer, parameter :: Dirac_write=0

! set 1 if you wanna check distribution of the simplicial complex
integer, parameter :: check_sub_sc=0

! bosonic and fermionic forces
double precision, save :: bosonic_force_Phi
double precision, save :: fermionic_force_Phi
double precision, save :: bosonic_force_A
double precision, save :: fermionic_force_A
integer, save :: b_phi_count, f_phi_count, b_A_count, f_A_count

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
integer, save :: NZF      ! number of nonzero elements of f(a,b,c)
! index of nonzero f(a,b,c)
integer, allocatable, save :: NZF_index(:,:)  ! NZF_index(1:3,1:NZF)
! value of nonzero f(a,b,c)
double precision, allocatable, save :: NZF_value(:) ! NZF_value(1:NZF)
! number of indices (b,c) with nonzero f(a,b,c) for given a
type(INDEXDAT), allocatable, save :: NZF_a(:) 
! number of nonzero elementes of f(a,b,e)*f(c,d,e)
integer, save :: NZFF

! index of nonzero f(a,b,e)*f(c,d,e)
integer, allocatable, save :: NZFF_index(:,:)  ! NZF_index(1:4,1:NZFF)
! value of nonzero f(a,b,e)*f(c,d,e)
double precision, allocatable, save :: NZFF_value(:) ! NZF_value(1:NZFF)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! data of the simplicial complex

!! global simplicial complex
integer, save :: global_num_sites ! number of sites
integer, save :: global_num_links ! number of links
integer, save :: global_num_faces ! number of faces
double precision, allocatable, save :: global_alpha_s(:)
double precision, allocatable, save :: global_alpha_l(:)
double precision, allocatable, save :: global_alpha_f(:)
double precision, allocatable, save :: global_beta_f(:)

integer, allocatable, save :: global_link_org(:) ! l番目のリンクのorigin
integer, allocatable, save :: global_link_tip(:) ! l番目のリンクのtop
type(SITE_LINKDAT), allocatable, save :: global_linktip_from_s(:) ! sから出るリンクのtip
type(SITE_LINKDAT), allocatable, save :: global_linkorg_to_s(:) ! sに入るリンクのorigin
type(LINKDAT), allocatable, save :: global_links_in_f(:) ! fに含まれるリンクのデータ
type(FACEDAT), allocatable, save :: global_face_in_l(:) ! lに含まれるface
type(SITEDAT), allocatable, save :: global_sites_in_f(:) ! face f に含まれるサイト



integer, save :: num_sites ! number of sites
integer, save :: num_links ! number of links
integer, save :: num_faces ! number of faces
type(SmpCom), save :: SC
integer, allocatable, save :: link_org(:) ! l番目のリンクのorigin
integer, allocatable, save :: link_tip(:) ! l番目のリンクのtop
!! <s,・>
type(SITE_LINKDAT), allocatable, save :: linktip_from_s(:) ! sから出るリンクのtip
!! <・,t>
type(SITE_LINKDAT), allocatable, save :: linkorg_to_s(:) ! sに入るリンクのorigin
!!
type(LINKDAT), allocatable, save :: links_in_f(:) ! fに含まれるリンクのデータ
!! 
type(FACEDAT), allocatable, save :: face_in_l(:) ! lに含まれるface
!!
type(SITEDAT), allocatable, save :: sites_in_f(:) ! face f に含まれるサイト

double precision, allocatable, save :: alpha_s(:)
double precision, allocatable, save :: alpha_l(:)
double precision, allocatable, save :: alpha_f(:)
double precision, allocatable, save :: beta_f(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! data of the LOCAL simplicial complex
integer :: num_sub_SC
integer, save :: num_local_sites ! number of sites for each sub_SC
integer, save :: num_local_links ! number of links for each sub_SC
integer, save :: num_local_faces ! number of faces for each sub_SC


integer, save :: num_necessary_sites ! number of sites needed for each rank
integer, save :: num_necessary_links ! number of links needed for each rank
integer, save :: num_necessary_faces ! number of faces needed for each rank


!integer, allocatable, save :: local_link_org(:) ! l番目のリンクのorigin
!integer, allocatable, save :: local_link_tip(:) ! l番目のリンクのtop
!!! <s,・>
!type(SITE_LINKDAT), allocatable :: local_linktip_from_s(:) ! sから出るリンクとそのtip
!!! <・,t>
!type(SITE_LINKDAT), allocatable :: local_linkorg_to_s(:) ! sに入るリンクとそのorigin
!!!
!type(LINKDAT), allocatable :: local_links_in_f(:) ! fに含まれるリンクのデータ
!!! 
!type(FACEDAT), allocatable :: local_face_in_l(:) ! lに含まれるface
!!!
!type(SITEDAT), allocatable :: local_sites_in_f(:) ! face f に含まれるサイト

!double precision, allocatable :: local_alpha_s(:)
!double precision, allocatable :: local_alpha_l(:)
!double precision, allocatable :: local_alpha_f(:)
!double precision, allocatable :: local_beta_f(:)


! (rank, local label) corresponding to a global label
#ifdef PARALLEL
type LOCAL_LABEL
  integer :: rank_ ! nodes
  integer :: label_ ! label
end type LOCAL_LABEL

integer, allocatable, save :: global_site_of_local(:) ! map from local site to global site
integer, allocatable, save:: global_link_of_local(:) ! map from local link to global link
integer, allocatable, save:: global_face_of_local(:) ! map from local face to global face

type(LOCAL_LABEL), allocatable, save :: local_site_of_global(:) ! map from global site to local site
type(LOCAL_LABEL), allocatable, save :: local_link_of_global(:) ! map from global link to local link
type(LOCAL_LABEL), allocatable, save :: local_face_of_global(:) ! map from global face to local face


integer, save :: num_send_sites, num_recv_sites
integer, save :: num_send_links, num_recv_links
integer, save :: num_send_faces, num_recv_faces

type(LOCAL_LABEL), allocatable, save :: send_sites(:) 
type(LOCAL_LABEL), allocatable, save :: send_links(:) 
type(LOCAL_LABEL), allocatable, save :: send_faces(:) 

type(LOCAL_LABEL), allocatable, save :: recv_sites(:) 
type(LOCAL_LABEL), allocatable, save :: recv_links(:) 
type(LOCAL_LABEL), allocatable, save :: recv_faces(:) 
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! controle paramaeters of simulation 
integer, save :: test_mode !! 0: Simulation mode, 1: Test mode, 
integer, save :: force_measurement !!  0:measure forces
integer, save :: fix_seed !! 0: use the value at the end of the previous simulation
                    !! 1: fix seed to the value in the parameter file
                    !! 2: determine by the system time
integer, save :: new_config !! 0: Old config, 1: New config
!integer :: read_alpha !! 0: alpha,beta=1 1: read ALPHA_BETA_FILE
integer, save :: save_med_step !! steps to write out configurations
integer, save :: obs_step !! steps to compute observables
integer, save :: reset_ite !! 0:iteration number continue, 1:iteration number reset 
integer, save :: save_config_step !! save config by this step. 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! paramaeters in simulation 
integer, save :: NMAT     ! rank of the gauge group
integer, save :: dimG     ! dimension of su(NMAT)=NMAT**2-1
integer, save :: m_omega ! the integer m to avoid degenerated vacua
integer, save :: num_ite ! number of iteration
integer, save :: total_ite ! present iterations
integer, save :: job_number ! 0:use previous job, (say)123: use config 122
double precision, save :: phys_mass_square_phi ! mass^2 of Phi in mass-dim 2
double precision, save :: mass_square_phi ! dimensionless mass^2 of Phi
double precision, save :: mass_f   ! mass of fermions
double precision, save :: LatticeSpacing ! lattice spacing
complex(kind(0d0)), save :: e_max


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Molecular evolutin
!integer :: Ntau_base ! Ntau
!integer :: Ntau 
double precision, save :: Tau ! trajectory length
double precision, save :: Dtau,Dtau_phi,Dtau_A   ! Dtau
integer, save :: Ntau

! For multi-step
integer, save :: FB_ratio ! Phiのfermion forceの計算回数をboson forceの回数の何倍にするか
integer, save :: Nfermion, Nboson
double precision, save :: Dtau_boson ! Dtau for boson
double precision, save :: Dtau_fermion ! Dtau for fermion

double precision, save :: overall_factor ! N/2g^2N=N/2a^2
double precision, save :: maximal_dist ! the largest allowed value of |U_f|

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! parameters for Remez algorithm
character(128), parameter :: Remez_1ovminus4="remez_Q24_1ov-4.dat"
character(128), parameter :: Remez_1ov8="remez_Q24_1ov8.dat" 
integer, save :: N_Remez4 ! number of the rational function in rational approximation 
integer, save :: N_Remez8 ! number of the rational function in rational approximation 
double precision, save :: Remez_min4
double precision, save :: Remez_max4
double precision, save :: Remez_min8
double precision, save :: Remez_max8
double precision, allocatable, save :: Remez_alpha4(:)
double precision, allocatable, save :: Remez_beta4(:)
double precision, allocatable, save :: Remez_alpha8(:)
double precision, allocatable, save :: Remez_beta8(:)
double precision, save :: Remez_factor4
double precision, save :: Remez_factor8
double precision, save :: epsilon ! range to stop CG solver
integer, save :: CG_max ! maximum number of CG iteration

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Files
character(128), save :: PAR_FILE_NAME="parameters.dat" ! parameter file
character(128), save :: INPUT_FILE_NAME="inputfile" ! input file
character(128), save :: SC_FILE_NAME ! structure of simplicial complex
!character(128), save :: ALPHA_BETA ! parameterfile of alpha_s, alpha_l, alpha_f, beta_f
character(128), save :: Fconfigin ! configuration input 
character(128), save :: Fconfigout ! configuration output in the end
character(128), save :: Foutput ! output file
character(128), save :: Fmedconf ! configurations during the simulation 
!character(128), save :: FsubSC ! file to describe sub simplicial complex
character(128), save :: Flatestconf="CONFIG/latestconf"


integer, parameter :: PAR_FILE = 10 ! parameter file
integer, parameter :: SC_FILE = 11 ! data of simplicial complex
!integer, parameter :: ALPHA_BETA_FILE = 12 ! factor of the action
integer, parameter :: Remez4_FILE = 13 ! remez parameter for -1/4
integer, parameter :: Remez8_FILE = 14 ! remez parameter for 1/8
integer, parameter :: IN_CONF_FILE = 20 ! input configuration file
integer, parameter :: OUTPUT_FILE = 21 ! output file
integer, parameter :: OUT_CONF_FILE = 22 ! output configuration file
integer, parameter :: MED_CONF_FILE = 23 ! configurations during the simulation 
!integer, parameter :: SUBSC_FILE = 24 ! file to describe sub simplicial complex
integer, parameter :: INPUT_FILE = 25 ! input file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Control Algorighms
!integer, parameter :: USE_CG = 0 ! 0: do not use CG, 1: use CG

contains 
#include "set_theory_parameters.f90"
#include "set_NZF.f90"
#include "set_global_simplicial_complex.f90"
#include "set_simulation_parameters.f90"
#include "set_local_data.f90"
#include "set_Remez_data.f90"

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




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine set_sc_alpha_beta_parallel
!use simplicial_complex
!implicit none
!integer l,origin,tip
!integer s,f
!integer :: fsize
!integer, allocatable :: fsites(:), faces_l(:), sites_f(:)
!integer :: FaceSize
!integer, allocatable :: sites(:)
!integer k,j
!character(128) tmp
!double precision :: alpha
!! open SC_FILE 
!if(MYRANK==0) then
!  open(SC_FILE, file=SC_FILE_NAME, status='old',action='READ')
!  read(SC_FILE,'()') ! skip 1 line
!  read(SC_FILE,*) num_sites
!  read(SC_FILE,*) num_links
!  read(SC_FILE,*) num_faces
!  read(SC_FILE,*) LatticeSpacing
!  read(SC_FILE,'()') ! skip 1 line
!  !write(*,*) num_sites,num_links,num_faces,LatticeSpacing
!endif
!  call MPI_BCAST(num_sites,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
!  call MPI_BCAST(num_links,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
!  call MPI_BCAST(num_faces,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
!  call MPI_BCAST(LatticeSpacing,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!  overall_factor=dble(NMAT)/(2d0*LatticeSpacing*LatticeSpacing)
!  mass_square_phi = phys_mass_square_phi * (LatticeSpacing*latticeSpacing)
!
!! initialize the simplicial complex
!call init_smpcom(SC,num_sites,num_links,num_faces)
!!ひとまずglobalなalpha,betaを設定してしまう。
!allocate( alpha_s(1:num_sites) )
!allocate( alpha_l(1:num_links) )
!allocate( alpha_f(1:num_faces) )
!allocate( beta_f(1:num_faces) )
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Set alpha_s
!if( MYRANK==0 ) then
!  do k=1,num_sites
!    read(SC_FILE,*) s,alpha
!    alpha_s(s)=alpha
!    !write(*,*) s, alpha_s(s)
!  enddo
!  read(SC_FILE,'()') ! skip 1 line
!endif
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Set alpha_l
!do k=1,num_links
!  if( MYRANK==0 ) then 
!    read(SC_FILE,*) l,origin,tip,alpha
!    alpha_l(l)=alpha
!  endif
!  call MPI_BCAST(l,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
!  call MPI_BCAST(origin,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
!  call MPI_BCAST(tip,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
!  !write(*,*) MYRANK,l,origin,tip
!  !write(*,*) l,origin,tip,alpha_l(l)
!  call put_link_sc(SC,l,origin,tip)
!enddo
!if( MYRANK==0 ) then 
!  read(SC_FILE,'()') ! skip 1 line
!endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Set alpha_f
!do k=1,num_faces
!  if(MYRANK==0) then
!    read(SC_FILE,"(I6,X,I6,X)",advance="no") f,fsize
!    !write(*,*) f,fsize
!  endif
!  call MPI_BCAST(f,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
!  call MPI_BCAST(fsize,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
!
!  allocate( fsites(1:fsize ) )
!
!! set sites constructing i'th face
!  if (MYRANK==0) then
!    do j=1,fsize
!      read(SC_FILE,'(I6,X)',advance='no') fsites(j) 
!    enddo
!    !write(*,*) fsites
!    read(SC_FILE,*) alpha_f(f),beta_f(f)
!  endif
!  call MPI_BCAST(fsites,fsize,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
!!write(*,*) MYRANK, f,fsites
!  call put_face_sc(SC,f,fsites)
!
!  deallocate( fsites )
!enddo
!if (MYRANK==0) then
!  close(SC_FILE)
!endif
!
!! set links
!allocate(link_org(1:num_links))
!allocate(link_tip(1:num_links))
!do l=1,num_links
!  call get_link_sc(sc,l,link_org(l),link_tip(l))
!enddo
!
!! tips of links from s
!allocate(linktip_from_s(1:num_sites))
!do s=1,num_sites
!call get_links_from_s_sc(sc,s,&
!  linktip_from_s(s)%labels_,&
!  linktip_from_s(s)%sites_)
!linktip_from_s(s)%num_=size(linktip_from_s(s)%labels_)
!enddo
!
!! origins of links to s
!allocate(linkorg_to_s(1:num_sites))
!do s=1,num_sites
!call get_links_to_s_sc(sc,s,&
!  linkorg_to_s(s)%labels_,&
!  linkorg_to_s(s)%sites_)
!linkorg_to_s(s)%num_=size(linkorg_to_s(s)%labels_)
!enddo
!
!! links included in a face f
!allocate(links_in_f(1:num_faces))
!do f=1,num_faces
!  call get_links_in_face_sc(sc,f,&
!    links_in_f(f)%num_,&
!    sites,&
!    links_in_f(f)%link_labels_,&
!    links_in_f(f)%link_dirs_)
!enddo
!
!! faces included in a link l
!allocate(face_in_l(1:num_links))
!do l=1,num_links
!  call get_faces_in_link_sc(sc,l,faces_l)
!  face_in_l(l)%num_=size(faces_l)
!  allocate(face_in_l(l)%label_(1:size(faces_l)))
!  face_in_l(l)%label_=faces_l
!enddo
!  
!! sites included in a face f
!allocate(sites_in_f(1:num_faces))
!do f=1,num_faces
!  call get_face_components_sc(sc,f,sites_f)
!  sites_in_f(f)%num_=size(sites_f)
!  allocate(sites_in_f(f)%label_(1:size(sites_f)))
!  sites_in_f(f)%label_=sites_f
!enddo
!
!end subroutine set_sc_alpha_beta_parallel
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine set_sc_alpha_beta
!use simplicial_complex
!implicit none
!integer l,origin,tip
!integer s,f
!integer :: fsize
!integer, allocatable :: fsites(:), faces_l(:), sites_f(:)
!integer :: FaceSize
!integer, allocatable :: sites(:)
!integer k,j
!character(128) tmp
!double precision :: alpha
!! open SC_FILE 
!open(SC_FILE, file=SC_FILE_NAME, status='old',action='READ')
!read(SC_FILE,'()') ! skip 1 line
!read(SC_FILE,*) num_sites
!read(SC_FILE,*) num_links
!read(SC_FILE,*) num_faces
!read(SC_FILE,*) LatticeSpacing
!read(SC_FILE,'()') ! skip 1 line
!
!! initialize the simplicial complex
!call init_smpcom(SC,num_sites,num_links,num_faces)
!!ひとまずglobalなalpha,betaを設定してしまう。
!allocate( alpha_s(1:num_sites) )
!allocate( alpha_l(1:num_links) )
!allocate( alpha_f(1:num_faces) )
!allocate( beta_f(1:num_faces) )
!
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Set alpha_s
!do k=1,num_sites
!  read(SC_FILE,*) s,alpha
!  alpha_s(s)=alpha
!enddo
!read(SC_FILE,'()') ! skip 1 line
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Set alpha_l
!do k=1,num_links
!  read(SC_FILE,*) l,origin,tip,alpha
!  alpha_l(l)=alpha
!  !write(*,*) MYRANK,l,origin,tip
!  call put_link_sc(SC,l,origin,tip)
!enddo
!read(SC_FILE,'()') ! skip 1 line
!
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Set alpha_f
!do k=1,num_faces
!  read(SC_FILE,"(I6,X,I6)",advance="no") f,fsize
!
!  allocate( fsites(1:fsize ) )
!
!! set sites constructing i'th face
!  do j=1,fsize
!    read(SC_FILE,'(I6)',advance='no') fsites(j) 
!  enddo
!  read(SC_FILE,*) alpha_f(f),beta_f(f)
!!write(*,*) MYRANK, f,fsites
!  call put_face_sc(SC,f,fsites)
!
!  deallocate( fsites )
!enddo
!close(SC_FILE)
!
!! set links
!allocate(link_org(1:num_links))
!allocate(link_tip(1:num_links))
!do l=1,num_links
!  call get_link_sc(sc,l,link_org(l),link_tip(l))
!enddo
!
!! tips of links from s
!allocate(linktip_from_s(1:num_sites))
!do s=1,num_sites
!call get_links_from_s_sc(sc,s,&
!  linktip_from_s(s)%labels_,&
!  linktip_from_s(s)%sites_)
!linktip_from_s(s)%num_=size(linktip_from_s(s)%labels_)
!enddo
!
!! origins of links to s
!allocate(linkorg_to_s(1:num_sites))
!do s=1,num_sites
!call get_links_to_s_sc(sc,s,&
!  linkorg_to_s(s)%labels_,&
!  linkorg_to_s(s)%sites_)
!linkorg_to_s(s)%num_=size(linkorg_to_s(s)%labels_)
!enddo
!
!! links included in a face f
!allocate(links_in_f(1:num_faces))
!do f=1,num_faces
!  call get_links_in_face_sc(sc,f,&
!    links_in_f(f)%num_,&
!    sites,&
!    links_in_f(f)%link_labels_,&
!    links_in_f(f)%link_dirs_)
!enddo
!
!! faces included in a link l
!allocate(face_in_l(1:num_links))
!do l=1,num_links
!  call get_faces_in_link_sc(sc,l,faces_l)
!  face_in_l(l)%num_=size(faces_l)
!  allocate(face_in_l(l)%label_(1:size(faces_l)))
!  face_in_l(l)%label_=faces_l
!enddo
!  
!! sites included in a face f
!allocate(sites_in_f(1:num_faces))
!do f=1,num_faces
!  call get_face_components_sc(sc,f,sites_f)
!  sites_in_f(f)%num_=size(sites_f)
!  allocate(sites_in_f(f)%label_(1:size(sites_f)))
!  sites_in_f(f)%label_=sites_f
!enddo
!
!end subroutine set_sc_alpha_beta
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine set_sc
!use simplicial_complex
!implicit none
!integer l,origin,tip
!integer s,f
!integer :: fsize
!integer, allocatable :: fsites(:), faces_l(:), sites_f(:)
!integer :: FaceSize
!integer, allocatable :: sites(:)
!integer k
!character(128) tmp
!! open SC_FILE 
!#ifdef PARALLEL
!if(MYRANK==0) then
!  open(SC_FILE, file=SC_FILE_NAME, status='old',action='READ')
!  read(SC_FILE,*) num_sites
!  read(SC_FILE,*) num_links
!  read(SC_FILE,*) num_faces
!endif
!  call MPI_BCAST(num_sites,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
!  call MPI_BCAST(num_links,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
!  call MPI_BCAST(num_faces,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
!#else
!open(SC_FILE, file=SC_FILE_NAME, status='old',action='READ')
!  read(SC_FILE,*) num_sites
!  read(SC_FILE,*) num_links
!  read(SC_FILE,*) num_faces
!#endif
!
!
!! initialize the simplicial complex
!call init_smpcom(SC,num_sites,num_links,num_faces)
!
!do l=1,num_links
!#ifdef PARALLEL
!if(MYRANK==0) then
!  read(SC_FILE,*) origin,tip
!endif
!  call MPI_BCAST(origin,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
!  call MPI_BCAST(tip,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
!  !write(*,*) MYRANK,l,origin,tip
!#else
!  read(SC_FILE,*) origin,tip
!#endif
!  call put_link_sc(SC,l,origin,tip)
!enddo
!
!do f=1,num_faces
!#ifdef PARALLEL
!if(MYRANK==0) then
!  read(SC_FILE,"(I1)",advance="no") fsize
!endif
!call MPI_BCAST(fsize,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
!#else
!  read(SC_FILE,"(I1)",advance="no") fsize
!#endif
!
!  allocate( fsites(1:fsize ) )
!
!  ! set sites constructing i'th face
!#ifdef PARALLEL
!  if (MYRANK==0) then
!    read(SC_FILE,*) (fsites(k),k=1,fsize)
!  endif
!  call MPI_BCAST(fsites,fsize,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
!  !write(*,*) MYRANK, f,fsites
!#else
!  read(SC_FILE,*) (fsites(k),k=1,fsize)
!#endif
!  call put_face_sc(SC,f,fsites)
!
!  deallocate( fsites )
!enddo
!#ifdef PARALLEL
!if (MYRANK==0) then
!  close(SC_FILE)
!endif
!#else 
!  close(SC_FILE)
!#endif
!
!
!
!!write(*,*) "==="
!!write(*,*) sc%FaceSet_(1)%FaceSize_
!!write(*,*) "==="
!
!! set links
!allocate(link_org(1:num_links))
!allocate(link_tip(1:num_links))
!do l=1,num_links
!  call get_link_sc(sc,l,link_org(l),link_tip(l))
!enddo
!
!! tips of links from s
!allocate(linktip_from_s(1:num_sites))
!do s=1,num_sites
!call get_links_from_s_sc(sc,s,&
!  linktip_from_s(s)%labels_,&
!  linktip_from_s(s)%sites_)
!linktip_from_s(s)%num_=size(linktip_from_s(s)%labels_)
!enddo
!
!! origins of links to s
!allocate(linkorg_to_s(1:num_sites))
!do s=1,num_sites
!call get_links_to_s_sc(sc,s,&
!  linkorg_to_s(s)%labels_,&
!  linkorg_to_s(s)%sites_)
!linkorg_to_s(s)%num_=size(linkorg_to_s(s)%labels_)
!enddo
!
!! links included in a face f
!allocate(links_in_f(1:num_faces))
!do f=1,num_faces
!  call get_links_in_face_sc(sc,f,&
!    links_in_f(f)%num_,&
!    sites,&
!    links_in_f(f)%link_labels_,&
!    links_in_f(f)%link_dirs_)
!enddo
!
!! faces included in a link l
!allocate(face_in_l(1:num_links))
!do l=1,num_links
!  call get_faces_in_link_sc(sc,l,faces_l)
!  face_in_l(l)%num_=size(faces_l)
!  allocate(face_in_l(l)%label_(1:size(faces_l)))
!  face_in_l(l)%label_=faces_l
!enddo
!  
!! sites included in a face f
!allocate(sites_in_f(1:num_faces))
!do f=1,num_faces
!  call get_face_components_sc(sc,f,sites_f)
!  sites_in_f(f)%num_=size(sites_f)
!  allocate(sites_in_f(f)%label_(1:size(sites_f)))
!  sites_in_f(f)%label_=sites_f
!enddo
!
!  
!! set size of Dirac matrix
!!sizeD=(dimG)*(num_sites+num_links+num_faces)
!
!end subroutine set_sc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! set Remez data
!subroutine set_Remez_data
!implicit none
!
!integer i
!
!#ifdef PARALLEL
!if (MYRANK==0) then
!#endif
!
!open(Remez4_FILE, file=Remez_1ovminus4, status='OLD',action='READ')
!read(Remez4_FILE,*) N_Remez4
!read(Remez4_FILE,*) Remez_min4
!read(Remez4_FILE,*) Remez_max4
!allocate( Remez_alpha4(0:N_Remez4) )
!allocate( Remez_beta4(1:N_Remez4) )
!do i=0,N_Remez4 
!  read(Remez4_FILE,*) Remez_alpha4(i)
!enddo
!do i=1,N_Remez4 
!  read(Remez4_FILE,*) Remez_beta4(i)
!enddo
!close(Remez4_FILE)
!
!#ifdef PARALLEL
!endif
!call MPI_BCAST(N_Remez4,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
!call MPI_BCAST(Remez_min4,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!call MPI_BCAST(Remez_max4,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!if( MYRANK .ne. 0) then
!  allocate( Remez_alpha4(0:N_Remez4) )
!  allocate( Remez_beta4(1:N_Remez4) )
!endif
!call MPI_BCAST(Remez_alpha4,N_Remez4+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!call MPI_BCAST(Remez_beta4,N_Remez4,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!#endif
!
!!! rescaling
!Remez_alpha4(0)=Remez_alpha4(0)*Remez_factor4**(-0.25d0)
!do i=1,N_Remez4
!  Remez_alpha4(i)=Remez_alpha4(i)*Remez_factor4**(0.75d0)
!  Remez_beta4(i)=Remez_beta4(i)*Remez_factor4
!enddo
!
!
!
!!!
!#ifdef PARALLEL
!if (MYRANK==0) then
!#endif
!open(Remez8_FILE, file=Remez_1ov8, status='OLD',action='READ')
!read(Remez8_FILE,*) N_Remez8
!read(Remez8_FILE,*) Remez_min8
!read(Remez8_FILE,*) Remez_max8
!allocate( Remez_alpha8(0:N_Remez8) )
!allocate( Remez_beta8(1:N_Remez8) )
!do i=0,N_Remez8 
!  read(Remez8_FILE,*) Remez_alpha8(i)
!enddo
!do i=1,N_Remez8 
!  read(Remez8_FILE,*) Remez_beta8(i)
!enddo
!close(Remez8_FILE)
!
!#ifdef PARALLEL
!endif
!call MPI_BCAST(N_Remez8,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
!call MPI_BCAST(Remez_min8,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!call MPI_BCAST(Remez_max8,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!if( MYRANK .ne. 0) then
!  allocate( Remez_alpha8(0:N_Remez8) )
!  allocate( Remez_beta8(1:N_Remez8) )
!endif
!call MPI_BCAST(Remez_alpha8,N_Remez8+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!call MPI_BCAST(Remez_beta8,N_Remez8,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!#endif
!
!!! rescaling
!Remez_alpha8(0)=Remez_alpha8(0)*Remez_factor8**(0.125d0)
!do i=1,N_Remez8
!  Remez_alpha8(i)=Remez_alpha8(i)*Remez_factor8**(1.125d0)
!  Remez_beta8(i)=Remez_beta8(i)*Remez_factor8
!enddo
!end subroutine set_Remez_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! set sub simplicial complex
!subroutine set_sub_SC
!#ifdef PARALLEL
!use parallel
!#endif
!implicit none
!integer :: s,l,f,i
!
!
!#ifdef PARALLEL
!  ! sub_SC file is read only by MYRANK==0. 
!  num_sub_SC=0
!  !! set (rank, label) for a global site/link/face label
!  allocate( local_site_of_global(1:num_sites) )
!  allocate( local_link_of_global(1:num_links) )
!  allocate( local_face_of_global(1:num_faces) )
!
!  !! SUBSC_FILEを読み取って、local_{site/link/face}_of_global(:)を設定
!  !! すなわち、各ノードが担当するsite/link/faceを設定
!  !!   num_local_sites/links/faces
!  !!   local_{site/link/face}_of_global(:)
!  call set_map_from_global_to_local_data
!
!  !call stop_for_test
!  !! 各ノードが計算に必要なsiteの情報と、site情報の送り先、受信元を設定する
!  !!   num_necessary_sites
!  !!   global_site_of_local(s)
!  !!   send_sites(s)
!  !!   recv_sites(s)
!  call set_local_sites
!  !! 各ノードが計算に必要なlinkの情報と、link情報の送り先、受信元を設定する
!  !!   num_necessary_links
!  !!   global_link_of_local(l)
!  !!   send_links(l)
!  !!   recv_links(l)
!  call set_local_links
!  !! 各ノードが計算に必要なfaceの情報を設定する
!  !!   num_necessary_faces
!  !!   global_face_of_local(l)
!  call set_local_faces
!
!  !write(*,*) MYRANK,":",global_face_of_local
!
!
!  !! localなlinkのorgとtipを設定
!  call set_local_link_org_and_tip
!  !! localなsiteから出入りするlocal linkのラベルを設定
!  call set_local_link_fromto_site
!  !! localなfaceに含まれるlinkのlocalラベルを設定
!  call set_local_links_in_f
!  !! localなlinkを共有するfaceのlocalラベルを設定
!  call set_local_face_in_l
!  !! localなfaceに含まれるsiteのlocalラベルを設定
!  !! ただし、使うのは sites_in_f(f)%label_(1)だけなので、
!  !! はじめから size=1 にしておく。
!  call set_local_sites_in_f
!
!
!#else
!  num_local_sites=num_sites
!  num_local_links=num_links
!  num_local_faces=num_faces
!  num_necessary_sites=num_sites
!  num_necessary_links=num_links
!  num_necessary_faces=num_faces
!  global_num_sites=num_sites
!  global_num_links=num_links
!  global_num_faces=num_faces
!  allocate( global_site_of_local(1:num_sites) )
!  allocate( global_link_of_local(1:num_links) )
!  allocate( global_face_of_local(1:num_faces) )
!
!  do s=1,num_sites
!    global_site_of_local(s)=s
!  enddo
!  do l=1,num_links
!    global_link_of_local(l)=l
!  enddo
!  do f=1,num_faces
!    global_face_of_local(f)=f
!  enddo
!#endif
!
!end subroutine set_sub_SC
!
!#ifdef PARALLEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine set_map_from_global_to_local_data
!use parallel
!implicit none
!
!integer :: tmp_num_local_sites,tmp_num_local_links,tmp_num_local_faces
!integer,allocatable :: tmp_global_site_of_local(:)
!integer :: s,l,f,rank,i,part
!
!if (MYRANK == 0) then
!  open(SUBSC_FILE, file=FsubSC, status='OLD',action='READ')
!  read(SUBSC_FILE,*) num_sub_SC
!  ! send num_sub_SC to all the other nodes
!  call MPI_BCAST(num_sub_SC,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
!
!  if (num_sub_SC .ne. NPROCS) then
!    write(*,*) "number of core is mismatch."
!    close(SUBSC_FILE)
!    call MPI_FINALIZE(IERR)
!    stop
!  endif
!  do part=1,num_sub_SC
!    ! sub_SCの構造は、
!    !   ・Num_Sub_SC
!    ! の後に、
!    !   ・RANK
!    !   ・number of sites in the node
!    !   ・そのnodeが担当するsiteのリスト
!    ! が順番に並んでいる。
!    
!    !! RANKとそれが担当するノードを読み込む
!    read(SUBSC_FILE,*) rank
!    read(SUBSC_FILE,*) tmp_num_local_sites
!    allocate( tmp_global_site_of_local(1:tmp_num_local_sites) )
!    read(SUBSC_FILE,*) (tmp_global_site_of_local(i),i=1,tmp_num_local_sites)
!    !  
!    !! set local_site_of_global(:) to a temporary variable
!    do i=1,tmp_num_local_sites
!      local_site_of_global(tmp_global_site_of_local(i))%rank_=rank
!      local_site_of_global(tmp_global_site_of_local(i))%label_=i
!    enddo
!    tmp_num_local_links=0
!    do l=1,num_links
!      do s=1,tmp_num_local_sites
!        if( tmp_global_site_of_local(s)==link_org(l) ) then
!          tmp_num_local_links=tmp_num_local_links+1
!          local_link_of_global(l)%rank_=rank
!          local_link_of_global(l)%label_=tmp_num_local_links
!        endif
!      enddo
!    enddo
!    tmp_num_local_faces=0
!    do f=1,num_faces
!      do s=1,tmp_num_local_sites
!        if( tmp_global_site_of_local(s) == sites_in_f(f)%label_(1)) then
!          tmp_num_local_faces=tmp_num_local_faces+1
!          i=i+1
!          local_face_of_global(f)%rank_=rank
!          local_face_of_global(f)%label_=tmp_num_local_faces
!        endif
!      enddo
!    enddo
!    !!
!    if( rank == 0 ) then
!      num_local_sites=tmp_num_local_sites
!      num_local_links=tmp_num_local_links
!      num_local_faces=tmp_num_local_faces
!    else
!      call MPI_SEND(tmp_num_local_sites,1,MPI_INTEGER,rank,1,MPI_COMM_WORLD,IERR)
!      call MPI_SEND(tmp_num_local_links,1,MPI_INTEGER,rank,2,MPI_COMM_WORLD,IERR)
!      call MPI_SEND(tmp_num_local_faces,1,MPI_INTEGER,rank,3,MPI_COMM_WORLD,IERR)
!    endif
!    deallocate( tmp_global_site_of_local ) 
!  enddo
!  close(SUBSC_FILE)
!  else
!    call MPI_BCAST(num_sub_SC,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
!    call MPI_RECV(num_local_sites,1,MPI_INTEGER,0,1,MPI_COMM_WORLD,ISTATUS,IERR)
!    call MPI_RECV(num_local_links,1,MPI_INTEGER,0,2,MPI_COMM_WORLD,ISTATUS,IERR)
!    call MPI_RECV(num_local_faces,1,MPI_INTEGER,0,3,MPI_COMM_WORLD,ISTATUS,IERR)
!endif
!call MPI_BCAST(local_site_of_global,num_sites,MPI_LOCAL_LABEL,0,MPI_COMM_WORLD,IERR)
!call MPI_BCAST(local_link_of_global,num_links,MPI_LOCAL_LABEL,0,MPI_COMM_WORLD,IERR)
!call MPI_BCAST(local_face_of_global,num_faces,MPI_LOCAL_LABEL,0,MPI_COMM_WORLD,IERR)
!
!end subroutine set_map_from_global_to_local_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!! 計算に必要なsiteの情報を設定
!subroutine set_local_sites
!use parallel
!implicit none
!
!integer s,l,i,nsend,nrecv,info,k
!integer :: tmp_global_site_of_local(1:num_sites)
!!integer :: tmp_send_sites(1:num_sites)
!!integer :: tmp_recv_sites(1:num_sites)
!type(LOCAL_LABEL) :: tmp_send_sites(1:num_sites)
!type(LOCAL_LABEL) :: tmp_recv_sites(1:num_sites)
!
!!! 計算を担当するsiteを tmp_global_site_of_local の最初の方に設定
!num_necessary_sites=0
!tmp_global_site_of_local=0
!do s=1,num_sites
!  if( local_site_of_global(s)%rank_ == MYRANK ) then
!    num_necessary_sites = num_necessary_sites + 1
!    tmp_global_site_of_local(num_necessary_sites)=s
!  endif
!enddo
!!! この段階で、num_necesarry_sites=num_local_sites
!!! linkの始点と終点が別のnodeに属していたら、
!!! tipの位置を、orgを担当するnodeに送る必要がある
!!! (1) 担当するlinkの始点と終点
!!! (2) 担当するサイトを始点とするlinkの終点
!!! (3) 担当するサイトを終点とするlinkの始点
!!! 重複があり得るので、逐一チェックしながら進める
!nsend=0
!nrecv=0
!!! (1) 担当するlinkの始点と終点
!do l=1,num_links
!  if( local_site_of_global(link_org(l))%rank_ &
!        /= local_site_of_global(link_tip(l))%rank_ ) then
!    !!!!!!!!!!!
!    if( local_site_of_global(link_tip(l))%rank_ == MYRANK ) then
!      !! 重複チェック
!      info=0
!      do i=1,nsend
!        if( tmp_send_sites(i)%rank_ == local_site_of_global(link_org(l))%rank_ &
!          .and. &
!            tmp_send_sites(i)%label_ == link_tip(l) ) then
!          info=1
!          exit
!        endif
!      enddo
!      if (info==0) then 
!        nsend=nsend+1
!        tmp_send_sites(nsend)%rank_ = local_site_of_global(link_org(l))%rank_ 
!        tmp_send_sites(nsend)%label_ = link_tip(l) 
!      endif
!    !!!!!!!!!!!
!    elseif( local_site_of_global(link_org(l))%rank_ == MYRANK ) then
!      !! 重複チェック
!      info=0
!      do i=1,nrecv
!        if( tmp_recv_sites(i)%rank_ == local_site_of_global(link_tip(l))%rank_ &
!          .and. &
!            tmp_recv_sites(i)%label_ == link_tip(l) ) then
!          info=1
!          exit
!        endif
!      enddo
!      if (info==0) then 
!        nrecv=nrecv+1
!        tmp_recv_sites(nrecv)%rank_ = local_site_of_global(link_tip(l))%rank_ 
!        tmp_recv_sites(nrecv)%label_ = link_tip(l) 
!        tmp_global_site_of_local(num_local_sites+nrecv)=link_tip(l)
!      endif
!    endif
!  endif
!enddo
!!! (2) 担当するサイトを終点とするlinkの始点
!do s=1,num_sites
!  do k=1,linkorg_to_s(s)%num_
!    if( local_site_of_global(linkorg_to_s(s)%sites_(k))%rank_ &
!        /= local_site_of_global(s)%rank_ ) then
!      if(local_site_of_global(linkorg_to_s(s)%sites_(k))%rank_ == MYRANK) then
!        !! 重複チェック
!        info=0
!        do i=1,nsend
!          if( tmp_send_sites(i)%rank_ &
!              == local_site_of_global(s)%rank_ &
!            .and. &
!              tmp_send_sites(i)%label_  &
!              == linkorg_to_s(s)%sites_(k) &
!              ) then
!            info=1
!            exit
!          endif
!        enddo
!        if (info==0) then 
!          !write(*,*) "(S)site",linkorg_to_s(s)%sites_(k),"from",MYRANK,"to",local_site_of_global(s)%rank_
!          nsend=nsend+1
!          tmp_send_sites(nsend)%rank_ = local_site_of_global(s)%rank_ 
!          tmp_send_sites(nsend)%label_ = linkorg_to_s(s)%sites_(k)
!        endif
!  !!!!!!!!!!!
!      elseif( local_site_of_global(s)%rank_ == MYRANK ) then
!        !! 重複チェック
!        info=0
!        do i=1,nrecv
!          if( tmp_recv_sites(i)%rank_ &
!            == local_site_of_global(linkorg_to_s(s)%sites_(k))%rank_ &
!            .and. &
!              tmp_recv_sites(i)%label_ == linkorg_to_s(s)%sites_(k) ) then
!            info=1
!            exit
!          endif
!        enddo
!        if (info==0) then 
!          !write(*,*) "(R)site",linkorg_to_s(s)%sites_(k),"from",local_site_of_global(linkorg_to_s(s)%sites_(k))%rank_,"to",MYRANK
!          nrecv=nrecv+1
!          tmp_recv_sites(nrecv)%rank_  &
!            = local_site_of_global(linkorg_to_s(s)%sites_(k))%rank_ 
!          tmp_recv_sites(nrecv)%label_ = linkorg_to_s(s)%sites_(k)
!          tmp_global_site_of_local(num_local_sites+nrecv) = linkorg_to_s(s)%sites_(k)
!        endif
!      endif
!    endif
!  enddo
!enddo
!
!!! (3) 担当するサイトを始点とするlinkの終点
!do s=1,num_sites
!  do k=1,linktip_from_s(s)%num_
!    if( local_site_of_global(linktip_from_s(s)%sites_(k))%rank_  &
!        /= local_site_of_global(s)%rank_ ) then
!      if(local_site_of_global(linktip_from_s(s)%sites_(k))%rank_ == MYRANK) then
!        !! 重複チェック
!        info=0
!        do i=1,nsend
!          if( tmp_send_sites(i)%rank_ &
!              == local_site_of_global(s)%rank_ &
!            .and. &
!              tmp_send_sites(i)%label_  &
!              == linktip_from_s(s)%sites_(k)) then
!            info=1
!            exit
!          endif
!        enddo
!        if (info==0) then 
!          nsend=nsend+1
!          tmp_send_sites(nsend)%rank_ = local_site_of_global(s)%rank_ 
!          tmp_send_sites(nsend)%label_ = linktip_from_s(s)%sites_(k)
!        endif
!  !!!!!!!!!!!
!      elseif( local_site_of_global(s)%rank_ == MYRANK ) then
!        !! 重複チェック
!        info=0
!        do i=1,nrecv
!          if( tmp_recv_sites(i)%rank_ &
!            == local_site_of_global(linktip_from_s(s)%sites_(k))%rank_ &
!            .and. &
!              tmp_recv_sites(i)%label_ == linktip_from_s(s)%sites_(k) ) then
!            info=1
!            exit
!          endif
!        enddo
!        if (info==0) then 
!          nrecv=nrecv+1
!          tmp_recv_sites(nrecv)%rank_  &
!            = local_site_of_global(linktip_from_s(s)%sites_(k))%rank_ 
!          tmp_recv_sites(nrecv)%label_ = linktip_from_s(s)%sites_(k)
!          tmp_global_site_of_local(num_local_sites+nrecv)=linktip_from_s(s)%sites_(k)
!        endif
!      endif
!    endif
!  enddo
!enddo
!
!!! 重複度を除いて配列を定義
!num_send_sites=nsend
!num_recv_sites=nrecv
!allocate( send_sites(num_send_sites) )
!allocate( recv_sites(num_recv_sites) )
!num_necessary_sites = num_local_sites + num_recv_sites
!allocate( global_site_of_local(1:num_necessary_sites) )
!
!!! 代入
!do s=1,num_send_sites
!  send_sites(s)%rank_ = tmp_send_sites(s)%rank_
!  send_sites(s)%label_ = local_site_of_global(tmp_send_sites(s)%label_)%label_
!enddo
!do s=1,num_recv_sites
!  recv_sites(s)%rank_ = tmp_recv_sites(s)%rank_
!  recv_sites(s)%label_ = num_local_sites+s
!enddo
!do s=1,num_necessary_sites
!  global_site_of_local(s) = tmp_global_site_of_local(s)
!enddo
!
!!! 通信する順番を定義
!
!
!
!end subroutine set_local_sites
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!! 計算に必要なlinkの情報を設定
!subroutine set_local_links
!use parallel
!implicit none
!
!integer s,l,f,i,j,k,ll,ll_label
!integer :: tmp_global_link_of_local(1:num_links)
!type(LOCAL_LABEL) :: tmp_send_links(1:num_links)
!type(LOCAL_LABEL) :: tmp_recv_links(1:num_links)
!integer :: nsend,nrecv,info
!
!
!do l=1,num_links
!  tmp_global_link_of_local(l)=0
!  tmp_send_links(l)%rank_=-1
!  tmp_send_links(l)%label_=-1
!  tmp_recv_links(l)%rank_=-1
!  tmp_recv_links(l)%label_=-1
!enddo
!
!!! 計算を担当するlinkを tmp_global_link_of_local の最初の方に設定
!num_necessary_links=0
!tmp_global_link_of_local=0
!do l=1,num_links
!  if( local_link_of_global(l)%rank_ == MYRANK ) then
!    num_necessary_links = num_necessary_links + 1
!    tmp_global_link_of_local(num_necessary_links)=l
!  endif
!enddo
!!! この段階で、num_necesarry_links=num_local_links
!!! 計算に必要なのは、
!!!   (1) 担当するsiteをtipとするlink
!!!   (2) 担当するfaceを構成するlink
!!!   (3) 担当するlinkを共有するfaceに含まれる全link
!!! 重複があり得るので、逐一チェックしながら進める
!
!!!   (1) 担当するsiteをtipとするlink
!nsend=0
!nrecv=0
!do l=1,num_links
!  if( local_site_of_global(link_tip(l))%rank_ /= &
!      local_site_of_global(link_org(l))%rank_ ) then
!    if( local_site_of_global(link_org(l))%rank_ == MYRANK ) then
!      !! 重複チェック
!      info=0
!      do i=1,nsend
!        if( tmp_send_links(i)%rank_ == local_site_of_global(link_tip(l))%rank_ &
!          .and. &
!            tmp_send_links(i)%label_ == l ) then
!          info=1
!          exit
!        endif
!      enddo
!      if( info == 0 ) then
!        nsend=nsend+1
!        tmp_send_links(nsend)%rank_ = local_site_of_global(link_tip(l))%rank_
!        tmp_send_links(nsend)%label_ = l
!      endif
!    !!!!!!!!!!!
!    elseif( local_site_of_global(link_tip(l))%rank_ == MYRANK ) then
!      !! 重複チェック
!      info=0
!      do i=1,nrecv
!        if( tmp_recv_links(i)%rank_ == local_site_of_global(link_org(l))%rank_ &
!          .and. &
!            tmp_recv_links(i)%label_ == l ) then
!          info=1
!          exit
!        endif
!      enddo
!      if (info==0) then 
!        nrecv=nrecv+1
!        tmp_recv_links(nrecv)%rank_ = local_site_of_global(link_org(l))%rank_ 
!        tmp_recv_links(nrecv)%label_ = l
!        tmp_global_link_of_local(num_local_links+nrecv)=l
!      endif
!    endif
!  endif
!enddo
!
!!!   (2) 担当するfaceを構成するlink
!do f=1,num_faces
!  do j=1,links_in_f(f)%num_
!    l=links_in_f(f)%link_labels_(j)
!    if( local_face_of_global(f)%rank_ /= local_link_of_global(l)%rank_ ) then
!      if( local_link_of_global(l)%rank_ == MYRANK ) then
!        !! 重複チェック
!        info=0
!        do i=1,nsend
!          if( tmp_send_links(i)%rank_ == local_face_of_global(f)%rank_ &
!            .and. &
!              tmp_send_links(i)%label_ == l ) then
!            info=1
!            exit
!          endif
!        enddo
!        if( info == 0 ) then
!          nsend=nsend+1
!          tmp_send_links(nsend)%rank_ = local_face_of_global(f)%rank_ 
!          tmp_send_links(nsend)%label_ = l
!        endif
!      !!!!!!!!!!!
!      elseif( local_face_of_global(f)%rank_ == MYRANK ) then
!        !! 重複チェック
!        info=0
!        do i=1,nrecv
!          if( tmp_recv_links(i)%rank_ == local_link_of_global(l)%rank_ &
!            .and. &
!              tmp_recv_links(i)%label_ == l ) then
!            info=1
!            exit
!          endif
!        enddo
!        if (info==0) then 
!          nrecv=nrecv+1
!          tmp_recv_links(nrecv)%rank_ = local_link_of_global(l)%rank_
!          tmp_recv_links(nrecv)%label_ = l
!          tmp_global_link_of_local(num_local_links+nrecv)=l
!        endif
!      endif
!    endif
!  enddo
!enddo
!
!!!   (3) 担当するlinkを共有するfaceに含まれる全link
!do l=1,num_links
!  do k=1,face_in_l(l)%num_
!    f=face_in_l(l)%label_(k)
!    do ll_label=1,links_in_f(f)%num_
!      ll=links_in_f(f)%link_labels_(ll_label)
!      if( local_link_of_global(l)%rank_ /= local_link_of_global(ll)%rank_ ) then 
!        if( local_link_of_global(ll)%rank_ == MYRANK ) then
!          !! 重複チェック
!          info=0
!          do i=1,nsend
!            if( tmp_send_links(i)%rank_ == local_link_of_global(l)%rank_ &
!                .and. &
!                tmp_send_links(i)%label_ == ll ) then 
!              info=1
!              exit
!            endif
!          enddo
!          if (info == 0 ) then
!            nsend=nsend+1
!            tmp_send_links(nsend)%rank_ = local_link_of_global(l)%rank_ 
!            tmp_send_links(nsend)%label_ = ll
!          endif
!        !!!!!!!!!!!
!        elseif( local_link_of_global(l)%rank_ == MYRANK ) then
!          !! 重複チェック
!          info=0
!          do i=1,nrecv
!            if( tmp_recv_links(i)%rank_ == local_link_of_global(ll)%rank_ &
!                .and. &
!                tmp_recv_links(i)%label_ == ll ) then 
!              info=1
!              exit
!            endif
!          enddo
!          if (info == 0 ) then
!            nrecv=nrecv+1
!            tmp_recv_links(nrecv)%rank_ = local_link_of_global(ll)%rank_ 
!            tmp_recv_links(nrecv)%label_ = ll
!            tmp_global_link_of_local(num_local_links+nrecv)=ll
!          endif
!        endif
!      endif
!    enddo
!  enddo
!enddo
!
!!! 重複度を除いて配列を定義
!num_send_links=nsend
!num_recv_links=nrecv
!allocate( send_links(num_send_links) )
!allocate( recv_links(num_recv_links) )
!num_necessary_links = num_local_links + num_recv_links
!allocate( global_link_of_local(1:num_necessary_links) )
!
!!! 代入
!do l=1,num_send_links
!  send_links(l)%rank_ = tmp_send_links(l)%rank_
!  send_links(l)%label_ = local_link_of_global(tmp_send_links(l)%label_)%label_
!enddo
!do l=1,num_recv_links
!  recv_links(l)%rank_ = tmp_recv_links(l)%rank_
!  recv_links(l)%label_ = num_local_links+l
!enddo
!do l=1,num_necessary_links
!  global_link_of_local(l) = tmp_global_link_of_local(l)
!enddo
!
!end subroutine set_local_links
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!! 計算に必要なfaceの情報を設定
!subroutine set_local_faces
!use parallel
!implicit none
!
!integer f,l,i,k
!type(LOCAL_LABEL) :: tmp_send_faces(1:num_faces)
!type(LOCAL_LABEL) :: tmp_recv_faces(1:num_faces)
!integer :: nsend,nrecv,info
!integer :: tmp_global_face_of_local(1:num_faces)
!
!!! 計算を担当するfaceを tmp_global_face_of_local の最初の方に設定
!num_necessary_faces=0
!do f=1,num_faces
!  if( local_face_of_global(f)%rank_ == MYRANK ) then
!    num_necessary_faces = num_necessary_faces + 1
!    tmp_global_face_of_local(num_necessary_faces)=f
!  endif
!enddo
!
!nsend=0
!nrecv=0
!!! 担当するlinkを共有するfaceが全て必要
!do l=1,num_links
!  do k=1,face_in_l(l)%num_
!    f=face_in_l(l)%label_(k)
!    if( local_link_of_global(l)%rank_ &
!        /= local_face_of_global(f)%rank_ ) then
!      if( local_face_of_global(f)%rank_ == MYRANK ) then
!        !重複チェック
!        info=0
!        do i=1,nsend
!          if( tmp_send_faces(i)%rank_ == local_link_of_global(l)%rank_ &
!              .and. &
!              tmp_send_faces(i)%label_ == f ) then 
!            info=1
!            exit
!          endif
!        enddo
!        if (info == 0 ) then
!          nsend=nsend+1
!          tmp_send_faces(nsend)%rank_ = local_link_of_global(l)%rank_ 
!          tmp_send_faces(nsend)%label_ = f
!        endif
!      !!!!!!!!!!!
!      elseif( local_link_of_global(l)%rank_ == MYRANK ) then
!        !! 重複チェック
!        info=0
!        do i=1,nrecv
!          if( tmp_recv_faces(i)%rank_ == local_face_of_global(f)%rank_ &
!              .and. &
!              tmp_recv_faces(i)%label_ == f ) then 
!            info=1
!            exit
!          endif
!        enddo
!        if (info == 0 ) then
!          nrecv=nrecv+1
!          tmp_recv_faces(nrecv)%rank_ = local_face_of_global(f)%rank_ 
!          tmp_recv_faces(nrecv)%label_ = f
!          tmp_global_face_of_local(num_local_faces+nrecv)=f
!        endif
!      endif
!    endif
!  enddo
!enddo
!
!
!!! 重複度を除いて配列を定義
!num_send_faces=nsend
!num_recv_faces=nrecv
!allocate( send_faces(num_send_faces) )
!allocate( recv_faces(num_recv_faces) )
!num_necessary_faces = num_local_faces + num_recv_faces
!allocate( global_face_of_local(1:num_necessary_faces) )
!
!!! 代入
!do i=1,num_send_faces
!  send_faces(i)%rank_ = tmp_send_faces(i)%rank_
!  send_faces(i)%label_ = local_face_of_global(tmp_send_faces(i)%label_)%label_
!enddo
!do i=1,num_recv_faces
!  recv_faces(i)%rank_ = tmp_recv_faces(i)%rank_
!  recv_faces(i)%label_ = num_local_faces+i
!enddo
!do f=1,num_necessary_faces
!  global_face_of_local(f) = tmp_global_face_of_local(f)
!enddo
!
!end subroutine set_local_faces
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!! それぞれのnode上でのlink_tipとlink_orgを設定
!subroutine set_local_link_org_and_tip
!use parallel
!implicit none
!
!integer info_o, info_t
!integer l,s
!
!allocate( local_link_org(1:num_necessary_links) )
!allocate( local_link_tip(1:num_necessary_links) )
!
!do l=1,num_necessary_links
!  info_o=0
!  info_t=0
!  do s=1,num_necessary_sites
!    if( global_site_of_local(s) == link_org(global_link_of_local(l)) ) then
!      local_link_org(l) = s
!      info_o=1
!    elseif( global_site_of_local(s) == link_tip(global_link_of_local(l)) ) then
!      local_link_tip(l) = s
!      info_t=1
!    endif
!    if( info_o * info_t == 1 ) exit
!  enddo
!  if( info_o==0 ) local_link_org(l)=0
!  if( info_t==0 ) local_link_tip(l)=0
!  !if( info_o * info_t == 0 ) then 
!  !  write(*,*) "something is wrong in local link", l, "in RANK",MYRANK
!  !  write(*,*) "  global_link:",global_link_of_local(l),&
!  !    "org:",(local_link_org(l)),&
!  !    "tip:",(local_link_tip(l))
!  !  call stop_for_test
!  !endif
!enddo
!
!
!end subroutine set_local_link_org_and_tip
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!! それぞれのnode上で、siteに出入りするlinkのラベルを設定
!subroutine set_local_link_fromto_site
!use parallel
!implicit none
!
!integer s,ss,l,num,info,i
!
!allocate(local_linktip_from_s(1:num_local_sites))
!allocate(local_linkorg_to_s(1:num_local_sites))
!do s=1,num_local_sites
!  !! local_linktip
!  num=linktip_from_s(global_site_of_local(s))%num_
!  local_linktip_from_s(s)%num_=num
!  allocate(local_linktip_from_s(s)%labels_(1:num) )
!  allocate(local_linktip_from_s(s)%sites_(1:num) )!使わないので設定はしない
!  do i=1,num
!    info=0
!    do l=1,num_necessary_links
!      if( linktip_from_s(global_site_of_local(s))%labels_(i) &
!          == global_link_of_local(l) ) then
!        local_linktip_from_s(s)%labels_(i)=l
!        info=1
!        exit
!      endif
!    enddo
!    if( info==0 ) then 
!      write(*,*) MYRANK,"Something happpend in set_local_linktip_from_s"
!      call stop_for_test
!    endif
!  enddo
!  !! local_linkorg
!  num=linkorg_to_s(global_site_of_local(s))%num_
!  local_linkorg_to_s(s)%num_= num
!  allocate(local_linkorg_to_s(s)%labels_(1:num) )
!  allocate(local_linkorg_to_s(s)%sites_(1:num) ) !使わないので設定はしない
!  do i=1,num
!    do l=1,num_necessary_links
!      info=0
!      if( linkorg_to_s(global_site_of_local(s))%labels_(i) &
!          == global_link_of_local(l) ) then
!        local_linkorg_to_s(s)%labels_(i)=l
!        info=1
!        exit
!      endif
!    enddo
!    if( info==0 ) then 
!      write(*,*) MYRANK,"Something happpend in set_local_linkorg_to_s 1"
!      call stop_for_test
!    endif
!  enddo
!enddo
!
!end subroutine set_local_link_fromto_site
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!! それぞれのnode上で、faceを構成するlinkのラベルを設定
!subroutine  set_local_links_in_f
!use parallel
!implicit none
!
!integer f,l,ll, num, g_link, dir
!
!allocate( local_links_in_f(1:num_necessary_faces) )
!do f=1,num_necessary_faces
!  num = links_in_f(global_face_of_local(f))%num_
!  local_links_in_f(f)%num_ = num
!  allocate( local_links_in_f(f)%link_labels_(1:num) )
!  allocate( local_links_in_f(f)%link_dirs_(1:num))
!
!  do l=1,num
!    g_link = links_in_f(global_face_of_local(f))%link_labels_(l)
!    dir = links_in_f(global_face_of_local(f))%link_dirs_(l)
!    do ll=1,num_necessary_links
!      if( global_link_of_local(ll) == g_link ) then
!        local_links_in_f(f)%link_labels_(l) = ll
!        local_links_in_f(f)%link_dirs_(l) = dir
!      endif
!    enddo
!  enddo
!enddo
!
!end subroutine  set_local_links_in_f
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!! それぞれのnode上で、linkに共有されるfaceのラベルを設定
!subroutine  set_local_face_in_l
!use parallel
!implicit none
!
!integer l,num ,i,f,g_face,info
!
!allocate( local_face_in_l(1:num_local_links) )
!do l=1,num_local_links
!  num = face_in_l(global_link_of_local(l))%num_
!  local_face_in_l(l)%num_ = num
!  allocate( local_face_in_l(l)%label_(1:num) )
!
!  do i=1,num
!    info=0
!    g_face = face_in_l(global_link_of_local(l))%label_(i)
!    do f=1,num_necessary_faces
!      if( global_face_of_local(f) == g_face ) then
!        local_face_in_l(l)%label_(i) = f
!        info=1
!      endif
!    enddo
!    if( info == 0 ) then 
!      write(*,*) "something happened in setting local_face_in_l."
!      write(*,*) MYRANK, f,i
!      call stop_for_test
!    endif
!  enddo
!enddo
!
!
!end subroutine  set_local_face_in_l
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!! それぞれのnode上で担当するfaceのlocalな代表点を設定
!subroutine set_local_sites_in_f
!use parallel
!implicit none
!
!integer f,s,global_s,info
!
!allocate( local_sites_in_f(1:num_necessary_faces) )
!do f=1,num_necessary_faces
!  local_sites_in_f(f)%num_=1 
!  allocate( local_sites_in_f(f)%label_(1:1) )
!
!  global_s=sites_in_f(global_face_of_local(f))%label_(1)
!  do s=1,num_necessary_sites
!    
!    if( global_s == global_site_of_local(s) ) then
!      local_sites_in_f(f)%label_(1) = s
!      info=1
!      exit
!    endif
!  enddo
!  if( info == 0 ) then
!    write(*,*) "something happened in set_local_sites_in_f"
!    call stop_for_test
!  endif
!enddo
!
!end subroutine set_local_sites_in_f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! set alpha and beta
!subroutine set_local_alpha_beta
!use parallel
!implicit none
!double precision tmp_alpha_s(1:global_num_sites) 
!double precision tmp_alpha_l(1:global_num_links) 
!double precision tmp_alpha_f(1:global_num_faces)
!double precision tmp_beta_f(1:global_num_faces)
!integer s,l,f
!integer :: local, rank,  tag
!
!
!if(MYRANK==0) then
!  tmp_alpha_s=alpha_s
!  tmp_alpha_l=alpha_l
!  tmp_alpha_f=alpha_f
!  tmp_beta_f=beta_f
!endif
!
!deallocate( alpha_s,alpha_l,alpha_f,beta_f)
!allocate ( alpha_s(1:num_necessary_sites) )
!allocate ( alpha_l(1:num_necessary_links) )
!allocate ( alpha_f(1:num_necessary_faces) )
!allocate ( beta_f(1:num_necessary_faces) )
! 
!!if(MYRANK==0) then
!!if (read_alpha==1) then
!!open(ALPHA_BETA_FILE, file=ALPHA_BETA, status='old',action='READ')
!!  read(ALPHA_BETA_FILE,'()') ! skip 1 line
!!  do s=1,global_num_sites
!!    read(ALPHA_BETA_FILE,*) tmp_alpha_s(s)
!!  enddo
!!  read(ALPHA_BETA_FILE,'()') !! skip 1 line 
!!  do l=1,global_num_links
!!    read(ALPHA_BETA_FILE,*) tmp_alpha_l(l)
!!  enddo
!!  read(ALPHA_BETA_FILE,'()')  !! skip 1 line 
!!  do f=1,global_num_faces
!!    read(ALPHA_BETA_FILE,*) tmp_alpha_f(f)
!!  enddo
!!  read(ALPHA_BETA_FILE,'()') !! skip 1 line 
!!  do f=1,global_num_faces
!!    read(ALPHA_BETA_FILE,*) tmp_beta_f(f)
!!  enddo
!!close(ALPHA_BETA_FILE)
!!else
!!    tmp_alpha_s=1d0
!!    tmp_alpha_l=1d0
!!    tmp_alpha_f=1d0
!!    tmp_beta_f=1d0
!!endif
!!#ifdef PARALLEL
!!endif
!
!do s=1,global_num_sites
!  local=local_site_of_global(s)%label_
!  rank=local_site_of_global(s)%rank_
!  tag=s
!  !if(MYRANK==0) then
!    !write(*,*) s,local,rank,tmp_alpha_s(s)
!  !endif
!  if( MYRANK == 0 ) then
!    if( rank == 0 ) then
!      alpha_s(local)=tmp_alpha_s(s)
!    else
!      call MPI_SEND(tmp_alpha_s(s),1,MPI_DOUBLE_PRECISION,rank,tag,MPI_COMM_WORLD,IERR)
!    endif
!  elseif( rank == MYRANK ) then
!    call MPI_RECV(alpha_s(local),1,MPI_DOUBLE_PRECISION, 0,tag,MPI_COMM_WORLD,ISTATUS,IERR)
!  endif
!enddo
!
!do l=1,global_num_links
!  local=local_link_of_global(l)%label_
!  rank=local_link_of_global(l)%rank_
!  tag=l+global_num_sites
!  if( MYRANK == 0 ) then
!    if( rank == 0 ) then
!      alpha_l(local)=tmp_alpha_l(l)
!    else
!      call MPI_SEND(tmp_alpha_l(l),1,MPI_DOUBLE_PRECISION,rank,tag,MPI_COMM_WORLD,IERR)
!    endif
!  elseif( rank == MYRANK ) then
!    call MPI_RECV(alpha_l(local),1,MPI_DOUBLE_PRECISION, 0,tag,MPI_COMM_WORLD,ISTATUS,IERR)
!  endif
!enddo
!
!do f=1,global_num_faces
!  local=local_face_of_global(f)%label_
!  rank=local_face_of_global(f)%rank_
!  tag=f+global_num_sites+global_num_links
!  if( MYRANK == 0 ) then
!    if( rank == 0 ) then
!      alpha_f(local)=tmp_alpha_f(f)
!    else
!      call MPI_SEND(tmp_alpha_f(f),1,MPI_DOUBLE_PRECISION,rank,tag,MPI_COMM_WORLD,IERR)
!    endif
!  elseif( rank == MYRANK ) then
!    call MPI_RECV(alpha_f(local),1,MPI_DOUBLE_PRECISION, 0,tag,MPI_COMM_WORLD,ISTATUS,IERR)
!  endif
!enddo
!
!do f=1,global_num_faces
!  local=local_face_of_global(f)%label_
!  rank=local_face_of_global(f)%rank_
!  tag=f+global_num_sites+global_num_links+global_num_faces
!  if( MYRANK == 0 ) then
!    if( rank == 0 ) then
!      beta_f(local)=tmp_beta_f(f)
!    else
!      call MPI_SEND(tmp_beta_f(f),1,MPI_DOUBLE_PRECISION,rank,tag,MPI_COMM_WORLD,IERR)
!    endif
!  elseif( rank == MYRANK ) then
!    call MPI_RECV(beta_f(local),1,MPI_DOUBLE_PRECISION, 0,tag,MPI_COMM_WORLD,ISTATUS,IERR)
!  endif
!enddo
!
!call syncronize_ab(alpha_s,'S')
!call syncronize_ab(alpha_l,'L')
!call syncronize_ab(alpha_f,'F')
!call syncronize_ab(beta_f,'F')
!
!
!end subroutine set_local_alpha_beta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! alpha, betaを通信するsubroutine
!#ifdef PARALLEL
!subroutine syncronize_ab(alpha,C)
!use parallel
!double precision, intent(inout) :: alpha(:)
!character, intent(in) :: C
!
!integer :: num_send, num_recv
!integer :: s_send
!integer :: s_recv
!integer :: local, rank, tag
!integer, allocatable :: ISEND(:), IRECV(:) ! for MPI_WAIT 
!
!!!!!!!!!
!
!if( C=='S' ) then 
!  num_send=num_send_sites
!  num_recv=num_recv_sites
!elseif( C=='L' ) then 
!  num_send=num_send_links
!  num_recv=num_recv_links
!elseif( C=='F' ) then 
!  num_send=num_send_faces
!  num_recv=num_recv_faces
!endif
!
!allocate(ISEND(1:num_send))
!allocate(IRECV(1:num_recv))
!do s_send=1,num_send
!  if( C=='S' ) then
!    local=send_sites(s_send)%label_
!    rank=send_sites(s_send)%rank_
!    tag=10000*rank + global_site_of_local(local)
!  elseif( C=='L' ) then
!    local=send_links(s_send)%label_
!    rank=send_links(s_send)%rank_
!    tag=10000*rank + global_link_of_local(local)
!  elseif( C=='F' ) then
!    local=send_faces(s_send)%label_
!    rank=send_faces(s_send)%rank_
!    tag=10000*rank + global_face_of_local(local)
!  endif
!
!  call MPI_ISEND(alpha(local),1,MPI_DOUBLE_PRECISION,rank,tag,MPI_COMM_WORLD,ISEND(s_send),IERR)
!enddo
!
!do s_recv=1,num_recv
!  if( C=='S' ) then
!    local=recv_sites(s_recv)%label_
!    rank=recv_sites(s_recv)%rank_
!    tag=10000*MYRANK + global_site_of_local(local)
!  elseif( C=='L' ) then
!    local=recv_links(s_recv)%label_
!    rank=recv_links(s_recv)%rank_
!    tag=10000*MYRANK + global_link_of_local(local)
!  elseif( C=='F' ) then
!    local=recv_faces(s_recv)%label_
!    rank=recv_faces(s_recv)%rank_
!    tag=10000*MYRANK + global_face_of_local(local)
!  endif
!
!  call MPI_IRECV(alpha(local),1,MPI_DOUBLE_PRECISION,rank,tag,MPI_COMM_WORLD,IRECV(s_recv),IERR)
!enddo
!
!do s_send=1,num_send
!  call MPI_WAIT(ISEND(s_send),ISTATUS,IERR)
!enddo
!do s_recv=1,num_recv
!  call MPI_WAIT(IRECV(s_recv),ISTATUS,IERR)
!enddo
!
!deallocate(ISEND, IRECV)
!end subroutine syncronize_ab
!#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!! それぞれのnode上で必要な alpha, beta を設定する
!subroutine set_local_alpha_beta
!use parallel
!implicit none
!
!integer s,l,f
!
!allocate(local_alpha_s(1:num_necessary_sites) )
!allocate(local_alpha_l(1:num_necessary_links) )
!allocate(local_alpha_f(1:num_necessary_faces) )
!allocate(local_beta_f(1:num_necessary_faces) )
!do s=1,num_necessary_sites
!  local_alpha_s(s) = alpha_s(global_site_of_local(s) )
!enddo
!
!do l=1,num_necessary_links
!  local_alpha_l(l) = alpha_l(global_link_of_local(l) )
!enddo
!
!do f=1,num_necessary_faces
!  local_alpha_f(f) = alpha_f(global_face_of_local(f) )
!  local_beta_f(f) = beta_f(global_face_of_local(f) )
!enddo
!
!
!end subroutine set_local_alpha_beta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! globalに定義していた諸々の数をlocalに置き換える
!subroutine switch_globalnum_to_localnum
!use parallel
!implicit none
!
!integer s,l,f,i,num
!
!global_num_sites = num_sites
!global_num_links = num_links
!global_num_faces = num_faces
!
!num_sites = num_local_sites
!num_links = num_local_links
!num_faces = num_local_faces
!
!deallocate( link_org, link_tip )
!allocate(link_org(1:num_necessary_links))
!allocate(link_tip(1:num_necessary_links))
!link_org = local_link_org
!link_tip = local_link_tip
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deallocate( linktip_from_s, linkorg_to_s )
!allocate(linktip_from_s(1:num_local_sites)) 
!allocate(linkorg_to_s(1:num_local_sites))
!do s=1,num_sites
!  num = local_linktip_from_s(s)%num_
!  linktip_from_s(s)%num_ = num
!  allocate(linktip_from_s(s)%labels_(1:num) )
!  allocate(linktip_from_s(s)%sites_(1:num) )
!  do i=1,num
!    linktip_from_s(s)%labels_(i) = local_linktip_from_s(s)%labels_(i)
!  enddo
!  !!!!!
!  num = local_linkorg_to_s(s)%num_
!  allocate(linkorg_to_s(s)%labels_(1:num) )
!  allocate(linkorg_to_s(s)%sites_(1:num) )
!  linkorg_to_s(s)%num_ = num
!  do i=1,num
!    linkorg_to_s(s)%labels_(i) = local_linkorg_to_s(s)%labels_(i)
!  enddo
!enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deallocate( links_in_f )
!allocate( links_in_f(1:num_necessary_faces) )
!do f=1,num_necessary_faces
!  num=local_links_in_f(f)%num_
!  links_in_f(f)%num_ = num
!  allocate( links_in_f(f)%link_labels_(1:num) )
!  allocate( links_in_f(f)%link_dirs_(1:num) )
!  links_in_f(f)%link_labels_=local_links_in_f(f)%link_labels_
!  links_in_f(f)%link_dirs_=local_links_in_f(f)%link_dirs_
!enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deallocate( face_in_l )
!allocate( face_in_l(1:num_links) )
!do l=1,num_links
!  num=local_face_in_l(l)%num_
!  face_in_l(l)%num_ = num
!  allocate( face_in_l(l)%label_(1:num) )
!  face_in_l(l)%label_=local_face_in_l(l)%label_
!enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deallocate( sites_in_f )
!allocate( sites_in_f(1:num_faces) )
!do f=1,num_faces
!  sites_in_f(f)%num_=1
!  allocate( sites_in_f(f)%label_(1:1) )
!  sites_in_f(f)%label_ = local_sites_in_f(f)%label_
!enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!deallocate( alpha_s )
!!allocate( alpha_s(1:num_necessary_sites) )
!!alpha_s = local_alpha_s
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!deallocate( alpha_l )
!!allocate( alpha_l(1:num_necessary_links) )
!!alpha_l = local_alpha_l
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!deallocate( alpha_f )
!!allocate( alpha_f(1:num_necessary_faces) )
!!alpha_f = local_alpha_f
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!deallocate( beta_f )
!!allocate( beta_f(1:num_necessary_faces) )
!!beta_f = local_beta_f
!
!end subroutine switch_globalnum_to_localnum
!
!#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!! local なラベルを出力する
!#ifdef PARALLEL
!subroutine check_local_sc
!use parallel
!implicit none
!
!integer :: rank, turn
!integer :: local,global
!
!turn=0
!if ( turn .ne. MYRANK ) then
!  call MPI_RECV(turn,1,MPI_INTEGER,MYRANK-1,MYRANK,MPI_COMM_WORLD,ISTATUS,IERR)
!endif
!write(*,*) "#####",MYRANK,"#####"
!
!write(*,*) "#!!!! local data to global data !!!!"
!write(*,*) "### site ###"
!do local=1,num_necessary_sites
!  write(*,*) "local site",local,"->","global",global_site_of_local(local)
!enddo
!write(*,*) "### link ###"
!do local=1,num_necessary_links
!  write(*,*) "local link",local,"->","global",global_link_of_local(local)
!enddo
!write(*,*) "### face ###"
!do local=1,num_necessary_faces
!  write(*,*) "local face",local,"->","global",global_face_of_local(local)
!enddo
!
!
!write(*,*) "#!!!! SEND and RECV !!!!"
!write(*,*) "### send site ###"
!do local=1,size(send_sites,1)
!  write(*,*) "send local site",send_sites(local)%label_,"to RANK",send_sites(local)%rank_
!enddo
!write(*,*) "### send link ###"
!do local=1,size(send_links,1)
!  write(*,*) "send local link",send_links(local)%label_,"to RANK",send_links(local)%rank_
!enddo
!write(*,*) "### send face ###"
!do local=1,size(send_faces,1)
!  write(*,*) "send local face",send_faces(local)%label_,"to RANK",send_faces(local)%rank_
!enddo
!write(*,*) "### recv site ###"
!do local=1,size(recv_sites,1)
!  write(*,*) "recv local site",recv_sites(local)%label_,"from RANK",recv_sites(local)%rank_
!enddo
!write(*,*) "### recv link ###"
!do local=1,size(recv_links,1)
!  write(*,*) "recv local link",recv_links(local)%label_,"from RANK",recv_links(local)%rank_
!enddo
!write(*,*) "### recv face ###"
!do local=1,size(recv_faces,1)
!  write(*,*) "recv local face",recv_faces(local)%label_,"from RANK",recv_faces(local)%rank_
!enddo
!
!write(*,*) "#!!!! SITE-LINK  !!!!"
!do local = 1,num_links
!  write(*,*) "local link",local,"is",link_org(local),link_tip(local)
!enddo
!
!write(*,*) "###"
!do local=1,num_sites
!  write(*,*) "local links from",local,";",linktip_from_s(local)%labels_
!enddo
!write(*,*) "###"
!do local=1,num_sites
!  write(*,*) "local links to",local,";",linkorg_to_s(local)%labels_
!enddo
!
!
!write(*,*)
!write(*,*)
!turn=MYRANK+1
!if( MYRANK .ne. NPROCS-1 ) then
!  call MPI_SEND(turn,1,MPI_INTEGER,MYRANK+1,MYRANK+1,MPI_COMM_WORLD,IERR)
!endif
!
!end subroutine check_local_sc
!#endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!! local なラベルを出力する
!#ifdef PARALLEL
!subroutine check_local_vals(PhiMat,UMat)
!use parallel
!implicit none
!
!complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
!
!integer :: rank, turn,i,j
!integer :: local,global
!complex(kind(0d0)) :: tmp
!
!turn=0
!if ( turn .ne. MYRANK ) then
!  call MPI_RECV(turn,1,MPI_INTEGER,MYRANK-1,MYRANK,MPI_COMM_WORLD,ISTATUS,IERR)
!endif
!write(*,*) "#####",MYRANK,"#####"
!
!write(*,*) "#!!!! local data to global data !!!!"
!write(*,*) "### Tr|Phi(s)|^2 ###"
!do local=1,num_necessary_sites
!  tmp=(0d0,0d0)
!  do i=1,NMAT
!    do j=1,NMAT
!      tmp=(PhiMat(i,j,local)*dconjg(PhiMat(i,j,local)))
!    enddo
!  enddo
!  write(*,*) "s=",global_site_of_local(local),":",dble(tmp)
!enddo
!write(*,*) "### sum_{ij} U(l)_{ij} ###"
!do local=1,num_necessary_links
!  tmp=(0d0,0d0)
!  do i=1,NMAT
!    do j=1,NMAT
!      tmp=tmp+UMat(i,j,local)
!    enddo
!  enddo
!  write(*,*) "l=",global_link_of_local(local),":",tmp
!enddo
!
!write(*,*)
!write(*,*)
!turn=MYRANK+1
!if( MYRANK .ne. NPROCS-1 ) then
!  call MPI_SEND(turn,1,MPI_INTEGER,MYRANK+1,MYRANK+1,MPI_COMM_WORLD,IERR)
!endif
!
!
!end subroutine check_local_vals
!#endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!! force の出力
!#ifdef PARALLEL
!subroutine check_force(force,C)
!use parallel
!implicit none
!
!integer, intent(in) :: C
!complex(kind(0d0)), intent(in) :: force(:,:,:)
!integer :: rank, turn,i,j,num
!integer :: local,global
!complex(kind(0d0)) :: tmp
!
!
!num=size(force,3)
!turn=0
!if ( turn .ne. MYRANK ) then
!  call MPI_RECV(turn,1,MPI_INTEGER,MYRANK-1,MYRANK,MPI_COMM_WORLD,ISTATUS,IERR)
!endif
!
!do local=1,num
!  tmp=(0d0,0d0)
!  do i=1,NMAT
!    do j=1,NMAT
!      tmp=(force(i,j,local)*dconjg(force(i,j,local)))
!    enddo
!  enddo
!  if(C==1) write(*,*) "s=",global_site_of_local(local),":",dble(tmp)
!  if(C==2) write(*,*) "l=",global_link_of_local(local),":",dble(tmp)
!  if(C==3) write(*,*) "f=",global_face_of_local(local),":",dble(tmp)
!enddo
!turn=MYRANK+1
!if( MYRANK .ne. NPROCS-1 ) then
!  call MPI_SEND(turn,1,MPI_INTEGER,MYRANK+1,MYRANK+1,MPI_COMM_WORLD,IERR)
!endif
!
!
!
!end subroutine check_force
!#endif
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!! テスト用に止めるルーチン
!subroutine stop_for_test
!#ifdef PARALLEL
!use parallel
!#endif
!implicit none
!
!#ifdef PARALLEL
!  call MPI_FINALIZE(IERR)
!#endif
!  stop
!end subroutine stop_for_test



end module global_parameters
