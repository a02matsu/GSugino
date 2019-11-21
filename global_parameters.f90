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
  integer, allocatable :: labels_(:) ! link label
  integer, allocatable :: sites_(:)  ! site label
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

type A_in_B
  integer :: num_ ! number of elements
  integer, allocatable :: label_(:) ! global link label
  complex(kind(0d0)), allocatable :: val_(:) ! 1~num_
end type A_in_B



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
integer, parameter :: p_mass=0

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

#ifdef COUNT_TIME
real :: t_site0,t_site1,t_site
real :: t_link0,t_link1,t_link
real :: t_face0,t_face1,t_face
#endif




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

type(FACEDAT), allocatable, save :: global_face_in_s(:) ! faces touching at site




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
!! 
integer, allocatable, save :: num_faces_in_s(:) ! number of faces touching on s

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
integer, save :: force_measurement !!  1:measure forces
integer, save :: eigen_measurement !!  1:measure min and max of DDdag
integer, save :: fix_seed !! 0: use the value at the end of the previous simulation
                    !! 1: fix seed to the value in the parameter file
                    !! 2: determine by the system time
integer, save :: new_config !! 0: Old config, 1: New config
integer, save :: branch_use !! in which branch simulation is performed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! branch
integer, save :: branch_mode !! 0: normal, 1: make branch
integer, save :: branch_root !! make branch from this config. default:0
integer, save :: branch_num !! number of branches to make (active iff branch_mode=1)


!integer :: read_alpha !! 0: alpha,beta=1 1: read ALPHA_BETA_FILE
integer, save :: save_med_step !! steps to write out configurations
integer, save :: obs_step !! steps to compute observables
integer, save :: reset_ite !! 0:iteration number continue, 1:iteration number reset 
integer, save :: save_config_step !! save config by this step. 
integer, save :: eval_eigen !! 0:SKIP the calculation of eigenvalues of Dirac


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
character(128), save :: Remez_1ovminus4
character(128), save :: Remez_1ov8
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
character(128), save :: INPUT_FILE_NAME !="inputfile" ! input file
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

end module global_parameters
