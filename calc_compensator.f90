!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! code for computing observables
!
! 2016/03/04 added Pfaffian using pfapack
! v01: separate the mass term from PCSC
program calccomp
use global_parameters
use initialization
use simplicial_complex
use observables
use Dirac_operator, only : make_Dirac
use matrix_functions, only : matrix_inverse, CalcPfaffian, matrix_eigenvalues, PfaffianLog
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! observables 
integer, parameter :: num_obs=10
double precision :: obs(1:num_obs)

character(30) :: obs_name(1:num_obs)
data obs_name/ &
"|C_naive|", &    ! obs(1)
"arg(C_naive)", &  ! obs(2)
"|C_trace|", &     ! obs(3)
"arg(C_trace)", &  ! obs(4)
"|C_det|", &      ! obs(5)
"arg(C_det)", &      ! obs(6)
"|C2_trace|", &   ! obs(7) [ sum tr(\Phi^2) ]^[...]
"arg(C2_trace)", & ! obs(8)
"|C2_det|", &     ! obs(9) [ sum det(\Phi) ]^[...]
"arg(C2_det)"/    ! obs(10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


integer :: ios, num_data, i, s, j
integer :: num, iargc

integer,parameter :: MED_FILE=10
integer :: OUT_FILE
character(128) :: medfile,outputfile

!! variables
complex(kind(0d0)), allocatable :: UMAT(:,:,:) ! unitary link variables
complex(kind(0d0)), allocatable :: PHI(:,:) ! complex scalar at sites
integer :: ite

complex(kind(0d0)), allocatable :: Dirac(:,:) 
complex(kind(0d0)), allocatable :: Dirac_inv(:,:) 
complex(kind(0d0)), allocatable :: eigenvalues(:)

!! for anomaly-phase quench
complex(kind(0d0)) :: compensator,pcsc_mass

!! for pseudo zeromodes
complex(kind(0d0)) :: eigen_phase1, eigen_phase2, ctmp
real(8) :: rand

!! for another pfaffian
real(8) logPfaffian
complex(kind(0d0)) :: phasePfaffian

num = iargc()
call getarg(1,medfile)
if( num == 1 ) then 
  OUT_FILE=6
else
  OUT_FILE=11
  call getarg(2,outputfile)
endif

open(unit=MED_FILE,file=medfile,action='read',form='unformatted')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! initialization
read(MED_FILE) NMAT
read(MED_FILE) LatticeSpacing
read(MED_FILE) SC_FILE_NAME
read(MED_FILE) ALPHA_BETA
read(MED_FILE) test_mode
read(MED_FILE) new_config
read(MED_FILE) fix_seed
read(MED_FILE) read_alpha
read(MED_FILE) config_step
read(MED_FILE) obs_step
read(MED_FILE) m_omega
read(MED_FILE) mass_square_phi
read(MED_FILE) mass_f
read(MED_FILE) Remez_factor4
read(MED_FILE) Remez_factor8
read(MED_FILE) epsilon
read(MED_FILE) CG_max
read(MED_FILE) num_ite
read(MED_FILE) Ntau
read(MED_FILE) Dtau
read(MED_FILE) R_phi
read(MED_FILE) R_A
read(MED_FILE) Fconfigin
read(MED_FILE) Foutput
read(MED_FILE) Fconfigout
read(MED_FILE) Fmedconf
dimG=NMAT*NMAT-1
Dtau_phi = R_phi * Dtau
Dtau_A = R_A * Dtau
!one_ov_2g2N=1d0/(2d0*LatticeSpacing*LatticeSpacing)
overall_factor=dble(NMAT)/(2d0*LatticeSpacing*LatticeSpacing)

call set_NZF
call set_sc
call set_alpha_beta
allocate( UMAT(1:NMAT,1:NMAT,1:num_links) )
allocate( PHI(1:dimG, 1:num_sites) )
allocate( Dirac(1:sizeD,1:sizeD) )
allocate( Dirac_inv(1:sizeD,1:sizeD) )
allocate( eigenvalues(1:sizeD) )

CALL RANDOM_SEED()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if( OUT_FILE /= 6 ) then
  open(unit=OUT_FILE,status='replace',file=outputfile,action='write')
endif
write(OUT_FILE,'(a)',advance='no') "# 1) iteration, " 
do i=1,num_obs
  write(OUT_FILE,'(I3,a1,a)',advance='no') i+1,")",trim(obs_name(i))
enddo
write(OUT_FILE,*) 

num_data=0
ios=0
do while (ios == 0)
  num_data=num_data+1
  read(MED_FILE,iostat=ios) ite
  read(MED_FILE,iostat=ios) UMAT
  read(MED_FILE,iostat=ios) PHI

  if( ios /= 0 ) exit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! prepare eigenvalues


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compute observables 
  call calc_compensator_naive(compensator,Phi)
  OBS(1)= abs(compensator)
  OBS(2)= arg(compensator)
  call calc_compensator_trace(compensator,Phi)
  OBS(3)= abs(compensator)
  OBS(4)= arg(compensator)
  call calc_compensator_det(compensator,Phi)
  OBS(5)= abs(compensator)
  OBS(6)= arg(compensator)
  call calc_compensator2_trace(compensator,Phi)
  OBS(7)= abs(compensator)
  OBS(8)= arg(compensator)
  call calc_compensator2_det(compensator,Phi)
  OBS(9)= abs(compensator)
  OBS(10)= arg(compensator)


  write(OUT_FILE,'(I6,2X)',advance='no') ite
  do i=1,num_obs
    write(OUT_FILE,'(E12.5,2X)',advance='no') OBS(i)
  enddo
  write(OUT_FILE,*) 
enddo

num_data=num_data-1
write(OUT_FILE,'(a,I10)') "# total number= ",num_data
if ( ios /= -1 ) then 
  write(OUT_FILE,'(a,a,a)') "# med file ",trim(medfile)," ends abnormally."
  write(OUT_FILE,'(a,I8)') "#   error code = ", ios
endif

end program calccomp


