!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! code for computing observables
!
! 2016/03/04 added Pfaffian using pfapack
! v01: separate the mass term from PCSC
! v02: include test mode : usage: ./calcobs**.exe [med-file] [0:usual, 1:
program calcobs
use global_parameters
use initialization
use simplicial_complex
use observables
use Dirac_operator, only : make_Dirac
use matrix_functions, only : matrix_inverse, CalcPfaffian, matrix_eigenvalues, PfaffianLog
use SUN_generators, only : make_traceless_matrix_from_modes
use hamiltonian, only : bosonic_action_site, bosonic_action_link, bosonic_action_face
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! observables 
integer, parameter :: num_obs=19
double precision :: obs(1:num_obs)

character(30) :: obs_name(1:num_obs)
data obs_name/ &
"Sb",&                  ! obs(1)
"TrX2",&                ! obs(2)
"min eigen of DDdag",&  ! obs(3)
"max eigen of DDdag",&  ! obs(4)
"|Pfaffian|", &         ! obs(5)
"arg(Pf)" ,&            ! obs(6)
!"Re(PCSC base)" ,&       ! obs(7)
!"Im(PCSC base)" ,&       ! obs(8)
"Re(mu^2 part PCSC)" ,&       ! obs(7)
"Im(mu^2 part pCSC)" ,&       ! obs(8)
"|C_naive|", &                ! obs(9)
"arg(C_naive)", &             ! obs(10)
"|C_trace|", &                ! obs(11)
"arg(C_trace)", &             ! obs(12)
"|C_det|", &                ! obs(13)
"arg(C_det)", &             ! obs(14)
"|C_trace_rand|", &             ! obs(15)
"arg(C_trace_rand)", &                ! obs(16)
"arg(Pf \ PZ(Re(eigen)>0))", &      ! obs(17) PZ:(N^2-1)*euler/2
"arg(Pf \ PZ(Re(eigen)<0))", &      ! obs(18) 
"arg(Pf \ PZ(Re(eigen):random))" \   ! obs(19) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


integer :: ios, num_data, i, j, s, l, f
integer :: num, iargc

integer,parameter :: MED_FILE=10
integer :: OUT_FILE
character(128) :: medfile,outputfile,char_testmode,char_numcon
integer :: testmode, numcon

!! variables
complex(kind(0d0)), allocatable :: UMAT(:,:,:) ! unitary link variables
complex(kind(0d0)), allocatable :: PHI(:,:) ! complex scalar at sites
integer :: ite
complex(kind(0d0)), allocatable ::PhiMat(:,:) ! Phi in a matrix form
real(8) :: SB

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
  testmode=0
elseif( num == 2 ) then
  call getarg(2,char_testmode) 
  read(char_testmode,*) testmode
  numcon=1
elseif( num >= 3 ) then
  call getarg(2,char_testmode) 
  read(char_testmode,*) testmode
  call getarg(3,char_numcon) 
  read(char_numcon,*) numcon
endif

OUT_FILE=6

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
allocate( PHIMAT(1:NMAT, 1:NMAT) )
allocate( Dirac(1:sizeD,1:sizeD) )
allocate( Dirac_inv(1:sizeD,1:sizeD) )
allocate( eigenvalues(1:sizeD) )

CALL RANDOM_SEED()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if( OUT_FILE /= 6 ) then
  open(unit=OUT_FILE,status='replace',file=outputfile,action='write')
endif
if ( testmode .ne. 1 ) then
  write(OUT_FILE,'(a)') "# mu^2 PCSC = mu^2/2g^2 \Sigma*Tr(Phi \eta)"
  write(OUT_FILE,'(a)') '# C1=Tr(Phi^{-(N^2-1)*\chi/2),'
  write(OUT_FILE,'(a)') '# C2=Tr(Phi^2)^{-(N^2-1)*\chi/4) ),'
  write(OUT_FILE,'(a)') '# C3=det(Phi)^{-(N^2-1)*\chi/2N}'
  write(OUT_FILE,'(a)',advance='no') "# 1) iteration, " 
  do i=1,num_obs
    write(OUT_FILE,'(I3,a1,a)',advance='no') i+1,")",trim(obs_name(i))
  enddo
  write(OUT_FILE,*) 
endif

num_data=0
ios=0
do while (ios == 0)
  num_data=num_data+1
  read(MED_FILE,iostat=ios) ite
  read(MED_FILE,iostat=ios) UMAT
  read(MED_FILE,iostat=ios) PHI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! TEST MODE 
if ( testmode == 1 .and. num_data == numcon) then
  write(*,'(a)') "# site, i, j, Re(Phi_s(i,j)), Im(Phi_s(i,j))"
  do s=1,num_sites
    call make_traceless_matrix_from_modes(phimat(:,:), NMAT, Phi(:,s))
    do i=1,NMAT
      do j=1,NMAT
        write(*,'(3I2,X,E12.5,X,E12.5)') s,i,j,dble(phimat(i,j)),&
          dble( (0d0,-1d0)*phimat(i,j) )
      enddo
    enddo
  enddo
  write(*,'(a)') "# link, i, j, Re(U_l(i,j)), Im(U_l(i,j))"
  do l=1,num_links
    do i=1,NMAT
      do j=1,NMAT
        write(*,'(3I2,X,E12.5,X,E12.5)') l,i,j,dble(UMAT(i,j,l)),&
          dble( (0d0,-1d0)*UMAT(i,j,l) )
      enddo
    enddo
  enddo
! parameters
  write(*,'(a)') "#######"
  write(*,'(a,E12.5)') "# a = ", LatticeSpacing
  write(*,'(a,E12.5)') "# mu^2 = ", mass_square_phi
  write(*,'(a,E12.5)') "# 1/2g^2 = ", overall_factor
! bosonic action
  write(*,'(a)') "#######"
  call bosonic_action_site(SB,Phi)
  write(*,'(a)',advance='no') "# Sb_s = "
  write(*,'(E12.5)') SB
  call bosonic_action_link(SB,UMAT,Phi)
  write(*,'(a)',advance='no') "# Sb_l = "
  write(*,'(E12.5)') SB
  call bosonic_action_face(SB,UMAT)
  write(*,'(a)',advance='no') "# Sb_f = "
  write(*,'(E12.5)') SB
! mu^2 part of PCSC
  call make_Dirac(Dirac_inv,UMAT,Phi) ! Dirac operator
  call matrix_inverse(Dirac_inv)      ! compute inverse
  !!
  call calc_pcsc_site(pcsc_mass,Phi,Dirac_inv)
  write(*,'(a)',advance='no') "# \Sigma_s QS_\mu = "
  write(*,'(E12.5,X,E12.5)') pcsc_mass
  !!
  call calc_pcsc_link(pcsc_mass,Phi,UMAT,Dirac_inv)
  write(*,'(a)',advance='no') "# \Sigma_l QS_\mu = "
  write(*,'(E12.5,X,E12.5)') pcsc_mass
  !!
  call calc_pcsc_face(pcsc_mass,Phi,UMAT,Dirac_inv)
  write(*,'(a)',advance='no') "# \Sigma_f QS_\mu = "
  write(*,'(E12.5,X,E12.5)') pcsc_mass
  !!
  call calc_compensator_trace(compensator,Phi)
  write(*,'(a)',advance='no') "# A_trace = "
  write(*,'(E12.5,X,E12.5)') compensator 
  call calc_compensator_det(compensator,Phi)
  write(*,'(a)',advance='no') "# A_det = "
  write(*,'(E12.5,X,E12.5)') compensator

  stop
endif
  if( ios /= 0 ) exit

  if( testmode .ne. 1) then 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! prepare eigenvalues
  call make_Dirac(Dirac,UMAT,Phi)
  call matrix_eigenvalues(eigenvalues,Dirac)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compute observables 
! bosonic action
  call calc_bosonic_action(OBS(1),UMAT,Phi)
! TrX^2
  call calc_TrX2(OBS(2),Phi)
! min eigen
  OBS(3)=dble(eigenvalues(1)*dconjg(eigenvalues(1)))
! max eigen
  OBS(4)=dble(eigenvalues(sizeD)*dconjg(eigenvalues(sizeD)))
! Pfaffian
  !call make_Dirac(Dirac,UMAT,Phi)
  call make_Dirac(Dirac,UMAT,Phi)
  call CalcPfaffian(OBS(5), OBS(6), Dirac)
! mu^2 part of PCSC
  call make_Dirac(Dirac_inv,UMAT,Phi) ! Dirac operator
  call matrix_inverse(Dirac_inv)      ! compute inverse
  call calc_pcsc_mass(pcsc_mass,Phi,UMAT,Dirac_inv)
  OBS(7)=dble(pcsc_mass)
  OBS(8)=dble((0d0,-1d0)*pcsc_mass)
  call calc_compensator_naive(compensator,Phi)
  OBS(9)= abs(compensator)
  OBS(10)= arg(compensator)
  call calc_compensator_trace(compensator,Phi)
  OBS(11)= abs(compensator)
  OBS(12)= arg(compensator)
  call calc_compensator_det(compensator,Phi)
  OBS(13)= abs(compensator)
  OBS(14)= arg(compensator)
  call calc_compensator25(compensator,Phi)
  OBS(15)= abs(compensator)
  OBS(16)= arg(compensator)
! Phase of pseudo zeromodes
!COMMENT: we choose Re(eigen) > 0
  eigen_phase1=(1d0,0d0)
  eigen_phase2=(1d0,0d0)
  do j=1, (NMAT*NMAT-1)*abs(num_sites-num_links+num_faces)/2
    if( dble(eigenvalues(2*j-1)) > 0 ) then
      eigen_phase1 = eigen_phase1 * (0d0,1d0)*eigenvalues(2*j-1) &
        / cmplx( abs( eigenvalues(2*j-1) ) )
      eigen_phase2 = eigen_phase2 * (0d0,1d0)*eigenvalues(2*j) &
        / cmplx( abs( eigenvalues(2*j) ) )
    else
      eigen_phase1 = eigen_phase1 * (0d0,1d0)*eigenvalues(2*j) &
        / cmplx( abs( eigenvalues(2*j) ) )
      eigen_phase2 = eigen_phase2 * (0d0,1d0)*eigenvalues(2*j-1) &
        / cmplx( abs( eigenvalues(2*j-1) ) )
    endif
  enddo
  OBS(17) = arg( exp( (0d0,1d0)*cmplx(OBS(6)) ) / eigen_phase1 )
  OBS(18) = arg( exp( (0d0,1d0)*cmplx(OBS(6)) ) / eigen_phase2 )
    
  eigen_phase1=(1d0,0d0)
  do j=1, (NMAT*NMAT-1)*abs(num_sites-num_links+num_faces)/2 
    if( dble(eigenvalues(2*j-1)) > 0 ) then 
        ctmp=eigenvalues(2*j-1)
    else
        ctmp=eigenvalues(2*j)
    endif
    call random_number(rand)
    if( rand < 0.5d0 ) then
      eigen_phase1 = eigen_phase1 * (0d0,1d0)*ctmp &
        / cmplx( abs( ctmp ) )
    else
      eigen_phase1 = eigen_phase1 * (-(0d0,1d0)*ctmp) &
        / cmplx( abs( ctmp ) )
    endif
  enddo
  OBS(19) = arg( exp( (0d0,1d0)*cmplx(OBS(6)) ) / eigen_phase1 )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  write(OUT_FILE,'(I6,2X)',advance='no') ite
  do i=1,num_obs
    write(OUT_FILE,'(E12.5,2X)',advance='no') OBS(i)
  enddo
  write(OUT_FILE,*) 

  endif
enddo

num_data=num_data-1
write(OUT_FILE,'(a,I10)') "# total number= ",num_data
if ( ios /= -1 ) then 
  write(OUT_FILE,'(a,a,a)') "# med file ",trim(medfile)," ends abnormally."
  write(OUT_FILE,'(a,I8)') "#   error code = ", ios
endif

end program calcobs


