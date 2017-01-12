!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! code for computing observables
!
! 2016/03/04 added Pfaffian using pfapack
! v01: separate the mass term from PCSC
! v02: include test mode : usage: ./calcobs**.exe [med-file] [0:usual, 1:
! v03: include deviation of WT in naive quench and partial Sf terms
! v04: include test-mode2: 
! v05: include SU(2) compensator and IZ compensator
! v06: include subtracted Pfaffian where PZ are identified as 
!      the smallest (dimG)*|Euler| eigen values with 
!        |Re(eigen)| < ¥epsilon and |Im(eigen)| < ¥epsilon
!      or
!        |Im(eigen)| > ¥epsilon
! usage: 
!   ./calcobs**.exe [med-file] 0 
!      :: usual
!   ./calcobs**.exe [med-file] 1 [config num] 
!      :: output configuration and observable for a specific configuration
!   ./calcobs**.exe [med-file] 2 
!      :: output U_f U_f^¥dagger - 1 
! v07: swhich from phi to phimat
program calcobs
use global_parameters
use initialization
use simplicial_complex
use observables
use Dirac_operator, only : make_Dirac
use matrix_functions, only : matrix_inverse, CalcPfaffian, matrix_eigenvalues, PfaffianLog, matrix_product
use SUN_generators, only : make_traceless_matrix_from_modes, trace_MTa
use hamiltonian, only : bosonic_action_site, bosonic_action_link, bosonic_action_face
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! observables 
integer, parameter :: num_obs=41
double precision :: obs(1:num_obs)

character(100) :: obs_name(1:num_obs)
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
"arg(Pf \ PZ(Re(eigen):random))", &   ! obs(19) 
"Re( 1/2 (QAtr*)/|Atr|^2 Sigma )", &  ! obs(20) 
"Im( 1/2 (QAtr*)/|Atr|^2 Sigma )", &  ! obs(21) 
"Re( Sf site )", & ! obs(22)
"Im( Sf site )", & ! obs(23)
"Re( Sf link )", & ! obs(24)
"Im( Sf link )", & ! obs(25)
"Re( Sf face )", & ! obs(26)
"Im( Sf face )", & ! obs(27)
"Re( PCSC_mass site )", & ! obs(28)
"Im( PCSC_mass site )", & ! obs(29)
"Re( PCSC_mass link )", & ! obs(30)
"Im( PCSC_mass link )", & ! obs(31)
"Re( PCSC_mass face )", & ! obs(32)
"Im( PCSC_mass face )", &! obs(33)
"|Ctr_SU2| (ONLY SU2)", &   ! obs(34)
"arg(Ctr_SU2) (ONLY SU2)", &! obs(35)
"|C_IZ1|", &   ! obs(36)
"arg(C_IZ1) ",&! obs(37)
"|C_IZ2|", &   ! obs(38)
"arg(C_IZ2)", &  ! obs(39)
"arg(Pf \ PZ(not real and origin, Re>0) )", &  ! obs(40)
"arg(Pf \ PZ(not real and origin, Re:random) )" /  ! obs(41)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


integer :: ios, num_data, i, j, s, l, f,a
integer :: num, iargc

integer,parameter :: MED_FILE=10
integer :: OUT_FILE
character(128) :: medfile,outputfile,char_testmode,char_numcon
integer :: testmode, numcon

!! variables
complex(kind(0d0)), allocatable :: UMAT(:,:,:) ! unitary link variables
complex(kind(0d0)), allocatable :: PHI(:,:) ! complex scalar at sites
integer :: ite
complex(kind(0d0)), allocatable ::PhiMat(:,:,:) ! Phi in a matrix form
real(8) :: SB

complex(kind(0d0)), allocatable :: Dirac(:,:) 
complex(kind(0d0)), allocatable :: Dirac_inv(:,:) 
complex(kind(0d0)), allocatable :: eigenvalues(:)

!! for anomaly-phase quench
complex(kind(0d0)) :: compensator,pcsc_mass,Sf_site, Sf_link, Sf_face

!! for pseudo zeromodes
complex(kind(0d0)) :: eigen_phase1, eigen_phase2, ctmp
real(8) :: rand

!! for another pfaffian
real(8) logPfaffian
complex(kind(0d0)) :: phasePfaffian

!! general observable
complex(kind(0d00)) :: observable

!! for testmode 2
complex(kind(0d0)), allocatable :: MAT(:,:)
real(8) :: rtmp

!! for complex pseudo zero modes
real :: minimal
integer :: numPZ


num = iargc()
call getarg(1,medfile)

if( num == 1 ) then 
  testmode=0
elseif( num == 2 ) then
  call getarg(2,char_testmode) 
  read(char_testmode,*) testmode
  if( testmode == 1 ) then
    numcon=1
  endif
elseif( num >= 3 ) then
  call getarg(2,char_testmode) 
  read(char_testmode,*) testmode
  if ( testmode == 1 ) then
    call getarg(3,char_numcon) 
    read(char_numcon,*) numcon
  endif
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
read(MED_FILE) save_med_step
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
allocate( PHIMAT(1:NMAT, 1:NMAT,1:num_sites) )
allocate( Dirac(1:sizeD,1:sizeD) )
allocate( Dirac_inv(1:sizeD,1:sizeD) )
allocate( eigenvalues(1:sizeD) )
allocate( MAT(1:NMAT,1:NMAT) )

CALL RANDOM_SEED()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if( OUT_FILE /= 6 ) then
  open(unit=OUT_FILE,status='replace',file=outputfile,action='write')
endif
if ( testmode == 0 ) then
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
  read(MED_FILE,iostat=ios) PHIMAT

  do s=1,num_sites
    do a=1,dimG
      call trace_MTa(Phi(a,s),PhiMat(:,:,s),a,NMAT)
    enddo
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! TEST MODE 1
if ( testmode == 1 .and. num_data == numcon) then
  write(*,'(a)') "# site, i, j, Re(Phi_s(i,j)), Im(Phi_s(i,j))"
  do s=1,num_sites
    !call make_traceless_matrix_from_modes(phimat(:,:,s), NMAT, Phi(:,s))
    do i=1,NMAT
      do j=1,NMAT
        write(*,'(3I2,X,E23.15,X,E23.15)') s,i,j,dble(phimat(i,j,s)),&
          dble( (0d0,-1d0)*phimat(i,j,s) )
      enddo
    enddo
  enddo
  write(*,'(a)') "# link, i, j, Re(U_l(i,j)), Im(U_l(i,j))"
  do l=1,num_links
    do i=1,NMAT
      do j=1,NMAT
        write(*,'(3I2,X,E23.15,X,E23.15)') l,i,j,dble(UMAT(i,j,l)),&
          dble( (0d0,-1d0)*UMAT(i,j,l) )
      enddo
    enddo
  enddo
! parameters
  write(*,'(a)') "#######"
  write(*,'(a,E23.15)') "# a = ", LatticeSpacing
  write(*,'(a,E23.15)') "# mu^2 = ", mass_square_phi
  write(*,'(a,E23.15)') "# 1/2g^2 = ", overall_factor
! bosonic action
  write(*,'(a)') "#######"
  call bosonic_action_site(SB,PhiMat)
  write(*,'(a)',advance='no') "# Sb_s = "
  write(*,'(E23.15)') SB
  call bosonic_action_link(SB,UMAT,PhiMat)
  write(*,'(a)',advance='no') "# Sb_l = "
  write(*,'(E23.15)') SB
  call bosonic_action_face(SB,UMAT)
  write(*,'(a)',advance='no') "# Sb_f = "
  write(*,'(E23.15)') SB
  call calc_bosonic_action(SB,UMAT,PhiMat)
  write(*,'(a)',advance='no') "### Sb_all = "
  write(*,'(E23.15)') SB
! mu^2 part of PCSC
  call make_Dirac(Dirac_inv,UMAT,PhiMat) ! Dirac operator
  call matrix_inverse(Dirac_inv)      ! compute inverse
  !!
  call calc_pcsc_site(pcsc_mass,Phi,Dirac_inv)
  write(*,'(a)',advance='no') "# \Sigma_s QS_\mu = "
  write(*,'(E23.15,X,E23.15)') pcsc_mass
  !!
  call calc_pcsc_link(pcsc_mass,Phi,UMAT,Dirac_inv)
  write(*,'(a)',advance='no') "# \Sigma_l QS_\mu = "
  write(*,'(E23.15,X,E23.15)') pcsc_mass
  !!
  call calc_pcsc_face(pcsc_mass,Phi,UMAT,Dirac_inv)
  write(*,'(a)',advance='no') "# \Sigma_f QS_\mu = "
  write(*,'(E23.15,X,E23.15)') pcsc_mass
  !!
  call calc_pcsc_mass(pcsc_mass,Phi,UMAT,Dirac_inv)
  write(*,'(a)',advance='no') "### \Sigma_all QS_\mu = "
  write(*,'(E23.15,X,E23.15)') pcsc_mass
  !!
  call calc_compensator_trace(compensator,Phi)
  write(*,'(a)',advance='no') "# A_trace = "
  write(*,'(E23.15,X,E23.15)') compensator 
  call calc_compensator_det(compensator,Phi)
  write(*,'(a)',advance='no') "# A_det = "
  write(*,'(E23.15,X,E23.15)') compensator

  stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! TEST MODE 2
elseif( testmode == 2 ) then
  write(*,'(I7)',advance='no') ite
  do l=1, num_links
    rtmp=0d0
    call matrix_product(MAT,UMAT(:,:,l),UMAT(:,:,l),'N','C')
    do i=1,NMAT
      do j=1,NMAT
        if( i == j ) then
          rtmp = rtmp + dble( (MAT(i,i)-(1d0,0d0))*conjg( MAT(i,i) - (1d0,0d0) ))
        else
          rtmp = rtmp + dble( MAT(i,j)*conjg( MAT(i,j) ) )
        endif
      enddo
    enddo
    write(*,'(E12.5)',advance='no') rtmp
  enddo
  write(*,*)
endif
  if( ios /= 0 ) exit

  if( testmode == 0) then 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! prepare eigenvalues
  call make_Dirac(Dirac,UMAT,PhiMat)
  call matrix_eigenvalues(eigenvalues,Dirac)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compute observables 
! bosonic action
  call calc_bosonic_action(OBS(1),UMAT,PhiMat)
! TrX^2
  call calc_TrX2(OBS(2),PhiMat)
! min eigen
  OBS(3)=dble(eigenvalues(1)*dconjg(eigenvalues(1)))
! max eigen
  OBS(4)=dble(eigenvalues(sizeD)*dconjg(eigenvalues(sizeD)))
! Pfaffian
  !call make_Dirac(Dirac,UMAT,Phi)
  call make_Dirac(Dirac,UMAT,PhiMat)
  call CalcPfaffian(OBS(5), OBS(6), Dirac)
! mu^2 part of PCSC
  call make_Dirac(Dirac_inv,UMAT,PhiMat) ! Dirac operator
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
  call QbarAtr_sigma(observable,Phi,UMAT,Dirac_inv)
  call calc_compensator_trace(compensator,Phi)
  OBS(20) = dble( observable &
    / ( (2d0,0d0) * conjg( compensator ) ))
  OBS(21) = dble( (0d0,-1d0) * observable &
    / ( (2d0,0d0) * conjg( compensator ) ))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call make_Dirac(Dirac_inv,UMAT,PhiMat) ! Dirac operator
  call matrix_inverse(Dirac_inv)      ! compute inverse
  call partial_Sf(Sf_site,Sf_link,Sf_face,Dirac,Dirac_inv)
  OBS(22) = dble(Sf_site)
  OBS(23) = dble( (0d0,-1d0)*Sf_site )
  !!
  OBS(24) = dble(Sf_link)
  OBS(25) = dble( (0d0,-1d0)*Sf_link )
  !!
  OBS(26) = dble(Sf_face)
  OBS(27) = dble( (0d0,-1d0)*Sf_face )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call calc_pcsc_site(pcsc_mass,Phi,Dirac_inv)
  OBS(28) = dble(pcsc_mass)
  OBS(29) = dble((0d0,-1d0)*pcsc_mass)
  !!
  call calc_pcsc_link(pcsc_mass,Phi,UMAT,Dirac_inv)
  OBS(30) = dble(pcsc_mass)
  OBS(31) = dble((0d0,-1d0)*pcsc_mass)
  !!
  call calc_pcsc_face(pcsc_mass,Phi,UMAT,Dirac_inv)
  OBS(32) = dble(pcsc_mass)
  OBS(33) = dble((0d0,-1d0)*pcsc_mass)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! SU(2) compensator
  call calc_compensator_trace_SU2(compensator,Phi)
  OBS(34)= abs(compensator)
  OBS(35)= arg(compensator)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Izikson-Zuber compensators
  call calc_compensator_IZ1(compensator,Phi,UMAT,Dirac_inv)
  OBS(36)= abs(compensator)
  OBS(37)= arg(compensator)
  call calc_compensator_IZ2(compensator,Phi,UMAT,Dirac_inv)
  OBS(38)= abs(compensator)
  OBS(39)= arg(compensator)

! Phase of pseudo zeromodes 
! PZ: smallest dimG*Euler eigenvalues except on the real axis. 
!COMMENT: we choose Re(eigen) > 0
  minimal=abs( eigenvalues(1) )*10d0
  eigen_phase1=(1d0,0d0)
  numPZ=0 
  j=0
  do while ( numPZ < abs((NMAT*NMAT-1)*abs(num_sites-num_links+num_faces)/2) )
    j=j+1
    if( dble(eigenvalues(2*j-1)) > 0 ) then
      if( (dble(eigenvalues(2*j-1)) > minimal .and. &
           abs(dble( (0d0,1d0)*eigenvalues(2*j-1) )) > minimal ) &
          .or. &
          (dble(eigenvalues(2*j-1)) < minimal) ) then
        eigen_phase1 = eigen_phase1 * (0d0,1d0)*eigenvalues(2*j-1) &
          / cmplx( abs( eigenvalues(2*j-1) ) )
        numPZ=numPZ+1
      endif
    else
      if( (dble(eigenvalues(2*j)) > minimal .and. &
           abs(dble( (0d0,1d0)*eigenvalues(2*j) )) > minimal ) &
          .or. &
          (dble(eigenvalues(2*j)) < minimal) ) then
        eigen_phase1 = eigen_phase1 * (0d0,1d0)*eigenvalues(2*j) &
          / cmplx( abs( eigenvalues(2*j) ) )
        numPZ=numPZ+1
      endif
    endif
  enddo
  OBS(40) = arg( exp( (0d0,1d0)*cmplx(OBS(6)) ) / eigen_phase1 )
! choose the sigunature of Re(PZ) randomely
  eigen_phase1=(1d0,0d0)
  numPZ=0 
  j=0
  do while ( numPZ < abs((NMAT*NMAT-1)*abs(num_sites-num_links+num_faces)/2) )
    j=j+1
    call random_number(rand)
    if( rand >= 0.5d0 ) then
      if( (dble(eigenvalues(2*j-1)) > minimal .and. &
           abs(dble( (0d0,1d0)*eigenvalues(2*j-1) )) > minimal ) &
          .or. &
          (dble(eigenvalues(2*j-1)) < minimal) ) then
        eigen_phase1 = eigen_phase1 * (0d0,1d0)*eigenvalues(2*j-1) &
          / cmplx( abs( eigenvalues(2*j-1) ) )
        numPZ=numPZ+1
      endif
    else
      if( (dble(eigenvalues(2*j)) > minimal .and. &
           abs(dble( (0d0,1d0)*eigenvalues(2*j) )) > minimal ) &
          .or. &
          (dble(eigenvalues(2*j)) < minimal) ) then
        eigen_phase1 = eigen_phase1 * (0d0,1d0)*eigenvalues(2*j) &
          / cmplx( abs( eigenvalues(2*j) ) )
        numPZ=numPZ+1
      endif
    endif
  enddo
  OBS(41) = arg( exp( (0d0,1d0)*cmplx(OBS(6)) ) / eigen_phase1 )

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


