!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! code for computing observables
!
! 2016/03/04 added Pfaffian using pfapack
program calcobs
use global_parameters
use initialization
use simplicial_complex
use observables
use Dirac_operator, only : make_Dirac
use matrix_functions, only : matrix_inverse, CalcPfaffian, matrix_eigenvalues, PfaffianLog
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! observables 
integer, parameter :: num_obs=24
double precision :: obs(1:num_obs)
character(30) :: obs_name(1:num_obs)
data obs_name/ &
"Sb",&                  ! obs(1)
"TrX2",&                ! obs(2)
"min eigen of DDdag",&  ! obs(3)
"max eigen of DDdag",&  ! obs(4)
"|Pfaffian|", &         ! obs(5)
"arg(Pf)" ,&            ! obs(6)
"Re(PCSC base)" ,&       ! obs(7)
"Im(PCSC base)" ,&       ! obs(8)
"|C_naive|", &                ! obs(9)
"arg(C_naive)", &             ! obs(10)
"|C_trace|", &                ! obs(11)
"arg(C_trace)", &             ! obs(12)
"|C_det|", &                ! obs(13)
"arg(C_det)", &             ! obs(14)
"arg(Pf*C_naive)", &              ! obs(15)
"arg(Pf*C_trace)", &              ! obs(16)
"arg(Pf*C_det)", &               ! obs(17)
"|C_trace_rand|", &             ! obs(18)
"arg(C_trace_rand)", &                ! obs(19)
"arg(Pf*C_trace_rand)", &                ! obs(20)
"arg(Pf \ PZ(Re(eigen)>0))", &      ! obs(21) PZ:(N^2-1)*euler/2
"arg(Pf \ PZ(Re(eigen)<0))", &      ! obs(22) 
"arg(Pf \ PZ(Re(eigen):random)" , &   ! obs(23) 
"arg(Pf \ PZ2(Re(eigen)>0)" /   ! obs(24) PZ2:(N^2-1)*euler
!"arg(prod_eigen Re>0)" /   ! obs(25) PZ2:(N^2-1)*euler
!"|Pfaffian|(org)", &         ! obs(24)
!"arg(Pf)(org)" /            ! obs(25)
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
complex(kind(0d0)) :: compensator,pcsc

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
write(OUT_FILE,'(a)') "# PCSC base = Sb - mu^2/2g^2 \Sigma*Tr(Phi \eta)"
write(OUT_FILE,'(a)') '# C1=Tr(Phi^{-(N^2-1)*\chi/2),'
write(OUT_FILE,'(a)') '# C2=Tr(Phi^2)^{-(N^2-1)*\chi/4) ),'
write(OUT_FILE,'(a)') '# C3=det(Phi)^{-(N^2-1)*\chi/2N}'
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
! PCSC
  call make_Dirac(Dirac_inv,UMAT,Phi) ! Dirac operator
  call matrix_inverse(Dirac_inv)      ! compute inverse
  call calc_pcsc(pcsc,Phi,UMAT,Dirac_inv)
  OBS(7)=dble(pcsc)
  OBS(8)=dble((0d0,-1d0)*pcsc)
  call calc_compensator_naive(compensator,Phi)
  OBS(9)= abs(compensator)
  OBS(10)= arg(compensator)
  OBS(15)= arg( exp( (0d0,1d0)*OBS(6) ) * compensator )
  call calc_compensator_trace(compensator,Phi)
  OBS(11)= abs(compensator)
  OBS(12)= arg(compensator)
  OBS(16)= arg( exp( (0d0,1d0)*OBS(6) ) * compensator )
  call calc_compensator_det(compensator,Phi)
  OBS(13)= abs(compensator)
  OBS(14)= arg(compensator)
  OBS(17)= arg( exp( (0d0,1d0)*OBS(6) ) * compensator )
  call calc_compensator25(compensator,Phi)
  OBS(18)= abs(compensator)
  OBS(19)= arg(compensator)
  OBS(20)= arg( exp( (0d0,1d0)*OBS(6) ) * compensator )
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
  OBS(21) = arg( exp( (0d0,1d0)*cmplx(OBS(6)) ) / eigen_phase1 )
  OBS(22) = arg( exp( (0d0,1d0)*cmplx(OBS(6)) ) / eigen_phase2 )
    
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
  OBS(23) = arg( exp( (0d0,1d0)*cmplx(OBS(6)) ) / eigen_phase1 )

  eigen_phase1=(1d0,0d0)
  do j=1, (NMAT*NMAT-1)*abs(num_sites-num_links+num_faces)
    if( dble(eigenvalues(2*j-1)) > 0 ) then
      eigen_phase1 = eigen_phase1 * (0d0,1d0)*eigenvalues(2*j-1) &
        / cmplx( abs( eigenvalues(2*j-1) ) )
    else
      eigen_phase1 = eigen_phase1 * (0d0,1d0)*eigenvalues(2*j) &
        / cmplx( abs( eigenvalues(2*j) ) )
    endif
  enddo
  OBS(24) = arg( exp( (0d0,1d0)*cmplx(OBS(6)) ) / eigen_phase1 )
  !call make_Dirac(Dirac,UMAT,Phi)
  !call PfaffianLog(Dirac,logPfaffian,phasePfaffian)
  !OBS(24)=exp( logPfaffian )
  !OBS(25)=arg(phasePfaffian)

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

end program calcobs


