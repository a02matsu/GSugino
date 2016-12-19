!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! code for computing observables
!
! 2016/03/04 added Pfaffian using pfapack
program calceigen
use global_parameters
use initialization
use simplicial_complex
use observables
use Dirac_operator, only : make_Dirac
use matrix_functions, only : matrix_inverse, CalcPfaffian, matrix_eigenvalues, PfaffianLog
implicit none


integer :: ios, num_data, i, s, j
integer :: num, iargc

integer,parameter :: MED_FILE=10
!integer,parameter :: EIGEN_FILE1=11
!integer,parameter :: EIGEN_FILE2=12
!integer,parameter :: EIGEN_FILE3=13
!integer,parameter :: EIGEN_FILE4=14
character(128) :: medfile! outputfile1,outputfile2,outputfile3,outputfile4

!! variables
complex(kind(0d0)), allocatable :: UMAT(:,:,:) ! unitary link variables
complex(kind(0d0)), allocatable :: PHI(:,:) ! complex scalar at sites
integer :: ite

complex(kind(0d0)), allocatable :: Dirac(:,:) 
complex(kind(0d0)), allocatable :: eigenvalues(:)
real(8) :: re_eig, im_eig

num = iargc()
call getarg(1,medfile)

!outputfile1="eigen_pp.dat"
!outputfile2="eigen_np.dat"
!outputfile3="eigen_nn.dat"
!outputfile4="eigen_pn.dat"
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
allocate( eigenvalues(1:sizeD) )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!open(unit=EIGEN_FILE1,status='replace',file=outputfile1,action='write')
!open(unit=EIGEN_FILE2,status='replace',file=outputfile2,action='write')
!open(unit=EIGEN_FILE3,status='replace',file=outputfile3,action='write')
!open(unit=EIGEN_FILE4,status='replace',file=outputfile4,action='write')

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

  write(*,'(I6,2X)',advance='no') ite
  !write(EIGEN_FILE1,'(I6,2X)',advance='no') ite
  !write(EIGEN_FILE2,'(I6,2X)',advance='no') ite
  !write(EIGEN_FILE3,'(I6,2X)',advance='no') ite
  !write(EIGEN_FILE4,'(I6,2X)',advance='no') ite
  do j=1,sizeD
    re_eig=dble(eigenvalues(j))
    im_eig=dble((0d0,-1d0)*eigenvalues(j))
      write(*,'(E12.5,2X)',advance='no') re_eig
      write(*,'(E12.5,2X)',advance='no') im_eig
    !if ( re_eig >= 0d0 .and. im_eig >= 0d0 ) then
      !write(EIGEN_FILE1,'(E12.5,2X)',advance='no') re_eig
      !write(EIGEN_FILE1,'(E12.5,2X)',advance='no') im_eig
    !elseif( re_eig < 0d0 .and. im_eig >= 0d0 ) then
      !write(EIGEN_FILE2,'(E12.5,2X)',advance='no') re_eig
      !write(EIGEN_FILE2,'(E12.5,2X)',advance='no') im_eig
    !elseif( re_eig < 0d0 .and. im_eig < 0d0 ) then
      !write(EIGEN_FILE3,'(E12.5,2X)',advance='no') re_eig
      !write(EIGEN_FILE3,'(E12.5,2X)',advance='no') im_eig
    !elseif( re_eig >= 0d0 .and. im_eig < 0d0 ) then
      !write(EIGEN_FILE4,'(E12.5,2X)',advance='no') re_eig
      !write(EIGEN_FILE4,'(E12.5,2X)',advance='no') im_eig
    !endif
  enddo
    write(*,*) 
    write(*,*) 
    write(*,*) 
    !write(EIGEN_FILE1,*) 
    !write(EIGEN_FILE1,*) 
    !write(EIGEN_FILE1,*) 
    !write(EIGEN_FILE2,*) 
    !write(EIGEN_FILE2,*) 
    !write(EIGEN_FILE2,*) 
    !write(EIGEN_FILE3,*) 
    !write(EIGEN_FILE3,*) 
    !write(EIGEN_FILE3,*) 
    !write(EIGEN_FILE4,*) 
    !write(EIGEN_FILE4,*) 
    !write(EIGEN_FILE4,*) 

enddo

end program calceigen
