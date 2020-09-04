!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! use as
!! % ./a.out [MEDFILE] [Ckind]
!!   Ckind=1 : A = sum_s ( Tr( (\phi + \epsilon)^2 )^{-r/4} )
!!   Ckind=2 : A = 1 / sum_s (Tr(\phi^2))^r/4 + \epsilon)
program main
use global_parameters
use initialization_calcobs
use simulation
use parallel
implicit none


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: iarg
character(128) :: config_file

integer, parameter :: num_operators=4

integer :: control
character(128) :: MEDFILE
character(10) :: C_Ckind
integer :: Ckind
integer, parameter :: N_MEDFILE=100

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Observables
complex(kind(0d0)), allocatable :: tmpmat(:,:) !tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Acomp, Acomp_reg
integer :: num_fermion, eular, ratio 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! misc
double precision, parameter :: max_epsilon=5.0d0
double precision, parameter :: interval=0.01d0
double precision :: regulator
integer :: ite
integer :: i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! initialization
iarg=iargc()
if( iarg < 1 ) then
  if (MYRANK==0) write(*,*) "use as a.out [MEDFILE] [kind of compensator]"
  stop
endif
call getarg(1,MEDFILE)
if( iarg >= 2 ) then
  call getarg(2,C_Ckind)
  read(C_Ckind,*) Ckind
else
  Ckind=1
endif
if( Ckind >= 3 .or. Ckind <= 0 ) then
  if (MYRANK==0) write(*,*) "Ckind = 1 or 2 at present."
  stop
endif

INPUT_FILE_NAME="inputfile"

call initialization 

eular=global_num_sites-global_num_links+global_num_faces 
ratio=(NMAT*NMAT-1)*eular/2
num_fermion=(global_num_sites+global_num_links+global_num_faces)*(NMAT*NMAT-1)

if( MYRANK==0 ) then
  open(N_MEDFILE, file=MEDFILE, status='OLD',action='READ',form='unformatted')
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! output measurements 
if( MYRANK == 0 ) then
  write(*,'(a,f6.2,a,f6.2)', advance='no') "# min=0.0, max=",max_epsilon,", interval=",interval
  if( Ckind == 1 ) then
    write(*,'(a)') "A = sum_s ( Tr( (\phi + \epsilon)^2 )^{-r/4} )"
  elseif( Ckind == 2 ) then 
    write(*,'(a)') "A = 1 / sum_s (Tr(\phi^2))^r/4 + \epsilon)"
  endif
endif
do 
  !! read configuration
  call read_config_from_medfile(Umat,PhiMat,ite,N_MEDFILE,control)
  call MPI_BCAST(control, 1, MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  if( control == 1 ) exit

  if( MYRANK == 0 ) then
    write(*,'(I7,2X)',advance='no') ite
  endif

  !"|Atr|", &
  call calc_trace_compensator(Acomp,PhiMat)

  if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') cdabs(Acomp)
  if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble(Acomp/cdabs(Acomp))
  if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble( (0d0,-1d0)*Acomp/cdabs(Acomp))


  !!!!!!!!!!!!!!!!
  regulator=interval 
  do while ( regulator < max_epsilon + interval )
    if( Ckind == 1 ) then
      call calc_regularized_compensator(Acomp_reg,PhiMat,regulator)
    elseif( Ckind == 2 ) then
      call calc_regularized_compensator_v2(Acomp_reg,PhiMat,regulator)
    endif
    if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') cdabs(Acomp_reg)
    if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble(Acomp_reg/cdabs(Acomp_reg))
    if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble( (0d0,-1d0)*Acomp_reg/cdabs(Acomp_reg))

    regulator = regulator + interval
  enddo


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if( MYRANK==0 ) then
    write(*,*)
  endif
enddo
if( MYRANK == 0 ) then
  close(N_MEDFILE)
endif

call stop_for_test

end program main


