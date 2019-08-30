!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
use global_parameters
use initialization_calcobs
use simulation
use parallel
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: iarg
character(128) :: config_file

integer :: control
character(128) :: MEDFILE
character(128) :: trphi2FILE
integer, parameter :: N_MEDFILE=100
integer, parameter :: N_trphi2FILE=101

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Observables
double precision, allocatable :: trphi2(:)
double precision :: mass_reweight
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! misc
integer :: ite
integer :: rank,tag
integer :: ls, gs
double precision :: rtmp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! initialization
iarg=iargc()
if( iarg <=1 ) then
  if (MYRANK==0) write(*,*) "use as a.out [MEDFILE] [trphi2FILE]"
  stop
endif
call getarg(1,MEDFILE)
call getarg(2,trphi2FILE)
INPUT_FILE_NAME="inputfile"
call initialization 

!allocate( UMAT(1:NMAT,1:NMAT,1:num_necessary_links) )
!allocate(  PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites) )
allocate( trphi2(1:num_sites) )

if( MYRANK==0 ) then
  open(N_MEDFILE, file=MEDFILE, status='OLD',action='READ',form='unformatted')
  open(N_trphi2FILE, file=trphi2FILE, status='REPLACE')
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! output measurements 
do 
  call read_config_from_medfile(Umat,PhiMat,ite,N_MEDFILE,control)
  if( control == 0 ) then 
    if( MYRANK == 0 ) then
      write(N_trphi2FILE,'(I7,2X)',advance='no') ite
    endif
    !!!!!!!!!!!!!!!!
    call calc_trphi2(trphi2, PhiMat)
    call calc_mass_reweight(mass_reweight,PhiMat)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! write trphi2
    do gs=1,global_num_sites
      ls=local_site_of_global(gs)%label_
      rank=local_site_of_global(gs)%rank_
      tag=gs
      if( MYRANK == rank ) then
        rtmp=trphi2(ls)
        if( MYRANK /= 0 ) then 
          call MPI_SEND(rtmp,1,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,IERR)
        endif
      endif
      if( MYRANK == 0 .and. rank /= 0 ) then
        call MPI_RECV(rtmp,1,MPI_DOUBLE_PRECISION,rank,tag,MPI_COMM_WORLD,ISTATUS,IERR)
      endif
      if( MYRANK==0 ) then
        write(N_trphi2FILE,'(E15.8,2X)',advance='no') rtmp
      endif
      call MPI_BARRIER(MPI_COMM_WORLD,IERR)
    enddo

    if( MYRANK==0 ) then
      write(N_trphi2FILE,'(E15.8,2X)') mass_reweight
    endif
  else
    if( MYRANK == 0 ) then
      close(N_MEDFILE)
      close(N_trphi2FILE)
    endif
    exit
  endif
enddo

call stop_for_test
end program main


