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
character(128) :: trF2FILE
integer, parameter :: N_MEDFILE=100
integer, parameter :: N_trphi2FILE=101
integer, parameter :: N_trF2FILE=102

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Observables
double precision, allocatable :: trF2(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! misc
integer :: ite
integer :: rank,tag
integer :: lf, gf
double precision :: rtmp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! initialization
iarg=iargc()
if( iarg <=1 ) then
  if (MYRANK==0) write(*,*) "use as a.out [MEDFILE] [trF2FILE]"
  stop
endif
call getarg(1,MEDFILE)
call getarg(2,trF2FILE)
INPUT_FILE_NAME="inputfile"
call initialization 

!allocate( UMAT(1:NMAT,1:NMAT,1:num_necessary_links) )
!allocate(  PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites) )
allocate( trF2(1:num_faces) )

open(N_MEDFILE, file=MEDFILE, status='OLD',action='READ',form='unformatted')
open(N_trF2FILE, file=trF2FILE, status='REPLACE')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! output measurements 
do 
  call read_config_from_medfile(Umat,PhiMat,ite,N_MEDFILE,control)
  if( control == 0 ) then 
    if( MYRANK == 0 ) then
      write(N_trF2FILE,'(I7,2X)',advance='no') ite
    endif
    !!!!!!!!!!!!!!!!
    call calc_trF2(trF2,Umat)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! write trF2
    do gf=1,global_num_faces
      lf=local_face_of_global(gf)%label_
      rank=local_face_of_global(gf)%rank_
      tag=gf
      if( MYRANK == rank ) then
        rtmp=trF2(lf)
        if( MYRANK /= 0 ) then 
          call MPI_SEND(rtmp,1,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,IERR)
        endif
      endif
      if( MYRANK == 0 .and. rank /= 0 ) then
        call MPI_RECV(rtmp,1,MPI_DOUBLE_PRECISION,rank,tag,MPI_COMM_WORLD,ISTATUS,IERR)
      endif
      if( MYRANK==0 ) then
        write(N_trF2FILE,'(E15.8,2X)',advance='no') rtmp
      endif
      call MPI_BARRIER(MPI_COMM_WORLD,IERR)
    enddo


    if( MYRANK==0 ) then
      write(N_trF2FILE,*)
    endif
  else
    if( MYRANK == 0 ) then
      close(N_MEDFILE)
      close(N_trF2FILE)
    endif
    exit
  endif
enddo

call stop_for_test
end program main


