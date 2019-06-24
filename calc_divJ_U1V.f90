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
character(128) :: DinvFILE
character(128) :: divJFILE
integer, parameter :: N_MEDFILE=100
integer, parameter :: N_DinvFILE=101
integer, parameter :: N_divJFILE=102

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Observables
complex(kind(0d0)), allocatable :: divJ1(:)
complex(kind(0d0)), allocatable :: divJ2(:)
complex(kind(0d0)), allocatable :: divJ(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! misc
integer :: ite
integer :: rank,tag
integer :: lf, gf
integer :: jj
double precision :: rtmp
complex(kind(0d0)) :: ctmp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! initialization
iarg=iargc()
if( iarg <=1 ) then
  if (MYRANK==0) write(*,*) "use as a.out [MEDFILE] [DinvFILE] [divJFILE]"
  stop
endif
call getarg(1,MEDFILE)
call getarg(2,DinvFILE)
call getarg(3,divJFILE)
INPUT_FILE_NAME="inputfile"
call initialization 

!allocate( UMAT(1:NMAT,1:NMAT,1:num_necessary_links) )
!allocate(  PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites) )
allocate( divJ(1:num_faces) )
allocate( divJ1(1:num_faces) )
allocate( divJ2(1:num_faces) )

if( MYRANK==0 ) then
  open(N_MEDFILE, file=MEDFILE, status='OLD',action='READ',form='unformatted')
  open(N_DinvFILE, file=DinvFILE, status='OLD',action='READ')
  open(N_divJFILE, file=divJFILE, status='REPLACE')
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! output measurements 
do 
  call read_config_from_medfile(Umat,PhiMat,ite,N_MEDFILE,control)
  call read_Dinv(ite, Geta_eta, Glambda_eta, Gchi_eta, &
                 Geta_lambda, Glambda_lambda, Gchi_lambda, &
                 Geta_chi, Glambda_chi, Gchi_chi, &
                 N_DinvFILE)

  if( control == 0 ) then 
    if( MYRANK == 0 ) then
      write(N_divJFILE,'(I7,2X)',advance='no') ite
    endif
    !!!!!!!!!!!!!!!!
    call calc_divJ_U1V(divJ,Glambda_eta,Gchi_lambda,UMAT)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! write divJ
    do gf=1,global_num_faces
      lf=local_face_of_global(gf)%label_
      rank=local_face_of_global(gf)%rank_
      tag=gf
      !do jj=1,2
      !  if( jj==1 ) then 
      !    divJ=divJ1
      !  else
      !    divJ=divJ2
      !  endif

      if( MYRANK == rank ) then
        ctmp=divJ(lf)
        if( MYRANK /= 0 ) then 
          call MPI_SEND(ctmp,1,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,IERR)
        endif
      endif
      if( MYRANK == 0 .and. rank /= 0 ) then
        call MPI_RECV(ctmp,1,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,ISTATUS,IERR)
      endif
      if( MYRANK==0 ) then
        write(N_divJFILE,'(E15.8,2X,E15.8,2X)',advance='no') &
          dble(ctmp), dble( (0d0,-1d0)*ctmp )
      endif
      call MPI_BARRIER(MPI_COMM_WORLD,IERR)
      !enddo
    enddo

    if( MYRANK==0 ) then
      write(N_divJFILE,*)
    endif
  else
    if( MYRANK == 0 ) then
      close(N_MEDFILE)
      close(N_divJFILE)
    endif
    exit
  endif
enddo

call stop_for_test
end program main


