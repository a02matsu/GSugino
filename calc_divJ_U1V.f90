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
character(128) :: F4FILE
integer, parameter :: N_MEDFILE=100
integer, parameter :: N_DinvFILE=101
integer, parameter :: N_divJFILE=102
integer, parameter :: N_F4FILE=103

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Observables
complex(kind(0d0)), allocatable :: divJ1(:)
complex(kind(0d0)), allocatable :: divJ2(:)
complex(kind(0d0)), allocatable :: divJ(:)
complex(kind(0d0)), allocatable :: Dinv(:,:)
complex(kind(0d0)), allocatable :: F4(:,:)
integer :: num_fermion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! misc
integer :: ite, ite2
integer :: rank,tag
integer :: ls, ll, lf
integer :: gs, gl, gf, gf1,gf2
integer :: jj,i,j,k,l,ios
double precision :: rtmp, itmp
complex(kind(0d0)) :: ctmp, tmp
complex(kind(0d0)) :: phase_pf

!complex(kind(0d0)),allocatable :: tttt(:,:)
!allocate(tttt(1:2,1:3))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! initialization
iarg=iargc()
if( iarg < 1 ) then
  if (MYRANK==0) write(*,*) "use as a.out [MEDFILE]"
  stop
endif
call getarg(1,MEDFILE)

!! DinvFILE
DinvFILE=trim(MEDFILE(1:index(MEDFILE,"/"))//"Dinv"//MEDFILE(index(MEDFILE,"_"):))

!! divJFILE
if( index(MEDFILE,"/") == 8 ) then ! MEDCONF/***
  divJFILE=trim("OBS/U1V"//MEDFILE(index(MEDFILE,"_"):))
else
  divJFILE=trim("OBS"//MEDFILE(8:index(MEDFILE,"/"))//"U1V"//MEDFILE(index(MEDFILE,"_"):))
endif

!! F4FILE ! 4-fermion part of OdJ
if( index(MEDFILE,"/") == 8 ) then ! MEDCONF/***
  F4FILE=trim("OBS/F4dJ"//MEDFILE(index(MEDFILE,"_"):))
else
  F4FILE=trim("OBS"//MEDFILE(8:index(MEDFILE,"/"))//"F4dJ"//MEDFILE(index(MEDFILE,"_"):))
endif

INPUT_FILE_NAME="inputfile"

call initialization 

num_fermion=(global_num_sites+global_num_links+global_num_faces)*(NMAT*NMAT-1)
allocate( Dinv(1:num_fermion, 1:num_fermion) )
allocate( divJ(1:num_faces) )
allocate( divJ1(1:num_faces) )
allocate( divJ2(1:num_faces) )
allocate( F4(1:global_num_faces,global_num_faces) )

if( MYRANK==0 ) then
  open(N_MEDFILE, file=MEDFILE, status='OLD',action='READ',form='unformatted')
  open(N_DinvFILE, file=DinvFILE, status='OLD',action='READ',form='unformatted')
  open(N_divJFILE, file=divJFILE, status='REPLACE')
  open(N_F4FILE, file=F4FILE, status='REPLACE')

  write(N_divJFILE,*) "# ite, Re(rot(J1)), Im(rot(J1), Re(div(J2)), Im(div(J2)), Re(DJ), Im(DJ)"
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! output measurements 
do 
  call read_config_from_medfile(Umat,PhiMat,ite,N_MEDFILE,control)
  
  if( MYRANK==0 ) then
    read(N_DinvFILE,iostat=ios) ite2
    if( ios == -1) then
      control=1
    else
      read(N_DinvFILE) Dinv
      read(N_DinvFILE) phase_pf
      read(N_DinvFILE) 
    endif
  endif
  call MPI_BCAST(control, 1, MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  if( control == 1 ) exit

  call make_fermion_correlation_from_Dinv(&
      Geta_eta, Glambda_eta, Gchi_eta, &
      Geta_lambda, Glambda_lambda, Gchi_lambda, &
      Geta_chi, Glambda_chi, Gchi_chi, &
      Dinv,num_fermion)

  if( control == 0 ) then 
    if( MYRANK == 0 ) then
      write(N_divJFILE,'(I7,2X)',advance='no') ite
      write(N_F4FILE,'(I7,2X)',advance='no') ite
    endif
    call calc_DJ_U1V(DivJ1,DivJ2,Glambda_eta,Gchi_lambda,Umat)

    call calc_OdJ_4F(F4, \
      Geta_eta,\
      Gchi_eta,\
      Geta_chi,\
      Geta_lambda,\
      Gchi_lambda,\
      Gchi_chi, \
      Phimat,Umat)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! write divJ
    do gf=1,global_num_faces
      lf=local_face_of_global(gf)%label_
      rank=local_face_of_global(gf)%rank_
      tag=gf

      do jj=1,3
        if( jj==1 ) then 
          divJ=divJ1
        elseif( jj==2 ) then
          divJ=divJ2
        else
          divJ=divJ1-divJ2
        endif

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
          write(N_divJFILE,'(E23.16,2X,E23.16,2X)',advance='no') &
            dble(ctmp), dble( (0d0,-1d0)*ctmp )
        endif
        call MPI_BARRIER(MPI_COMM_WORLD,IERR)

      enddo
    enddo
    if( MYRANK==0 ) then
      write(N_divJFILE,*)
    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! write OdJ
    if( MYRANK==0 ) then
      do gf1=1,global_num_faces
        do gf2=1,global_num_faces
          write(N_F4FILE,'(E23.16,2X,E23.16,2X)',advance='no') &
            dble(F4(gf1,gf2)), dble( (0d0,-1d0)*F4(gf1,gf2) )
        enddo
      enddo
      write(N_F4FILE,*)
    endif
  else
    exit
  endif
enddo
if( MYRANK == 0 ) then
  close(N_MEDFILE)
  close(N_divJFILE)
  close(N_DinvFILE)
endif

!write(*,*) MYRANK, control
!stop
call stop_for_test

end program main

#include  "Measurement/FermionCorrelation_from_Dinv.f90"


