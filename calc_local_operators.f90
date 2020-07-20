!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
use global_parameters
use initialization_calcobs
use simulation
use parallel
implicit none

!! 1) tr(\phibar^2)^r(f)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: iarg
character(128) :: config_file

integer, parameter :: num_operators=1
integer :: control
character(128) :: MEDFILE
character(128) :: DinvFILE
character(128) :: operatorFILE(1:num_operators)
integer, parameter :: N_MEDFILE=100
integer, parameter :: N_DinvFILE=101
integer :: N_operatorFILE(1:num_operators) !=102

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Observables
complex(kind(0d0)), allocatable :: phibar(:), phibar_site(:)
complex(kind(0d0)), allocatable :: Dinv(:,:)
integer :: num_fermion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! misc
integer :: ite, ite2
integer :: rank,tag
integer :: lf, gf, ls
integer :: jj,j,i,k,ios
double precision :: rtmp, itmp
complex(kind(0d0)) :: ctmp
complex(kind(0d0)) :: phase_pf


do i=1,num_operators
  N_operatorFILE(i)=101+i
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! initialization
iarg=iargc()
if( iarg < 1 ) then
  if (MYRANK==0) write(*,*) "use as a.out [MEDFILE]"
  stop
endif
call getarg(1,MEDFILE)
DinvFILE=trim("MEDCONF/Dinv"//MEDFILE(18:))
INPUT_FILE_NAME="inputfile"

!! op1=tr(phibar^2)^r
operatorFILE(1)=trim("OBS/phibar"//MEDFILE(18:))
allocate( phibar(1:num_faces) )
allocate( phibar_site(1:num_necessary_sites) )

call initialization 

num_fermion=(global_num_sites+global_num_links+global_num_faces)*(NMAT*NMAT-1)
allocate( Dinv(1:num_fermion, 1:num_fermion) )

if( MYRANK==0 ) then
  open(N_MEDFILE, file=MEDFILE, status='OLD',action='READ',form='unformatted')
  open(N_DinvFILE, file=DinvFILE, status='OLD',action='READ')
  open(N_operatorFILE(1), file=operatorFILE(1), status='REPLACE')
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! output measurements 
do 
  call read_config_from_medfile(Umat,PhiMat,ite,N_MEDFILE,control)
  
  if( MYRANK==0 ) then
    read(N_DinvFILE,'(I10,2X)',advance='no',iostat=ios) ite2
    if( ios == -1) control=1
    if( control==0 ) then 
      do j=1,num_fermion
        do i=1,num_fermion
          read(N_DinvFILE,'(E23.15,2X,E23.15,2X)',advance='no') &
            rtmp,itmp
            Dinv(i,j)=dcmplx(rtmp)+(0d0,1d0)*itmp
        enddo
      enddo
      read(N_DinvFILE,'(E23.15,2X,E23.15,2X)') rtmp, itmp
      phase_pf=dcmplx(rtmp)+(0d0,1d0)*dcmplx(itmp)
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
      write(N_operatorFILE(1),'(I7,2X)',advance='no') ite
    endif
    !!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! write phibar
    call calc_phibar_compensator_site(phibar_site,PhiMat)
    phibar=(0d0,0d0)
    do lf=1,num_faces
      do k=1,sites_in_f(lf)%num_
        ls=sites_in_f(lf)%label_(k)
        phibar(lf)=phibar(lf)+phibar_site(ls)
      enddo
      phibar(lf)=phibar(lf)/dcmplx( sites_in_f(lf)%num_ )
    enddo
    do gf=1,global_num_faces
      lf=local_face_of_global(gf)%label_
      rank=local_face_of_global(gf)%rank_
      tag=gf

      if( MYRANK == rank ) then
        ctmp=phibar(lf)
        if( MYRANK /= 0 ) then 
          call MPI_SEND(ctmp,1,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,IERR)
        endif
      endif
      if( MYRANK == 0 .and. rank /= 0 ) then
        call MPI_RECV(ctmp,1,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,ISTATUS,IERR)
      endif
      if( MYRANK==0 ) then
        write(N_operatorFILE(1),'(E23.16,2X,E23.16,2X)',advance='no') &
          dble(ctmp), dble( (0d0,-1d0)*ctmp )
      endif
      call MPI_BARRIER(MPI_COMM_WORLD,IERR)
    enddo

    if( MYRANK==0 ) then
      do i=1,num_operators
        write(N_operatorFILE(i),*)
      enddo
    endif
  else
    exit
  endif
enddo
if( MYRANK == 0 ) then
  close(N_MEDFILE)
  close(N_DinvFILE)
  do i=1,num_operators
    close(N_operatorFILE(i))
  enddo
endif

!write(*,*) MYRANK, control
!stop
call stop_for_test
end program main

#include  "Measurement/FermionCorrelation_from_Dinv.f90"


