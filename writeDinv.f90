!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Measurement 
!! % a.out [MEDFILE]
!! 
!! [MEDFILE] must be the name 
!!   "MEDCONF/medconfig_..."
!! If there is not Dinv_... in the same directory, this code create it. 
!! Otherwize, this code read it. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
use global_parameters
use initialization_calcobs
use simulation
use parallel
use matrix_functions, only : matrix_inverse, calcpfaffian, matrix_eigenvalues
implicit none

character(128) :: MEDFILE
character(128) :: DinvFILE
character(128) :: EigenFILE
integer, parameter :: N_MEDFILE=100
integer, parameter :: N_DINVFILE=101
integer, parameter :: N_EIGENFILE=102
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: l,ll,s,ls,lf,gf,tag,rank,i,j, ite, ite2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: iarg
character(128) :: config_file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: control
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
character(128) :: FMT1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer :: access !! check the existence of file
integer :: exist_dinv

complex(kind(0d0)), allocatable :: Dinv(:,:)
double precision :: rtmp,ctmp
complex(kind(0d0)), allocatable :: Dirac(:,:)
complex(kind(0d0)), allocatable :: Dirac2(:,:)
double precision :: abs_pf!, arg_pf,
double precision :: re_phase, im_phase
!double precision :: re_tmp1, im_tmp1, re_tmp2, im_tmp2
complex(kind(0d0)) :: phase_pf
complex(kind(0d0)), allocatable :: eigenvals(:)
complex(kind(0d0)) :: tmp,tmp1,tmp2
integer :: k,s1,s2,l1,l2,f1,f2,a,b
integer :: ios,info,num_fermion

!!!!!


iarg=iargc()
if( iarg ==0 ) then
  if (MYRANK==0) write(*,*) "use as a.out [MEDFILE]"
  stop
endif
call getarg(1,MEDFILE)
!call getarg(2,DinvFILE)
DinvFILE=trim("MEDCONF/Dinv"//MEDFILE(18:))
EigenFILE=trim("MEDCONF/Eigen"//MEDFILE(18:))
INPUT_FILE_NAME="inputfile"
  
call initialization 
num_fermion=(global_num_sites+global_num_links+global_num_faces)*(NMAT*NMAT-1)


allocate( Dinv(1:num_fermion, 1:num_fermion) )
allocate( Dirac(1:num_fermion, 1:num_fermion) )
allocate( Dirac2(1:num_fermion, 1:num_fermion) )
allocate( Eigenvals(1:num_fermion) )


! check DinvFile
exist_dinv=access( trim(adjustl(DinvFILE)), ' ' )
!exist_dinv=1

if( exist_dinv==0 ) then 
  stop
endif


if( MYRANK == 0 ) then
  open(N_MEDFILE, file=MEDFILE, status='OLD',action='READ',form='unformatted')
  open(N_DinvFILE, file=DinvFILE, status='REPLACE')
  open(N_EigenFILE, file=EigenFILE, status='REPLACE')
endif 

      
do    
  call read_config_from_medfile(Umat,PhiMat,ite,N_MEDFILE,control)
  call MPI_BCAST(control, 1, MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  if( control == 1 ) exit

  !! Dinvをつくる
  call construct_Dirac(Dirac,Umat,PhiMat) !correct!
  !!
  Dirac2=Dirac
  if( MYRANK==0 ) call CalcPfaffian(abs_pf,phase_pf,Dirac2)
  !!
  Dinv=Dirac
  if( MYRANK==0 ) call matrix_inverse(Dinv)
  !!
  if( MYRANK==0 ) call matrix_eigenvalues(eigenvals,Dirac)
  !! Dinvを書き出す
  if( MYRANK==0 ) then
    write(N_DinvFILE,'(I10,2X)',advance='no') ite
    do j=1,num_fermion
      do i=1,num_fermion
        write(N_DinvFILE,'(E23.15,2X,E23.15,2X)',advance='no') &
          dble(Dinv(i,j)), dble(Dinv(i,j)*(0d0,-1d0))
      enddo
    enddo
    write(N_DinvFILE,'(E23.15,2X,E23.15,2X)') &
      dble(phase_pf), dble((0d0,-1d0)*phase_pf)
    !! eigenvalues
    do i=1,num_fermion
      write(N_EigenFILE,'(E23.15,2X,E23.15,2X)',advance='no') &
        dble(eigenvals(i)), dble(eigenvals(i)*(0d0,-1d0))
    enddo
    write(N_EigenFILE,*)
  endif
enddo

if( MYRANK == 0 ) then
  close(N_MEDFILE)
  close(N_DinvFILE)
  close(N_EigenFILE)
endif

contains
#include "Measurement/construct_Dirac.f90"

end program main

