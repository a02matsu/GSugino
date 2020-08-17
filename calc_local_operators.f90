!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
use global_parameters
use initialization_calcobs
use simulation
use matrix_functions, only : trace_mm, hermitian_conjugate, matrix_power, matrix_inverse
use parallel
implicit none


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: iarg
character(128) :: config_file

integer, parameter :: num_operators=4
!! 1) tr(\phibar^2)^{r/2}(f)
!! 2) tr(\phi^2)^{-r/2}(f)
!! 3) tr(Y\phibar^r)
!! 4) tr(Y\phi^{-r})

integer :: control
character(128) :: MEDFILE
character(128) :: DinvFILE
character(128) :: operatorFILE(1:num_operators)
integer, parameter :: N_MEDFILE=100
integer, parameter :: N_DinvFILE=101
integer :: N_operatorFILE(1:num_operators) !=102

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Observables
complex(kind(0d0)), allocatable :: phibar(:), phibar_site(:) !! for (1)
complex(kind(0d0)), allocatable :: phi_face(:), phi_site(:) !! for (2) 
complex(kind(0d0)), allocatable :: Yphibar(:) !! for (3) 
complex(kind(0d0)), allocatable :: Yphi(:) !! for (4) 
complex(kind(0d0)), allocatable :: Dinv(:,:)
complex(kind(0d0)), allocatable :: tmpmat(:,:) !tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)), allocatable :: tmpmat2(:,:) !tmpmat2(1:NMAT,1:NMAT)
complex(kind(0d0)), allocatable :: Uf(:,:) !Uf(1:NMAT,1:NMAT), Ymat(1:NMAT,1:NMAT)
complex(kind(0d0)), allocatable :: Ymat(:,:) !Uf(1:NMAT,1:NMAT), Ymat(1:NMAT,1:NMAT)
integer :: num_fermion
integer :: eular, ratio
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! misc
integer :: ite, ite2
integer :: rank,tag
integer :: lf, gf, ls
integer :: jj,j,i,k,ios
double precision :: rtmp, itmp
complex(kind(0d0)) :: ctmp
complex(kind(0d0)) :: phase_pf
double precision :: radius, phase



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

allocate(tmpmat(1:NMAT,1:NMAT))
allocate(tmpmat2(1:NMAT,1:NMAT))
allocate(Uf(1:NMAT,1:NMAT))
allocate(Ymat(1:NMAT,1:NMAT))
!! for 1) tr(\phibar^2)^{r/2}(f)
operatorFILE(1)=trim("OBS/trphibar"//MEDFILE(18:))
allocate( phibar(1:num_faces) )
allocate( phibar_site(1:num_necessary_sites) )
!! for 2) tr(\phi^2)^{-r/2}(f)
operatorFILE(2)=trim("OBS/trphi"//MEDFILE(18:))
allocate( phi_face(1:num_faces) )
allocate( phi_site(1:num_necessary_sites) )
!! for 3) tr(Y\phibar^r)
operatorFILE(3)=trim("OBS/Yphibar"//MEDFILE(18:))
allocate( Yphibar(1:num_faces) )
!! for 4) tr(Y\phi^{-r})
operatorFILE(4)=trim("OBS/Yphi"//MEDFILE(18:))
allocate( Yphi(1:num_faces) )

call initialization 

eular=global_num_sites-global_num_links+global_num_faces 
ratio=(NMAT*NMAT-1)*eular/2
num_fermion=(global_num_sites+global_num_links+global_num_faces)*(NMAT*NMAT-1)
allocate( Dinv(1:num_fermion, 1:num_fermion) )

if( MYRANK==0 ) then
  open(N_MEDFILE, file=MEDFILE, status='OLD',action='READ',form='unformatted')
  open(N_DinvFILE, file=DinvFILE, status='OLD',action='READ')
  do i=1,num_operators
    open(N_operatorFILE(i), file=operatorFILE(i), status='REPLACE')
  enddo
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! output measurements 
do 
  !! read configuration
  call read_config_from_medfile(Umat,PhiMat,ite,N_MEDFILE,control)
  call MPI_BCAST(control, 1, MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  if( control == 1 ) exit

  !! read Dirac inverse
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
  !! make Dinv
  call make_fermion_correlation_from_Dinv(&
      Geta_eta, Glambda_eta, Gchi_eta, &
      Geta_lambda, Glambda_lambda, Gchi_lambda, &
      Geta_chi, Glambda_chi, Gchi_chi, &
      Dinv,num_fermion)

  if( control == 0 ) then 
    if( MYRANK == 0 ) then
      do i=1,num_operators
        write(N_operatorFILE(i),'(I7,2X)',advance='no') ite
      enddo
    endif
    !!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! 1) tr(\phibar^2)^{r/2}(f)
    phibar_site=(0d0,0d0)
    call calc_phibar_compensator_site(phibar_site,PhiMat)
    phibar=(0d0,0d0)
    do lf=1,num_faces
      do k=1,sites_in_f(lf)%num_
        ls=sites_in_f(lf)%label_(k)
        phibar(lf)=phibar(lf)+phibar_site(ls)
      enddo
      phibar(lf)=phibar(lf)/dcmplx( sites_in_f(lf)%num_ )
    enddo
    call write_operator(phibar, N_operatorFILE(1))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! 2) tr(\phi^2)^{-r/2}(f)
    phi_site=(0d0,0d0)
    do ls=1,num_necessary_sites
      ctmp=(0d0,0d0)
      call trace_MM(ctmp, phimat(:,:,ls),phimat(:,:,ls))
      ctmp=ctmp/dcmplx(NMAT)
      radius=cdabs(ctmp)
      phase=atan2(dble(ctmp),dble(ctmp*(0d0,-1d0)))

      phi_site(ls)=dcmplx(radius**(-dble(ratio)/2d0)) &
        * cdexp( (0d0,1d0)*dcmplx(phase*dble(-ratio/2d0)) ) 
    enddo
    phi_face=(0d0,0d0)
    do lf=1,num_faces
      do k=1,sites_in_f(lf)%num_
        ls=sites_in_f(lf)%label_(k)
        phi_face(lf)=phi_face(lf)+phi_site(ls)
      enddo
      phi_face(lf)=phi_face(lf)/dcmplx( sites_in_f(lf)%num_ )
    enddo
    call write_operator(phi_face, N_operatorFILE(2))

write(*,*) "test5"
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! 3) tr(Y\phibar^r)
    !! 4) tr(Y\phi^{-r})
    Yphibar=(0d0,0d0)
    Yphi=(0d0,0d0)
    do lf=1,num_faces
      !! Omega
      call Make_face_variable(Uf,lf,UMAT)
      call Make_moment_map_adm(Ymat,Uf)
      Ymat = Ymat * (0d0,0.5d0)*beta_f(lf)*Ymat
      !! Yphibar
      ls=sites_in_f(lf)%label_(1)
      call hermitian_conjugate(tmpmat2,phimat(:,:,ls))
      call matrix_power(tmpmat,tmpmat2,ratio)
      call trace_mm(Yphibar(lf), Ymat, tmpmat)
      !! Yphi
      tmpmat2=phimat(:,:,ls)
      call matrix_inverse(tmpmat2)
      call matrix_power(tmpmat,tmpmat2,ratio)
      call trace_mm(Yphi(lf), Ymat, tmpmat)
    enddo
    Yphibar=Yphibar/dcmplx(NMAT)
    Yphi=Yphi/dcmplx(NMAT)
    
    call write_operator(Yphibar, N_operatorFILE(3))
    call write_operator(Yphi, N_operatorFILE(4))
write(*,*) "test6"


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

contains
  subroutine write_operator(obs, NFILE)
  use global_parameters
  use parallel
  implicit none
  
  complex(kind(0d0)), intent(in) :: obs(1:num_faces)
  integer, intent(in) :: NFILE
  
  integer :: gf, lf, rank, tag
  complex(kind(0d0)) :: ctmp

  do gf=1,global_num_faces
    lf=local_face_of_global(gf)%label_
    rank=local_face_of_global(gf)%rank_
    tag=gf

    if( MYRANK == rank ) then
      ctmp=obs(lf)
      if( MYRANK /= 0 ) then 
        call MPI_SEND(ctmp,1,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,IERR)
      endif
    endif
    if( MYRANK == 0 .and. rank /= 0 ) then
      call MPI_RECV(ctmp,1,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,ISTATUS,IERR)
    endif
    if( MYRANK==0 ) then
      write(NFILE,'(E23.16,2X,E23.16,2X)',advance='no') &
        dble(ctmp), dble( (0d0,-1d0)*ctmp )
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,IERR)
  enddo
  end subroutine write_operator
end program main

#include  "Measurement/FermionCorrelation_from_Dinv.f90"


