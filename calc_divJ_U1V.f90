!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
use global_parameters
use initialization_calcobs
use simulation
use matrix_functions, only :matrix_product, hermitian_conjugate, make_unit_matrix, matrix_exp, matrix_trace, trace_mm
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
character(128) :: WTFILE
integer, parameter :: N_MEDFILE=100
integer, parameter :: N_DinvFILE=101
integer, parameter :: N_divJFILE=102
integer, parameter :: N_F4FILE=103
integer, parameter :: N_WTFILE=104

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Observables
complex(kind(0d0)), allocatable :: divJ1(:)
complex(kind(0d0)), allocatable :: divJ2(:)
complex(kind(0d0)), allocatable :: divJ(:)
complex(kind(0d0)), allocatable :: Dinv(:,:)
complex(kind(0d0)), allocatable :: F4(:,:)
complex(kind(0d0)), allocatable :: localFC(:) !! for (7) 
complex(kind(0d0)), allocatable :: Uf(:,:) !Uf(1:NMAT,1:NMAT)
complex(kind(0d0)), allocatable :: Ymat(:,:) !Ymat(1:NMAT,1:NMAT)
complex(kind(0d0)), allocatable :: phibar_p(:,:,:)
complex(kind(0d0)), allocatable :: tmpmat(:,:)

integer :: num_fermion
integer :: ratio
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! misc
integer :: ite, ite2
integer :: rank,tag
integer :: ls, ll, lf
integer :: gs, gl, gf, gf1,gf2
integer :: jj,i,j,k,l,ios,p
double precision :: rtmp, itmp
complex(kind(0d0)) :: ctmp, tmp,ctmp2
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


!! WT identity
if( index(MEDFILE,"/") == 8 ) then ! MEDCONF/***
  WTFILE=trim("OBS/WTU1V"//MEDFILE(index(MEDFILE,"_"):))
else
  WTFILE=trim("OBS"//MEDFILE(8:index(MEDFILE,"/"))//"WTU1V"//MEDFILE(index(MEDFILE,"_"):))
endif

INPUT_FILE_NAME="inputfile"

call initialization 

num_fermion=(global_num_sites+global_num_links+global_num_faces)*(NMAT*NMAT-1)
ratio = dimG*(global_num_sites-global_num_links+global_num_faces)/2
allocate( Dinv(1:num_fermion, 1:num_fermion) )
allocate( divJ(1:num_faces) )
allocate( divJ1(1:num_faces) )
allocate( divJ2(1:num_faces) )
allocate( F4(1:global_num_faces,global_num_faces) )
allocate( phibar_p(1:NMAT,1:NMAT,0:ratio) )
allocate(Uf(1:NMAT,1:NMAT))
allocate(Ymat(1:NMAT,1:NMAT))
allocate(tmpmat(1:NMAT,1:NMAT))
allocate( localFC(1:global_num_faces) )


if( MYRANK==0 ) then
  open(N_MEDFILE, file=MEDFILE, status='OLD',action='READ',form='unformatted')
  open(N_DinvFILE, file=DinvFILE, status='OLD',action='READ',form='unformatted')
  open(N_divJFILE, file=divJFILE, status='REPLACE')
  open(N_F4FILE, file=F4FILE, status='REPLACE')
  open(N_WTFILE, file=WTFILE, status='REPLACE')

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
      write(N_WTFILE,'(I7,2X)',advance='no') ite
    endif


    !! divJ
    call calc_DJ_U1V(DivJ1,DivJ2,Glambda_eta,Gchi_lambda,Umat)
    !! Qtr(\chi\phibar^r)
    localFC=(0d0,0d0)
    do lf=1, num_faces
      ls=sites_in_f(lf)%label_(1)
      gf=global_face_of_local(lf)
      gs=global_sites_in_f(gf)%label_(1)
    
      !! phibar_p = \bar(\PhiMat)^p
      call make_unit_matrix(phibar_p(:,:,0))
      call hermitian_conjugate(phibar_p(:,:,1), PhiMat(:,:,ls))
      do k=2,ratio
        call matrix_product(phibar_p(:,:,k),phibar_p(:,:,k-1),phibar_p(:,:,1))
      enddo
      !! Omega
      call Make_face_variable(Uf,lf,UMAT)
      call Make_moment_map_adm(Ymat,Uf)
      Ymat = Ymat * (0d0,0.5d0)*beta_f(lf)

      !! bosonic part
      call trace_mm(ctmp, Ymat, tmpmat)
      localFC(gf) = localFC(gf) + ctmp / dcmplx(NMAT)
      !! fermionic part
      do p=0,ratio-1
        do l=1,NMAT
          do k=1,NMAT
            do j=1,NMAT
              do i=1,NMAT
                localFC(gf) = localFC(gf) &
                + phibar_p(i,j,ratio-1-p)*phibar_p(k,l,p)&
                  * Geta_chi(j,k,l,i,gs,lf) / dcmplx(NMAT)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    do gf=1,global_num_faces
      rank = local_face_of_global(gf)%rank_
      call MPI_BCAST(localFC(gf),1, MPI_DOUBLE_COMPLEX,rank,MPI_COMM_WORLD,IERR)
    enddo
    localFC = localFC / dcmplx( LatticeSpacing**(ratio+2) )


    !! 4-fermi terms
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

    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! write WT identity
    do gf2=1,global_num_faces
      lf=local_face_of_global(gf2)%label_
      rank=local_face_of_global(gf2)%rank_
      tag=gf2
      !! ctmp = dJ(gf2)
      if( MYRANK == rank ) then
        ctmp = divJ1(lf) - divJ2(lf)
        if( MYRANK /= 0 ) then 
          call MPI_SEND(ctmp,1,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,IERR)
        endif
      endif
      if( MYRANK == 0 .and. rank /= 0 ) then
        call MPI_RECV(ctmp,1,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,ISTATUS,IERR)
      endif
      !!
      do gf1=1,global_num_faces
        if( MYRANK==0 ) then
          ctmp2 = localFC(gf1) * ctmp + F4(gf1,gf2)
          write(N_WTFILE,'(E23.16,2X,E23.16,2X)',advance='no') &
            dble(ctmp2), dble( (0d0,-1d0)*ctmp2 )
        endif
      enddo
    enddo
    if( MYRANK==0 ) then
      write(N_WTFILE,*)
    endif

  else
    exit
  endif
enddo
if( MYRANK == 0 ) then
  close(N_MEDFILE)
  close(N_divJFILE)
  close(N_DinvFILE)
  close(N_F4FILE)
  close(N_WTFILE)
endif

!write(*,*) MYRANK, control
!stop
call stop_for_test

end program main

#include  "Measurement/FermionCorrelation_from_Dinv.f90"


