module global_calcobs
implicit none

character(128), parameter :: PARAFILE="parameters_calcobs.dat"
character(128), allocatable :: MEDFILE(:)
integer, parameter :: num_calcobs=4 ! 考えているobservableの数
integer :: trig_obs(1:num_calcobs)
integer :: sizeM,sizeN

double precision :: Sb, TrX2
complex(kind(0d0)) :: Acomp_tr ! trace compensator
complex(kind(0d0)) :: Acomp_VM ! van der monde compensator
complex(kind(0d0)) :: APQ_phase ! A*/|A|
complex(kind(0d0)) :: min_eigen
complex(kind(0d0)), allocatable :: WT1(:)
complex(kind(0d0)), allocatable :: WT2(:)
integer :: num_fermion ! total fermion number
integer :: num_sitelink ! total fermion number

integer, parameter :: N_MEDFILE=100
integer, parameter :: N_PARAFILE=101

integer :: Sb_computed ! if Sb_computed=1, Sb has been already computed

end module global_calcobs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
use global_parameters
use global_calcobs
use simulation
use parallel
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! variables
complex(kind(0d0)), allocatable :: UMAT(:,:,:) ! unitary link variables
complex(kind(0d0)), allocatable :: PHIMAT(:,:,:) ! complex scalar at sites
complex(kind(0d0)), allocatable :: tmpMAT(:,:) 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! FOR PARALLEL
type(SITE_DIST), allocatable,save :: local_site_list(:) !(0:NPROCS-1)

! for type genrand_srepr
!integer :: MPI_GENRAND_SREPR, MPI_LOCAL_LABEL, MPI_GENRAND_STATE
integer IBLOCK1(1),IDISP1(1),ITYPE1(1)
integer IBLOCK2(1),IDISP2(1),ITYPE2(1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
complex(kind(0d0)) tmpobs1, tmpobs2
complex(kind(0d0)) XiPhiEta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: seed
integer :: l,ll,s,ls,tag,rank,i,j, ite
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: iarg
character(128) :: config_file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: control

iarg=iargc()
allocate( MEDFILE(1:iarg) )
if( iarg ==0 ) then
  if (MYRANK==0) write(*,*) "use as a.out [MEDFILE-1] [MEDFILE-2]..."
  stop
endif
do i=1,iarg
  call getarg(i,MEDFILE(i))
enddo

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call MPI_INIT(IERR)
call MPI_COMM_SIZE(MPI_COMM_WORLD,NPROCS,IERR)
call MPI_COMM_RANK(MPI_COMM_WORLD,MYRANK,IERR)

! setup MPI_GENRAND_SREPR
IBLOCK1(1)=1
IDISP1(1)=0
ITYPE1(1)=MPI_CHARACTER
call MPI_TYPE_STRUCT(1,IBLOCK1,IDISP1,ITYPE1,MPI_GENRAND_SREPR,IERR)
call MPI_TYPE_COMMIT(MPI_GENRAND_SREPR,IERR)

! setup MPI_LOCAL_LABEL
IBLOCK2(1)=2
IDISP2(1)=0
ITYPE2(1)=MPI_INTEGER
call MPI_TYPE_STRUCT(1,IBLOCK2,IDISP2,ITYPE2,MPI_LOCAL_LABEL,IERR)
call MPI_TYPE_COMMIT(MPI_LOCAL_LABEL,IERR)

allocate(local_site_list(0:NPROCS-1) )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call set_theory_parameters(seed)
! set the structure constant of SU(NMAT)
call set_NZF
! set global data of the simplicial complex
call set_global_simplicial_complex 
! set global parameters
overall_factor=dble(NMAT)/(2d0*LatticeSpacing*LatticeSpacing)
mass_square_phi = phys_mass_square_phi * (LatticeSpacing*latticeSpacing)
! set simulation data
!   ここで、inputfileを読み込み、各コアに必要な情報を設定して、
!   Rank=0にたまっているalphaを振り分ける。
!   この段階で、
!      num_sites,num_links,num_faces
!      local_{site,link,face}_of_globalの(1:num_{sites,links,faces})が確定する。
call set_simulation_parameters(local_site_list)

! set local data
call set_local_data(local_site_list)

allocate( UMAT(1:NMAT,1:NMAT,1:num_necessary_links) )
allocate( PHIMAT(1:NMAT,1:NMAT, 1:num_necessary_sites) )
allocate( tmpMAT(1:NMAT,1:NMAT) )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! read parameters_calcobs.dat
if( MYRANK == 0 ) then
  open(N_PARAFILE, file=PARAFILE, status='OLD',action='READ')
  !read(N_PARAFILE,*) MEDFILE
  read(N_PARAFILE,*) sizeM
  read(N_PARAFILE,*) sizeN
  do i=1,num_calcobs
    read(N_PARAFILE,*) trig_obs(i)
  enddo
endif
call MPI_BCAST(trig_obs, num_calcobs, MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(sizeM, 1, MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(sizeN, 1, MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  
num_fermion=(global_num_sites+global_num_links+global_num_faces)*(NMAT*NMAT-1)
num_sitelink=(global_num_sites+global_num_links)*(NMAT*NMAT-1)

allocate( WT1(1:sizeN/2-1) )
allocate( WT2(1:sizeN/2-1) )

do i=1, iarg

  if( MYRANK == 0 ) then
    write(*,*) "##", trim(MEDFILE(i))
    open(N_MEDFILE, file=MEDFILE(i), status='OLD',action='READ',form='unformatted')
  endif
  
  do
    call read_config_from_medfile(Umat,PhiMat,ite,N_MEDFILE,control)
    if( control == 0 ) then 
      call calc_trace_compensator(Acomp_tr,PhiMat)
      !call calc_VM_compensator(Acomp_VM,PhiMat)
      APQ_phase = dconjg(Acomp_tr) / cdabs(Acomp_tr) 
      Sb_computed=0
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! calculate observables
      if( MYRANK == 0 ) then
        write(*,'(I7,2X)',advance='no') ite
      endif
      !! trace X^2
      if( trig_obs(1) == 0 ) then 
        call calc_TrX2(TrX2,PhiMat)
        if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') TrX2
      endif

      !! bosonic action
      if( trig_obs(2) == 0 ) then 
        Sb_computed=1
        call calc_bosonic_action(Sb,Umat,PhiMat)
        if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') Sb
      endif
    
      !! trivial WT
      if( trig_obs(3) == 0 ) then 
        !do i=0,1
        if ( Sb_computed==0 ) call calc_bosonic_action(Sb,Umat,PhiMat)
        call calc_XiPhiEta(XiPhiEta,Umat,PhiMat,1)
        if( MYRANK == 0 ) then
          tmpobs1= cdabs(Acomp_tr) * &
            ( dcmplx(Sb) &
            + dcmplx(0.5d0*mass_square_phi)*XiPhiEta &
            - dcmplx(0.5d0*dble(num_sitelink)) )  
          !tmpobs2= cdabs(Acomp_VM) * &
            !( dcmplx(Sb) &
            !+ dcmplx(0.5d0*mass_square_phi)*XiPhiEta &
            !- dcmplx(0.5d0*dble(num_sitelink)) )  
          write(*,'(E15.8,2X,E15.8,2X)',advance='no') dble(tmpobs1), dble((0d0,-1d0)*tmpobs1)
        endif
        !enddo
      endif
    
      !! WT identity 1
      if( trig_obs(4) == 0 ) then
        !call calc_WT12(WT1,WT2,Umat,PhiMat,1,21)
        call calc_WT12_average(WT1,WT2,Umat,PhiMat,20,sizeM,sizeN)
        if( MYRANK==0 ) then
          do j=1,sizeN/2-1
            tmpobs1 = APQ_phase * WT1(j)
            tmpobs2 = APQ_phase * WT2(j)
            write(*,'(E15.8,2X,E15.8,2X,E15.8,2X,E15.8,2X)',advance='no') &
              dble(tmpobs1), dble((0d0,-1d0)*tmpobs1), &
              dble(tmpobs2), dble((0d0,-1d0)*tmpobs2) 
          enddo
        endif
      endif
      !! minimal eigenvalu of DD^\dagger
!      if( trig_obs(4) == 0 ) then
!        call min_eigen_DdagD(min_eigen,Umat,PhiMat)
!        if( MYRANK == 0 ) then
!          write(*,'(E15.8,2X)',advance='no') dble(min_eigen)
!        endif
!      endif

      if(MYRANK==0) write(*,*)
    else
      exit
    endif
  enddo

  if( MYRANK == 0 ) then
    close(N_MEDFILE)
  endif

enddo


end program main




