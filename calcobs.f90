!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Measurement 
!! % a.out [MEDFILE]
!! 
!! [MEDFILE] must be the name 
!!   "MEDCONF/medconfig_..."
!! If there is not Dinv_... in the same directory, this code create it. 
!! Otherwize, this code read it. 
module global_calcobs
implicit none

character(128), parameter :: PARAFILE="parameters_calcobs.dat"
character(128) :: MEDFILE
character(128) :: DinvFILE
character(128) :: EigenFILE
integer, parameter :: num_calcobs=87 ! 考えているobservableの数
character(128) :: name_obs(1:num_calcobs) = (/ &
  "Re(phase Pf)",  &
  "Im(phase Pf)", &
  "|Atr_pre|", &
  "Re(phase Atr_pre)", &
  "Im(phase Atr_pre)", &
  "SbS", &
  "SbL", &
  "SbF", &
  "dimG*NF/2", &
  "Re(Sf_site)", &
  "Im(Sf_site)", &
  "Re(Sf_link1)", &
  "Im(Sf_link1)", &
  "Re(Sf_link2)", &
  "Im(Sf_link2)", &
  "Re(Sf_face1)", &
  "Im(Sf_face1)", &
  "Re(Sf_face2)", &
  "Im(Sf_face2)", &
  "all Sf", &
  "Re(WTmass_site)", &
  "Im(WTmass_site)", &
  "Re(WTmass_link)", &
  "Im(WTmass_link)", &
  "Re(WTmass_face)", &
  "Im(WTmass_face)", &
  "Re(WT site)", &
  "Im(WT site)", &
  "Re(WT link)", &
  "Im(WT link)", &
  "Re(WT face)", &
  "Im(WT face)", &
  "|Atr|", &
  "Re(phase Atr)", &
  "Im(phase Atr)", &
  "|Atr2|", &
  "Re(phase Atr2)", &
  "Im(phase Atr2)", &
  "|Areg(0.001)|", &
  "Re(phase Areg(0.01))", &
  "Im(phase Areg(0.01))", &
  "|Aface|", &
  "Re(phase Aface)", &
  "Im(phase Aface)", &
  "|Af-SF4 site|", &
  "Re(phase Af-SF4 site)", &
  "Im(phase Af-SF4 site)", &
  "|Af-SF4 link|", &
  "Re(phase Af-SF4 link)", &
  "Im(phase Af-SF4 link)", &
  "|Af-SF4 face|", &
  "Re(phase Af-SF4 face)", &
  "Im(phase Af-SF4 face)", &
  "|AYphi|", &
  "Re(phase AYphi)", &
  "Im(phase AYphi)", &
  "|AYphibar|", &
  "Re(phase AYphibar)", &
  "Im(phase AYphibar)", &
  "Re(Q(AYphibar)XiS)", &
  "Im(Q(AYphibar)XiS)", &
  "Re(Q(AYphibar)XiL)", &
  "Im(Q(AYphibar)XiL)", &
  "Re(Q(AYphibar)XiF)", &
  "Im(Q(AYphibar)XiF)", &
  "|Aphibar|", &
  "Re(phase Aphibar)", &
  "Im(phase Aphibar)", &
  "Re(Q(Aphibar)XiS)", &
  "Im(Q(Aphibar)XiS)", &
  "Re(Q(Aphibar)XiL)", &
  "Im(Q(Aphibar)XiL)", &
  "Re(Q(Aphibar)XiF)", &
  "Im(Q(Aphibar)XiF)", &
  "Tr|phi^2|", &
  "Re(WT Atr site)", &
  "Im(WT Atr site)", &
  "Re(WT Atr link)", &
  "Im(WT Atr link)", &
  "Re(WT Atr face)", &
  "Im(WT Atr face)", &
  "Re(WT Aface site)", &
  "Im(WT Aface site)", &
  "Re(WT Aface link)", &
  "Im(WT Aface link)", &
  "Re(WT Aface face)", &
  "Im(WT Aface face)" &
  /)
  !"Re(exactSf_link1)", &
  !"Im(exactSf_link1)", &
  !"Re(exactSf_link2)", &
  !"Im(exactSf_link2)", &
  !"Re(exactSf_face2)", &
  !"Im(exactSf_face2)"  &
  !"Re(XiQS_site)", &
  !"Im(XiQS_site)", &
  !"Re(XiQS_link)", &
  !"Im(XiQS_link)", &
  !"Re(XiQS_face)", &
  !"Im(XiQS_face)", &

!integer :: trig_obs(1:num_calcobs)
integer :: sizeM,sizeN

double precision :: Sb, SbS, SbL, SbF
complex(kind(0d0)) :: SfL2
complex(kind(0d0)) :: Acomp_trpre ! trace compensator
complex(kind(0d0)) :: Acomp_tr ! trace compensator
complex(kind(0d0)) :: Acomp_tr2 ! trace compensator2
complex(kind(0d0)) :: Acomp_face ! face compensator
complex(kind(0d0)) :: Acomp_Yphi ! face compensator
complex(kind(0d0)) :: Acomp_reg3 ! regularized compensator
complex(kind(0d0)) :: Acomp_reg2 ! regularized compensator
complex(kind(0d0)) :: Acomp_reg1 ! regularized compensator
complex(kind(0d0)) :: Acomp_reg05 ! regularized compensator
complex(kind(0d0)) :: Acomp_reg01 ! regularized compensator
complex(kind(0d0)) :: CSF_site ! 4-fermi term in Aface-(SFsite+mass)
complex(kind(0d0)) :: CSF_link ! 4-fermi term in Aface-(SFlink+mass)
complex(kind(0d0)) :: CSF_face ! 4-fermi term in Aface-(SFface+mass)
complex(kind(0d0)) :: Acomp_phibar ! phia compensator
!complex(kind(0d0)) :: Atr_phase ! A*/|A|
!complex(kind(0d0)) :: Aface_phase ! A*/|A|
!complex(kind(0d0)) :: Aphibar_phase ! A*/|A|
complex(kind(0d0)) :: QC_Xi_site ! Q(Aphibar).\Xi 
complex(kind(0d0)) :: QC_Xi_link ! Q(Aphibar).\Xi 
complex(kind(0d0)) :: QC_Xi_face ! Q(Aphibar).\Xi 
complex(kind(0d0)) :: min_eigen
complex(kind(0d0)) :: mass_cont_site
complex(kind(0d0)) :: mass_cont_link
complex(kind(0d0)) :: mass_cont_face
complex(kind(0d0)) :: WT_site
complex(kind(0d0)) :: WT_link
complex(kind(0d0)) :: WT_face
complex(kind(0d0)) :: Sf1, Sf2, Sf3, Sf4, Sf5
complex(kind(0d0)) :: XiQS
complex(kind(0d0)) :: tmp_obs
!complex(kind(0d0)), allocatable :: WT1(:)
!complex(kind(0d0)), allocatable :: WT2(:)
complex(kind(0d0)) :: WT1, WT2
double precision :: trphi2
double precision :: pf_arg
integer :: num_fermion ! total fermion number
integer :: num_sitelink ! total fermion number

integer, parameter :: N_MEDFILE=100
integer, parameter :: N_PARAFILE=101
integer, parameter :: N_DinvFILE=102
integer, parameter :: N_EigenFILE=103

integer :: Sb_computed ! if Sb_computed=1, Sb has been already computed

contains
!!!
subroutine make_format(FMT1,num)
implicit none

character(128), intent(out) :: FMT1
integer, intent(in) :: num
character(128), parameter :: FMT0="E15.8,2X"

character(128) :: tmp
integer :: i

tmp = trim("(")
do i=1,num-1
  FMT1 = tmp // trim(FMT0) // ","
  tmp = FMT1
  write(*,*)  tmp // trim(FMT0) // ","
enddo
FMT1 = tmp // trim(FMT0) // ")"
write(*,*) FMT1
end subroutine make_format

#include "Measurement/FermionCorrelation_from_Dinv.f90"
#include "Measurement/construct_Dirac.f90"
#include "Measurement/writeout_Dirac.f90"



end module global_calcobs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
use global_parameters
use global_calcobs
use initialization_calcobs
use simulation
use parallel
use matrix_functions, only : matrix_inverse, calcpfaffian, matrix_eigenvalues
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: l,ll,s,ls,lf,gf,tag,rank,i,j, ite, ite2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: iarg
character(128) :: config_file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: control
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
character(128) :: FMT1
integer :: pos_current
integer :: pos_operator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
complex(kind(0d0)) tmpobs1, tmpobs2
complex(kind(0d0)) XiPhiEta

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
integer :: ios

!!!!!
complex(kind(0d0)), allocatable :: DJ1(:), DJ2(:)


iarg=iargc()
if( iarg ==0 ) then
  !if (MYRANK==0) write(*,*) "use as a.out [MEDFILE] [DinvFILE]"
  if (MYRANK==0) write(*,*) "use as a.out [MEDFILE]"
  stop
endif
call getarg(1,MEDFILE)
!call getarg(2,DinvFILE)
if( iarg == 1 ) then
  DinvFILE=trim("MEDCONF/Dinv"//MEDFILE(18:))
else
  call getarg(2,DinvFILE)
endif
if( iarg <= 2 ) then
  EigenFILE=trim("MEDCONF/Eigen"//MEDFILE(18:))
else
  call getarg(3,EigenFILE)
endif
INPUT_FILE_NAME="inputfile"
  
call initialization 
num_fermion=(global_num_sites+global_num_links+global_num_faces)*(NMAT*NMAT-1)
num_sitelink=(global_num_sites+global_num_links)*(NMAT*NMAT-1)


allocate( Dinv(1:num_fermion, 1:num_fermion) )
allocate( Dirac(1:num_fermion, 1:num_fermion) )
allocate( Dirac2(1:num_fermion, 1:num_fermion) )
allocate( eigenvals(1:num_fermion) )

allocate( DJ1(1:num_faces) )
allocate( DJ2(1:num_faces) )

! check DinvFile
exist_dinv=access( trim(adjustl(DinvFILE)), ' ' )
!exist_dinv=1

! write contents
if( MYRANK == 0 ) then
  write(*,'(a)',advance='no') "#"
    write(*,'(a)',advance='no') "1) ite, "
  do i=1, num_calcobs
    write(*,'(I3,a,a,",")',advance='no') i+1, ") ", trim(name_obs(i))
  enddo
  !!!!!!!!
  write(*,*) "##", trim(MEDFILE)
  open(N_MEDFILE, file=MEDFILE, status='OLD',action='READ',form='unformatted')
  if( exist_dinv==0 ) then 
    !open(N_DinvFILE, file=DinvFILE, status='OLD',action='READ')
    open(N_DinvFILE, file=DinvFILE, status='OLD',action='READ',form='unformatted')
  else
    !open(N_DinvFILE, file=DinvFILE, status='REPLACE')
    open(N_DinvFILE, file=DinvFILE, status='REPLACE',form='unformatted')
    open(N_EigenFILE, file=EigenFILE, status='REPLACE')
  endif
endif 

      
do    
  call read_config_from_medfile(Umat,PhiMat,ite,N_MEDFILE,control)
  call MPI_BCAST(control, 1, MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  if( control == 1 ) exit

  !! Dinvが存在するときは読み込む
  if( exist_dinv==0 ) then 
    if( MYRANK==0 ) then 
      !read(N_DinvFILE,'(I10,2X)',advance='no',iostat=ios) ite2
      read(N_DinvFILE,iostat=ios) ite2
      if( ios == -1) then 
        control = 1
      else
        !do j=1,num_fermion
        !  do i=1,num_fermion
        !    read(N_DinvFILE) rtmp,ctmp
        !      Dinv(i,j)=dcmplx(rtmp)+(0d0,1d0)*ctmp
        !  enddo
        !enddo
        !read(N_DinvFILE,'(E23.15,2X,E23.15,2X)') re_phase, im_phase
        !phase_pf=dcmplx(re_phase)+(0d0,1d0)*dcmplx(im_phase)
        read(N_DinvFILE) Dinv
        read(N_DinvFILE) phase_pf
        read(N_DinvFILE) 
      endif
    endif
    call MPI_BCAST(control, 1, MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
    if( control == 1 ) exit
  
    call MPI_BCAST(ite,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ite2,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
    if( ite .ne. ite2 ) then 
      write(*,*) "MEDFILE and Dinv do not match."
      control=1
      !call stop_for_test
    endif
    call MPI_BCAST(control, 1, MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
    if( control == 1 ) exit
  !! さもなくば作る
  else
    call construct_Dirac(Dirac,Umat,PhiMat) !correct!
    !!
    Dirac2=Dirac
    if( MYRANK==0 ) call CalcPfaffian(abs_pf,phase_pf,Dirac2)
    !!
    Dinv=Dirac
    if( MYRANK==0 ) call matrix_inverse(Dinv)
    !!
    if( MYRANK==0 ) call matrix_eigenvalues(eigenvals,Dirac)

    !write(N_DinvFILE,'(I10,2X)',advance='no') ite
    write(N_DinvFILE) ite
    !do j=1,num_fermion
    !  do i=1,num_fermion
    !    write(N_DinvFILE,'(E23.15,2X,E23.15,2X)',advance='no') &
    !      dble(Dinv(i,j)), dble(Dinv(i,j)*(0d0,-1d0))
    !  enddo
    !enddo
    write(N_DinvFILE) Dinv
    !write(N_DinvFILE,'(E23.15,2X,E23.15,2X)') &
      !dble(phase_pf), dble((0d0,-1d0)*phase_pf)
    write(N_DinvFILE) phase_pf
    write(N_DinvFILE) !! necessary?? 
    !! eigenvalues
    do i=1,num_fermion
      write(N_EigenFILE,'(E23.15,2X,E23.15,2X)',advance='no') &
        dble(eigenvals(i)), dble(eigenvals(i)*(0d0,-1d0))
    enddo
    write(N_EigenFILE,*)
  endif

  call make_fermion_correlation_from_Dinv(&
    Geta_eta, Glambda_eta, Gchi_eta, &
    Geta_lambda, Glambda_lambda, Gchi_lambda, &
    Geta_chi, Glambda_chi, Gchi_chi, &
    Dinv,num_fermion)

  !!!!!!!!!!!!!!!!!!!!!!!
  !!! CHECK ROUTINES !!!!
   call check_DinvPF(&
    Geta_eta, Glambda_eta, Gchi_eta, &
    Geta_lambda, Glambda_lambda, Gchi_lambda, &
    Geta_chi, Glambda_chi, Gchi_chi, &
    Umat,PhiMat,1)
  !!!!!!!!!!!!!!!!!!!!!!!

  if( control == 0 ) then 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! calculate observables
    if( MYRANK == 0 ) then
      write(*,'(I7,2X)',advance='no') ite
    endif

    !write(*,*) ite
    !do ll=1,global_num_links
    !  if(MYRANK==local_link_of_global(ll)%rank_) then
    !    write(*,*) Umat(:,:,local_link_of_global(ll)%label_)
    !  endif
    !enddo

    !! Pfaffian phase
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble(phase_pf)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble( (0d0,-1d0)*phase_pf)
      !if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble((0d0,-1d0)*Sf)

    !"|Atr_pre|", &
    Acomp_trpre=(0d0,0d0)
      call calc_trace_compensator_pre(Acomp_trpre,PhiMat)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') cdabs(Acomp_trpre)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble(Acomp_trpre/cdabs(Acomp_trpre))
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble( (0d0,-1d0)*Acomp_trpre ) /cdabs(Acomp_trpre)


    Sb_computed=0
    !! SbS
      call calc_bosonic_action_site(SbS,PhiMat)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  SbS

    !! SbL
      call calc_bosonic_action_link(SbL,Umat,PhiMat)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  SbL

    !! SbF
      call calc_bosonic_action_face(SbF,Umat)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  SbF

    !! fermion number in (5.8)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  &
        dble(global_num_faces*dimG)*0.5d0 

    !! Sf_site
      call calc_Sf_site(Sf1,Geta_eta,PhiMat)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble(Sf1)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble((0d0,-1d0)*Sf1)

    !! Sf_link1
      call calc_Sf_link1(Sf2,Geta_lambda,Umat,PhiMat)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble(Sf2)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble((0d0,-1d0)*Sf2)

    !! Sf_link2
      call calc_Sf_link2(Sf3, PhiMat, Umat, Glambda_lambda)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble(Sf3)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble((0d0,-1d0)*Sf3)

    !! Sf_face1
      call calc_Sf_face1(Sf4,Gchi_chi,PhiMat)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble(Sf4)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble((0d0,-1d0)*Sf4)

    !! Sf_face2
      call calc_Sf_face2(Sf5,Glambda_chi,Umat)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble(Sf5)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble((0d0,-1d0)*Sf5)
      !if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble((0d0,-1d0)*Sf)

    !! Total Sf
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') &
      dble(Sf1+Sf2+Sf3+Sf4+Sf5)

    !! mass contribution in WT_site
      call mass_contribution_site(mass_cont_site,Geta_eta,PhiMat)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble(mass_cont_site)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble((0d0,-1d0)*mass_cont_site)

    !! mass contribution in WT_link
      call  mass_contribution_link(mass_cont_link,Glambda_eta,Umat,PhiMat)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble(mass_cont_link)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble((0d0,-1d0)*mass_cont_link)

    !! mass contribution in WT_face
      call mass_contribution_face(mass_cont_face,Gchi_eta,Umat,PhiMat)
       if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble(mass_cont_face)
       if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble((0d0,-1d0)*mass_cont_face)

    !! WT identity site
       WT_site=dcmplx(SbS)+Sf1+mass_cont_site
       if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble(WT_site)
       if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble((0d0,-1d0)*WT_site)

    !! WT identity link
       WT_link=dcmplx(SbL)+Sf2+Sf3+mass_cont_link
       if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble(WT_link)
       if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble((0d0,-1d0)*WT_link)

    !! WT identity face
       WT_face=dcmplx(SbF)+Sf4+Sf5+mass_cont_face+(0.5d0,0d0)*dcmplx( (NMAT*NMAT-1)*(global_num_faces) )
       if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble(WT_face)
       if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble((0d0,-1d0)*WT_face)

    !"|Atr|", &
    Acomp_tr=(0d0,0d0)
      call calc_trace_compensator(Acomp_tr,PhiMat)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') cdabs(Acomp_tr)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble(Acomp_tr/cdabs(Acomp_tr))
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble( (0d0,-1d0)*Acomp_tr/cdabs(Acomp_tr))

    !"|Atr|", &
    Acomp_tr=(0d0,0d0)
      call calc_trace_compensator2(Acomp_tr2,PhiMat)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') cdabs(Acomp_tr2)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble(Acomp_tr2/cdabs(Acomp_tr2))
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble( (0d0,-1d0)*Acomp_tr2/cdabs(Acomp_tr2))

    !"|Areg(0.1)|", &
      call calc_regularized_compensator(Acomp_reg01,PhiMat,1d-3)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') cdabs(Acomp_reg01)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble(Acomp_reg01/cdabs(Acomp_reg01))
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble( (0d0,-1d0)*Acomp_reg01/cdabs(Acomp_reg01))

    !! Face compensator
      call calc_face_compensator(&
        Acomp_face,CSF_site,CSF_link,CSF_face,&
        Umat,PhiMat,&
        Geta_eta, Geta_lambda, Geta_chi, Gchi_eta, Gchi_lambda, Gchi_chi) 
     !call calc_face_compensator(Acomp_face,Umat,PhiMat,Geta_chi)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') cdabs(Acomp_face)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble(Acomp_face/cdabs(Acomp_face))
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble( (0d0,-1d0)*Acomp_face/cdabs(Acomp_face))

    !"|Af-SF4_site|", &
      !call calc_4fermi_in_CSFsite(CSF_site, Umat, Phimat, Geta_eta, Gchi_eta )
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') cdabs(CSF_site)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble(CSF_site/cdabs(CSF_site))
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble( (0d0,-1d0)*CSF_site/cdabs(CSF_site))

    !"|Af-SF4_link|", &
      !call calc_4fermi_in_CSFlink(CSF_link, Umat, Phimat, Geta_eta, Gchi_eta, Geta_lambda, Gchi_lambda )
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') cdabs(CSF_link)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble(CSF_link/cdabs(CSF_link))
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble( (0d0,-1d0)*CSF_link/cdabs(CSF_link))

    !"|Af-SF4_face|", &
      !call calc_4fermi_in_CSFface(CSF_face, Umat, Phimat, Geta_eta, Gchi_eta, Geta_chi, Gchi_chi, Geta_lambda, Gchi_lambda )
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') cdabs(CSF_face)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble(CSF_face/cdabs(CSF_face))
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble( (0d0,-1d0)*CSF_face/cdabs(CSF_face))


    !"|AYphi|", &
      call calc_Yphi_compensator(Acomp_Yphi,PhiMat,Umat)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') cdabs(Acomp_Yphi)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble(Acomp_Yphi/cdabs(Acomp_Yphi))
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble( (0d0,-1d0)*Acomp_Yphi/cdabs(Acomp_Yphi))

    !"|AYphibar|", &
      call calc_Yphibar_compensator(Acomp_phibar,PhiMat,Umat)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') cdabs(Acomp_phibar)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble(Acomp_phibar/cdabs(Acomp_phibar))
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble( (0d0,-1d0)*Acomp_phibar/cdabs(Acomp_phibar))

    !"QCYphibar_Xi", &
      call calc_QCYphibar_Xi(QC_Xi_site,QC_Xi_link,QC_Xi_face,&
        Geta_eta,Geta_lambda,Geta_chi,&
        Gchi_eta,Gchi_lambda,Gchi_chi,&
        Umat,PhiMat)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  &
        dble(QC_Xi_site / (Acomp_phibar/dcmplx(cdabs(Acomp_phibar))))
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  &
        dble((0d0,-1d0) * QC_Xi_site / (Acomp_phibar/dcmplx(cdabs(Acomp_phibar))))
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  &
        dble(QC_Xi_link / (Acomp_phibar/dcmplx(cdabs(Acomp_phibar))))
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  &
        dble((0d0,-1d0) * QC_Xi_link / (Acomp_phibar/dcmplx(cdabs(Acomp_phibar))))
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  &
        dble(QC_Xi_face / (Acomp_phibar/dcmplx(cdabs(Acomp_phibar))))
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  &
        dble((0d0,-1d0) * QC_Xi_face / (Acomp_phibar/dcmplx(cdabs(Acomp_phibar))))


    !"|Aphibar|", &
      call calc_phibar_compensator(Acomp_phibar,PhiMat)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') cdabs(Acomp_phibar)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble(Acomp_phibar/cdabs(Acomp_phibar))
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble( (0d0,-1d0)*Acomp_phibar/cdabs(Acomp_phibar))

    !"QCphibar_Xi", &
      call calc_QCphibar_Xi(QC_Xi_site,QC_Xi_link,QC_Xi_face,Geta_eta,Geta_lambda,Geta_chi,Umat,PhiMat)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  &
        dble(QC_Xi_site / (Acomp_phibar/dcmplx(cdabs(Acomp_phibar))))
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  &
        dble((0d0,-1d0) * QC_Xi_site / (Acomp_phibar/dcmplx(cdabs(Acomp_phibar))))
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  &
        dble(QC_Xi_link / (Acomp_phibar/dcmplx(cdabs(Acomp_phibar))))
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  &
        dble((0d0,-1d0) * QC_Xi_link / (Acomp_phibar/dcmplx(cdabs(Acomp_phibar))))
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  &
        dble(QC_Xi_face / (Acomp_phibar/dcmplx(cdabs(Acomp_phibar))))
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  &
        dble((0d0,-1d0) * QC_Xi_face / (Acomp_phibar/dcmplx(cdabs(Acomp_phibar))))


    !"tr|\phi^2|", &
      call calc_trphi2(trphi2,PhiMat)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') trphi2

!    !! Sf_link1
!      call calc_Qexact_Sf_link1(Sf2,Geta_lambda,Umat,PhiMat)
!      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble(Sf2)
!      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble((0d0,-1d0)*Sf2)
!
!    !! Sf_link2
!      call calc_Qexact_Sf_link2(Sf3, PhiMat, Umat, Glambda_lambda)
!      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble(Sf3)
!      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble((0d0,-1d0)*Sf3)
!
!    !! Sf_face2
!      call calc_Qexact_Sf_face2(Sf5,Glambda_chi,Umat)
!      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble(Sf5)
!      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble((0d0,-1d0)*Sf5)
      !if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble((0d0,-1d0)*Sf)

    !! XiQS_site
      !call calc_Xisite_QS(XiQS,Geta_eta,Glambda_eta,Gchi_eta,PhiMat,Umat)
      !if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble(XiQS)
      !if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble((0d0,-1d0)*XiQS)

    !! XiQS_link
      !call calc_Xilink_QS(XiQS,Geta_lambda,Glambda_lambda,Gchi_lambda,PhiMat,Umat)
      !if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble(XiQS)
      !if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble((0d0,-1d0)*XiQS)

    !! XiQS_face
      !call calc_Xiface_QS(XiQS,Geta_chi,Glambda_chi,Gchi_chi,PhiMat,Umat)
      !if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble(XiQS)
      !if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble((0d0,-1d0)*XiQS)

!      !! divergence of U(1)V current
!      call calc_DJ_U1V(DJ1,DJ2,Glambda_eta,Glambda_chi,Umat)
!      tmp=(0d0,0d0)
!      do j=1,NMAT
!        do i=1,NMAT
!          tmp=tmp+PhiMat(i,j,1)*dconjg(PhiMat(i,j,1))
!        enddo
!      enddo
!      do gf=1,global_num_faces
!        rank=local_face_of_global(gf)%rank_
!        lf=local_face_of_global(gf)%label_
!        tag=gf
!        if( MYRANK == 0 ) then 
!          if( MYRANK == rank ) then
!            tmp1 = DJ1(lf)
!            tmp2 = DJ2(lf)
!          else
!            call MPI_RECV(tmp1,1,MPI_DOUBLE_COMPLEX,rank,2*tag-1,MPI_COMM_WORLD,ISTATUS,IERR)
!            call MPI_RECV(tmp2,1,MPI_DOUBLE_COMPLEX,rank,2*tag,MPI_COMM_WORLD,ISTATUS,IERR)
!          endif
!        else
!          if( MYRANK == rank ) then
!            call MPI_SEND(DJ1(lf),1,MPI_DOUBLE_COMPLEX,0,2*tag-1,MPI_COMM_WORLD,IERR)
!            call MPI_SEND(DJ2(lf),1,MPI_DOUBLE_COMPLEX,0,2*tag,MPI_COMM_WORLD,IERR)
!          endif
!        endif
!        tmp=(1d0,0d0)
!        if( MYRANK == 0 ) then 
!          write(*,'(E15.8,2X)',advance='no')  dble(tmp*tmp1)
!          write(*,'(E15.8,2X)',advance='no')  dble((0d0,-1d0)*tmp*tmp1)
!          write(*,'(E15.8,2X)',advance='no')  dble(tmp*tmp2)
!          write(*,'(E15.8,2X)',advance='no')  dble((0d0,-1d0)*tmp*tmp2)
!        endif
!      enddo

    ! WT with Atr site
      tmp_obs=WT_site*cdabs(Acomp_tr)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble(tmp_obs)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble((0d0,-1d0)*tmp_obs)

    ! WT with Atr link
      tmp_obs=WT_link*cdabs(Acomp_tr)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble(tmp_obs)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble((0d0,-1d0)*tmp_obs)

    ! WT with Atr face
      tmp_obs=WT_face*cdabs(Acomp_tr)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble(tmp_obs)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble((0d0,-1d0)*tmp_obs)

    ! WT with Aface site
      tmp_obs=WT_site*cdabs(Acomp_face) + CSF_site*dconjg(Acomp_face)/cdabs(Acomp_face)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble(tmp_obs)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble((0d0,-1d0)*tmp_obs)

    ! WT with Aface link
      tmp_obs=WT_link*cdabs(Acomp_face) + CSF_link*dconjg(Acomp_face)/cdabs(Acomp_face)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble(tmp_obs)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble((0d0,-1d0)*tmp_obs)

    ! WT with Aface face
      tmp_obs=WT_face*cdabs(Acomp_face) + CSF_link*dconjg(Acomp_face)/cdabs(Acomp_face)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble(tmp_obs)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble((0d0,-1d0)*tmp_obs)

    if(MYRANK==0) write(*,*)
  else
    exit
  endif
enddo

if( MYRANK == 0 ) then
  close(N_MEDFILE)
  close(N_DinvFILE)
endif
end program main


