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
integer, parameter :: num_calcobs=42 ! 考えているobservableの数
character(128) :: name_obs(1:num_calcobs) = (/ &
  "Re(phase Pf)",  &
  "Im(phase Pf)", &
  "|Atr|", &
  "Re(phase Atr)", &
  "Im(phase Atr)", &
  "|Aface|", &
  "Re(phase Aface)", &
  "Im(phase Aface)", &
  "|Aphibar|", &
  "Re(phase Aphibar)", &
  "Im(phase Aphibar)", &
  "Re(Q(Aphibar)Xi)", &
  "Im(Q(Aphibar)Xi)", &
  "SbS", &
  "SbL", &
  "SbF", &
  "dimG*(NS+NL)/2", &
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
  "Tr|phi^2|", &
  "Re(exactSf_link1)", &
  "Im(exactSf_link1)", &
  "Re(exactSf_link2)", &
  "Im(exactSf_link2)", &
  "Re(exactSf_face2)", &
  "Im(exactSf_face2)"  &
  /)
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
complex(kind(0d0)) :: Acomp_tr ! trace compensator
complex(kind(0d0)) :: Acomp_face ! face compensator
complex(kind(0d0)) :: Acomp_phibar ! phia compensator
!complex(kind(0d0)) :: Atr_phase ! A*/|A|
!complex(kind(0d0)) :: Aface_phase ! A*/|A|
!complex(kind(0d0)) :: Aphibar_phase ! A*/|A|
complex(kind(0d0)) :: QC_Xi ! Q(Aphibar).\Xi 
complex(kind(0d0)) :: min_eigen
complex(kind(0d0)) :: mass_cont
complex(kind(0d0)) :: WT_site
complex(kind(0d0)) :: WT_link
complex(kind(0d0)) :: WT_face
complex(kind(0d0)) :: Sf1, Sf2, Sf3, Sf4, Sf5
complex(kind(0d0)) :: XiQS
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
integer :: l,ll,s,ls,tag,rank,i,j, ite, ite2
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
complex(kind(0d0)) :: phase_pf
complex(kind(0d0)), allocatable :: eigenvals(:)
complex(kind(0d0)) :: tmp
integer :: k,s1,s2,l1,l2,f1,f2,a,b
integer :: ios

iarg=iargc()
if( iarg ==0 ) then
  !if (MYRANK==0) write(*,*) "use as a.out [MEDFILE] [DinvFILE]"
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
num_sitelink=(global_num_sites+global_num_links)*(NMAT*NMAT-1)


allocate( Dinv(1:num_fermion, 1:num_fermion) )
allocate( Dirac(1:num_fermion, 1:num_fermion) )
allocate( Dirac2(1:num_fermion, 1:num_fermion) )
allocate( eigenvals(1:num_fermion) )

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
    open(N_DinvFILE, file=DinvFILE, status='OLD',action='READ')
  else
    open(N_DinvFILE, file=DinvFILE, status='REPLACE')
    open(N_EigenFILE, file=EigenFILE, status='REPLACE')
  endif
endif 

      
do    
  call read_config_from_medfile(Umat,PhiMat,ite,N_MEDFILE,control)

  if( exist_dinv==0 ) then 
    read(N_DinvFILE,'(I10,2X)',advance='no',iostat=ios) ite2
    if( ios == -1) exit
    do j=1,num_fermion
      do i=1,num_fermion
        read(N_DinvFILE,'(E23.15,2X,E23.15,2X)',advance='no') &
          rtmp,ctmp
          Dinv(i,j)=rtmp+(0d0,1d0)*ctmp
      enddo
    enddo
    read(N_DinvFILE,'(E23.15,2X,E23.15,2X)') re_phase, im_phase
    phase_pf=dcmplx(re_phase)+(0d0,1d0)*dcmplx(im_phase)





    call MPI_BCAST(ite,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ite2,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
    if( ite .ne. ite2 ) then 
      write(*,*) "MEDFILE and Dinv do not match."
      call stop_for_test
    endif
  else
    call construct_Dirac(Dirac,Umat,PhiMat) !correct!

    !call writeout_Dirac(Dirac)
!    if(MYRANK==0) then
!      write(*,'(i9)') ite
!      do i=1,size(Dirac,1)
!        do j=1,size(Dirac,2)
!          if( cdabs(Dirac(i,j)) > 1d-8) write(*,'(i5,2x,i5,2x,E15.8,2X,E15.8,2X)')i,j, dble(Dirac(i,j)),dble((0d0,-1d0)*Dirac(i,j))
!        enddo
!      enddo
!    endif


    
    !!
    Dirac2=Dirac
    if( MYRANK==0 ) call CalcPfaffian(abs_pf,phase_pf,Dirac2)
    !!
    Dinv=Dirac
    if( MYRANK==0 ) call matrix_inverse(Dinv)
    !!
    if( MYRANK==0 ) call matrix_eigenvalues(eigenvals,Dirac)

    write(N_DinvFILE,'(I10,2X)',advance='no') ite
    do j=1,num_fermion
      do i=1,num_fermion
        write(N_DinvFILE,'(E23.15,2X,E23.15,2X)',advance='no') &
          dble(Dinv(i,j)), dble(Dinv(i,j)*(0d0,-1d0))
      enddo
    enddo
    write(N_DinvFILE,'(E23.15,2X,E23.15,2X)',advance='no') &
      dble(phase_pf), dble((0d0,-1d0)*phase_pf)
    !! eigenvalues
    do i=1,num_fermion
      write(N_EigenFILE,'(E23.15,2X,E23.15,2X)',advance='no') &
        dble(eigenvals(i)), dble(eigenvals(i)*(0d0,-1d0))
    enddo
  endif

  call make_fermion_correlation_from_Dinv(&
    Geta_eta, Glambda_eta, Gchi_eta, &
    Geta_lambda, Glambda_lambda, Gchi_lambda, &
    Geta_chi, Glambda_chi, Gchi_chi, &
    Dinv)

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

    Sb_computed=0

    !! Pfaffian phase
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble(phase_pf)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble( (0d0,-1d0)*phase_pf)
      !if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble((0d0,-1d0)*Sf)

    !"|Atr|", &
      call calc_trace_compensator(Acomp_tr,PhiMat)
      !call calc_VM_compensator(Acomp_VM,PhiMat)
      !APQ_phase = dconjg(Acomp_tr) / cdabs(Acomp_tr) 
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') cdabs(Acomp_tr)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble(Acomp_tr/cdabs(Acomp_tr))
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble( (0d0,-1d0)*Acomp_tr/cdabs(Acomp_tr))

    !"|Aface|", &
      call calc_face_compensator(Acomp_face,Umat,PhiMat,Geta_chi)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') cdabs(Acomp_face)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble(Acomp_face/cdabs(Acomp_face))
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble( (0d0,-1d0)*Acomp_face/cdabs(Acomp_face))

    !"|Aphibar|", &
      call calc_phibar_compensator(Acomp_phibar,PhiMat)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') cdabs(Acomp_phibar)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble(Acomp_phibar/cdabs(Acomp_phibar))
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') dble( (0d0,-1d0)*Acomp_phibar/cdabs(Acomp_phibar))

    !"|Aphibar|", &
      call calc_QC_Xi(QC_Xi,Geta_eta,Geta_lambda,Geta_chi,Umat,PhiMat)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  &
        dble(QC_Xi / (Acomp_phibar/dcmplx(cdabs(Acomp_phibar))))
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  &
        dble((0d0,-1d0) * QC_Xi / (Acomp_phibar/dcmplx(cdabs(Acomp_phibar))))

    !! SbS
      call calc_bosonic_action_site(SbS,PhiMat)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  SbS

    !! SbL
      call calc_bosonic_action_link(SbL,Umat,PhiMat)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  SbL

    !! SbF
      call calc_bosonic_action_face(SbF,Umat)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  SbF

    !! mass_cont
!      call mass_contribution(mass_cont,Geta_eta,Geta_lambda,Geta_chi,Umat,PhiMat)
!      if( MYRANK == 0 ) write(*,'(E15.8,2X,E15.8,2X)',advance='no')  &
!        dble(mass_cont), &
!        dble((0d0,-1d0)*mass_cont)

    !! SfL2
!      call calc_Sf_link2(SfL2,PhiMat,Umat,Glambda_lambda)
!      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') &
!        dble(SfL2)
!      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') &
!        dble( (0d0,-1d0)*SfL2 )

    !! fermion number in trivial WT
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') &
        0.5d0*dble( (NMAT*NMAT-1)*(global_num_sites+global_num_links) ) 

    !! fermion number in (5.8)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  &
        dble(global_num_faces*dimG)*0.5d0 

!    !! trivial WT for site
!      call calc_siteWT(WT_site, Geta_eta,PhiMat)
!      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble(WT_site)
!
!    !! trivial WT for link
!      call calc_linkWT(WT_link,Glambda_eta,Geta_lambda,Glambda_lambda,PhiMat,Umat)
!      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble(WT_link)

!    !! trivial WT for face
!      call calc_faceWT(WT_face,Glambda_chi,Gchi_eta,Gchi_chi,PhiMat,Umat)
!      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble(WT_face)

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
      call mass_contribution_site(mass_cont,Geta_eta,PhiMat)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble(mass_cont)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble((0d0,-1d0)*mass_cont)

    !! mass contribution in WT_link
      call  mass_contribution_link(mass_cont,Glambda_eta,Umat,PhiMat)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble(mass_cont)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble((0d0,-1d0)*mass_cont)

    !! mass contribution in WT_face
      call mass_contribution_face(mass_cont,Gchi_eta,Umat,PhiMat)
       if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble(mass_cont)
       if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble((0d0,-1d0)*mass_cont)

    !"tr|\phi^2|", &
      call calc_trphi2(trphi2,PhiMat)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') trphi2

    !! Sf_link1
      call calc_Qexact_Sf_link1(Sf2,Geta_lambda,Umat,PhiMat)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble(Sf2)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble((0d0,-1d0)*Sf2)

    !! Sf_link2
      call calc_Qexact_Sf_link2(Sf3, PhiMat, Umat, Glambda_lambda)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble(Sf3)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble((0d0,-1d0)*Sf3)

    !! Sf_face2
      call calc_Qexact_Sf_face2(Sf5,Glambda_chi,Umat)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble(Sf5)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  dble((0d0,-1d0)*Sf5)
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


