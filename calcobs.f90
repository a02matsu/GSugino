module global_calcobs
implicit none

character(128), parameter :: PARAFILE="parameters_calcobs.dat"
character(128) :: MEDFILE
character(128) :: DinvFILE
integer, parameter :: num_calcobs=6 ! 考えているobservableの数
character(128) :: name_obs(1:num_calcobs) = (/ &
  "|Atr|", &
  "Re(triv.WT)", &
  "Im(triv.WT)", &
  "SbS+SbL-rhs", &
  "SbL-2SbF-Re(SfL2)", &
  "Im(SfL2)" &
  /)
!integer :: trig_obs(1:num_calcobs)
integer :: sizeM,sizeN

double precision :: Sb, SbS, SbL, SbF
complex(kind(0d0)) :: SfL2
complex(kind(0d0)) :: Acomp_tr ! trace compensator
complex(kind(0d0)) :: Acomp_VM ! van der monde compensator
complex(kind(0d0)) :: APQ_phase ! A*/|A|
complex(kind(0d0)) :: min_eigen
!complex(kind(0d0)), allocatable :: WT1(:)
!complex(kind(0d0)), allocatable :: WT2(:)
complex(kind(0d0)) :: WT1, WT2
integer :: num_fermion ! total fermion number
integer :: num_sitelink ! total fermion number

integer, parameter :: N_MEDFILE=100
integer, parameter :: N_PARAFILE=101
integer, parameter :: N_DinvFILE=102

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

end module global_calcobs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
use global_parameters
use global_calcobs
use initialization_calcobs
use simulation
use parallel
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: l,ll,s,ls,tag,rank,i,j, ite
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

iarg=iargc()
if( iarg <=1 ) then
  if (MYRANK==0) write(*,*) "use as a.out [MEDFILE] [DinvFILE]"
  stop
endif
call getarg(1,MEDFILE)
call getarg(2,DinvFILE)
INPUT_FILE_NAME="inputfile"
  
call initialization 
num_fermion=(global_num_sites+global_num_links+global_num_faces)*(NMAT*NMAT-1)
num_sitelink=(global_num_sites+global_num_links)*(NMAT*NMAT-1)


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
  open(N_DinvFILE, file=DinvFILE, status='OLD',action='READ')
endif 
      
do    
  call read_config_from_medfile(Umat,PhiMat,ite,N_MEDFILE,control)
  call read_Dinv(ite, Geta_eta, Glambda_eta, Gchi_eta, &
         Geta_lambda, Glambda_lambda, Gchi_lambda, &
         Geta_chi, Glambda_chi, Gchi_chi, &
         N_DinvFILE)
  if( control == 0 ) then 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! calculate observables
    if( MYRANK == 0 ) then
      write(*,'(I7,2X)',advance='no') ite
    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! test for fermionic operators 
    !call check_DinvPF(&
    !    Geta_eta, Glambda_eta, Gchi_eta, &
    !    Geta_lambda, Glambda_lambda, Gchi_lambda, &
    !    Geta_chi, Glambda_chi, Gchi_chi, &
    !    Umat,PhiMat)

    !"|Atr|", &
      call calc_trace_compensator(Acomp_tr,PhiMat)
      !call calc_VM_compensator(Acomp_VM,PhiMat)
      APQ_phase = dconjg(Acomp_tr) / cdabs(Acomp_tr) 
      Sb_computed=0
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') cdabs(Acomp_tr)
    !"Sb", &
      !call calc_bosonic_action(Sb,Umat,PhiMat)
      call calc_bosonic_action_site(SbS,PhiMat)
      call calc_bosonic_action_link(SbL,Umat,PhiMat)
      call calc_bosonic_action_face(SbF,Umat)
      if( MYRANK == 0 ) Sb=SbS+SbL+SbF
      !if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') Sb
    !! trivial WT
      !call calc_XiPhiEta(XiPhiEta,&
      !  Geta_eta, Glambda_eta, Gchi_eta, &
      !  Geta_lambda, Glambda_lambda, Gchi_lambda, &
      !  Geta_chi, Glambda_chi, Gchi_chi, &
      !  Umat,PhiMat)
      !tmpobs1= &
      !  cdabs(Acomp_tr) * &
      !  (&
      !  dcmplx(Sb) &
      !  + dcmplx(0.5d0*mass_square_phi)*XiPhiEta &
      !  - dcmplx(0.5d0*dble(num_sitelink)) &
      !  )  
      call calc_trivialWT(tmpobs1,Geta_eta,Geta_lambda,Geta_chi,Umat,PhiMat)
      if( MYRANK == 0 ) write(*,'(E15.8,2X,E15.8,2X)',advance='no') &
        dble(tmpobs1), dble((0d0,-1d0)*tmpobs1)

    !! relations in bosonic action 1
    !!  Sb_site + Sb_face = #(face)/2*dimG
    !! we need also compensator
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no')  &
        SbS + SbF - dble(global_num_faces*dimG)*0.5d0 

    !! relations in bosonic action 2
    !!  Sb_link - 2 Sb_face = lambda.lambda(U.\bar{phi}.U^-1 + \bar{phi})
    !! we need also compensator
      call calc_Sf_link2(SfL2,PhiMat,Umat,Glambda_lambda)
      !tmpobs1 = &
        !dcmplx(cdabs(Acomp_tr)) * (dcmplx(SbL - 2d0*SbF) - SfL2)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') &
        SbL - 2d0*SbF - dble(SfL2)
      if( MYRANK == 0 ) write(*,'(E15.8,2X)',advance='no') &
        dble( (0d0,-1d0)*SfL2 )

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


