!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ver.01: correct the normalization of the PCSC relation 
!! ver.04: bug fix for mass part of PCSC
!! ver.05: include WT id. in naive quench
!! ver.06: added compensator for SU(2) (triple cover version)
!module observables
!use global_parameters
!use global_subroutines
!implicit none
!
!contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_bosonic_action(Sb,UMAT,PhiMat)
!use hamiltonian
implicit none

double precision, intent(out) :: Sb
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

double precision :: SB_S,SB_L,SB_F,SB_local

SB_S=0d0
SB_L=0d0
SB_F=0d0
SB=0d0
call bosonic_action_site(SB_S,PhiMat)
call bosonic_action_link(SB_L,UMAT,PhiMat)
call bosonic_action_face(SB_F,UMAT)

Sb_local=SB_S+SB_L+SB_F

#ifdef PARALLEL
call MPI_REDUCE(SB_local,SB,1,MPI_DOUBLE_PRECISION, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)
#else
SB=SB_local
#endif

end subroutine calc_bosonic_action

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_TrX2(TrX2, PhiMat)
implicit none

double precision, intent(out) :: TrX2
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
double precision :: tmp

integer :: s,i,j

tmp=0d0
do s=1,num_sites
  do i=1,NMAT
    do j=1,NMAT
      tmp=tmp+ dble( PhiMat(i,j,s)*conjg( PhiMat(i,j,s) ) )
    enddo
  enddo
enddo
#ifdef PARALLEL
call MPI_REDUCE(tmp,TrX2,1,MPI_DOUBLE_PRECISION, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)
#else
TrX2=tmp
#endif

TrX2 = TrX2 / (2d0 * dble(global_num_sites) )

end subroutine calc_TrX2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_trace_compensator(Acomp,PhiMat)
use parallel
implicit none

complex(kind(0d0)), intent(out) :: Acomp
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) tmp,A_tmp
double precision :: radius, phase, ratio
integer :: s,i,j,eular

eular=global_num_sites-global_num_links+global_num_faces
ratio=dble(-(NMAT*NMAT-1)*eular)/4d0
Acomp=(0d0,0d0)
A_tmp=(0d0,0d0)
do s=1,num_sites
  tmp=(0d0,0d0)
  do i=1,NMAT
    do j=1,NMAT
      tmp=tmp+PhiMat(i,j,s)*dconjg(PhiMat(i,j,s))
    enddo
  enddo
  tmp=(tmp/dcmplx(dble(NMAT)))
  radius=cdabs(tmp)
  phase=atan2(dble(tmp),dble(tmp*(0d0,-1d0)))

  A_tmp=A_tmp + dcmplx(radius**ratio) * cdexp( (0d0,1d0)*dcmplx(phase*ratio) )
enddo

call MPI_REDUCE(A_tmp,Acomp,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)
  
Acomp=Acomp/dcmplx(dble(global_num_sites))

end subroutine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_eigenvalues_Dirac(eigenvalues,UMAT,PhiMat)
use Dirac_operator, only : make_Dirac
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
!complex(kind(0d0)), intent(inout) :: eigenvalues(1:sizeD)
complex(kind(0d0)), intent(inout) :: eigenvalues(:)

!complex(kind(0d0)) :: Dirac(1:sizeD,1:sizeD)
!complex(kind(0d0)) VL(1,1:sizeD), VR(1,1:sizeD)
!double precision RWORK(1:2*sizeD)
!complex(kind(0d0)) WORK(2*(sizeD)),tako1
complex(kind(0d0)), allocatable :: Dirac(:,:)
complex(kind(0d0)), allocatable :: VL(:,:), VR(:,:)
double precision, allocatable :: RWORK(:)
complex(kind(0d0)), allocatable :: WORK(:)
complex(kind(0d0)) tako1
character JOBVL,JOBVR
integer info,lwork
integer i,j
integer :: sizeD

sizeD=size(eigenvalues,1)

allocate( Dirac(1:sizeD,1:sizeD) )
allocate( VL(1,1:sizeD), VR(1,1:sizeD) )
allocate( RWORK(1:2*sizeD) )
allocate( WORK(2*(sizeD)) )

lwork=2*sizeD
JOBVL='N'
JOBVR='N'

call make_Dirac(Dirac,UMAT,PhiMat)
!if( MYRANK==0 ) then 
!do i=1,sizeD
!  do j=1,sizeD
!    write(*,*) i,j,Dirac(i,j)
!  enddo
!enddo
!endif

#ifdef PARALLEL
if( MYRANK==0 ) then 
#endif
call ZGEEV(JOBVL,JOBVR,sizeD,&
     DIRAC,sizeD,eigenvalues,VL,1,VR,1,WORK,lwork,RWORK,info)
#ifdef PARALLEL
endif
#endif

! sort the eigenvalues
do i=1,sizeD
 do j=i+1,sizeD
  tako1 = eigenvalues(i)
  if(abs(eigenvalues(j)).LT.abs(eigenvalues(i))) then 
    eigenvalues(i) = eigenvalues(j)
    eigenvalues(j) = tako1
  endif
 enddo
enddo

end subroutine calc_eigenvalues_Dirac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_eigenvalues_DdagD(eigenvalues,UMAT,PhiMat)
use Dirac_operator, only : make_DdagD
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
!complex(kind(0d0)), intent(inout) :: eigenvalues(1:sizeD)
complex(kind(0d0)), intent(inout) :: eigenvalues(:)

!complex(kind(0d0)) :: Dirac(1:sizeD,1:sizeD)
!complex(kind(0d0)) VL(1,1:sizeD), VR(1,1:sizeD)
!double precision RWORK(1:2*sizeD)
!complex(kind(0d0)) WORK(2*(sizeD)),tako1
complex(kind(0d0)), allocatable :: DdagD(:,:)
complex(kind(0d0)), allocatable :: VL(:,:), VR(:,:)
double precision, allocatable :: RWORK(:)
complex(kind(0d0)), allocatable :: WORK(:)
complex(kind(0d0)) tako1
character JOBVL,JOBVR
integer info,lwork
integer i,j
integer :: sizeD

sizeD=size(eigenvalues,1)

allocate( DdagD(1:sizeD,1:sizeD) )
allocate( VL(1,1:sizeD), VR(1,1:sizeD) )
allocate( RWORK(1:2*sizeD) )
allocate( WORK(2*(sizeD)) )

lwork=2*sizeD
JOBVL='N'
JOBVR='N'

call make_DdagD(DdagD,UMAT,PhiMat)
!if( MYRANK==0 ) then 
!do i=1,sizeD
!  do j=1,sizeD
!    write(*,*) i,j,DdagD(i,j)
!  enddo
!enddo
!endif

#ifdef PARALLEL
if( MYRANK==0 ) then 
#endif
call ZGEEV(JOBVL,JOBVR,sizeD,&
     DdagD,sizeD,eigenvalues,VL,1,VR,1,WORK,lwork,RWORK,info)
#ifdef PARALLEL
endif
#endif

! sort the eigenvalues
do i=1,sizeD
 do j=i+1,sizeD
  tako1 = eigenvalues(i)
  if(abs(eigenvalues(j)).LT.abs(eigenvalues(i))) then 
    eigenvalues(i) = eigenvalues(j)
    eigenvalues(j) = tako1
  endif
 enddo
enddo

end subroutine calc_eigenvalues_DdagD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine calc_eigenvalues_Dirac2(eigenvalues,Dirac)
!implicit none
!
!complex(kind(0d0)), intent(inout) :: eigenvalues(1:sizeD)
!complex(kind(0d0)), intent(inout) :: Dirac(1:sizeD,1:sizeD)
!
!complex(kind(0d0)) VL(1,1:sizeD), VR(1,1:sizeD)
!double precision RWORK(1:2*sizeD)
!complex(kind(0d0)) WORK(2*(sizeD)),tako1
!character JOBVL,JOBVR
!integer info,lwork
!integer i,j
!
!lwork=2*sizeD
!JOBVL='N'
!JOBVR='N'
!
!call ZGEEV(JOBVL,JOBVR,sizeD,&
!     DIRAC,sizeD,eigenvalues,VL,1,VR,1,WORK,lwork,RWORK,info)
!
!! sort the eigenvalues
!do i=1,sizeD
! do j=i+1,sizeD
!  tako1 = eigenvalues(i)
!  if(abs(eigenvalues(j)).LT.abs(eigenvalues(i))) then 
!    eigenvalues(i) = eigenvalues(j)
!    eigenvalues(j) = tako1
!  endif
! enddo
!enddo
!
!end subroutine calc_eigenvalues_Dirac2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_min_and_max_of_eigenvalues_Dirac(minimal,maximal,UMAT,PhiMat)
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(inout) :: minimal, maximal

!complex(kind(0d0)) :: eigenvalues(1:sizeD)
complex(kind(0d0)), allocatable :: eigenvalues(:)
integer :: sizeD

#ifdef PARALLEL
sizeD=dimG*(global_num_sites+global_num_links+global_num_faces)
#else
sizeD=dimG*(num_sites+num_links+num_faces)
#endif
allocate(eigenvalues(1:sizeD))

call calc_eigenvalues_Dirac(eigenvalues,UMAT,PhiMat)

minimal=eigenvalues(1)
maximal=eigenvalues(sizeD)

end subroutine calc_min_and_max_of_eigenvalues_Dirac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compensator = 1/(NMAT*Nsite) * \sum_s Tr( \Phi^{ -(N^2-1)*\chi_h / 2 } )
!subroutine calc_compensator_naive(C,Phi)
!use SUN_generators, only : make_traceless_matrix_from_modes
!use matrix_functions, only : matrix_power, matrix_inverse
!implicit none
!
!complex(kind(0d0)), intent(out) :: C
!complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
!complex(kind(0d0)) :: tmpPhi(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: pmat(1:NMAT,1:NMAT)
!!! data of simplicial complex
!integer :: Euler,s,i
!
!Euler=num_sites - num_links + num_faces
!C=(0d0,0d0)
!do s=1,num_sites
!  call make_traceless_matrix_from_modes(tmpPhi,NMAT,Phi(:,s))
!  if (Euler > 0) then 
!    call matrix_inverse(tmpPhi)
!  endif
!  call matrix_power(pmat,tmpPhi,abs((NMAT*NMAT-1)*Euler)/2)
!  do i=1,NMAT
!    C=C+pmat(i,i)
!  enddo
!enddo
!C=C/ dble(num_sites*NMAT)
!
!end subroutine calc_compensator_naive
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! compensator_trace = 1/(NMAT*Nsite) * \sum_s Tr( \Phi^2 )^{ -(N^2-1)*\chi_h/4 } ) )
!subroutine calc_compensator_trace(C,PhiMat)
!#ifdef PARALLEL
!use parallel
!#endif
!implicit none
!
!complex(kind(0d0)), intent(out) :: C
!complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
!complex(kind(0d0)) :: trace
!real(8) :: trace_abs, trace_arg, num
!!! data of simplicial complex
!integer :: s,i,j
!complex(kind(0d0)) :: tmp
!
!num=dble( (NMAT*NMAT-1)*(global_num_sites-global_num_links+global_num_faces) )/4d0
!
!C=(0d0,0d0)
!tmp=(0d0,0d0)
!do s=1,num_sites
!  trace=(0d0,0d0)
!  do i=1,NMAT
!    do j=1,NMAT
!      trace=trace+PhiMat(i,j,s)*PhiMat(j,i,s)
!    enddo
!  enddo
!  trace=trace/cmplx(dble(NMAT))
!  !trace_abs = sqrt( trace * conjg(trace) ) 
!  trace_abs = cdabs(trace)
!  trace_arg = arg(trace)
!  tmp=tmp+cmplx(trace_abs**(-num)) &
!      * exp( (0d0,-1d0)*trace_arg*num )
!      !* ( cmplx( cos( trace_arg * num ) ) - (0d0,1d0)*cmplx( sin( trace_arg * num ) ) )
!enddo
!#ifdef PARALLEL
!call MPI_REDUCE(tmp,C,1,MPI_DOUBLE_COMPLEX, &
!  MPI_SUM,0,MPI_COMM_WORLD,IERR)
!#else
!C=tmp
!#endif
!C=C/ cmplx(dble(global_num_sites))
!
!end subroutine calc_compensator_trace
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! compensator_trace 
!!! = 1/(NMAT*Nsite) 
!!!   * [ \sum_s Tr( \Phi^2 ) ]^{ -(N^2-1)*\chi_h/4 } ) 
!subroutine calc_compensator2_trace(C,Phi)
!use SUN_generators, only : make_traceless_matrix_from_modes
!use matrix_functions, only : matrix_power, matrix_inverse
!implicit none
!
!complex(kind(0d0)), intent(out) :: C
!complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
!complex(kind(0d0)) :: trace
!real(8) :: C_abs, C_arg, num
!!! data of simplicial complex
!integer :: s,i,a
!
!num=dble( (NMAT*NMAT-1)*(num_sites-num_links+num_faces) )/4d0
!C=(0d0,0d0)
!do s=1,num_sites
!  trace=(0d0,0d0)
!  do a=1,dimG
!    trace=trace+Phi(a,s)*Phi(a,s)
!  enddo
!  C=C+trace/cmplx(dble(NMAT))
!  !trace_abs = abs(trace)
!  !trace_arg = arg(trace)
!  !C=C+cmplx(trace_abs**(-num)) &
!      !* exp( (0d0,-1d0)*trace_arg*num )
!      !* ( cmplx( cos( trace_arg * num ) ) - (0d0,1d0)*cmplx( sin( trace_arg * num ) ) )
!enddo
!C = C/cmplx(dble(num_sites))
!C_abs = abs(C)
!C_arg = arg(C)
!!write(*,*) C_abs, num, C_abs**(-num)
!C = cmplx( C_abs**(-num) ) &
!    * exp( (0d0,-1d0) * cmplx( C_arg * num ) ) 
!
!end subroutine calc_compensator2_trace
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! compensator25 = 1/(NMAT*Nsite) * \sum_s Tr( \Phi^2 )^{ -(N^2-1)*\chi_h/4 } ) )
!!!  where the branch of sqrt is chosen by randum number 
!subroutine calc_compensator25(C,Phi)
!use SUN_generators, only : make_traceless_matrix_from_modes
!use matrix_functions, only : matrix_power, matrix_inverse
!implicit none
!
!complex(kind(0d0)), intent(out) :: C
!complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
!complex(kind(0d0)) :: Phimat(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: pmat(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: trace
!real(8) :: trace_abs, trace_arg, num
!real(8) :: rand, argument
!!! data of simplicial complex
!integer :: s,i
!
!num=dble( (NMAT*NMAT-1)*(num_sites-num_links+num_faces) )/4d0
!C=(0d0,0d0)
!do s=1,num_sites
!  trace=(0d0,0d0)
!  call make_traceless_matrix_from_modes(Phimat,NMAT,Phi(:,s))
!  call matrix_power(pmat,Phimat,2)
!  do i=1,NMAT
!    trace=trace+pmat(i,i)
!  enddo
!  trace_abs = abs(trace)
!  trace_arg = arg(trace)
!  call random_number(rand)
!  if (rand < 0.5d0) then 
!      argument = trace_arg*num
!  else
!      argument = (2d0*acos(-1d0) + trace_arg)*num
!  endif
!  C=C+cmplx(trace_abs**(-num)) &
!      * ( cmplx( cos( argument ) ) - (0d0,1d0)*cmplx( sin( argument ) ) )
!enddo
!C=C/ cmplx(dble(num_sites*NMAT))
!
!end subroutine calc_compensator25
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! compensator_det 
!!!  = 1/Nsite ¥sum_s ( Det(¥Phi_s) )^{ -(N^2-1)*¥chi_h / 2N }
!subroutine calc_compensator_det(C,Phi)
!use SUN_generators, only : make_traceless_matrix_from_modes
!use matrix_functions, only : matrix_determinant
!implicit none
!
!complex(kind(0d0)), intent(out) :: C
!complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
!complex(kind(0d0)) :: matPhi(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: det_arg, det
!real(8) :: det_abs, num,logdet
!!! data of simplicial complex
!integer :: s,i
!
!num = dble( (NMAT*NMAT-1) * (num_sites - num_links + num_faces) ) / dble( 2*NMAT )
!C=(0d0,0d0)
!do s=1,num_sites
!  call make_traceless_matrix_from_modes(matPhi,NMAT,Phi(:,s))
!  call matrix_determinant(logdet,det_arg,matPhi)
!  det_abs=dexp(logdet)
!  det_arg=arg(det)
!  C = C + cmplx( det_abs**(-num) ) &
!      * ( cmplx( cos( det_arg*num ) ) - (0d0,1d0)*cmplx( sin( det_arg*num ) ) )
!enddo
!C=C / dble(num_sites)
!end subroutine calc_compensator_det
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! compensator_det 
!!!  = [ 1/Nsite ¥sum_s Det(¥Phi_s) ]^{ -(N^2-1)*¥chi_h / 2N }
!subroutine calc_compensator2_det(C,Phi)
!use SUN_generators, only : make_traceless_matrix_from_modes
!use matrix_functions, only : matrix_determinant
!implicit none
!
!complex(kind(0d0)), intent(out) :: C
!complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
!complex(kind(0d0)) :: matPhi(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: argdet
!real(8) :: C_abs, C_arg, num, logdet
!!! data of simplicial complex
!integer :: s,i
!
!num = dble( (NMAT*NMAT-1) * (num_sites - num_links + num_faces) ) / dble( 2*NMAT )
!C=(0d0,0d0)
!do s=1,num_sites
!  call make_traceless_matrix_from_modes(matPhi,NMAT,Phi(:,s))
!  call matrix_determinant(logdet,argdet,matPhi)
!  C=C+dexp(logdet)*argdet
!enddo
!C = C / cmplx( dble( num_sites ) )
!C_abs=abs(C)
!C_arg=arg(C)
!C = cmplx( C_abs ** (-num) ) &
!    * exp( (0d0,-1d0) * cmplx( C_arg * num ) )
!end subroutine calc_compensator2_det


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! site part of the PCSC
!subroutine calc_pcsc_site(obs,Phi,Dirac_inv)
!use SUN_generators, only : make_traceless_matrix_from_modes, trace_MTa
!use matrix_functions, only : matrix_commutator
!implicit none
!
!complex(kind(0d0)), intent(out) :: obs
!complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
!complex(kind(0d0)), intent(in) :: Dirac_inv(:,:)
!complex(kind(0d0)) :: phi_mat(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: comm(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: comm_a
!
!integer :: s,t,a,b
!integer :: sizeD
!
!sizeD=size(Dirac_inv,1)
!
!!obs=(0d0,0d0)
!!do s=1,sizeD
!  !do t=1,sizeD
!    !obs=obs+(Dirac_inv(s,t)+Dirac_inv(t,s))&
!      !*conjg( Dirac_inv(s,t)+Dirac_inv(t,s)) 
!  !enddo
!!enddo
!!write(*,*) obs
!
!
!obs=(0d0,0d0)
!do s=1,num_sites
!  call make_traceless_matrix_from_modes(phi_mat, NMAT, Phi(:,s))
!  call matrix_commutator(comm, phi_mat, phi_mat, 'N', 'C')
!  do a=1,dimG
!    call trace_MTa(comm_a, comm, a, NMAT)
!    do t=1,num_sites
!      do b=1,dimG
!        obs = obs &
!          + (0.25d0, 0d0) * cmplx(alpha_s(s) * overall_factor**2 * mass_square_phi) &
!            * comm_a * Phi(b,t) &
!            * Dirac_inv(site_index(a,s,NMAT), site_index(b,t,NMAT)) 
!      enddo
!    enddo
!  enddo
!enddo
!
!end subroutine calc_pcsc_site
!        
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! link part of the PCSC
!subroutine calc_pcsc_link(obs,Phi,UMAT,Dirac_inv)
!use SUN_generators, only : make_traceless_matrix_from_modes, trace_MTa
!!use matrix_functions, only : matrix_product
!!use Dirac_operator, only : make_Dirac
!implicit none
!
!complex(kind(0d0)), intent(out) :: obs
!complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
!complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)), intent(in) :: Dirac_inv(:,:)
!complex(kind(0d0)) :: bPhi(1:dimG,1:num_sites)
!complex(kind(0d0)) :: dphi_mat(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: dphi_a
!
!!complex(kind(0d0)) :: Dirac(1:sizeD,1:sizeD), prod(1:sizeD,1:sizeD)
!
!integer :: s,l,a,b
!integer :: sizeD
!
!sizeD=size(Dirac_inv,1)
!
!!call make_Dirac(Dirac,UMAT,Phi)
!!call matrix_product(prod,Dirac,Dirac_inv,"N","N")
!!do s=1,sizeD
!  !write(*,*) abs(prod(s,s))
!!enddo
!
!bPhi=conjg(Phi)
!
!obs=(0d0,0d0)
!do l=1,num_links
!  !tmp=(0d0,0d0)
!  call make_diff_phi(dphi_mat,l,UMAT,bPhi)
!  do a=1,dimG
!    call trace_MTa(dphi_a,dphi_mat,a,NMAT)
!    do s=1,num_sites
!      do b=1,dimG
!        obs=obs &
!          + (0d0, -1d0) * cmplx(alpha_l(l) * overall_factor**2 * mass_square_phi) &
!            * dphi_a * Phi(b,s) &
!            * Dirac_inv(link_index(a,l,NMAT,num_sites),site_index(b,s,NMAT)) 
!      enddo
!    enddo
!  enddo
!enddo
!
!end subroutine calc_pcsc_link
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! face part of the PCSC
!subroutine calc_pcsc_face(obs,Phi,UMAT,Dirac_inv)
!use SUN_generators, only : trace_MTa
!use matrix_functions, only : matrix_power
!implicit none
!
!complex(kind(0d0)), intent(out) :: obs
!complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
!complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)), intent(in) :: Dirac_inv(:,:)
!complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT), Ufm(1:NMAT,1:NMAT),Omega(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: Omega_a
!
!integer :: s,f,a,b
!integer :: sizeD
!
!sizeD=size(Dirac_inv,1)
!
!obs=(0d0,0d0)
!do f=1,num_faces
!  call Make_face_variable(Uf,f,UMAT)
!  call matrix_power(Ufm,Uf,m_omega)
!  call make_moment_map(Omega,Ufm)
!  do a=1,dimG
!    call trace_MTa(Omega_a,Omega,a,NMAT)
!    do s=1,num_sites
!      do b=1,dimG
!        obs=obs &
!          + (0d0, -0.5d0) * cmplx(alpha_f(f) * beta_f(f) &
!                                * overall_factor**2 * mass_square_phi) &
!            * Omega_a * Phi(b,s) &
!            * Dirac_inv(face_index(a,f,NMAT,num_sites,num_links),site_index(b,s,NMAT)) 
!      enddo
!    enddo
!  enddo
!enddo
!
!end subroutine calc_pcsc_face
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! PCSC mass part
!!!  pcsc_mass = "mu^2/2g^2 \Sigma*Tr(Phi \eta)"
!subroutine calc_pcsc_mass(pcsc_mass,Phi,UMAT,Dirac_inv)
!implicit none
!
!complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
!complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)), intent(in) :: Dirac_inv(:,:)
!complex(kind(0d0)), intent(out) :: pcsc_mass
!
!complex(kind(0d0)) :: obs_s
!complex(kind(0d0)) :: obs_l
!complex(kind(0d0)) :: obs_f
!real(8) :: Sb
!
!call calc_pcsc_site(obs_s,Phi,Dirac_inv)
!call calc_pcsc_link(obs_l,Phi,UMAT,Dirac_inv)
!call calc_pcsc_face(obs_f,Phi,UMAT,Dirac_inv)
!
!pcsc_mass= obs_s + obs_l + obs_f
!end subroutine calc_pcsc_mass
! 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! PCSC
!!! Let C be a compensator and obs be the one we compute here.
!!! Then
!!! <obs*|C|> = (Nc^2-1)/2*(Nsite+Nlink) * <|C|>
!subroutine calc_pcsc(obs,Phi,UMAT,Dirac_inv)
!implicit none
!
!complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
!complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)), intent(in) :: Dirac_inv(:,:)
!!complex(kind(0d0)), intent(in) :: C
!complex(kind(0d0)), intent(out) :: obs
!
!complex(kind(0d0)) :: obs_s
!complex(kind(0d0)) :: obs_l
!complex(kind(0d0)) :: obs_f
!real(8) :: Sb
!
!call calc_pcsc_site(obs_s,Phi,Dirac_inv)
!call calc_pcsc_link(obs_l,Phi,UMAT,Dirac_inv)
!call calc_pcsc_face(obs_f,Phi,UMAT,Dirac_inv)
!
!call calc_bosonic_action(Sb,UMAT,Phi)
!
!
!
!obs=cmplx(Sb) + obs_s + obs_l + obs_f
!end subroutine calc_pcsc
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Naive WT identity
!! deviation of the WT identity in naive quench 
!subroutine QbarAtr_sigma(obs,Phi,UMAT,Dirac_inv)
!use SUN_generators, only : make_traceless_matrix_from_modes, trace_MTa
!use matrix_functions, only : matrix_power,  matrix_commutator
!implicit none
!
!complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
!complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)), intent(in) :: Dirac_inv(:,:)
!complex(kind(0d0)), intent(out) :: obs
!
!complex(kind(0d0)) :: bPhi(1:dimG,1:num_sites)
!complex(kind(0d0)) :: tmp, trbPHi, bPhi_eta_Sigma
!complex(kind(0d0)) :: phi_mat(1:NMAT,1:NMAT), comm(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: dphi_mat(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: comm_b, dphi_b
!real(8) :: abs_z, arg_z
!
!complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT), Ufm(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: Omega(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: Omega_b
!
!integer :: euler
!integer :: a,b,s,l,f,t
!
!euler = num_sites - num_links + num_faces
!bPhi=conjg( Phi )
!
!obs=(0d0,0d0)
!do s=1,num_sites
!  !! (1/Ns tr( \bar\Phi_s^2 ) )^{ - dimG \chi / 4 - 1 }
!  trbPhi=(0d0,0d0)
!  do a=1,dimG
!    trbPhi = trbPhi + bPhi(a,s) * bPhi(a,s)
!  enddo
!  trbPhi = trbPhi / cmplx(dble(NMAT))
!  abs_z = abs( trbPhi ) 
!  arg_z = arg( trbPhi ) 
!
!  trbPhi = abs_z**( -dble(dimG*euler)/4d0 - 1d0 ) &
!    * exp( (0d0,1d0) * arg_z * ( -dble(dimG*euler)/4d0 - 1d0 ) )
!
!  !! bPhi_eta_Sigma = 1/NMAT* \sum_a( \bar\Phi_s^a \eta_s^a \Sigma )
!  bPhi_eta_Sigma=(0d0,0d0)
!  do a=1,dimG
!    ! \Sigma_s
!    do t=1,num_sites
!      call make_traceless_matrix_from_modes(phi_mat, NMAT, Phi(:,t))
!      call matrix_commutator(comm, phi_mat, phi_mat, 'N', 'C')
!      do b=1,dimG
!        call trace_MTa(comm_b, comm, b, NMAT)
!          bPhi_eta_Sigma = bPhi_eta_Sigma &
!            + cmplx(alpha_s(t) * overall_factor * 0.25d0 ) * comm_b & 
!              * bPhi(a,s) &
!              * Dirac_inv(site_index(a,s,NMAT), site_index(b,t,NMAT)) 
!      enddo
!    enddo
!    ! \Sigma_l
!    do l=1,num_links
!      call make_diff_phi(dphi_mat,l,UMAT,bPhi)
!      do b=1,dimG
!        call trace_MTa(dphi_b,dphi_mat,b,NMAT)
!        bPhi_eta_Sigma = bPhi_eta_Sigma &
!          + cmplx(alpha_l(l) * overall_factor) * (0d0,-1d0) * dphi_b &
!            * bPhi(a,s) &
!            * Dirac_inv(site_index(a,s,NMAT), link_index(b,l,NMAT,num_sites)) 
!      enddo
!    enddo
!    ! ¥Sigma_f
!    do f=1,num_faces
!      call Make_face_variable(Uf,f,UMAT)
!      call matrix_power(Ufm,Uf,m_omega)
!      call make_moment_map(Omega,Ufm)
!      do b=1,dimG
!        call trace_MTa(Omega_b,Omega,b,NMAT)
!        bPhi_eta_Sigma = bPhi_eta_Sigma &
!          + cmplx(alpha_f(f) * beta_f(f) * overall_factor) * (0d0,-0.5d0) * Omega_b &
!            * bPhi(a,s) &
!            * Dirac_inv(site_index(a,s,NMAT), face_index(b,f,NMAT,num_sites,num_links)) 
!      enddo
!    enddo
!  enddo
!  obs = obs + trbPhi * bPhi_eta_Sigma / cmplx(dble( NMAT ))
!enddo
!obs = obs * cmplx( dble( -dimG*euler ) / dble( 2*num_sites ) )
!end subroutine QbarAtr_sigma


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Sf part 
subroutine partial_Sf(obs_site, obs_link, obs_face,Dirac,Dirac_inv)
implicit none

complex(kind(0d0)), intent(out) :: obs_site, obs_link, obs_face
complex(kind(0d0)), intent(in) :: Dirac(:,:), Dirac_inv(:,:)
integer :: s, t, l, l2, f, f2, a, b

!!!
obs_site=(0d0,0d0)
do s=1,num_sites
  do t=1,num_sites
    do a=1,dimG
      do b=1,dimG
        obs_site = obs_site + (0.5d0,0d0) &
          * Dirac( site_index(a,s,NMAT), site_index(b,t,NMAT) ) &
          * Dirac_inv( site_index(b,t,NMAT), site_index(a,s,NMAT) )
      enddo
    enddo
  enddo
enddo
!!!
obs_link=(0d0,0d0)
do l=1,num_links
  do t=1,num_sites
    do a=1,dimG
      do b=1,dimG
        obs_link = obs_link + (0.5d0,0d0) &
          * Dirac( link_index(a,l,NMAT,num_sites), site_index(b,t,NMAT) ) &
          * Dirac_inv( site_index(b,t,NMAT), link_index(a,l,NMAT,num_sites) )
        !!
        obs_link = obs_link + (0.5d0,0d0) &
          * Dirac( site_index(a,t,NMAT), link_index(b,l,NMAT,num_sites) ) &
          * Dirac_inv( link_index(b,l,NMAT,num_sites), site_index(a,t,NMAT) )
      enddo
    enddo
  enddo
enddo
do l=1,num_links
  do l2=1,num_links
    do a=1,dimG
      do b=1,dimG
        obs_link = obs_link + (0.5d0,0d0) &
          * Dirac( link_index(a,l,NMAT,num_sites), link_index(b,l2,NMAT,num_sites) ) &
          * Dirac_inv( link_index(b,l2,NMAT,num_sites), link_index(a,l,NMAT,num_sites) )
      enddo
    enddo
  enddo
enddo
!!!
obs_face=(0d0,0d0)
do f=1,num_faces
  do l=1,num_links
    do a=1,dimG
      do b=1,dimG
        obs_face = obs_face + (0.5d0,0d0) &
          * Dirac( face_index(a,f,NMAT,num_sites,num_links), link_index(b,l,NMAT,num_sites) ) &
          * Dirac_inv( link_index(b,l,NMAT,num_sites), face_index(a,f,NMAT,num_sites,num_links) )
        !!
        obs_face = obs_face + (0.5d0,0d0) &
          * Dirac( link_index(a,l,NMAT,num_sites), face_index(b,f,NMAT,num_sites,num_links) ) &
          * Dirac_inv( face_index(b,f,NMAT,num_sites,num_links), link_index(a,l,NMAT,num_sites) )
      enddo
    enddo
  enddo
enddo
do f=1,num_faces
  do f2=1,num_faces
    do a=1,dimG
      do b=1,dimG
        obs_face = obs_face + (0.5d0,0d0) &
          * Dirac( face_index(a,f,NMAT,num_sites,num_links), face_index(b,f2,NMAT,num_sites,num_links) ) &
          * Dirac_inv( face_index(b,f2,NMAT,num_sites,num_links), face_index(a,f,NMAT,num_sites,num_links) )
      enddo
    enddo
  enddo
enddo
end subroutine partial_Sf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compensator_trace_SU2 (using eigenvalue)
!! VARID ONLY FOR SU(2)
!subroutine calc_compensator_trace_SU2(C,Phi)
!use SUN_generators, only : make_traceless_matrix_from_modes
!use matrix_functions, only : matrix_eigenvalues
!implicit none
!
!complex(kind(0d0)), intent(out) :: C
!complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
!complex(kind(0d0)) :: phi_mat(1:NMAT,1:NMAT), eigen(1:NMAT)
!real(8) :: r, theta
!!! data of simplicial complex
!integer :: s,a,num
!
!num = -3*(num_sites-num_links+num_faces) / 2
!
!
!C=(0d0,0d0)
!do s=1,num_sites
!  call make_traceless_matrix_from_modes(phi_mat, NMAT, Phi(:,s))
!  call matrix_eigenvalues(eigen,phi_mat)
!  r=abs(eigen(1))
!  theta=arg(eigen(1))
!  C=C+cmplx( r**num )* exp( (0d0,1d0) * cmplx( dble(num) * theta ) )
!enddo
!C=C/ cmplx(dble(num_sites))
!
!end subroutine calc_compensator_trace_SU2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compensator_IZ1 
!subroutine calc_compensator_IZ1(compensator,Phi,UMAT,Dirac_inv)
!use SUN_generators, only : make_traceless_matrix_from_modes, TtimesM
!use matrix_functions, only : matrix_product, matrix_inverse, matrix_power
!implicit none
!
!complex(kind(0d0)), intent(out) :: compensator
!complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
!complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)), intent(in) :: Dirac_inv(:,:)
!complex(kind(0d0)) :: phi_tip(1:NMAT,1:NMAT), phi_org(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: UphiUdag(1:NMAT,1:NMAT),tmpmat(1:NMAT,1:NMAT), base_mat(1:NMAT,1:NMAT)
!!! data of simplicial complex
!complex(kind(0d0)) :: f_abc
!integer :: l,n,a,b,c,k, i,j
!
!k = -dimG*(num_sites-num_links+num_faces) / 2
!
!compensator=(0d0,0d0)
!do l=1,num_links
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  !! construct 2 ¥Phi.U.¥Phi.U^¥dag + ¥lambda_l¥lambda_l (U.¥Phi.U^dag + ¥Phi)
!  call make_traceless_matrix_from_modes(phi_org,NMAT,Phi(:,link_org(l)))
!  call make_traceless_matrix_from_modes(phi_tip,NMAT,Phi(:,link_tip(l)))
!
!  call matrix_product(tmpmat,UMAT(:,:,l),phi_tip,'N','N')
!  call matrix_product(UphiUdag,tmpmat,UMAT(:,:,l),'N','C')
!
!  ! 2 ¥Phi U ¥Phi U^¥dagger
!  call matrix_product(base_mat, (2d0,0d0)*phi_org, UphiUdag,'N','N')
!  
!  do n=1,NZF
!    a=NZF_index(1,n)
!    b=NZF_index(2,n)
!    c=NZF_index(3,n)
!    f_abc=cmplx(NZF_value(n))
!
!    call TtimesM(tmpmat, UphiUdag+phi_org, c, NMAT)
!
!    base_mat=base_mat &
!      + (0d0,0.5d0) * f_abc * Dirac_inv( link_index(a,l,NMAT,num_sites), link_index(b,l,NMAT,num_sites) ) &
!        * tmpmat
!  enddo
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ! tmpmat = ¥PHi_org^{k-2}
!  if( k-2 < 0 ) then
!    call matrix_inverse(phi_org)
!  endif
!  call matrix_power(tmpmat,phi_org,abs(k-2))
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  do i=1,NMAT
!    do j=1,NMAT
!      compensator=compensator+tmpmat(i,j)*base_mat(j,i)
!    enddo
!  enddo
!
!enddo
!
!compensator = compensator / cmplx( dble( num_links * NMAT ) )
!end subroutine calc_compensator_IZ1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compensator_IZ2 
!subroutine calc_compensator_IZ2(compensator,Phi,UMAT,Dirac_inv)
!use SUN_generators, only : make_traceless_matrix_from_modes, TtimesM
!use matrix_functions, only : matrix_product, matrix_inverse, matrix_power
!implicit none
!
!complex(kind(0d0)), intent(out) :: compensator
!complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
!complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)), intent(in) :: Dirac_inv(:,:)
!complex(kind(0d0)) :: phi_tip(1:NMAT,1:NMAT), phi_org(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: UphiUdag(1:NMAT,1:NMAT),tmpmat(1:NMAT,1:NMAT), base_mat(1:NMAT,1:NMAT)
!!! data of simplicial complex
!complex(kind(0d0)) :: f_abc, ctmp
!real(8) :: abs_c, arg_c
!integer :: l,n,a,b,c,k, i,j
!
!k = -dimG*(num_sites-num_links+num_faces) / 2
!
!compensator=(0d0,0d0)
!do l=1,num_links
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  !! construct 2 ¥Phi.U.¥Phi.U^¥dag + ¥lambda_l¥lambda_l (U.¥Phi.U^dag + ¥Phi)
!  call make_traceless_matrix_from_modes(phi_org,NMAT,Phi(:,link_org(l)))
!  call make_traceless_matrix_from_modes(phi_tip,NMAT,Phi(:,link_tip(l)))
!
!  call matrix_product(tmpmat,UMAT(:,:,l),phi_tip,'N','N')
!  call matrix_product(UphiUdag,tmpmat,UMAT(:,:,l),'N','C')
!
!  ! 2 ¥Phi U ¥Phi U^¥dagger
!  call matrix_product(base_mat, (2d0,0d0)*phi_org, UphiUdag,'N','N')
!  
!  do n=1,NZF
!    a=NZF_index(1,n)
!    b=NZF_index(2,n)
!    c=NZF_index(3,n)
!    f_abc=cmplx(NZF_value(n))
!
!    call TtimesM(tmpmat, UphiUdag+phi_org, c, NMAT)
!
!    base_mat=base_mat &
!      + (0d0,0.5d0) * f_abc * Dirac_inv( link_index(a,l,NMAT,num_sites), link_index(b,l,NMAT,num_sites) ) &
!        * tmpmat
!  enddo
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  ctmp=(0d0,0d0)
!  do i=1,NMAT
!    ctmp=ctmp+base_mat(i,i)
!  enddo
!  ctmp=ctmp / cmplx(dble(NMAT))
!  abs_c=abs(ctmp)
!  arg_c=arg(ctmp)
!  ctmp=cmplx( abs_c**( dble(k) / 2d0 ) ) &
!         *exp( (0d0,1d0) * cmplx(arg_c* dble(k) / 2d0) )
!
!  compensator = compensator + ctmp / cmplx(dble( num_links ))
!enddo
!
!end subroutine calc_compensator_IZ2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! evaluate PCSC with an compensator 
!subroutine calc_PCSC(PCSC_RE,PCSC_IM,UMAT,PhiMat,C_type)
!#ifdef PARALLEL
!use parallel
!#endif
!implicit none
!
!double precision, intent(out) :: PCSC_RE,PCSC_IM
!complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
!integer, intent(in) :: C_type 
!
!complex(kind(0d0)) :: PCSC
!complex(kind(0d0)) :: A_comp
!double precision :: Sb
!complex(kind(0d0)) :: XiPhiEta
!complex(kind(0d0)) :: tmp
!
!PCSC=(0d0,0d0)
!if( C_type==1 ) call calc_compensator_trace(A_comp,PhiMat) !RANK=0 has the value
!call calc_bosonic_action(Sb,Umat,PhiMat)
!call calc_XiPhiEta(XiPhiEta,Umat,PhiMat)
!
!PCSC= dcmplx(cdabs(A_comp)) &
!  *( dcmplx(Sb) &
!     - dcmplx( dble((NMAT*NMAT-1)*(global_num_sites+global_num_links))/2d0 ) &
!     + dcmplx( mass_square_phi*0.5d0 ) * XiPhiEta ) 
!  
!
!PCSC_RE=dble(PCSC)
!PCSC_IM=dble((0d0,-1d0)*PCSC)
!
!end subroutine calc_PCSC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compute \Xi Tr(Phi eta)
subroutine calc_XiPhiEta(XiPhiEta,Umat,PhiMat)
#ifdef PARALLEL
use parallel
#endif
implicit none

complex(kind(0d0)), intent(out) :: XiPhiEta
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) Xi_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) Xi_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) Xi_chi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)) DinvXi_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) DinvXi_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) DinvXi_chi(1:NMAT,1:NMAT,1:num_faces)

complex(kind(0d0)) tmp
integer :: i,j,s
integer :: info

call make_XiVec(Xi_eta,Xi_lambda,Xi_chi,Umat,PhiMat)
call calc_DinvF(DinvXi_eta,DinvXi_lambda,DinvXi_chi,Xi_eta,Xi_lambda,Xi_chi,Umat,PhiMat,info)

XiPhiEta=(0d0,0d0)
if( info==1 ) then
#ifdef PARALLEL
  if(MYRANK==0) then
#endif
  write(*,*) "# Dirac is singular"
#ifdef PARALLEL
  endif
#endif
  return
endif

tmp=(0d0,0d0)
do s=1,num_sites
  do j=1,NMAT
    do i=1,NMAT
      tmp=tmp + DinvXi_eta(i,j,s)*PhiMat(j,i,s)
    enddo
  enddo
enddo
#ifdef PARALLEL
call MPI_REDUCE(tmp,XiPhiEta,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)
#else
XiPhiEta=tmp
#endif

end subroutine calc_XiPhiEta


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Maximal eigenvalue of (D^dag D)
subroutine max_eigen_DdagD(eigen,Umat,PhiMat)
use Dirac_operator, only :prod_DdagD
implicit none

complex(kind(0d0)), intent(out) :: eigen
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) :: eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: chi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)) :: Xeta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: Xlambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: Xchi(1:NMAT,1:NMAT,1:num_faces)

complex(kind(0d0)) :: gauss(1:(NMAT*NMAT-1)*(num_sites+num_links+num_faces) ) 
double precision :: gauss2(1:2*(NMAT*NMAT-1)*(num_sites+num_links+num_faces) ) 

complex(kind(0d0)) :: val,tmp_val, num,denom,eigen_pre
integer :: i,j,s,l,f,info,ite

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! random eta/lambda/chi
call BoxMuller2(gauss2, (NMAT*NMAT-1)*(num_sites+num_links+num_faces) )
val=(0d0,0d0)
do i=1,(NMAT*NMAT-1)*(num_sites+num_links+num_faces)
  gauss(i)=(gauss2(2*i-1) + gauss2(2*i)*(0d0,1d0)) * dcmplx(dsqrt(0.5d0))
  val=val+gauss(i)*dconjg(gauss(i))
enddo
val=dcmplx(dsqrt( dble(val) ))
do i=1,(NMAT*NMAT-1)*(num_sites+num_links+num_faces)
  gauss(i)=gauss(i)/val
enddo
call vec_to_mat(eta(:,:,1:num_sites),lambda(:,:,1:num_links),chi(:,:,1:num_faces),gauss)
#ifdef PARALLEL
  call syncronize_sites(eta)
  call syncronize_links(lambda)
  call syncronize_faces(chi)
#endif

eigen=(0d0,0d0)
do ite=1,10000
  !!
  eigen_pre=eigen
  !!
  call prod_DdagD(Xeta,Xlambda,Xchi,eta,lambda,chi,Umat,PhiMat)
  call InnerProd(num, &
    Xeta, Xlambda, Xchi, &
    eta, lambda, chi)
  call InnerProd(denom, &
    eta, lambda, chi, &
    eta, lambda, chi)
  eigen=num/denom
  !!
  if( cdabs(eigen-eigen_pre) < epsilon) exit
  eta(:,:,1:num_sites)=Xeta
  lambda(:,:,1:num_links)=Xlambda
  chi(:,:,1:num_faces)=Xchi
  !! normalization
  tmp_val=(0d0,0d0)
  do s=1,num_sites
    do j=1,NMAT
      do i=1,NMAT
        tmp_val=tmp_val+eta(i,j,s)*dconjg( eta(i,j,s) )
      enddo
    enddo
  enddo
  do l=1,num_links
    do j=1,NMAT
      do i=1,NMAT
        tmp_val=tmp_val+lambda(i,j,l)*dconjg( lambda(i,j,l) )
      enddo
    enddo
  enddo
  do f=1,num_faces
    do j=1,NMAT
      do i=1,NMAT
        tmp_val=tmp_val+chi(i,j,f)*dconjg( chi(i,j,f) )
      enddo
    enddo
  enddo
  val=(0d0,0d0)
#ifdef PARALLEL
  call MPI_ALLREDUCE(tmp_val,val,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
#else
  val=tmp_val
#endif
  val=dcmplx(dsqrt(dble(val)))
  do s=1,num_sites
    do j=1,NMAT
      do i=1,NMAT
        eta(i,j,s)=eta(i,j,s)/val
      enddo
    enddo
  enddo
  do l=1,num_links
    do j=1,NMAT
      do i=1,NMAT
        lambda(i,j,l)=lambda(i,j,l)/val
      enddo
    enddo
  enddo
  do f=1,num_faces
    do j=1,NMAT
      do i=1,NMAT
        chi(i,j,f)=chi(i,j,f)/val
      enddo
    enddo
  enddo
#ifdef PARALLEL
  call syncronize_sites(eta)
  call syncronize_links(lambda)
  call syncronize_faces(chi)
#endif
enddo

end subroutine max_eigen_DdagD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Minimal eigenvalue of (D^dag D)
subroutine min_eigen_DdagD(eigen,Umat,PhiMat)
use Dirac_operator, only :prod_DdagD
implicit none

complex(kind(0d0)), intent(out) :: eigen
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) :: eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: chi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)) :: Xeta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: Xlambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: Xchi(1:NMAT,1:NMAT,1:num_faces)

complex(kind(0d0)) :: gauss(1:(NMAT*NMAT-1)*(num_sites+num_links+num_faces) ) 
double precision :: gauss2(1:2*(NMAT*NMAT-1)*(num_sites+num_links+num_faces) ) 

complex(kind(0d0)) :: val,tmp_val, num,denom,eigen_pre
integer :: i,j,s,l,f,info,ite

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! random eta/lambda/chi
call BoxMuller2(gauss2, (NMAT*NMAT-1)*(num_sites+num_links+num_faces) )
val=(0d0,0d0)
do i=1,(NMAT*NMAT-1)*(num_sites+num_links+num_faces)
  gauss(i)=(gauss2(2*i-1) + gauss2(2*i)*(0d0,1d0)) * dcmplx(dsqrt(0.5d0))
  val=val+gauss(i)*dconjg(gauss(i))
enddo
val=dcmplx(dsqrt( dble(val) ))
do i=1,(NMAT*NMAT-1)*(num_sites+num_links+num_faces)
  gauss(i)=gauss(i)/val
enddo
call vec_to_mat(eta(:,:,1:num_sites),lambda(:,:,1:num_links),chi(:,:,1:num_faces),gauss)
#ifdef PARALLEL
  call syncronize_sites(eta)
  call syncronize_links(lambda)
  call syncronize_faces(chi)
#endif

eigen=(0d0,0d0)
do ite=1,10000
  !!
  eigen_pre=eigen
  !!
  call calc_DdagDinvF(Xeta,Xlambda,Xchi,eta,lambda,chi,Umat,Phimat,info)
  call InnerProd(num, &
    Xeta, Xlambda, Xchi, &
    eta, lambda, chi)
  call InnerProd(denom, &
    eta, lambda, chi, &
    eta, lambda, chi)
  eigen=num/denom
  !!
  if( cdabs(eigen-eigen_pre) < epsilon) exit
  eta(:,:,1:num_sites)=Xeta
  lambda(:,:,1:num_links)=Xlambda
  chi(:,:,1:num_faces)=Xchi
  !! normalization
  tmp_val=(0d0,0d0)
  do s=1,num_sites
    do j=1,NMAT
      do i=1,NMAT
        tmp_val=tmp_val+eta(i,j,s)*dconjg( eta(i,j,s) )
      enddo
    enddo
  enddo
  do l=1,num_links
    do j=1,NMAT
      do i=1,NMAT
        tmp_val=tmp_val+lambda(i,j,l)*dconjg( lambda(i,j,l) )
      enddo
    enddo
  enddo
  do f=1,num_faces
    do j=1,NMAT
      do i=1,NMAT
        tmp_val=tmp_val+chi(i,j,f)*dconjg( chi(i,j,f) )
      enddo
    enddo
  enddo
  val=(0d0,0d0)
#ifdef PARALLEL
  call MPI_ALLREDUCE(tmp_val,val,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
#else
  val=tmp_val
#endif
  val=dcmplx(dsqrt(dble(val)))
  do s=1,num_sites
    do j=1,NMAT
      do i=1,NMAT
        eta(i,j,s)=eta(i,j,s)/val
      enddo
    enddo
  enddo
  do l=1,num_links
    do j=1,NMAT
      do i=1,NMAT
        lambda(i,j,l)=lambda(i,j,l)/val
      enddo
    enddo
  enddo
  do f=1,num_faces
    do j=1,NMAT
      do i=1,NMAT
        chi(i,j,f)=chi(i,j,f)/val
      enddo
    enddo
  enddo
#ifdef PARALLEL
  call syncronize_sites(eta)
  call syncronize_links(lambda)
  call syncronize_faces(chi)
#endif
enddo
eigen=(1d0,0d0)/eigen


end subroutine min_eigen_DdagD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! (D D†~)^{-1} F
subroutine calc_DdagDinvF(X_eta,X_lambda,X_chi,eta,lambda,chi,Umat,Phimat,info)
use Dirac_operator, only :prod_DdagD
#ifdef PARALLEL
use parallel
use global_subroutines, only : syncronize_sites, syncronize_links, syncronize_faces
#endif
implicit none 

complex(kind(0d0)), intent(out) :: X_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(out) :: X_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(out) :: X_chi(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)), intent(in) :: eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: chi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
integer, intent(out) :: info

complex(kind(0d0)) :: r_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: r_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: r_chi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)) :: p_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: p_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: p_chi(1:NMAT,1:NMAT,1:num_necessary_faces)
!complex(kind(0d0)) :: s_eta(1:NMAT,1:NMAT,1:num_sites)
!complex(kind(0d0)) :: s_lambda(1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)) :: s_chi(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: tmp_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: tmp_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: tmp_chi(1:NMAT,1:NMAT,1:num_faces)

complex(kind(0d0)) :: rkrk,alpha_k,tmp_c1,beta_k
integer :: MaxIte=10000

integer :: ite,s,l,f


!! 初期配置
!Xvec=(0d0,0d0)
X_eta=(0d0,0d0)
X_lambda=(0d0,0d0)
X_chi=(0d0,0d0)
!r_vec = Bvec
r_eta=eta
r_lambda=lambda
r_chi=chi

p_eta=r_eta
p_lambda=r_lambda
p_chi=r_chi

info=0
!! iteration start
do ite=1,MaxIte
    ! (1) construct \alpha_k
  ! rkrk = (r_k, r_k)
  call InnerProd(rkrk, &
    r_eta, r_lambda, r_chi, &
    r_eta, r_lambda, r_chi)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!== DEPEND ON ProdMat YOUR MADE ==
  call Prod_DdagD(&
    tmp_eta,tmp_lambda,tmp_chi, &
    p_eta,p_lambda,p_chi, &
    Umat,PhiMat)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call InnerProd(tmp_c1, &
    p_eta,p_lambda,p_chi, &
    tmp_eta, tmp_lambda, tmp_chi)
  alpha_k = rkrk / dcmplx(dble(tmp_c1))
 
 ! (2) update Xvec and r_vec
  X_eta = X_eta + dcmplx(alpha_k)*p_eta(:,:,1:num_sites)
  X_lambda = X_lambda + dcmplx(alpha_k)*p_lambda(:,:,1:num_links)
  X_chi = X_chi + dcmplx(alpha_k)*p_chi(:,:,1:num_faces)

  r_eta(:,:,1:num_sites) = r_eta(:,:,1:num_sites) - dcmplx(alpha_k)*tmp_eta
  r_lambda(:,:,1:num_links) = r_lambda(:,:,1:num_links) - dcmplx(alpha_k)*tmp_lambda
  r_chi(:,:,1:num_faces) = r_chi(:,:,1:num_faces) - dcmplx(alpha_k)*tmp_chi
#ifdef PARALLEL
  call syncronize_sites(r_eta)
  call syncronize_faces(r_chi)
  call syncronize_links(r_lambda)
#endif
    
 ! (3) conservation check
  call InnerProd( tmp_c1, &
    r_eta(:,:,1:num_sites), r_lambda(:,:,1:num_links), r_chi(:,:,1:num_faces), &
    r_eta(:,:,1:num_sites), r_lambda(:,:,1:num_links), r_chi(:,:,1:num_faces))
  if ( dsqrt( dble(tmp_c1) ) < epsilon ) then
    return
  endif

! (4) update p_k --> p_{k+1}
!construct beta_k
  !call InnerProd(tmp_c1, r_vec, r_vec)
  beta_k = tmp_c1 / rkrk
  p_eta(:,:,1:num_sites) = r_eta(:,:,1:num_sites) + dcmplx(beta_k) * p_eta(:,:,1:num_sites)
  p_lambda(:,:,1:num_links) = r_lambda(:,:,1:num_links) + dcmplx(beta_k) * p_lambda(:,:,1:num_links)
  p_chi(:,:,1:num_faces) = r_chi(:,:,1:num_faces) + dcmplx(beta_k) * p_chi(:,:,1:num_faces)
#ifdef PARALLEL
  call syncronize_sites(p_eta)
  call syncronize_faces(p_chi)
  call syncronize_links(p_lambda)
#endif
enddo

info = 1
return

end subroutine calc_DdagDinvF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Dirac inverse by CG
subroutine calc_DinvF(X_eta,X_lambda,X_chi,eta,lambda,chi,Umat,Phimat,info)
use Dirac_operator, only : prod_Dirac,prod_DiracDag
#ifdef PARALLEL
use parallel
use global_subroutines, only : syncronize_sites, syncronize_links, syncronize_faces
#endif
implicit none 

complex(kind(0d0)), intent(out) :: X_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(out) :: X_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(out) :: X_chi(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)), intent(in) :: eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: chi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
integer, intent(out) :: info

complex(kind(0d0)) :: r_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: r_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: r_chi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)) :: p_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: p_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: p_chi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)) :: s_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: s_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: s_chi(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: Dp_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: Dp_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: Dp_chi(1:NMAT,1:NMAT,1:num_faces)

complex(kind(0d0)) :: rkrk,sksk,sksk2,DpDp,lambda_k,mu_k
integer :: MaxIte=10000

integer :: ite,s,l,f

!! 初期配置
!Xvec=(0d0,0d0)
X_eta=(0d0,0d0)
X_lambda=(0d0,0d0)
X_chi=(0d0,0d0)
!r_vec = Bvec
r_eta=eta
r_lambda=lambda
r_chi=chi
! p_vec = s_vec = D^\dag r_vec
call prod_DiracDag(p_eta(:,:,1:num_sites),p_lambda(:,:,1:num_links),p_chi(:,:,1:num_faces),r_eta,r_lambda,r_chi,Umat,Phimat)
#ifdef PARALLEL
  call syncronize_sites(p_eta)
  call syncronize_faces(p_chi)
  call syncronize_links(p_lambda)
#endif
s_eta=p_eta(:,:,1:num_sites)
s_lambda=p_lambda(:,:,1:num_links)
s_chi=p_chi(:,:,1:num_faces)

!! preparation
call InnerProd(sksk, &
 s_eta, s_lambda, s_chi, &
 s_eta, s_lambda, s_chi)

info=0
!! iteration start
do ite=1,MaxIte
  ! (1) construct lambda_k 
  ! rkrk = (r_k, r_k)
  !call InnerProd(rkrk, r_vec, r_vec)
  call Prod_Dirac(Dp_eta,Dp_lambda,Dp_chi, &
    p_eta,p_lambda,p_chi,Umat,PhiMat)
  call InnerProd(DpDp, &
    Dp_eta, Dp_lambda, Dp_chi, &
    Dp_eta, Dp_lambda, Dp_chi)
  lambda_k = dexp( dlog(dble(sksk)) - dlog(dble(DpDp)) ) 
   
  X_eta=X_eta+dcmplx(lambda_k)*p_eta(:,:,1:num_sites)
  X_lambda=X_lambda+dcmplx(lambda_k)*p_lambda(:,:,1:num_links)
  X_chi=X_chi+dcmplx(lambda_k)*p_chi(:,:,1:num_faces)

 
 ! (2) update r_vec
  !r_vec = r_vec - dcmplx(alpha_k)*tmp_vec
  r_eta(:,:,1:num_sites)=r_eta(:,:,1:num_sites)-dcmplx(lambda_k)*Dp_eta
  r_lambda(:,:,1:num_links)=r_lambda(:,:,1:num_links)-dcmplx(lambda_k)*Dp_lambda
  r_chi(:,:,1:num_faces)=r_chi(:,:,1:num_faces)-dcmplx(lambda_k)*Dp_chi
#ifdef PARALLEL
  call syncronize_sites(r_eta)
  call syncronize_faces(r_chi)
  call syncronize_links(r_lambda)
#endif
    
 ! (3) conservation check
  !call InnerProd( tmp_c1, r_vec, r_vec)
  call InnerProd(rkrk, &
    r_eta(:,:,1:num_sites),r_lambda(:,:,1:num_links),r_chi(:,:,1:num_faces),&
    r_eta(:,:,1:num_sites),r_lambda(:,:,1:num_links),r_chi(:,:,1:num_faces))

  if ( dsqrt( dble(rkrk) ) < epsilon ) then
    return
  endif

! (4) update p_k --> p_{k+1}
!construct beta_k
  !call InnerProd(tmp_c1, r_vec, r_vec)
  !beta_k = dble(tmp_c1) / dble(rkrk)
  !p_vec = r_vec + dcmplx(beta_k) * p_vec

  call prod_DiracDag(s_eta,s_lambda,s_chi,r_eta,r_lambda,r_chi,Umat,Phimat)

  call InnerProd(sksk2, &
   s_eta, s_lambda, s_chi, &
   s_eta, s_lambda, s_chi)
  mu_k = dcmplx( dexp( dlog(dble(sksk2)) - dlog(dble(sksk)) ) ) 

  !! skskを更新
  sksk=sksk2
  !! p_vecを更新
  p_eta(:,:,1:num_sites) = s_eta(:,:,1:num_sites) + dcmplx(mu_k) * p_eta(:,:,1:num_sites)
  p_lambda(:,:,1:num_links) = s_lambda(:,:,1:num_links) + dcmplx(mu_k) * p_lambda(:,:,1:num_links)
  p_chi(:,:,1:num_faces) = s_chi(:,:,1:num_faces) + dcmplx(mu_k) * p_chi(:,:,1:num_faces)
#ifdef PARALLEL
  call syncronize_sites(p_eta)
  call syncronize_faces(p_chi)
  call syncronize_links(p_lambda)
#endif
enddo

info = 1
return
end subroutine calc_DinvF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to take inner product for matrix value
!!  (V1,V2) = V1^\dagger_i V2_i
subroutine InnerProd(v1v2, &
   vec1_eta, vec1_lambda, vec1_chi, &
   vec2_eta, vec2_lambda, vec2_chi)
#ifdef PARALLEL
use parallel
#endif
implicit none

complex(kind(0d0)), intent(in) :: vec1_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(in) :: vec1_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: vec1_chi(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)), intent(in) :: vec2_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(in) :: vec2_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: vec2_chi(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)), intent(out) :: v1v2

integer :: s,l,f,i,j
complex(kind(0d0)) :: v1v2_tmp

v1v2=(0d0,0d0)
v1v2_tmp=(0d0,0d0)
do s=1,num_sites
  do j=1,NMAT
    do i=1,NMAT
      v1v2_tmp=v1v2_tmp+dconjg( vec1_eta(i,j,s) ) *vec2_eta(i,j,s) 
    enddo
  enddo
enddo
do l=1,num_links
  do j=1,NMAT
    do i=1,NMAT
      v1v2_tmp=v1v2_tmp+dconjg( vec1_lambda(i,j,l) ) *vec2_lambda(i,j,l) 
    enddo
  enddo
enddo
do f=1,num_faces
  do j=1,NMAT
    do i=1,NMAT
      v1v2_tmp=v1v2_tmp+dconjg( vec1_chi(i,j,f) ) *vec2_chi(i,j,f)
    enddo
  enddo
enddo
#ifdef PARALLEL
call MPI_REDUCE(v1v2_tmp,v1v2,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(v1v2,1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
#else
v1v2=v1v2_tmp
#endif 

end subroutine InnerProd


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to make S_eta, S_lambda, S_chi of
!!  \Zia = Tr(eta S_eta) + Tr(lambda S_lambda) + Tr(chi S_chi)
subroutine make_XiVec(Xi_eta,Xi_lambda,Xi_chi,Umat,PhiMat)
use matrix_functions, only : matrix_commutator,matrix_3_product
#ifdef PARALLEL
use parallel
#endif
implicit none

complex(kind(0d0)), intent(out) :: Xi_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(out) :: Xi_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(out) :: Xi_chi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT)
integer :: s,l,f
integer :: i,j

do s=1,num_sites
  call matrix_commutator(tmpmat,PhiMat(:,:,s),Phimat(:,:,s),'N','C')
  Xi_eta(:,:,s)=dcmplx(alpha_s(s)*0.25d0*overall_factor)*tmpmat
enddo

do l=1,num_links
  do j=1,NMAT
    do i=1,NMAT
      tmpmat(i,j)=-dconjg( PhiMat(j,i,link_org(l)) )
    enddo
  enddo
  call matrix_3_product(tmpmat,&
    Umat(:,:,l),PhiMat(:,:,link_tip(l)),Umat(:,:,l),&
    'N','C','C',(1d0,0d0),'ADD') 
  Xi_lambda(:,:,l)=(0d0,-1d0)*dcmplx(alpha_l(l)*overall_factor)*tmpmat
enddo

do f=1,num_faces
  call make_face_variable(Uf,f,Umat)
  if(m_omega==0) call make_moment_map0(tmpmat,Uf)
  if(m_omega==-1) call make_moment_map_adm(tmpmat,Uf)
  Xi_chi(:,:,f)=(0d0,-0.5d0)*dcmplx(alpha_f(f)*beta_f(f)*overall_factor)*tmpmat
enddo

#ifdef PARALLEL
  call syncronize_sites(Xi_eta)
  call syncronize_links(Xi_lambda)
  call syncronize_faces(Xi_chi)
#endif
end subroutine make_XiVec


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to test Q-transformation of S
subroutine check_QS(Umat,PhiMat)
use matrix_functions, only : matrix_commutator, matrix_3_product
use Dirac_operator, only : Prod_Dirac
#ifdef PARALLEL
use parallel
#endif
implicit none

complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) :: Qeta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: Qlambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: Qchi(1:NMAT,1:NMAT,1:num_necessary_faces)

complex(kind(0d0)) :: tmp_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: tmp_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: tmp_chi(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: Omega(1:NMAT,1:NMAT,1:num_faces)

complex(kind(0d0)) :: Bforce_s(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: Bforce_l(1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)) :: Fforce_s(1:NMAT,1:NMAT,1:num_sites)
!complex(kind(0d0)) :: Fforce_l(1:NMAT,1:NMAT,1:num_links)

complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
double precision :: tmp,QS
integer :: info,s,l,f,i,j,triger

call make_Qfermion(Qeta,Qlambda,Qchi,Omega,Umat,PhiMat)
call make_bosonic_force_nomass(Bforce_s,Bforce_l,Umat,PhiMat)

call Prod_Dirac(tmp_eta,tmp_lambda,tmp_chi,Qeta,Qlambda,Qchi,UMAT,Phimat)
! Q^2 \Omega を care する
!do f=1,num_faces
!  call matrix_commutator(tmpmat,PhiMat(:,:,sites_in_f(f)%label_(1)),Omega(:,:,f))
!  tmp_chi(:,:,f)=tmp_chi(:,:,f)+(0d0,1d0)*dcmplx(beta_f(f))*tmpmat
!enddo

if(MYRANK==0) write(*,*) "# QS = 0 ?"
!write(*,*) tmp_chi
QS=0d0
tmp=0d0
do s=1,num_sites
  !tmp=0d0
  do i=1,NMAT
    do j=1,NMAT
      tmp_eta(i,j,s)=-tmp_eta(i,j,s)+dconjg(Bforce_s(j,i,s))
      tmp=tmp+dble( tmp_eta(i,j,s)*dconjg(tmp_eta(i,j,s)) )
    enddo
  enddo
  !write(*,*) "# site",global_site_of_local(s),tmp
enddo
call MPI_REDUCE(tmp,QS,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
if(MYRANK==0) write(*,*) "#   site:",QS

QS=0d0
tmp=0d0
do l=1,num_links
  !tmp=0d0
  do i=1,NMAT
    do j=1,NMAT
      tmp_lambda(i,j,l) = -tmp_lambda(i,j,l)+Bforce_l(i,j,l)
      tmp=tmp+dble( tmp_lambda(i,j,l)*dconjg(tmp_lambda(i,j,l)) )
    enddo
  enddo
  !write(*,*) "# link",global_link_of_local(l),tmp
enddo
call MPI_REDUCE(tmp,QS,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
if(MYRANK==0) write(*,*) "#   link:",QS

QS=0d0
tmp=0d0
do f=1,num_faces
  !tmp=0d0
  do i=1,NMAT
    do j=1,NMAT
      tmp_chi(i,j,f)=-tmp_chi(j,i,f)
      tmp=tmp+dble( tmp_chi(i,j,f)*dconjg(tmp_chi(i,j,f)) )
    enddo
  enddo
  !write(*,*) "# face",global_face_of_local(f),tmp
enddo
call MPI_REDUCE(tmp,QS,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
if(MYRANK==0) write(*,*) "#   face:",QS

end subroutine check_QS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to make Q\Psi
subroutine make_Qfermion(Qeta,Qlambda,Qchi,Omega,Umat,PhiMat)
use matrix_functions, only : matrix_commutator, matrix_3_product
#ifdef PARALLEL
use parallel
#endif
implicit none

complex(kind(0d0)), intent(out) :: Qeta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(out) :: Qlambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(out) :: Qchi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), intent(out) :: Omega(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
!complex(kind(0d0)), intent(in) :: PF_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)) :: comm(1:NMAT,1:NMAT)

integer :: s,l,f,i,j
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT)

do s=1,num_sites
  call matrix_commutator(Qeta(:,:,s),PhiMat(:,:,s),PhiMat(:,:,s),'N','C')
enddo

do l=1,num_links
  Qlambda(:,:,l)=(0d0,-1d0)*PhiMat(:,:,link_org(l))
  call matrix_3_product(Qlambda(:,:,l),&
    Umat(:,:,l),PhiMat(:,:,link_tip(l)),Umat(:,:,l),&
    'N','N','C',(0d0,1d0),'ADD')
enddo

do f=1,num_faces
  call Make_face_variable(Uf,f,UMAT)
  if(m_omega == 0) then 
    call Make_moment_map0(Omega(:,:,f),Uf)
  elseif(m_omega == -1) then
    call Make_moment_map_adm(Omega(:,:,f),Uf)
  endif
  Qchi(:,:,f)=(0d0,-0.5d0)*dcmplx(beta_f(f))*Omega(:,:,f)
enddo


#ifdef PARALLEL
call syncronize_sites(Qeta)
call syncronize_links(Qlambda)
call syncronize_faces(Qchi)
#endif

end subroutine make_Qfermion



!end module observables
