!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ver.01: correct the normalization of the PCSC relation 
module observables
use global_parameters
use global_subroutines
implicit none

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_bosonic_action(Sb,UMAT,Phi)
use hamiltonian
implicit none

double precision, intent(out) :: Sb
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)

double precision :: SB_S,SB_L,SB_F

SB_S=0d0
SB_L=0d0
SB_F=0d0
call bosonic_action_site(SB_S,Phi)
call bosonic_action_link(SB_L,UMAT,Phi)
call bosonic_action_face(SB_F,UMAT)

Sb=SB_S+SB_L+SB_F

end subroutine calc_bosonic_action


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_TrX2(TrX2, Phi)
implicit none

double precision, intent(out) :: TrX2
complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)

integer :: a,s

TrX2=0d0
do s=1,num_sites
  do a=1,dimG
    TrX2=TrX2+ dble( Phi(a,s)*conjg( Phi(a,s) ) )
  enddo
enddo

TrX2 = TrX2 / (2d0 * dble(num_sites) )

end subroutine calc_TrX2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_eigenvalues_Dirac(eigenvalues,UMAT,Phi)
use Dirac_operator, only : make_Dirac
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)), intent(inout) :: eigenvalues(1:sizeD)

complex(kind(0d0)) :: Dirac(1:sizeD,1:sizeD)
complex(kind(0d0)) VL(1,1:sizeD), VR(1,1:sizeD)
double precision RWORK(1:2*sizeD)
complex(kind(0d0)) WORK(2*(sizeD)),tako1
character JOBVL,JOBVR
integer info,lwork
integer i,j

lwork=2*sizeD
JOBVL='N'
JOBVR='N'

call make_Dirac(Dirac,UMAT,Phi)
!do i=1,sizeD
!  do j=1,sizeD
!    write(*,*) i,j,Dirac(i,j)
!  enddo
!enddo


call ZGEEV(JOBVL,JOBVR,sizeD,&
     DIRAC,sizeD,eigenvalues,VL,1,VR,1,WORK,lwork,RWORK,info)

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
subroutine calc_min_and_max_of_eigenvalues_Dirac(minimal,maximal,UMAT,Phi)
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)), intent(inout) :: minimal, maximal

complex(kind(0d0)) :: eigenvalues(1:sizeD)

call calc_eigenvalues_Dirac(eigenvalues,UMAT,Phi)

minimal=eigenvalues(1)
maximal=eigenvalues(sizeD)

end subroutine calc_min_and_max_of_eigenvalues_Dirac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compensator = 1/(NMAT*Nsite) * \sum_s Tr( \Phi^{ -(N^2-1)*\chi_h / 2 } )
subroutine calc_compensator_naive(C,Phi)
use SUN_generators, only : make_traceless_matrix_from_modes
use matrix_functions, only : matrix_power, matrix_inverse
implicit none

complex(kind(0d0)), intent(out) :: C
complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)) :: tmpPhi(1:NMAT,1:NMAT)
complex(kind(0d0)) :: pmat(1:NMAT,1:NMAT)
!! data of simplicial complex
integer :: Euler,s,i

Euler=num_sites - num_links + num_faces
C=(0d0,0d0)
do s=1,num_sites
  call make_traceless_matrix_from_modes(tmpPhi,NMAT,Phi(:,s))
  if (Euler > 0) then 
    call matrix_inverse(tmpPhi)
  endif
  call matrix_power(pmat,tmpPhi,abs((NMAT*NMAT-1)*Euler)/2)
  do i=1,NMAT
    C=C+pmat(i,i)
  enddo
enddo
C=C/ dble(num_sites*NMAT)

end subroutine calc_compensator_naive

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compensator_trace = 1/(NMAT*Nsite) * \sum_s Tr( \Phi^2 )^{ -(N^2-1)*\chi_h/4 } ) )
subroutine calc_compensator_trace(C,Phi)
implicit none

complex(kind(0d0)), intent(out) :: C
complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)) :: trace
real(8) :: trace_abs, trace_arg, num
!! data of simplicial complex
integer :: s,i,a

num=dble( (NMAT*NMAT-1)*(num_sites-num_links+num_faces) )/4d0

C=(0d0,0d0)
do s=1,num_sites
  trace=(0d0,0d0)
  do a=1,dimG
    trace=trace+Phi(a,s)*Phi(a,s)
  enddo
  trace=trace/cmplx(dble(NMAT))
  trace_abs = sqrt( trace * conjg(trace) ) 
  trace_arg = arg(trace)
  C=C+cmplx(trace_abs**(-num)) &
      * exp( (0d0,-1d0)*trace_arg*num )
      !* ( cmplx( cos( trace_arg * num ) ) - (0d0,1d0)*cmplx( sin( trace_arg * num ) ) )
enddo
C=C/ cmplx(dble(num_sites))

end subroutine calc_compensator_trace

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compensator_trace 
!! = 1/(NMAT*Nsite) 
!!   * [ \sum_s Tr( \Phi^2 ) ]^{ -(N^2-1)*\chi_h/4 } ) 
subroutine calc_compensator2_trace(C,Phi)
use SUN_generators, only : make_traceless_matrix_from_modes
use matrix_functions, only : matrix_power, matrix_inverse
implicit none

complex(kind(0d0)), intent(out) :: C
complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)) :: trace
real(8) :: C_abs, C_arg, num
!! data of simplicial complex
integer :: s,i,a

num=dble( (NMAT*NMAT-1)*(num_sites-num_links+num_faces) )/4d0
C=(0d0,0d0)
do s=1,num_sites
  trace=(0d0,0d0)
  do a=1,dimG
    trace=trace+Phi(a,s)*Phi(a,s)
  enddo
  C=C+trace/cmplx(dble(NMAT))
  !trace_abs = abs(trace)
  !trace_arg = arg(trace)
  !C=C+cmplx(trace_abs**(-num)) &
      !* exp( (0d0,-1d0)*trace_arg*num )
      !* ( cmplx( cos( trace_arg * num ) ) - (0d0,1d0)*cmplx( sin( trace_arg * num ) ) )
enddo
C = C/cmplx(dble(num_sites))
C_abs = abs(C)
C_arg = arg(C)
!write(*,*) C_abs, num, C_abs**(-num)
C = cmplx( C_abs**(-num) ) &
    * exp( (0d0,-1d0) * cmplx( C_arg * num ) ) 

end subroutine calc_compensator2_trace

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compensator25 = 1/(NMAT*Nsite) * \sum_s Tr( \Phi^2 )^{ -(N^2-1)*\chi_h/4 } ) )
!!  where the branch of sqrt is chosen by randum number 
subroutine calc_compensator25(C,Phi)
use SUN_generators, only : make_traceless_matrix_from_modes
use matrix_functions, only : matrix_power, matrix_inverse
implicit none

complex(kind(0d0)), intent(out) :: C
complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)) :: Phimat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: pmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: trace
real(8) :: trace_abs, trace_arg, num
real(8) :: rand, argument
!! data of simplicial complex
integer :: s,i

num=dble( (NMAT*NMAT-1)*(num_sites-num_links+num_faces) )/4d0
C=(0d0,0d0)
do s=1,num_sites
  trace=(0d0,0d0)
  call make_traceless_matrix_from_modes(Phimat,NMAT,Phi(:,s))
  call matrix_power(pmat,Phimat,2)
  do i=1,NMAT
    trace=trace+pmat(i,i)
  enddo
  trace_abs = abs(trace)
  trace_arg = arg(trace)
  call random_number(rand)
  if (rand < 0.5d0) then 
      argument = trace_arg*num
  else
      argument = (2d0*acos(-1d0) + trace_arg)*num
  endif
  C=C+cmplx(trace_abs**(-num)) &
      * ( cmplx( cos( argument ) ) - (0d0,1d0)*cmplx( sin( argument ) ) )
enddo
C=C/ cmplx(dble(num_sites*NMAT))

end subroutine calc_compensator25


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compensator_det 
!!  = 1/Nsite ¥sum_s ( Det(¥Phi_s) )^{ -(N^2-1)*¥chi_h / 2N }
subroutine calc_compensator_det(C,Phi)
use SUN_generators, only : make_traceless_matrix_from_modes
use matrix_functions, only : matrix_determinant
implicit none

complex(kind(0d0)), intent(out) :: C
complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)) :: matPhi(1:NMAT,1:NMAT)
complex(kind(0d0)) :: det
real(8) :: det_abs, det_arg, num
!! data of simplicial complex
integer :: s,i

num = dble( (NMAT*NMAT-1) * (num_sites - num_links + num_faces) ) / dble( 2*NMAT )
C=(0d0,0d0)
do s=1,num_sites
  call make_traceless_matrix_from_modes(matPhi,NMAT,Phi(:,s))
  call matrix_determinant(det,matPhi)
  det_abs=abs(det)
  det_arg=arg(det)
  C = C + cmplx( det_abs**(-num) ) &
      * ( cmplx( cos( det_arg*num ) ) - (0d0,1d0)*cmplx( sin( det_arg*num ) ) )
enddo
C=C / dble(num_sites)
end subroutine calc_compensator_det

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compensator_det 
!!  = [ 1/Nsite ¥sum_s Det(¥Phi_s) ]^{ -(N^2-1)*¥chi_h / 2N }
subroutine calc_compensator2_det(C,Phi)
use SUN_generators, only : make_traceless_matrix_from_modes
use matrix_functions, only : matrix_determinant
implicit none

complex(kind(0d0)), intent(out) :: C
complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)) :: matPhi(1:NMAT,1:NMAT)
complex(kind(0d0)) :: det
real(8) :: C_abs, C_arg, num
!! data of simplicial complex
integer :: s,i

num = dble( (NMAT*NMAT-1) * (num_sites - num_links + num_faces) ) / dble( 2*NMAT )
C=(0d0,0d0)
do s=1,num_sites
  call make_traceless_matrix_from_modes(matPhi,NMAT,Phi(:,s))
  call matrix_determinant(det,matPhi)
  C=C+det
enddo
C = C / cmplx( dble( num_sites ) )
C_abs=abs(C)
C_arg=arg(C)
C = cmplx( C_abs ** (-num) ) &
    * exp( (0d0,-1d0) * cmplx( C_arg * num ) )
end subroutine calc_compensator2_det


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! site part of the PCSC
subroutine calc_pcsc_site(obs,Phi,Dirac_inv)
use SUN_generators, only : make_traceless_matrix_from_modes, trace_MTa
use matrix_functions, only : matrix_commutator
implicit none

complex(kind(0d0)), intent(out) :: obs
complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)), intent(in) :: Dirac_inv(1:sizeD,1:sizeD)
!complex(kind(0d0)), intent(in) :: C
complex(kind(0d0)) :: phi_mat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: comm(1:NMAT,1:NMAT)
complex(kind(0d0)) :: comm_a
complex(kind(0d0)) :: tmp

integer :: s,t,a,b

obs=(0d0,0d0)

do s=1,num_sites
  tmp=(0d0,0d0)
  call make_traceless_matrix_from_modes(phi_mat,NMAT,Phi(:,s))
  call matrix_commutator(comm,phi_mat,phi_mat,'N','C')
  do a=1,dimG
    call trace_MTa(comm_a,comm,a,NMAT)
    do t=1,num_sites
      do b=1,dimG
        tmp=tmp &
          +Dirac_inv(site_index(a,s),site_index(b,t)) & 
          * comm_a * Phi(b,t) 
      enddo
    enddo
  enddo
  obs=obs+tmp*alpha_s(s)
enddo

!obs = obs * overall_factor*overall_factor*0.5d0*mass_square_phi!*C
obs = obs * overall_factor*overall_factor*mass_square_phi!*C

end subroutine calc_pcsc_site
        

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! link part of the PCSC
subroutine calc_pcsc_link(obs,Phi,UMAT,Dirac_inv)
use SUN_generators, only : make_traceless_matrix_from_modes, trace_MTa
implicit none

complex(kind(0d0)), intent(out) :: obs
complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: Dirac_inv(1:sizeD,1:sizeD)
!complex(kind(0d0)), intent(in) :: C
complex(kind(0d0)) :: bPhi(1:dimG,1:num_sites)
complex(kind(0d0)) :: dphi_mat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: dphi_a
complex(kind(0d0)) :: tmp

integer :: s,l,a,b

bPhi=conjg(Phi)

obs=(0d0,0d0)
do l=1,num_links
  tmp=(0d0,0d0)
  call make_diff_phi(dphi_mat,l,UMAT,bPhi)
  do a=1,dimG
    call trace_MTa(dphi_a,dphi_mat,a,NMAT)
    do s=1,num_sites
      do b=1,dimG
        tmp=tmp &
          +Dirac_inv(link_index(a,l),site_index(b,s)) &
          * dphi_a * Phi(b,s)
      enddo
    enddo
  enddo
  obs=obs+tmp*alpha_l(l)
enddo

!obs = obs * (0d0,-1d0) * overall_factor*overall_factor*0.5d0*mass_square_phi !* C
obs = obs * (0d0,-1d0) * overall_factor*overall_factor*mass_square_phi !* C

end subroutine calc_pcsc_link


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! face part of the PCSC
subroutine calc_pcsc_face(obs,Phi,UMAT,Dirac_inv)
use SUN_generators, only : trace_MTa
use matrix_functions, only : matrix_power
implicit none

complex(kind(0d0)), intent(out) :: obs
complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: Dirac_inv(1:sizeD,1:sizeD)
!complex(kind(0d0)), intent(in) :: C
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT), Ufm(1:NMAT,1:NMAT),Omega(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Omega_a
complex(kind(0d0)) :: tmp

integer :: s,f,a,b

obs=(0d0,0d0)
do f=1,num_faces
  tmp=(0d0,0d0)
  call Make_face_variable(Uf,f,UMAT)
  call matrix_power(Ufm,Uf,m_omega)
  call make_moment_map(Omega,Ufm)
  do a=1,dimG
    call trace_MTa(Omega_a,Omega,a,NMAT)
    do s=1,num_sites
      do b=1,dimG
        tmp=tmp &
          +Dirac_inv(face_index(a,f),site_index(b,s)) &
          * Omega_a * Phi(b,s)
      enddo
    enddo
  enddo
  obs=obs+tmp*alpha_f(f)*beta_f(f)
enddo

obs = obs * (0d0,-0.5d0) * overall_factor*overall_factor*mass_square_phi! * C

end subroutine calc_pcsc_face


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! PCSC mass part
!!  pcsc_mass = "mu^2/2g^2 \Sigma*Tr(Phi \eta)"
subroutine calc_pcsc_mass(pcsc_mass,Phi,UMAT,Dirac_inv)
implicit none

complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: Dirac_inv(1:sizeD,1:sizeD)
!complex(kind(0d0)), intent(in) :: C
complex(kind(0d0)), intent(out) :: pcsc_mass

complex(kind(0d0)) :: obs_s
complex(kind(0d0)) :: obs_l
complex(kind(0d0)) :: obs_f
real(8) :: Sb

call calc_pcsc_site(obs_s,Phi,Dirac_inv)
call calc_pcsc_link(obs_l,Phi,UMAT,Dirac_inv)
call calc_pcsc_face(obs_f,Phi,UMAT,Dirac_inv)

pcsc_mass= obs_s + obs_l + obs_f
end subroutine calc_pcsc_mass
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! PCSC
!! Let C be a compensator and obs be the one we compute here.
!! Then
!! <obs*|C|> = (Nc^2-1)/2*(Nsite+Nlink) * <|C|>
subroutine calc_pcsc(obs,Phi,UMAT,Dirac_inv)
implicit none

complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: Dirac_inv(1:sizeD,1:sizeD)
!complex(kind(0d0)), intent(in) :: C
complex(kind(0d0)), intent(out) :: obs

complex(kind(0d0)) :: obs_s
complex(kind(0d0)) :: obs_l
complex(kind(0d0)) :: obs_f
real(8) :: Sb

call calc_pcsc_site(obs_s,Phi,Dirac_inv)
call calc_pcsc_link(obs_l,Phi,UMAT,Dirac_inv)
call calc_pcsc_face(obs_f,Phi,UMAT,Dirac_inv)

call calc_bosonic_action(Sb,UMAT,Phi)



obs=cmplx(Sb) + obs_s + obs_l + obs_f
end subroutine calc_pcsc
 

end module observables
