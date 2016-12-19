!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ver.01: correct the normalization of the PCSC relation 
!! ver.04: bug fix for mass part of PCSC
!! ver.05: include WT id. in naive quench
!! ver.06: added compensator for SU(2) (triple cover version)
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
  !trace_abs = sqrt( trace * conjg(trace) ) 
  trace_abs = abs(trace)
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
complex(kind(0d0)) :: phi_mat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: comm(1:NMAT,1:NMAT)
complex(kind(0d0)) :: comm_a

integer :: s,t,a,b

!obs=(0d0,0d0)
!do s=1,sizeD
  !do t=1,sizeD
    !obs=obs+(Dirac_inv(s,t)+Dirac_inv(t,s))&
      !*conjg( Dirac_inv(s,t)+Dirac_inv(t,s)) 
  !enddo
!enddo
!write(*,*) obs


obs=(0d0,0d0)
do s=1,num_sites
  call make_traceless_matrix_from_modes(phi_mat, NMAT, Phi(:,s))
  call matrix_commutator(comm, phi_mat, phi_mat, 'N', 'C')
  do a=1,dimG
    call trace_MTa(comm_a, comm, a, NMAT)
    do t=1,num_sites
      do b=1,dimG
        obs = obs &
          + (0.25d0, 0d0) * cmplx(alpha_s(s) * overall_factor**2 * mass_square_phi) &
            * comm_a * Phi(b,t) &
            * Dirac_inv(site_index(a,s), site_index(b,t)) 
      enddo
    enddo
  enddo
enddo

end subroutine calc_pcsc_site
        

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! link part of the PCSC
subroutine calc_pcsc_link(obs,Phi,UMAT,Dirac_inv)
use SUN_generators, only : make_traceless_matrix_from_modes, trace_MTa
!use matrix_functions, only : matrix_product
!use Dirac_operator, only : make_Dirac
implicit none

complex(kind(0d0)), intent(out) :: obs
complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: Dirac_inv(1:sizeD,1:sizeD)
complex(kind(0d0)) :: bPhi(1:dimG,1:num_sites)
complex(kind(0d0)) :: dphi_mat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: dphi_a

!complex(kind(0d0)) :: Dirac(1:sizeD,1:sizeD), prod(1:sizeD,1:sizeD)

integer :: s,l,a,b

!call make_Dirac(Dirac,UMAT,Phi)
!call matrix_product(prod,Dirac,Dirac_inv,"N","N")
!do s=1,sizeD
  !write(*,*) abs(prod(s,s))
!enddo

bPhi=conjg(Phi)

obs=(0d0,0d0)
do l=1,num_links
  !tmp=(0d0,0d0)
  call make_diff_phi(dphi_mat,l,UMAT,bPhi)
  do a=1,dimG
    call trace_MTa(dphi_a,dphi_mat,a,NMAT)
    do s=1,num_sites
      do b=1,dimG
        obs=obs &
          + (0d0, -1d0) * cmplx(alpha_l(l) * overall_factor**2 * mass_square_phi) &
            * dphi_a * Phi(b,s) &
            * Dirac_inv(link_index(a,l),site_index(b,s)) 
      enddo
    enddo
  enddo
enddo

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
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT), Ufm(1:NMAT,1:NMAT),Omega(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Omega_a

integer :: s,f,a,b

obs=(0d0,0d0)
do f=1,num_faces
  call Make_face_variable(Uf,f,UMAT)
  call matrix_power(Ufm,Uf,m_omega)
  call make_moment_map(Omega,Ufm)
  do a=1,dimG
    call trace_MTa(Omega_a,Omega,a,NMAT)
    do s=1,num_sites
      do b=1,dimG
        obs=obs &
          + (0d0, -0.5d0) * cmplx(alpha_f(f) * beta_f(f) &
                                * overall_factor**2 * mass_square_phi) &
            * Omega_a * Phi(b,s) &
            * Dirac_inv(face_index(a,f),site_index(b,s)) 
      enddo
    enddo
  enddo
enddo

end subroutine calc_pcsc_face


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! PCSC mass part
!!  pcsc_mass = "mu^2/2g^2 \Sigma*Tr(Phi \eta)"
subroutine calc_pcsc_mass(pcsc_mass,Phi,UMAT,Dirac_inv)
implicit none

complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: Dirac_inv(1:sizeD,1:sizeD)
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
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Naive WT identity
!! deviation of the WT identity in naive quench 
subroutine QbarAtr_sigma(obs,Phi,UMAT,Dirac_inv)
use SUN_generators, only : make_traceless_matrix_from_modes, trace_MTa
use matrix_functions, only : matrix_power,  matrix_commutator
implicit none

complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: Dirac_inv(1:sizeD,1:sizeD)
complex(kind(0d0)), intent(out) :: obs

complex(kind(0d0)) :: bPhi(1:dimG,1:num_sites)
complex(kind(0d0)) :: tmp, trbPHi, bPhi_eta_Sigma
complex(kind(0d0)) :: phi_mat(1:NMAT,1:NMAT), comm(1:NMAT,1:NMAT)
complex(kind(0d0)) :: dphi_mat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: comm_b, dphi_b
real(8) :: abs_z, arg_z

complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT), Ufm(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Omega(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Omega_b

integer :: euler
integer :: a,b,s,l,f,t

euler = num_sites - num_links + num_faces
bPhi=conjg( Phi )

obs=(0d0,0d0)
do s=1,num_sites
  !! (1/Ns tr( \bar\Phi_s^2 ) )^{ - dimG \chi / 4 - 1 }
  trbPhi=(0d0,0d0)
  do a=1,dimG
    trbPhi = trbPhi + bPhi(a,s) * bPhi(a,s)
  enddo
  trbPhi = trbPhi / cmplx(dble(NMAT))
  abs_z = abs( trbPhi ) 
  arg_z = arg( trbPhi ) 

  trbPhi = abs_z**( -dble(dimG*euler)/4d0 - 1d0 ) &
    * exp( (0d0,1d0) * arg_z * ( -dble(dimG*euler)/4d0 - 1d0 ) )

  !! bPhi_eta_Sigma = 1/NMAT* \sum_a( \bar\Phi_s^a \eta_s^a \Sigma )
  bPhi_eta_Sigma=(0d0,0d0)
  do a=1,dimG
    ! \Sigma_s
    do t=1,num_sites
      call make_traceless_matrix_from_modes(phi_mat, NMAT, Phi(:,t))
      call matrix_commutator(comm, phi_mat, phi_mat, 'N', 'C')
      do b=1,dimG
        call trace_MTa(comm_b, comm, b, NMAT)
          bPhi_eta_Sigma = bPhi_eta_Sigma &
            + cmplx(alpha_s(t) * overall_factor * 0.25d0 ) * comm_b & 
              * bPhi(a,s) &
              * Dirac_inv(site_index(a,s), site_index(b,t)) 
      enddo
    enddo
    ! \Sigma_l
    do l=1,num_links
      call make_diff_phi(dphi_mat,l,UMAT,bPhi)
      do b=1,dimG
        call trace_MTa(dphi_b,dphi_mat,b,NMAT)
        bPhi_eta_Sigma = bPhi_eta_Sigma &
          + cmplx(alpha_l(l) * overall_factor) * (0d0,-1d0) * dphi_b &
            * bPhi(a,s) &
            * Dirac_inv(site_index(a,s), link_index(b,l)) 
      enddo
    enddo
    ! ¥Sigma_f
    do f=1,num_faces
      call Make_face_variable(Uf,f,UMAT)
      call matrix_power(Ufm,Uf,m_omega)
      call make_moment_map(Omega,Ufm)
      do b=1,dimG
        call trace_MTa(Omega_b,Omega,b,NMAT)
        bPhi_eta_Sigma = bPhi_eta_Sigma &
          + cmplx(alpha_f(f) * beta_f(f) * overall_factor) * (0d0,-0.5d0) * Omega_b &
            * bPhi(a,s) &
            * Dirac_inv(site_index(a,s), face_index(b,f)) 
      enddo
    enddo
  enddo
  obs = obs + trbPhi * bPhi_eta_Sigma / cmplx(dble( NMAT ))
enddo
obs = obs * cmplx( dble( -dimG*euler ) / dble( 2*num_sites ) )
end subroutine QbarAtr_sigma


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Sf part 
subroutine partial_Sf(obs_site, obs_link, obs_face,Dirac,Dirac_inv)
implicit none

complex(kind(0d0)), intent(out) :: obs_site, obs_link, obs_face
complex(kind(0d0)), intent(in) :: Dirac(1:sizeD,1:sizeD), Dirac_inv(1:sizeD,1:sizeD)
integer :: s, t, l, l2, f, f2, a, b

!!!
obs_site=(0d0,0d0)
do s=1,num_sites
  do t=1,num_sites
    do a=1,dimG
      do b=1,dimG
        obs_site = obs_site + (0.5d0,0d0) &
          * Dirac( site_index(a,s), site_index(b,t) ) &
          * Dirac_inv( site_index(b,t), site_index(a,s) )
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
          * Dirac( link_index(a,l), site_index(b,t) ) &
          * Dirac_inv( site_index(b,t), link_index(a,l) )
        !!
        obs_link = obs_link + (0.5d0,0d0) &
          * Dirac( site_index(a,t), link_index(b,l) ) &
          * Dirac_inv( link_index(b,l), site_index(a,t) )
      enddo
    enddo
  enddo
enddo
do l=1,num_links
  do l2=1,num_links
    do a=1,dimG
      do b=1,dimG
        obs_link = obs_link + (0.5d0,0d0) &
          * Dirac( link_index(a,l), link_index(b,l2) ) &
          * Dirac_inv( link_index(b,l2), link_index(a,l) )
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
          * Dirac( face_index(a,f), link_index(b,l) ) &
          * Dirac_inv( link_index(b,l), face_index(a,f) )
        !!
        obs_face = obs_face + (0.5d0,0d0) &
          * Dirac( link_index(a,l), face_index(b,f) ) &
          * Dirac_inv( face_index(b,f), link_index(a,l) )
      enddo
    enddo
  enddo
enddo
do f=1,num_faces
  do f2=1,num_faces
    do a=1,dimG
      do b=1,dimG
        obs_face = obs_face + (0.5d0,0d0) &
          * Dirac( face_index(a,f), face_index(b,f2) ) &
          * Dirac_inv( face_index(b,f2), face_index(a,f) )
      enddo
    enddo
  enddo
enddo
end subroutine partial_Sf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compensator_trace_SU2 (using eigenvalue)
!! VARID ONLY FOR SU(2)
subroutine calc_compensator_trace_SU2(C,Phi)
use SUN_generators, only : make_traceless_matrix_from_modes
use matrix_functions, only : matrix_eigenvalues
implicit none

complex(kind(0d0)), intent(out) :: C
complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)) :: phi_mat(1:NMAT,1:NMAT), eigen(1:NMAT)
real(8) :: r, theta
!! data of simplicial complex
integer :: s,a,num

num = -3*(num_sites-num_links+num_faces) / 2


C=(0d0,0d0)
do s=1,num_sites
  call make_traceless_matrix_from_modes(phi_mat, NMAT, Phi(:,s))
  call matrix_eigenvalues(eigen,phi_mat)
  r=abs(eigen(1))
  theta=arg(eigen(1))
  C=C+cmplx( r**num )* exp( (0d0,1d0) * cmplx( dble(num) * theta ) )
enddo
C=C/ cmplx(dble(num_sites))

end subroutine calc_compensator_trace_SU2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compensator_IZ1 
subroutine calc_compensator_IZ1(compensator,Phi,UMAT,Dirac_inv)
use SUN_generators, only : make_traceless_matrix_from_modes, TtimesM
use matrix_functions, only : matrix_product, matrix_inverse, matrix_power
implicit none

complex(kind(0d0)), intent(out) :: compensator
complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: Dirac_inv(1:sizeD,1:sizeD)
complex(kind(0d0)) :: phi_tip(1:NMAT,1:NMAT), phi_org(1:NMAT,1:NMAT)
complex(kind(0d0)) :: UphiUdag(1:NMAT,1:NMAT),tmpmat(1:NMAT,1:NMAT), base_mat(1:NMAT,1:NMAT)
!! data of simplicial complex
complex(kind(0d0)) :: f_abc
integer :: l,n,a,b,c,k, i,j

k = -dimG*(num_sites-num_links+num_faces) / 2

compensator=(0d0,0d0)
do l=1,num_links
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! construct 2 ¥Phi.U.¥Phi.U^¥dag + ¥lambda_l¥lambda_l (U.¥Phi.U^dag + ¥Phi)
  call make_traceless_matrix_from_modes(phi_org,NMAT,Phi(:,link_org(l)))
  call make_traceless_matrix_from_modes(phi_tip,NMAT,Phi(:,link_tip(l)))

  call matrix_product(tmpmat,UMAT(:,:,l),phi_tip,'N','N')
  call matrix_product(UphiUdag,tmpmat,UMAT(:,:,l),'N','C')

  ! 2 ¥Phi U ¥Phi U^¥dagger
  call matrix_product(base_mat, (2d0,0d0)*phi_org, UphiUdag,'N','N')
  
  do n=1,NZF
    a=NZF_index(1,n)
    b=NZF_index(2,n)
    c=NZF_index(3,n)
    f_abc=cmplx(NZF_value(n))

    call TtimesM(tmpmat, UphiUdag+phi_org, c, NMAT)

    base_mat=base_mat &
      + (0d0,0.5d0) * f_abc * Dirac_inv( link_index(a,l), link_index(b,l) ) &
        * tmpmat
  enddo
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! tmpmat = ¥PHi_org^{k-2}
  if( k-2 < 0 ) then
    call matrix_inverse(phi_org)
  endif
  call matrix_power(tmpmat,phi_org,abs(k-2))
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do i=1,NMAT
    do j=1,NMAT
      compensator=compensator+tmpmat(i,j)*base_mat(j,i)
    enddo
  enddo

enddo

compensator = compensator / cmplx( dble( num_links * NMAT ) )
end subroutine calc_compensator_IZ1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compensator_IZ2 
subroutine calc_compensator_IZ2(compensator,Phi,UMAT,Dirac_inv)
use SUN_generators, only : make_traceless_matrix_from_modes, TtimesM
use matrix_functions, only : matrix_product, matrix_inverse, matrix_power
implicit none

complex(kind(0d0)), intent(out) :: compensator
complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: Dirac_inv(1:sizeD,1:sizeD)
complex(kind(0d0)) :: phi_tip(1:NMAT,1:NMAT), phi_org(1:NMAT,1:NMAT)
complex(kind(0d0)) :: UphiUdag(1:NMAT,1:NMAT),tmpmat(1:NMAT,1:NMAT), base_mat(1:NMAT,1:NMAT)
!! data of simplicial complex
complex(kind(0d0)) :: f_abc, ctmp
real(8) :: abs_c, arg_c
integer :: l,n,a,b,c,k, i,j

k = -dimG*(num_sites-num_links+num_faces) / 2

compensator=(0d0,0d0)
do l=1,num_links
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! construct 2 ¥Phi.U.¥Phi.U^¥dag + ¥lambda_l¥lambda_l (U.¥Phi.U^dag + ¥Phi)
  call make_traceless_matrix_from_modes(phi_org,NMAT,Phi(:,link_org(l)))
  call make_traceless_matrix_from_modes(phi_tip,NMAT,Phi(:,link_tip(l)))

  call matrix_product(tmpmat,UMAT(:,:,l),phi_tip,'N','N')
  call matrix_product(UphiUdag,tmpmat,UMAT(:,:,l),'N','C')

  ! 2 ¥Phi U ¥Phi U^¥dagger
  call matrix_product(base_mat, (2d0,0d0)*phi_org, UphiUdag,'N','N')
  
  do n=1,NZF
    a=NZF_index(1,n)
    b=NZF_index(2,n)
    c=NZF_index(3,n)
    f_abc=cmplx(NZF_value(n))

    call TtimesM(tmpmat, UphiUdag+phi_org, c, NMAT)

    base_mat=base_mat &
      + (0d0,0.5d0) * f_abc * Dirac_inv( link_index(a,l), link_index(b,l) ) &
        * tmpmat
  enddo
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ctmp=(0d0,0d0)
  do i=1,NMAT
    ctmp=ctmp+base_mat(i,i)
  enddo
  ctmp=ctmp / cmplx(dble(NMAT))
  abs_c=abs(ctmp)
  arg_c=arg(ctmp)
  ctmp=cmplx( abs_c**( dble(k) / 2d0 ) ) &
         *exp( (0d0,1d0) * cmplx(arg_c* dble(k) / 2d0) )

  compensator = compensator + ctmp / cmplx(dble( num_links ))
enddo

end subroutine calc_compensator_IZ2


end module observables
