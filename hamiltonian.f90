!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! module to compute Hamiltonian
module hamiltonian
use global_parameters
use global_subroutines
implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Hamiltonian
!! Hamiltonian is defined by 
!!   S_b + S_f + P^2
!! where P^2 is the conjugate momenta. 
!! S_f is connected with Dirac matrix D as
!!  S_f = 1/2 \Psi_A D_{AB} \Psi_B
!! with 
!!  D_{AB}=-D_{BA}
!! Make sure that the fermionic action includes 
!! the prefactor 1/2.
subroutine Make_Hamiltonian(Htotal,CGite,info,UMAT,PhiMat,PF,P_A,P_Phi)
!use SUN_generators, only : make_traceless_matrix_from_modes
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: PhiMAT(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)), intent(in) :: PF(1:sizeD)
complex(kind(0d0)), intent(in) :: P_Phi(1:dimG,1:num_sites)
double precision, intent(in) :: P_A(1:dimG,1:num_links)
integer, intent(inout) :: CGite,info

double precision, intent(out) :: Htotal
integer a,s,l
double precision :: SB_S,SB_L,SB_F,SB_M, SF !,SB_T

!complex(kind(0d0)):: PhiMat(1:NMAT,1:NMAT,1:num_sites)
!do s=1,num_sites
!call make_traceless_matrix_from_modes(PhiMat(:,:,s),NMAT,Phi(:,s))
!enddo

CGite=0
info=0
!SB_T=0d0
SB_M=0d0
SB_S=0d0
SB_L=0d0
SB_F=0d0
SF=0d0
call bosonic_action_mass(SB_M,PhiMat)
call bosonic_action_site(SB_S,PhiMat)
call bosonic_action_link(SB_L,UMAT,PhiMat)
call bosonic_action_face(SB_F,UMAT)
call fermionic_action(SF,CGite,info,UMAT,PhiMat,PF)
!call bosonic_action_test(SB_T,UMAT)

Htotal = SB_M
Htotal = Htotal+SB_S
Htotal = Htotal+SB_L
Htotal = Htotal+SB_F
Htotal = Htotal+SF 
!Htotal=Htotal+SB_T

! over all factor
Htotal = Htotal 

do s=1,num_sites
  do a=1,dimG
    Htotal = Htotal &
      + dble(P_phi(a,s)*dconjg(P_phi(a,s)))
  enddo
enddo

do l=1,num_links
  do a=1,dimG
    Htotal = Htotal &
      + P_A(a,l)*P_A(a,l)*0.5d0
  enddo
enddo

end subroutine Make_Hamiltonian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Parts
!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! test
!subroutine bosonic_action_test(SB_T,UMAT)
!use matrix_functions, only : matrix_power
!implicit none
!
!double precision, intent(out) :: SB_T
!complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
!
!complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT,1:num_faces)
!complex(kind(0d0)) :: Ufm(1:NMAT,1:NMAT,1:num_faces)
!complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
!integer :: f,i,j
!integer :: m_tmp
!
!m_tmp=m_omega
!do f=1,num_faces
!  call Make_face_variable(Uf(:,:,f),f,UMAT) 
!  call matrix_power(NMAT,Uf(:,:,f),m_tmp,Ufm(:,:,f))
!enddo
!
!tmpmat=(0d0,0d0)
!SB_T=0d0
!do f=1,num_faces
!  do i=1,NMAT
!    do j=1,NMAT
!      tmpmat(i,j)=(0d0,-1d0)*(Ufm(i,j,f)-dconjg( Ufm(j,i,f) ))
!    enddo
!  enddo
!  do i=1,NMAT
!    do j=1,NMAT
!      SB_T=SB_T+0.5d0*dble(tmpmat(i,j)*tmpmat(j,i))
!    enddo
!  enddo
!enddo
!
!end subroutine bosonic_action_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! mass action
subroutine bosonic_action_mass(SB_M,PhiMat)
implicit none

complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
double precision, intent(out) :: SB_M
integer :: s,i,j


SB_M=0d0
do s=1,num_sites
  do i=1,NMAT
    do j=1,NMAT
    SB_M=SB_M+dble(PhiMat(i,j,s)*dconjg(PhiMat(i,j,s)))
    enddo
  enddo
enddo
SB_M=SB_M*mass_square_phi & 
  * overall_factor
  !* one_ov_2g2N
end subroutine bosonic_action_mass

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! site action
subroutine bosonic_action_site(SB_S,PhiMat)
use matrix_functions, only : matrix_commutator
implicit none

complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
double precision, intent(out) :: SB_S

complex(kind(0d0)), allocatable :: tmpmat(:,:)
complex(kind(0d0)) :: phi_mat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: comm(1:NMAT,1:NMAT)
complex(kind(0d0)) :: ctmp
double precision :: tmp
integer s
integer a,b,c
integer i,j

SB_S=0d0

do s=1,num_sites
  !call make_traceless_matrix_from_modes(phi_mat,NMAT,Phi(:,s))
  call matrix_commutator(comm,phimat(:,:,s),phimat(:,:,s),'N','C')
  ctmp=(0d0,0d0)
  do i=1,NMAT
    do j=1,NMAT
      ctmp=ctmp+0.25d0*comm(i,j)*comm(j,i)
    enddo
  enddo
  SB_S=SB_S+alpha_s(s)*dble(ctmp)
enddo

SB_S=SB_S * overall_factor !* one_ov_2g2N

!! compute [\phi_s,\bar^phi_s]^a
!allocate(tmpmat(1:dimG,1:num_sites) )
!tmpmat=(0d0,0d0)
!do s=1,num_sites
!do i=1,NZF
!  a=NZF_index(1,i)
!  b=NZF_index(2,i)
!  c=NZF_index(3,i)
!  tmpmat(a,s)=tmpmat(a,s)+im_unit*NZF_value(i)*Phi(b,s)*dconjg(Phi(c,s))
!enddo
!!do a=1,dimG
!!write(*,*) tmpmat(a,s)
!!enddo
!enddo
!
!! S_b^{(s)}=\alpha_s tr( 1/4 [\phi_s, \bar\phi_s]^2 )
!do s=1,num_sites
!  tmp=0d0
!  do a=1,dimG
!    tmp=tmp+dble(tmpmat(a,s)*tmpmat(a,s))
!  enddo
!  SB_S=SB_S+alpha_s(s)*0.25d0*tmp
!enddo

end subroutine bosonic_action_site

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! link action
subroutine bosonic_action_link(SB_L,UMAT,PhiMat)
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
double precision, intent(out) :: SB_L
integer :: i,j
integer :: l
double precision :: tmp
!complex(kind(0d0)) :: Phi_tip(1:NMAT,1:NMAT), Phi_org(1:NMAT,1:NMAT)
complex(kind(0d0)) :: dPhi(1:NMAT,1:NMAT)


SB_L=0d0
do l=1,num_links
  call Make_diff_PhiMat(dPhi, l,UMAT,PhiMat)

  tmp=0d0
  do i=1,NMAT
  do j=1,NMAT
    tmp=tmp+dble(dPhi(j,i)*dconjg(dPhi(j,i)))
  enddo
  enddo

  SB_L=SB_L+alpha_l(l)*tmp
enddo

SB_L=SB_L*overall_factor !*one_ov_2g2N

end subroutine bosonic_action_link


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! face action
subroutine bosonic_action_face(SB_F,UMAT)
use matrix_functions, only : matrix_power
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
double precision, intent(out) :: SB_F

integer :: f
complex(kind(0d0)) :: Omega(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT), Ufm(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmp
integer :: i,j

SB_F=0d0
do f=1,num_faces
  call Make_face_variable(Uf,f,UMAT)
  call matrix_power(Ufm,Uf,m_omega)
  call Make_moment_map(Omega,Ufm)
  tmp=(0d0,0d0)
  do j=1,NMAT
    do i=1,NMAT
      tmp=tmp+Omega(i,j)*dconjg(Omega(i,j))
    enddo
  enddo
  !write(*,*) tmp
  SB_F=SB_F+0.25d0*alpha_f(f)*beta_f(f)*beta_f(f)*dble(tmp)
enddo
SB_F=SB_F * overall_factor !*one_ov_2g2N
end subroutine bosonic_action_face

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! fermionic action 
subroutine fermionic_action(SF,CGite,info,UMAT,PhiMat,PF)
use rational_algorithm
use Dirac_operator
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
!complex(kind(0d0)):: Phi(1:dimG,1:num_sites)
complex(kind(0d0)), intent(in) :: PF(1:sizeD)
double precision, intent(out) :: SF
integer, intent(inout) :: info,CGite

complex(kind(0d0)) :: SF_CG
complex(kind(0d0)) :: chi(1:sizeD,1:N_Remez4)
complex(kind(0d0)) :: DdagD4_PF(1:sizeD)
integer :: i,r
!! for test
integer, parameter :: test=0 ! set 1 in testing
complex(kind(0d0)) :: chi_direct(1:sizeD,1:N_Remez4)
complex(kind(0d0)) :: DdagD(1:sizeD,1:sizeD),tmpmat(1:sizeD,1:sizeD)
complex(kind(0d0)) :: DdagD_m1ov4(1:sizeD,1:sizeD)
complex(kind(0d0)) :: DdagD4_PF_direct(1:sizeD)
complex(kind(0d0)) :: DdagD4_PF_direct2(1:sizeD)
complex(kind(0d0)) :: SF_direct
complex(kind(0d0)) :: tmpvec(1:sizeD),tmp,tmp2
double precision :: distance
integer :: j
!integer :: s,a
!do s=1,num_sites
  !do a=1,dimG
    !call trace_MTa(Phi(a,s),PhiMat(:,:,s),a,NMAT)
  !enddo
!enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! computation using CG
!  (D\dag D)^{-1/4}.PF by CG
!call calc_matrix_rational_power(&
!  DdagD4_PF, PF, sizeD, N_Remez4, epsilon, CG_max, info, CGite, &
!  Remez_alpha4, Remez_beta4, prod_DdagD)
call mmBiCG( chi, PF, Remez_beta4, sizeD, N_Remez4, epsilon, CG_max, info, &
             CGite, UMAT, PhiMat, Prod_DdagD )
DdagD4_PF=Remez_alpha4(0)*PF
do r=1,N_Remez4
  DdagD4_PF=DdagD4_PF + Remez_alpha4(r)*chi(:,r)
enddo
! SF
SF_CG=(0d0,0d0)
do i=1,sizeD
  SF_CG=SF_CG + dconjg(PF(i)) * DdagD4_PF(i) 
enddo 
! changed 2015/11/08
!SF=dble(SF_CG) / overall_factor ! / one_ov_2g2N
SF=dble(SF_CG) 


!if (test==1) then 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! direct computation 
!!  (D\dag D)^(-1/4) 
!  call make_DdagD(DdagD,UMAT,Phi)
!  call matrix_rational_power(DdagD_m1ov4,DdagD,-0.25d0,sizeD)
!  ! (D\dag D)^(-1/4).PF
!  DdagD4_PF_direct=(0d0,0d0)
!  do i=1,sizeD
!    do j=1,sizeD
!      DdagD4_PF_direct(i)=DdagD4_PF_direct(i)+DdagD_m1ov4(i,j)*PF(j)
!    enddo
!  enddo
!  ! SF
!  SF_direct=(0d0,0d0)
!  do i=1,sizeD
!    SF_direct=SF_direct + dconjg(PF(i)) * DdagD4_PF_direct(i) 
!  enddo 
!  write(*,*) "rational CG (D\dag D)^(-1/4):",SF_CG
!  write(*,*) "direct (D\dagD)^(-1/4): )",SF_direct
!  !!!!!!!!!!!!!!!!!!!!!!
!  
!  !!!!!!!!!!!!!!!!!!!!!
!  ! check (D\dag D + \beta_r)^(-1) . PF
!  chi_direct=(0d0,0d0)
!  call make_DdagD(DdagD,UMAT,Phi)
!  DdagD4_PF_direct2=Remez_alpha4(0)*PF
!  do r=1,N_Remez4
!    tmpmat=DdagD
!    do i=1,sizeD
!      tmpmat(i,i)=tmpmat(i,i)+Remez_beta4(r)
!    enddo
!    call matrix_inverse(sizeD,tmpmat)
!    tmpvec=(0d0,0d0)
!    do i=1,sizeD
!      do j=1,sizeD
!        chi_direct(i,r)=chi_direct(i,r)+tmpmat(i,j)*PF(j)
!      enddo
!    enddo
!  enddo
!  do r=1,N_Remez4
!    DdagD4_PF_direct2=DdagD4_PF_direct2+Remez_alpha4(r)*chi_direct(:,r)
!  enddo
!  
!  distance=0d0
!  do r=1,N_Remez4
!    do i=1,sizeD
!      tmp=chi(i,r)-chi_direct(i,r)
!      distance=distance+dble( tmp * dconjg(tmp) )
!  enddo
!  enddo
!  write(*,*) "chi vs chi_direct:", distance
!  
!  distance=0d0
!  tmp2=(0d0,0d0)
!  do i=1,sizeD
!    tmp=DdagD4_PF(i) - DdagD4_PF_direct2(i)
!    distance=distance+dble( tmp * dconjg(tmp))
!  enddo
!  write(*,*) "CG vs direct of rational (D\dag D)-(-1/4):", dsqrt(distance)
!
!endif


end subroutine fermionic_action


end module hamiltonian

