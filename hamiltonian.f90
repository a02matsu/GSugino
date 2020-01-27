!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!complex(kind(0d0)) Phi_BAK(1:dimG,1:num_sites)
!! module to compute Hamiltonian
!module hamiltonian
!use global_parameters
!use global_subroutines
!implicit none

!contains

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
subroutine Make_Hamiltonian(Htotal,CGite,info,UMAT,PhiMat,PF_eta,PF_lambda,PF_chi,P_AMat,P_PhiMat)
!use SUN_generators, only : make_traceless_matrix_from_modes
#ifdef PARALLEL
use parallel
#endif
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMAT(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: PF_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: PF_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PF_chi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), intent(in) :: P_PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(in) :: P_AMat(1:NMAT,1:NMAT,1:num_links)
integer, intent(inout) :: CGite,info

double precision, intent(out) :: Htotal
double precision :: Hlocal,tmp
integer a,s,l,i,j
double precision :: SB_S,SB_L,SB_F,SB_M, SF !,SB_T

!Hlocal=site_abs(P_PhiMat(:,:,1:num_sites))
!if(MYRANK==0) write(*,*) "P_Phi:",Hlocal
!Hlocal=link_abs(P_AMat(:,:,1:num_links))
!if(MYRANK==0) write(*,*) "P_Amat:",Hlocal

CGite=0
info=0
!SB_T=0d0
SF=0d0
Hlocal=0d0
Htotal=0d0
if(pb_mass==0) call bosonic_action_mass(SB_M,PhiMat)
if(pb_site==0) call bosonic_action_site(SB_S,PhiMat)
if(pb_link==0) call bosonic_action_link(SB_L,UMAT,PhiMat)
if(pb_face==0) call bosonic_action_face(SB_F,UMAT)
if(pf==0) call fermionic_action(SF,CGite,info,UMAT,PhiMat,PF_eta,PF_lambda,PF_chi)

Hlocal = SB_M+SB_S+SB_L+SB_F+SF 

do s=1,num_sites
  !do a=1,dimG
  do j=1,NMAT
    do i=1,NMAT
      Hlocal = Hlocal &
        + dble( dconjg(P_PhiMat(i,j,s)) * P_PhiMat(i,j,s) )
    enddo
  enddo
enddo

do l=1,num_links
  do j=1,NMAT
    do i=1,NMAT
      Hlocal = Hlocal + dble( P_AMat(i,j,l) * dcmplx( P_AMat(j,i,l) ))*0.5d0
    enddo
  enddo
enddo

#ifdef PARALLEL
call MPI_REDUCE(Hlocal,Htotal,1,MPI_DOUBLE_PRECISION, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)
#else
Htotal=Hlocal
#endif


end subroutine Make_Hamiltonian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! mass action
subroutine bosonic_action_mass(SB_M,PhiMat)
implicit none

complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
double precision, intent(out) :: SB_M
integer :: s,i,j


SB_M=0d0
do s=1,num_sites
  do i=1,NMAT
    do j=1,NMAT
    SB_M=SB_M+dble(PhiMat(i,j,s)*dconjg(PhiMat(i,j,s)))!*alpha_s(s)
    enddo
  enddo
enddo
SB_M=SB_M*mass_square_phi*0.5d0 & 
  * overall_factor
end subroutine bosonic_action_mass

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! site action
subroutine bosonic_action_site(SB_S,PhiMat)
use matrix_functions, only : matrix_commutator
implicit none

complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
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
  call matrix_commutator(comm,phimat(:,:,s),phimat(:,:,s),'N','C')
  ctmp=(0d0,0d0)
  do i=1,NMAT
    do j=1,NMAT
      ctmp=ctmp+(0.25d0,0d0)*comm(i,j)*comm(j,i)
    enddo
  enddo
  SB_S=SB_S+alpha_s(s)*dble(ctmp)
enddo

SB_S=SB_S * overall_factor !* one_ov_2g2N

end subroutine bosonic_action_site

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! link action
subroutine bosonic_action_link(SB_L,UMAT,PhiMat)
use matrix_functions, only : matrix_3_product
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
double precision, intent(out) :: SB_L
integer :: i,j
integer :: l
double precision :: tmp
complex(kind(0d0)) :: dPhi(1:NMAT,1:NMAT)


SB_L=0d0
do l=1,num_links
  !write(*,*) l,link_org(l),link_tip(l)
  !call matrix_3_product(dPhi,Umat(:,:,l),PhiMat(:,:,link_tip(l)),Umat(:,:,l),'N','N','C',(1d0,0d0),'ADD')
  call matrix_3_product(dPhi,Umat(:,:,l),PhiMat(:,:,link_tip(l)),Umat(:,:,l),'N','N','C')
  dPhi=dPhi*U1Rfactor(l)**2
  dPhi=dPhi-PhiMat(:,:,link_org(l))

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

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
double precision, intent(out) :: SB_F

integer :: f
complex(kind(0d0)) :: Omega(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT), Ufm(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmp
integer :: i,j

SB_F=0d0
do f=1,num_faces
  call Make_face_variable(Uf,f,UMAT)
  if(m_omega == 0) then 
    call Make_moment_map0(Omega,Uf)
  elseif(m_omega == -1) then
    call Make_moment_map_adm(Omega,Uf)
  else
    call matrix_power(Ufm,Uf,m_omega)
    call Make_moment_map(Omega,Ufm)
  endif 
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
subroutine fermionic_action(SF,CGite,info,UMAT,PhiMat,PF_eta,PF_lambda,PF_chi)
use rational_algorithm
use Dirac_operator
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)), intent(in) :: PF_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: PF_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PF_chi(1:NMAT,1:NMAT,1:num_necessary_faces)

complex(kind(0d0)) :: chi_eta(1:NMAT,1:NMAT,1:num_sites,N_Remez4)
complex(kind(0d0)) :: chi_lambda(1:NMAT,1:NMAT,1:num_links,N_Remez4)
complex(kind(0d0)) :: chi_chi(1:NMAT,1:NMAT,1:num_faces,N_Remez4)

complex(kind(0d0)) :: DdagD4_PF_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: DdagD4_PF_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: DdagD4_PF_chi(1:NMAT,1:NMAT,1:num_faces)

double precision, intent(out) :: SF
integer, intent(inout) :: info,CGite
integer :: r,s,l,f,i,j

complex(kind(0d0)) :: SF_CG
!complex(kind(0d0)) :: chi(1:sizeD,1:N_Remez4)
!complex(kind(0d0)) :: DdagD4_PF(1:sizeD)
!! for test
!integer, parameter :: test=0 ! set 1 in testing
!complex(kind(0d0)) :: chi_direct(1:sizeD,1:N_Remez4)
!complex(kind(0d0)) :: DdagD(1:sizeD,1:sizeD),tmpmat(1:sizeD,1:sizeD)
!complex(kind(0d0)) :: DdagD_m1ov4(1:sizeD,1:sizeD)
!complex(kind(0d0)) :: DdagD4_PF_direct(1:sizeD)
!complex(kind(0d0)) :: DdagD4_PF_direct2(1:sizeD)
!complex(kind(0d0)) :: SF_direct
!complex(kind(0d0)) :: tmpvec(1:sizeD),tmp,tmp2
!double precision :: distance


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! computation using CG
!  (D\dag D)^{-1/4}.PF by CG
!call calc_matrix_rational_power(&
!  DdagD4_PF, PF, sizeD, N_Remez4, epsilon, CG_max, info, CGite, &
!  Remez_alpha4, Remez_beta4, prod_DdagD)
!call vec_to_mat(PF_eta,PF_lambda,PF_chi,PF)
call mmBiCG( &
  chi_eta,chi_lambda,chi_chi,&
  PF_eta,PF_lambda,PF_chi, &
  Remez_beta4, epsilon, CG_max, info, CGite, UMAT, PhiMat, Prod_DdagD)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! if CG is failed, we should reject
if( info==1 ) return
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DdagD4_PF_eta=Remez_alpha4(0)*PF_eta
DdagD4_PF_lambda=Remez_alpha4(0)*PF_lambda
DdagD4_PF_chi=Remez_alpha4(0)*PF_chi
do r=1,N_Remez4
  DdagD4_PF_eta=DdagD4_PF_eta + Remez_alpha4(r)*chi_eta(:,:,:,r)
  DdagD4_PF_lambda=DdagD4_PF_lambda + Remez_alpha4(r)*chi_lambda(:,:,:,r)
  DdagD4_PF_chi=DdagD4_PF_chi + Remez_alpha4(r)*chi_chi(:,:,:,r)
enddo
! SF
SF_CG=(0d0,0d0)
do s=1,num_sites
  do j=1,NMAT
    do i=1,NMAT
      SF_CG=SF_CG + dconjg( PF_eta(i,j,s) ) * DdagD4_PF_eta(i,j,s)
    enddo
  enddo
enddo
do l=1,num_links
  do j=1,NMAT
    do i=1,NMAT
      SF_CG=SF_CG + dconjg( PF_lambda(i,j,l) ) * DdagD4_PF_lambda(i,j,l)
    enddo
  enddo
enddo
do f=1,num_faces
  do j=1,NMAT
    do i=1,NMAT
      SF_CG=SF_CG + dconjg( PF_chi(i,j,f) ) * DdagD4_PF_chi(i,j,f)
    enddo
  enddo
enddo

! changed 2015/11/08
!SF=dble(SF_CG) / overall_factor ! / one_ov_2g2N
SF=dble(SF_CG) 



end subroutine fermionic_action


!end module hamiltonian

