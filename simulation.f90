!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! module for molecular dynamics
module simulation
use global_parameters
use global_subroutines
use hamiltonian
use forces
use mt95
implicit none

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine test_hamiltonian(UMAT,PhiMat)
use mt95
use Dirac_operator
implicit none

complex(kind(0d0)), intent(inout) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)), intent(inout) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)

complex(kind(0d0)) :: P_Phi(1:dimG,1:num_sites)
double precision :: P_A(1:dimG,1:num_links)
complex(kind(0d0)) :: PF(1:sizeD)
!complex(kind(0d0)) :: Phi_BAK(1:dimG,1:num_sites)
complex(kind(0d0)) :: PhiMat_BAK(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: UMAT_BAK(1:NMAT,1:NMAT,1:num_links)
double precision :: Hold,Hnew
integer :: n
integer :: seed,CGite,info
type(genrand_state) :: state
type(genrand_srepr) :: srepr

!do s=1,num_sites
!call make_traceless_matrix_from_modes(PhiMat(:,:,s),NMAT,Phi(:,s))
!enddokk

call check_Dirac(UMAT,PhiMat)
write(*,*) "# test hamiltonian"
write(*,*) "# Ntau*Dtau=",Ntau*Dtau

!! backup
PhiMat_BAK=PhiMat
!Phi_BAK=Phi
UMAT_BAK=UMAT

call genrand_init( get=state )
n=1
do while ( n > 0 ) 
Ntau=Ntau*2
Dtau_phi=Dtau_phi/2d0
Dtau_A=Dtau_A/2d0

!! fix the random seed
!if( fix_seed == 1 ) then
  !seed=12345
call genrand_init( put=state )
!endif
!! set random momentum
call set_randomP(P_A,P_Phi)
! produce pseudo-fermion
call make_pseudo_fermion(PF,UMAT,PhiMat)
!! calculate Hamiltonian 
call Make_Hamiltonian(Hold,CGite,info,UMAT,PhiMat,PF,P_A,P_Phi)
!! molecular evolution
!call molecular_evolution(UMAT,Phi,PF,P_A,P_Phi,info)
call molecular_evolution_Omelyan(UMAT,PhiMat,PF,P_A,P_Phi,info)
  !write(*,*) UMAT-UMAT_BAK
  !write(*,*) "!!!"
  !write(*,*) Phi-Phi_BAK
!! calculate Hamiltonian 

!do s=1,num_sites
!call make_traceless_matrix_from_modes(PhiMat(:,:,s),NMAT,Phi(:,s))
!enddo
call Make_Hamiltonian(Hnew,CGite,info,UMAT,PhiMat,PF,P_A,P_Phi)
!! metropolice
write(*,'(I5,e20.10)') Ntau, dabs(Hnew-Hold)!, Hold, Hnew
!! return to the original values
PhiMat=PhiMat_bak
!Phi=Phi_Bak
UMAT=UMAT_bak

!srepr=state ! mt95では"="がassignmentされている
!write(*,*) srepr%repr
enddo


end subroutine test_hamiltonian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine HybridMonteCarlo(UMAT,PhiMat,seed,total_ite)
use output
use SUN_generators, only : make_traceless_matrix_from_modes,trace_MTa
use SUN_generators, only : trace_MTa
implicit none

complex(kind(0d0)), intent(inout) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)), intent(inout) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
integer, intent(inout) :: seed
integer, intent(inout) :: total_ite

!complex(kind(0d0)) :: P_A(1:dimG,1:num_links)
double precision :: P_A(1:dimG,1:num_links)
complex(kind(0d0)) :: P_Phi(1:dimG,1:num_sites)
complex(kind(0d0)) :: P_PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: PF(1:sizeD)

!complex(kind(0d0)) Phi_BAK(1:dimG,1:num_sites)
complex(kind(0d0)) PhiMat_BAK(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) UMAT_BAK(1:NMAT,1:NMAT,1:num_links)
double precision Hold, Hnew
integer :: ite
integer :: accept
type(genrand_state) :: state
type(genrand_srepr) :: srepr

integer :: t_start, t_end, t_rate, t_max
integer :: CGite1, CGite2, info1, info2, info
double precision :: diff

integer s,a

!! prepare intermediate file
open(unit=MED_CONF_FILE,status='replace',file=Fmedconf,action='write',form='unformatted')
call write_basic_info_to_medfile(MED_CONF_FILE)

!! prepare output file
open(unit=OUTPUT_FILE,status='replace',file=Foutput,action='write')

write(*,*) "# start Monte Carlo simulation"

accept=0
  do s=1,num_sites
    do a=1,dimG
      call trace_MTa(Phi(a,s),PhiMat(:,:,s),a,NMAT)
    enddo
  enddo
call write_header(seed,UMAT,PhiMat)
!! measure the time used in this simulation
call system_clock(t_start)
do ite=total_ite+1,total_ite+num_ite
  !! set random momentuet
  call set_randomP(P_A,P_Phi)
  do s=1,num_sites
    call make_traceless_matrix_from_modes(P_PhiMat(:,:,s),NMAT,P_Phi(:,s))
  enddo
  !! produce pseudo-fermion
  call make_pseudo_fermion(PF,UMAT,PhiMat)
  !! calculate Hamiltonian 
  call Make_Hamiltonian(Hold,CGite1,info1,UMAT,PhiMat,PF,P_A,P_PhiMat)
  !! backup
  PhiMat_BAK=PhiMat
  !Phi_BAK=Phi
  UMAT_BAK=UMAT
  !! molecular evolution
  call molecular_evolution_Omelyan(UMAT,PhiMat,PF,P_A,P_Phi,info)
  !do s=1,num_sites
  !call make_traceless_matrix_from_modes(PhiMat(:,:,s),NMAT,Phi(:,s))
  !enddo
  if( info == 1 ) then
    PhiMat=PhiMat_BAK
    !Phi=Phi_BAK
    UMAT=UMAT_BAK
    write(*,*) "### CAUTION: CG iterations reaches to the maximal."
  else
    !! check distance of Uf from the origin
    !call check_vacuum(UMAT)
    !! calculate Hamiltonian 
  do s=1,num_sites
    call make_traceless_matrix_from_modes(P_PhiMat(:,:,s),NMAT,P_Phi(:,s))
  enddo
    call Make_Hamiltonian(Hnew,CGite2,info2,UMAT,PhiMat,PF,P_A,P_PhiMat)
    !! metropolice
    call Metropolice_test(Hnew-Hold,PhiMat_Bak,UMAT_BAK,PhiMat,UMAT,accept)
  endif
  !! write out the configuration
  !do s=1,num_sites
    !do a=1,dimG
      !call trace_MTa(Phi(a,s),PhiMat(:,:,s),a,NMAT)
    !enddo
  !enddo
   if ( mod(ite,config_step) == 0 ) then
       write(MED_CONF_FILE) ite
       write(MED_CONF_FILE) UMAT
       write(MED_CONF_FILE) PHIMAT
       endif
   !! write out the observables 
   if ( mod(ite,obs_step) == 0 ) then
       call write_observables(PhiMat,UMAT,ite,accept,Hnew-Hold,total_ite,CGite1)
   endif
enddo

!! write the simulation time
call system_clock(t_end, t_rate, t_max)
if ( t_end < t_start ) then
  diff = dble((t_max - t_start) + t_end + 1) / dble(t_rate)
else
  diff = dble(t_end - t_start) / dble(t_rate)
endif
write(OUTPUT_FILE,*) "# Total time: ",diff,"[s]"
write(*,*) "# Total time: ",diff,"[s]"


!! write the final configuration 
  !do s=1,num_sites
    !do a=1,dimG
      !call trace_MTa(Phi(a,s),PhiMat(:,:,s),a,NMAT)
    !enddo
  !enddo
open(unit=OUT_CONF_FILE,status='replace',file=Fconfigout,action='write',form='unformatted')
write(OUT_CONF_FILE) ite-1
write(OUT_CONF_FILE) UMAT
write(OUT_CONF_FILE) PhiMat
call genrand_init( get=state )
srepr=state ! mt95では"="がassignmentされている
write(OUT_CONF_FILE) srepr%repr


close( MED_CONF_FILE )
close( OUTPUT_FILE )

end subroutine HybridMonteCarlo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine writedown_config_action_and_fores(UMAT,PhiMat,seed)
use hamiltonian
use SUN_generators, only : make_traceless_matrix_from_modes,trace_MTa
implicit none

complex(kind(0d0)), intent(inout) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(inout) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
integer, intent(inout) :: seed

complex(kind(0d0)) :: PF(1:sizeD)
double precision :: SB_S,SB_L,SB_F,SB_M, SF,SB
complex(kind(0d0)) :: dSdPhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: dSdPhiMat_final(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: dSdPhi(1:dimG,1:num_sites)
double precision :: dSdA(1:dimG,1:num_links)
double precision :: dSdA_final(1:dimG,1:num_links)
complex(kind(0d0)) :: tmp
integer :: info,CGIte
integer :: i,j,a,s,l

dSdPhiMat_final=(0d0,0d0)
dSdPhi=(0d0,0d0)
dSdA_final=0d0
!! write down UMAT
write(*,*) "# link, i, j, UMAT(i,j,link)"
do l=1,num_links
  do i=1,NMAT
    do j=1,NMAT
      write(*,*) l,i,j,UMAT(i,j,l)
    enddo
  enddo
enddo
write(*,*)
write(*,*) "###################"
write(*,*) "# site, a, Phi(a,site) (a=1..dim(G))"
do s=1,num_sites
  do a=1,dimG
    call trace_MTa(tmp,PhiMat(:,:,s),a,NMAT)
    write(*,*) s,a,tmp
  enddo
enddo
!! produce pseudo-fermion
call make_pseudo_fermion(PF,UMAT,PhiMat)
write(*,*) "###################"
write(*,*) "# I, PseudoFermion (eta_{s,a},lambda_{l,a},chi_{f,a}の順番, a=1..dim(G))"
do a=1,sizeD
    write(*,*) a, PF(a)
enddo


 
CGite=0
info=0
!SB_T=0d0
SB_M=0d0
SB_S=0d0
SB_L=0d0
SB_F=0d0
SB=0d0
SF=0d0

!! actions
call bosonic_action_mass(SB_M,PhiMat)
call bosonic_action_site(SB_S,PhiMat)
call bosonic_action_link(SB_L,UMAT,PhiMat)
call bosonic_action_face(SB_F,UMAT)
SB=SB_M+SB_S+SB_F
call fermionic_action(SF,CGite,info,UMAT,PhiMat,PF)

write(*,'(a)',advance='no') "SB_mass = "
write(*,*) SB_M
write(*,'(a)',advance='no') "SB_site = "
write(*,*) SB_S
write(*,'(a)',advance='no') "SB_link = "
write(*,*) SB_L
write(*,'(a)',advance='no') "SB_face = "
write(*,*) SB_F
write(*,'(a)',advance='no') "S_fermion = "
write(*,*) SF


!call Make_force(dSdPhi,dSdA,UMAT,Phi,PF,info)

dSdPhiMat=(0d0,0d0)
call Make_bosonic_force_Phi_mass(dSdPhiMat,PhiMat)
dSdPhiMat_final=dSdPhiMat_final+dSdPhiMat
dSdPhiMat=(0d0,0d0)
call Make_bosonic_force_Phi_site(dSdPhiMat,PhiMat)
dSdPhiMat_final=dSdPhiMat_final+dSdPhiMat
dSdPhiMat=(0d0,0d0)
call Make_bosonic_force_Phi_link(dSdPhiMat,UMAT,PhiMat)
dSdPhiMat_final=dSdPhiMat_final+dSdPhiMat
write(*,*) "##############################"
write(*,*) "#site,a,dSdPhi from SB"
do s=1,num_sites
  do a=1,dimG
    call trace_MTa(tmp,dSdPhiMat_final(:,:,s),a,NMAT)
    write(*,*) s,a,2d0*conjg(tmp)
  enddo
enddo

dSdA=0d0
call Make_bosonic_force_A_link(dSdA,UMAT,PhiMat)
dSdA_final=dSdA_final+dSdA
dSdA=0d0
call Make_bosonic_force_A_face(dSdA,UMAT)
dSdA_final=dSdA_final+dSdA
write(*,*) "##############################"
write(*,*) "#link,a,dSdA from SB"
do l=1,num_links
  do a=1,dimG
    write(*,*) l,a,dSdA_final(a,l)
  enddo
enddo

dSdPhi=(0d0,0d0)
dSdA=0d0
call Make_fermionic_force(dSdPhi,dSdA,UMAT,PhiMat,PF,info)
write(*,*) "##############################"
write(*,*) "#site,a,dSdPhi from SF"
do s=1,num_sites
  do a=1,dimG
    write(*,*) s,a,2d0*conjg(dSdPhi(a,s))
  enddo
enddo
write(*,*) "##############################"
write(*,*) "#link,a,dSdA from SF"
do l=1,num_links
  do a=1,dimG
    write(*,*) l,a,dSdA(a,l)
  enddo
enddo



end subroutine writedown_config_action_and_fores

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_randomP(P_A,P_Phi)
implicit none

double precision, intent(inout) :: P_A(1:dimG,1:num_links)
complex(kind(0d0)), intent(inout) :: P_Phi(1:dimG,1:num_sites)
double precision, allocatable :: gauss(:)
integer s,l,f,a
integer i

!! set random momentum
call BoxMuller(gauss,(dimG)*num_sites)
do s=1,num_sites
  do a=1,dimG
    i=dimG*(s-1)+a
    P_Phi(a,s)=dcmplx(dsqrt(0.5d0)*gauss(2*i-1)) &
        + im_unit*dcmplx(dsqrt(0.5d0)*gauss(2*i))
  enddo
enddo
call BoxMuller( gauss,(dimG*num_links+1)/2 )
do l=1,num_links
  do a=1,dimG
    i=dimG*(l-1)+a
    P_A(a,l)=gauss(i)
  enddo
enddo


end subroutine set_randomP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Molecular Evolution by Omelyan integrator
subroutine molecular_evolution_Omelyan(UMAT,PhiMat,PF,P_A,P_Phi,info)
use SUN_generators, only : make_traceless_matrix_from_modes, trace_mta
!use matrix_functions, only : MATRIX_EXP
!use observables, only : calc_min_and_max_of_eigenvalues_Dirac
implicit none

double precision, parameter :: lambda=0.1931833275d0
complex(kind(0d0)), intent(inout) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)), intent(inout) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(in) :: PF(1:sizeD)
double precision, intent(inout) :: P_A(1:dimG,1:num_links)
complex(kind(0d0)), intent(inout) :: P_Phi(1:dimG,1:num_sites)
integer, intent(inout) :: info
complex(kind(0d0)) :: dSdPhi(1:dimG,1:num_sites)
double precision :: dSdA(1:dimG,1:num_links)

integer :: s,l,f,a
integer :: i
integer :: j,k,ii

complex(kind(0d0)) :: minimal, maximal

do s=1,num_sites
  do a=1,dimG
    call trace_MTa(Phi(a,s),PhiMat(:,:,s),a,NMAT)
  enddo
enddo

!! first step
! momentum
call update_momentum(P_Phi,P_A,PhiMat,UMAT,PF,info,Dtau_phi*lambda,Dtau_A*lambda)
if( info == 1 ) return

!! main steps 
!!        Val:Dtau/2 
!!    --> Mom:Dtau*(1-2\lambda)*Dtau
!!    --> Val:Dtau/2
!!    --> Mom:Dtau*2\lambda*Dtau
do i=1, Ntau-1
  !write(*,*) "i= ",i
! variables
  call update_PhiMat(PhiMat,Phi,P_phi,Dtau_phi*0.5d0)
  call update_UMAT(UMAT,P_A,0.5d0*Dtau_A)
! momentum
  call update_momentum(P_Phi,P_A,PhiMat,UMAT,PF,info,Dtau_phi*(1d0-2d0*lambda),Dtau_A*(1d0-2d0*lambda))
  if( info == 1 ) return
! variables
  call update_PhiMat(PhiMat,Phi,P_phi,Dtau_phi*0.5d0)
  call update_UMAT(UMAT,P_A,0.5d0*Dtau_A)
! momentum
  call update_momentum(P_Phi,P_A,PhiMat,UMAT,PF,info,Dtau_phi*2d0*lambda,Dtau_A*2d0*lambda)
  if ( info == 1 ) return
enddo

!! final step
! variables
call update_PhiMat(PhiMat,Phi,P_phi,Dtau_phi*0.5d0)
call update_UMAT(UMAT,P_A,Dtau_A*0.5d0)
! momentum
call update_momentum(P_Phi,P_A,PhiMat,UMAT,PF,info,Dtau_phi*(1d0-2d0*lambda),Dtau_A*(1d0-2d0*lambda))
if ( info == 1 ) return 
! variables
call update_PhiMat(PhiMat,Phi,P_phi,Dtau_phi*0.5d0)
call update_UMAT(UMAT,P_A,Dtau_A*0.5d0)
! momentum
call update_momentum(P_Phi,P_A,PhiMat,UMAT,PF,info,Dtau_phi*(lambda),Dtau_A*(lambda))
if ( info == 1 ) return

end subroutine molecular_evolution_Omelyan


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Molecular Evolution by Leap Frog integrator
subroutine molecular_evolution(UMAT,Phi,PF,P_A,P_Phi,info)
!use matrix_functions, only : MATRIX_EXP
implicit none

complex(kind(0d0)), intent(inout) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(inout) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)), intent(in) :: PF(1:sizeD)
double precision, intent(inout) :: P_A(1:dimG,1:num_links)
complex(kind(0d0)), intent(inout) :: P_Phi(1:dimG,1:num_sites)
integer, intent(inout) :: info

complex(kind(0d0)) :: dSdPhi(1:dimG,1:num_sites)
double precision :: dSdA(1:dimG,1:num_links)
integer :: s,l,f,a
integer :: i
integer :: j,k,ii



!! freeze UMAT 
!P_A=0d0   !<== これを入れるとダメになる。なぜ？？


!! first step
! momentum
call Make_force(dSdPhi,dSdA,UMAT,Phi,PF,info)
if( info == 1 ) return
do s=1,num_sites
  do a=1,dimG
    P_Phi(a,s)=P_Phi(a,s) - (0.5d0,0d0) * Dtau_phi * dSdPhi(a,s)
  enddo
enddo
do l=1,num_links
  do a=1,dimG
    P_A(a,l)=P_A(a,l) - 0.5d0 * Dtau_A * dSdA(a,l)
  enddo
enddo

!! main steps
do i=1, Ntau-1
! variables
  do s=1,num_sites
    do a=1,dimG
      Phi(a,s)=Phi(a,s) + Dtau_phi * dconjg( P_phi(a,s) )
    enddo
  enddo
  call update_UMAT(UMAT,P_A,Dtau_A)
! momentum
  call Make_force(dSdPhi,dSdA,UMAT,Phi,PF,info)
  if( info == 1 ) return
  do s=1,num_sites
    do a=1,dimG
      P_Phi(a,s)=P_Phi(a,s) -  Dtau_phi * dSdPhi(a,s)
    enddo
  enddo
  do l=1,num_links
    do a=1,dimG
      P_A(a,l)=P_A(a,l) - Dtau_A * dSdA(a,l)
    enddo
  enddo
enddo

!! final step
! variables
do s=1,num_sites
  do a=1,dimG
    Phi(a,s)=Phi(a,s) + Dtau_phi * dconjg( P_phi(a,s) )
  enddo
enddo
call update_UMAT(UMAT,P_A,Dtau_A)
! momentum
call Make_force(dSdPhi,dSdA,UMAT,Phi,PF,info)
if( info == 1 ) return
do s=1,num_sites
  do a=1,dimG
    P_Phi(a,s)=P_Phi(a,s) - 0.5d0 * Dtau_phi * dSdPhi(a,s)
  enddo
enddo
do l=1,num_links
  do a=1,dimG
    P_A(a,l)=P_A(a,l) - 0.5d0 * Dtau_A * dSdA(a,l)
  enddo
enddo
end subroutine molecular_evolution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! update UMAT
subroutine update_UMAT(UMAT,P_A,Dtau)
use matrix_functions, only : MATRIX_EXP
use SUN_generators, only : make_traceless_matrix_from_modes
implicit none

complex(kind(0d0)),intent(inout) :: UMAT(1:NMAT,1:NMAT,1:num_links)
double precision,intent(in) :: P_A(1:dimG,1:num_links)
double precision, intent(in) :: Dtau
complex(kind(0d0)) :: Plmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: dU(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
integer :: l

do l=1,num_links
  !write(*,*) "test000"
  call Make_traceless_matrix_from_modes(tmpmat, NMAT, dcmplx(P_A(:,l)))
  Plmat=(0d0,1d0)*Dtau*tmpmat
  !write(*,*) "test001"
  !call MATRIX_EXP(NMAT,Plmat,dU)
  call MATRIX_EXP(dU,Plmat)
  tmpmat=UMAT(:,:,l)
  !write(*,*) "test002"
  call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
    dU, NMAT, &
    tmpmat, NMAT, &
    (0d0,0d0), UMAT(:,:,l), NMAT)
  !write(*,*) "test003"
enddo


end subroutine update_UMAT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! update PhiMat
subroutine update_PhiMat(PhiMat,Phi,P_phi,Dtau)
use SUN_generators, only : make_traceless_matrix_from_modes
implicit none

complex(kind(0d0)),intent(inout) :: PhiMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)),intent(inout) :: Phi(1:dimG,1:num_links)
complex(kind(0d0)),intent(in) :: P_Phi(1:dimG,1:num_links)
double precision, intent(in) :: Dtau
integer :: s,a

do s=1,num_sites
  do a=1,dimG
    Phi(a,s)=Phi(a,s) + Dtau * dconjg( P_phi(a,s) )
  enddo
enddo
do s=1,num_sites
  call make_traceless_matrix_from_modes(PhiMat(:,:,s),NMAT,Phi(:,s))
enddo

end subroutine update_PhiMat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! update momentum
subroutine update_momentum(P_Phi,P_A,PhiMat,UMAT,PF,info,deltaPPhi,deltaA)
implicit none

double precision, intent(inout) :: P_A(1:dimG,1:num_links)
complex(kind(0d0)), intent(inout) :: P_Phi(1:dimG,1:num_sites)
complex(kind(0d0)), intent(inout) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(inout) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(in) :: PF(1:sizeD)
double precision, intent(in) :: deltaPPhi,deltaA
integer, intent(out) :: info

complex(kind(0d0)) :: dSdPhi(1:dimG,1:num_sites)
double precision :: dSdA(1:dimG,1:num_links)
integer :: s,a,l

call Make_force(dSdPhi,dSdA,UMAT,PhiMat,PF,info)
if( info == 1) return
do s=1,num_sites
  do a=1,dimG
    P_Phi(a,s)=P_Phi(a,s) - dSdPhi(a,s) * deltaPPhi 
  enddo
enddo
do l=1,num_links
  do a=1,dimG
    P_A(a,l)=P_A(a,l) - dSdA(a,l) * deltaA
  enddo
enddo

end subroutine update_momentum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! make pseudo fermion
subroutine make_pseudo_fermion(PF,UMAT,PhiMat)
!use matrix_functions,  only : BoxMuller2
use SUN_generators, only : trace_MTa
use Dirac_operator
use rational_algorithm
implicit none


complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)), intent(out) :: PF(1:sizeD)

complex(kind(0d0)) :: gauss(1:sizeD)
double precision :: gauss2(1:2*sizeD),rtmp
integer :: i,j,info,CGite

!integer :: s,a
!do s=1,num_sites
  !do a=1,dimG
    !call trace_MTa(Phi(a,s),PhiMat(:,:,s),a,NMAT)
  !enddo
!enddo

call BoxMuller2(gauss2,sizeD)
do i=1,sizeD
  gauss(i)=(gauss2(2*i-1) + gauss2(2*i)*(0d0,1d0)) * dcmplx(dsqrt(0.5d0))
enddo

call calc_matrix_rational_power(&
  PF, gauss, sizeD, N_Remez8, epsilon, CG_max, info, CGite, &
  Remez_alpha8, Remez_beta8, UMAT,PhiMat, prod_DdagD)

end subroutine make_pseudo_fermion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Metropolice_test(delta_Ham,PhiMat_Bak,UMAT_BAK,PhiMat,UMAT,accept)
implicit none

double precision, intent(in) :: delta_Ham
!complex(kind(0d0)), intent(in) :: Phi_BAK(1:dimG,1:num_sites)
complex(kind(0d0)), intent(in) :: PhiMat_BAK(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(in) :: UMAT_BAK(1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)), intent(inout) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)), intent(inout) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(inout) :: UMAT(1:NMAT,1:NMAT,1:num_links)
integer, intent(inout) :: accept

double precision :: random_num

if( delta_Ham <= 0d0 ) then 
  accept=accept+1
else
  call genrand_real3(random_num)

  if( exp( -delta_Ham ) > random_num ) then
    accept=accept+1
  else
    !Phi=Phi_BAK
    PhiMat=PhiMat_BAK
    UMAT=UMAT_BAK
  endif
endif
return
end subroutine Metropolice_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine check_vacuum(UMAT)
use matrix_functions, only : matrix_norm, make_unit_matrix
!use global_subroutines
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
integer l,i,j,f
double precision norm
complex(kind(0d0)) UNITM(1:NMAT,1:NMAT)
complex(kind(0d0)) tmp(1:NMAT,1:NMAT)
complex(kind(0d0)) Uf(1:NMAT,1:NMAT)
double precision min_dist
double precision dist(1:NMAT-1)

!write(*,*) "===== check distance from 1 ==========="
if( NMAT <= 4 ) then 
  min_dist=2d0*sqrt(2d0)
else
  min_dist=dsin(PI/dble(NMAT))*2d0*sqrt(2d0)
endif


call make_unit_matrix(UNITM)
do f=1,num_faces
  call Make_face_variable(Uf,f,UMAT)
  tmp=UNITM-Uf
  !call matrix_norm(tmp,NMAT,norm)
  call matrix_norm(norm,tmp)
  write(*,*) f,min_dist-norm
enddo


end subroutine check_vacuum


end module simulation
