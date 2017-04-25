!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! module for molecular dynamics
module simulation
use global_parameters
use mt95
implicit none

!!!!!!!!!!!!!!
!! for output
integer, parameter :: num_obs=2
character(10) :: obs_name(1:num_obs)
data obs_name/ "Sb","TrX2" /
double precision :: OBS(1:num_obs) ! 1) bosonic action SB

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine test_hamiltonian(UMAT,PhiMat)
use SUN_generators, only : make_traceless_matrix_from_modes
use mt95
!use Dirac_operator
implicit none

complex(kind(0d0)), intent(inout) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)), intent(inout) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)

complex(kind(0d0)) :: P_Phi(1:dimG,1:num_sites)
complex(kind(0d0)) :: P_PhiMat(1:NMAT,1:NMAT,1:num_sites)
double precision :: P_A(1:dimG,1:num_links)
complex(kind(0d0)) :: P_AMat(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: PF(1:sizeD)
!complex(kind(0d0)) :: Phi_BAK(1:dimG,1:num_sites)
complex(kind(0d0)) :: PhiMat_BAK(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: UMAT_BAK(1:NMAT,1:NMAT,1:num_links)
double precision :: Hold,Hnew
integer :: n
integer :: seed,CGite,info,s,l,a
complex(kind(0d0)) :: tmp
type(genrand_state) :: state
type(genrand_srepr) :: srepr

!do s=1,num_sites
!call make_traceless_matrix_from_modes(PhiMat(:,:,s),NMAT,Phi(:,s))
!enddokk

call check_Dirac(UMAT,PhiMat)
write(*,'(a)') " # test hamiltonian"
write(*,'(a,I8)') " # NMAT=",NMAT
write(*,'(a,I4)') " # fix_seed=",fix_seed
write(*,'(a,I2)') " # read_alpha=",read_alpha
write(*,'(a,I3)') " # m_omega=",m_omega
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
call set_randomP(P_AMat,P_PhiMat)
! produce pseudo-fermion
call make_pseudo_fermion(PF,UMAT,PhiMat)
!! calculate Hamiltonian 
call Make_Hamiltonian(Hold,CGite,info,UMAT,PhiMat,PF,P_AMat,P_PhiMat)
!! molecular evolution
!call molecular_evolution(UMAT,Phi,PF,P_A,P_Phi,info)
call molecular_evolution_Omelyan(UMAT,PhiMat,PF,P_AMat,P_PhiMat,info)

call Make_Hamiltonian(Hnew,CGite,info,UMAT,PhiMat,PF,P_AMat,P_PhiMat)
!! metropolice
write(*,'(I5,e20.10)') Ntau, dabs(Hnew-Hold)!, Hold, Hnew
!! return to the original values
PhiMat=PhiMat_bak
UMAT=UMAT_bak

!srepr=state ! mt95では"="がassignmentされている
!write(*,*) srepr%repr
enddo


end subroutine test_hamiltonian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine HybridMonteCarlo(UMAT,PhiMat,seed,total_ite)
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
complex(kind(0d0)) :: P_AMat(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: P_Phi(1:dimG,1:num_sites)
complex(kind(0d0)) :: P_PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: PF(1:sizeD)

!complex(kind(0d0)) Phi_BAK(1:dimG,1:num_sites)
complex(kind(0d0)) PhiMat_BAK(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) UMAT_BAK(1:NMAT,1:NMAT,1:num_links)
double precision Hold, Hnew, ratio, ave_dHam
integer :: ite
integer :: accept
type(genrand_state) :: state
type(genrand_srepr) :: srepr

integer :: t_start, t_end, t_rate, t_max
integer :: CGite1, CGite2, info1, info2, info
double precision :: diff

integer s,a,l
complex(kind(0d0)) :: tmp

!! prepare intermediate file
open(unit=MED_CONF_FILE,status='replace',file=Fmedconf,action='write',form='unformatted')
call write_basic_info_to_medfile(MED_CONF_FILE)

!! prepare output file
open(unit=OUTPUT_FILE,status='replace',file=Foutput,action='write')

write(*,*) "# start Monte Carlo simulation"

ave_dHam=0d0
accept=0
call write_header(seed,UMAT,PhiMat)
!! measure the time used in this simulation
call system_clock(t_start)
do ite=total_ite+1,total_ite+num_ite
  !! set random momentuet
  call set_randomP(P_AMat,P_PhiMat)
  !! produce pseudo-fermion
  call make_pseudo_fermion(PF,UMAT,PhiMat)
  !do l=1,num_links
    !call Make_traceless_matrix_from_modes(P_Amat(:,:,l), NMAT, dcmplx(P_A(:,l)))
  !enddo
  !! calculate Hamiltonian 
  call Make_Hamiltonian(Hold,CGite1,info1,UMAT,PhiMat,PF,P_AMat,P_PhiMat)
  !! backup
  PhiMat_BAK=PhiMat
  !Phi_BAK=Phi
  UMAT_BAK=UMAT
  !! molecular evolution
  call molecular_evolution_Omelyan(UMAT,PhiMat,PF,P_AMat,P_PhiMat,info)
  do s=1,num_sites
    call traceless_projection(PhiMat(:,:,s))
  enddo
  if( info == 1 ) then
    PhiMat=PhiMat_BAK
    !Phi=Phi_BAK
    UMAT=UMAT_BAK
    write(*,*) "### CAUTION: CG iterations reaches to the maximal."
  else
    !call check_vacuum(UMAT)
    !! calculate Hamiltonian 
    call Make_Hamiltonian(Hnew,CGite2,info2,UMAT,PhiMat,PF,P_AMat,P_PhiMat)
    !! metropolice
    ave_dHam=ave_dHam+abs(Hnew-Hold)
    call Metropolice_test(Hnew-Hold,PhiMat_Bak,UMAT_BAK,PhiMat,UMAT,accept)
    !!
    !! check distance of Uf from the origin
    call check_distance(info,ratio,UMAT)
    if ( info .ne. 0 ) then 
      write(*,*) "!!! CAUTION: plaquette variables are out of proper region !!!!"
    endif
  endif
  !! write out the configuration
   if ( mod(ite,save_med_step) == 0 ) then
       write(MED_CONF_FILE) ite
       write(MED_CONF_FILE) UMAT
       write(MED_CONF_FILE) PHIMAT
       endif
   !! write out the observables 
   if ( mod(ite,obs_step) == 0 ) then
       call write_observables(PhiMat,UMAT,ite,accept,Hnew-Hold,total_ite,CGite1,ratio)
   endif
   !! write out config file
   if ( mod(ite,save_config_step) == 0 ) then
     call write_config_file(ite,UMAT,PhiMat,state,srepr)
   endif 
enddo

!! write the simulation time
call system_clock(t_end, t_rate, t_max)
if ( t_end < t_start ) then
  diff = dble((t_max - t_start) + t_end + 1) / dble(t_rate)
else
  diff = dble(t_end - t_start) / dble(t_rate)
endif
ave_dHam = ave_dHam/dble(num_ite)
write(OUTPUT_FILE,*) "# Total time: ",diff,"[s]"
write(OUTPUT_FILE,*) "# Time/ite: ",diff/dble(num_ite),"[s]"
write(OUTPUT_FILE,*) "# Time/ite/Ntau: ",diff/dble(num_ite*Ntau),"[s]"
write(OUTPUT_FILE,*) "# ave_dHam=",ave_dHam
write(OUTPUT_FILE,*) "# acceptance=",dble(accept)/dble(num_ite)
write(*,*) "# Total time: ",diff,"[s]"
write(*,*) "# Step time: ",diff/dble(num_ite),"[s]"
write(*,*) "# Time/ite/Ntau: ",diff/dble(num_ite*Ntau),"[s]"
write(*,*) "# ave_dHam=",ave_dHam
write(*,*) "# acceptance=",dble(accept)/dble(num_ite)


call write_config_file(ite,UMAT,PhiMat,state,srepr)

close( MED_CONF_FILE )
close( OUTPUT_FILE )

end subroutine HybridMonteCarlo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_config_file(ite,UMAT,PhiMat,state,srepr)
implicit none

integer, intent(in) :: ite
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
type(genrand_state), intent(inout) :: state
type(genrand_srepr), intent(inout) :: srepr

open(unit=OUT_CONF_FILE,status='replace',file=Fconfigout,action='write',form='unformatted')
write(OUT_CONF_FILE) ite-1
write(OUT_CONF_FILE) UMAT
write(OUT_CONF_FILE) PhiMat
call genrand_init( get=state )
srepr=state ! mt95では"="がassignmentされている
write(OUT_CONF_FILE) srepr%repr
close( OUT_CONF_FILE)

end subroutine write_config_file



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine traceless_projection(Mat)
implicit none

complex(kind(0d0)), intent(inout) :: MAT(:,:)
complex(kind(0d0)) trace
integer :: N,i

N=size(MAT(:,1))

trace=(0d0,0d0)
do i=1,N
  trace=trace+MAT(i,i)
enddo
do i=1,N
  MAT(i,i)=MAT(i,i)-trace/cmplx(dble( N ))
enddo

end subroutine traceless_projection


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_randomP(P_AMat,P_PhiMat)
use SUN_generators, only : make_traceless_matrix_from_modes
implicit none

double precision :: P_A(1:dimG)
complex(kind(0d0)), intent(inout) :: P_AMat(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(inout) :: P_PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: P_Phi(1:dimG)
double precision, allocatable :: gauss(:)
integer s,l,f,a
integer i

!! set random momentum
call BoxMuller(gauss,(dimG)*num_sites)
do s=1,num_sites
  do a=1,dimG
    i=dimG*(s-1)+a
    P_Phi(a)=dcmplx(dsqrt(0.5d0)*gauss(2*i-1)) &
        + im_unit*dcmplx(dsqrt(0.5d0)*gauss(2*i))
  enddo
  call make_traceless_matrix_from_modes(P_PhiMat(:,:,s),NMAT,P_Phi)
enddo

call BoxMuller( gauss,(dimG*num_links+1)/2 )
do l=1,num_links
  do a=1,dimG
    i=dimG*(l-1)+a
    P_A(a)=gauss(i)
  enddo
  call make_traceless_matrix_from_modes(P_AMat(:,:,l),NMAT,dcmplx( P_A(:) ))
enddo


end subroutine set_randomP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Molecular Evolution by Omelyan integrator
subroutine molecular_evolution_Omelyan(UMAT,PhiMat,PF,P_AMat,P_PhiMat,info)
use SUN_generators, only : make_traceless_matrix_from_modes, trace_mta
!use matrix_functions, only : MATRIX_EXP
!use observables, only : calc_min_and_max_of_eigenvalues_Dirac
implicit none

double precision, parameter :: lambda=0.1931833275d0
complex(kind(0d0)), intent(inout) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(inout) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(in) :: PF(1:sizeD)
complex(kind(0d0)), intent(inout) :: P_PhiMat(1:NMAT,1:NMAT,1:num_sites)
!double precision, intent(inout) :: P_A(1:dimG,1:num_links)
complex(kind(0d0)), intent(inout) :: P_AMat(1:NMAT,1:NMAT,1:num_links)
integer, intent(inout) :: info

integer :: s,l,f,a
integer :: i
integer :: j,k,ii

complex(kind(0d0)) :: minimal, maximal,tmp
!do l=1,num_links
  !call Make_traceless_matrix_from_modes(P_Amat(:,:,l), NMAT, dcmplx(P_A(:,l)))
!enddo

!! first step
! momentum
call update_momentum(P_PhiMat,P_AMat,PhiMat,UMAT,PF,info,Dtau_phi*lambda,Dtau_A*lambda)
if( info == 1 ) return

!! main steps 
!!        Val:Dtau/2 
!!    --> Mom:Dtau*(1-2\lambda)*Dtau
!!    --> Val:Dtau/2
!!    --> Mom:Dtau*2\lambda*Dtau
do i=1, Ntau-1
! variables
  call update_PhiMat(PhiMat,P_phiMat,Dtau_phi*0.5d0)
  call update_UMAT(UMAT,P_AMat,0.5d0*Dtau_A)
! momentum
  call update_momentum(P_PhiMat,P_AMat,PhiMat,UMAT,PF,info,Dtau_phi*(1d0-2d0*lambda),Dtau_A*(1d0-2d0*lambda))
  if( info == 1 ) return
! variables
  call update_PhiMat(PhiMat,P_phimat,Dtau_phi*0.5d0)
  call update_UMAT(UMAT,P_AMat,0.5d0*Dtau_A)
! momentum
  call update_momentum(P_PhiMat,P_AMat,PhiMat,UMAT,PF,info,Dtau_phi*2d0*lambda,Dtau_A*2d0*lambda)
  if ( info == 1 ) return
enddo

!! final step
! variables
call update_PhiMat(PhiMat,P_phimat,Dtau_phi*0.5d0)
call update_UMAT(UMAT,P_AMat,Dtau_A*0.5d0)
! momentum
call update_momentum(P_PhiMat,P_AMat,PhiMat,UMAT,PF,info,Dtau_phi*(1d0-2d0*lambda),Dtau_A*(1d0-2d0*lambda))
if ( info == 1 ) return 
! variables
call update_PhiMat(PhiMat,P_phimat,Dtau_phi*0.5d0)
call update_UMAT(UMAT,P_AMat,Dtau_A*0.5d0)
! momentum
call update_momentum(P_PhiMat,P_AMat,PhiMat,UMAT,PF,info,Dtau_phi*(lambda),Dtau_A*(lambda))
if ( info == 1 ) return

!do l=1,num_links
  !do a=1,dimG
    !call Trace_MTa(tmp,P_AMat(:,:,l),a,NMAT)
    !P_A(a,l)=dble(tmp)
  !enddo
!enddo

end subroutine molecular_evolution_Omelyan


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! update UMAT
subroutine update_UMAT(UMAT,P_AMat,Dtau)
use matrix_functions, only : MATRIX_EXP
use SUN_generators, only : make_traceless_matrix_from_modes
implicit none

complex(kind(0d0)),intent(inout) :: UMAT(1:NMAT,1:NMAT,1:num_links)
!double precision,intent(in) :: P_A(1:dimG,1:num_links)
complex(kind(0d0)), intent(in) :: P_AMat(1:NMAT,1:NMAT,1:num_links)
double precision, intent(in) :: Dtau
complex(kind(0d0)) :: Plmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: dU(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
integer :: l

!do l=1,num_links
  !call Make_traceless_matrix_from_modes(P_Amat(:,:,l), NMAT, dcmplx(P_A(:,l)))
!enddo
do l=1,num_links
  Plmat=(0d0,1d0)*Dtau*P_Amat(:,:,l)
  !call MATRIX_EXP(NMAT,Plmat,dU)
  call traceless_projection(Plmat)
  call MATRIX_EXP(dU,Plmat)
  tmpmat=UMAT(:,:,l)
  call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
    dU, NMAT, &
    tmpmat, NMAT, &
    (0d0,0d0), UMAT(:,:,l), NMAT)
enddo


end subroutine update_UMAT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! update PhiMat
subroutine update_PhiMat(PhiMat,P_phimat,Dtau)
use SUN_generators, only : make_traceless_matrix_from_modes
implicit none

complex(kind(0d0)),intent(inout) :: PhiMAT(1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)),intent(inout) :: Phi(1:dimG,1:num_links)
complex(kind(0d0)),intent(in) :: P_PhiMat(1:NMAT,1:NMAT,1:num_links)
double precision, intent(in) :: Dtau
integer :: s,i,j

do s=1,num_sites
  !do a=1,dimG
  do j=1,NMAT
    do i=1,NMAT
      PhiMat(i,j,s)=PhiMat(i,j,s)+Dtau*dconjg( P_PhiMat(j,i,s) )
      !Phi(a,s)=Phi(a,s) + Dtau * dconjg( P_phi(a,s) )
    enddo
  enddo
enddo
!do s=1,num_sites
  !call make_traceless_matrix_from_modes(PhiMat(:,:,s),NMAT,Phi(:,s))
!enddo

end subroutine update_PhiMat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! update momentum
subroutine update_momentum(P_PhiMat,P_AMat,PhiMat,UMAT,PF,info,deltaPPhi,deltaA)
use SUN_generators, only : make_traceless_matrix_from_modes, trace_mta
implicit none

double precision :: P_A(1:dimG,1:num_links)
complex(kind(0d0)), intent(inout) :: P_AMat(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(inout) :: P_PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(inout) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(inout) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(in) :: PF(1:sizeD)
complex(kind(0d0)) :: tmp
double precision, intent(in) :: deltaPPhi,deltaA
integer, intent(out) :: info

complex(kind(0d0)) :: dSdPhi(1:NMAT,1:NMAT,1:num_sites)
double precision :: dSdA_org(1:dimG,1:num_links)
complex(kind(0d0)) :: dSdA(1:NMAT,1:NMAT,1:num_links)
integer :: s,a,l

call Make_force(dSdPhi,dSdA,UMAT,PhiMat,PF,info)
if( info == 1) return
do s=1,num_sites
  P_PhiMat(:,:,s)=P_PhiMat(:,:,s) - dSdPhi(:,:,s) * deltaPPhi
enddo

!do l=1,num_links
  !call make_traceless_matrix_from_modes(dSdA(:,:,l),NMAT,dcmplx( dSdA_org(:,l) ))
!enddo
do l=1,num_links
    P_AMat(:,:,l) = P_AMat(:,:,l) - dSdA(:,:,l) * dcmplx(deltaA)
enddo

end subroutine update_momentum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! make pseudo fermion
subroutine make_pseudo_fermion(PF,UMAT,PhiMat)
!use matrix_functions,  only : BoxMuller2
use SUN_generators, only : trace_MTa
!use Dirac_operator
use rational_algorithm
implicit none


complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)), intent(out) :: PF(1:sizeD)

complex(kind(0d0)) :: gauss(1:sizeD)
double precision :: gauss2(1:2*sizeD),rtmp
integer :: i,j,info,CGite

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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Hamiltonian
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
subroutine Make_Hamiltonian(Htotal,CGite,info,UMAT,PhiMat,PF,P_AMat,P_PhiMat)
!use SUN_generators, only : make_traceless_matrix_from_modes
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: PhiMAT(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(in) :: PF(1:sizeD)
!complex(kind(0d0)), intent(in) :: P_Phi(1:dimG,1:num_sites)
complex(kind(0d0)), intent(in) :: P_PhiMat(1:NMAT,1:NMAT,1:num_sites)
!double precision, intent(in) :: P_A(1:dimG,1:num_links)
complex(kind(0d0)), intent(in) :: P_AMat(1:NMAT,1:NMAT,1:num_links)
integer, intent(inout) :: CGite,info

double precision, intent(out) :: Htotal
integer a,s,l,i,j
double precision :: SB_S,SB_L,SB_F,SB_M, SF !,SB_T


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
  !do a=1,dimG
  do j=1,NMAT
    do i=1,NMAT
      Htotal = Htotal &
        + dble( conjg(P_PhiMat(i,j,s)) * P_PhiMat(i,j,s) )
        !+ dble(P_phi(a,s)*dconjg(P_phi(a,s)))
    enddo
  enddo
enddo

do l=1,num_links
  do j=1,NMAT
    do i=1,NMAT
      Htotal = Htotal + dble( P_AMat(i,j,l) * dcmplx( P_AMat(j,i,l) ))*0.5d0
    enddo
  enddo
enddo

end subroutine Make_Hamiltonian

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
  if(m_omega == 0) then 
    call Make_moment_map0(Omega,Uf)
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
subroutine fermionic_action(SF,CGite,info,UMAT,PhiMat,PF)
use rational_algorithm
!use Dirac_operator
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
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


end subroutine fermionic_action


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! END Hamiltonian
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Begin Force
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! dS/Dphi 
subroutine Make_force(dSdPhi,dSdA,UMAT,PhiMat,PF,info)
use SUN_generators, only : trace_mta, make_traceless_matrix_from_modes
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(in) :: PF(1:sizeD)
complex(kind(0d0)), intent(out) :: dSdPhi(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(out) :: dSdA(1:NMAT,1:NMAT,1:num_links)
integer, intent(inout) :: info

complex(kind(0d0)) :: dSdPhi_boson_mass(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: dSdPhi_boson_site(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: dSdPhi_boson_link(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: dSdPhi_fermion(1:NMAT,1:NMAT,1:num_sites)
!!!
complex(kind(0d0)) :: dSdA_boson_link(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: dSdA_boson_face(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: dSdA_fermion(1:NMAT,1:NMAT,1:num_links)

integer :: s,a,ii,jj,l
complex(kind(0d0)) :: tmp

dSdPhi=(0d0,0d0)
dSdPhi_boson_mass=(0d0,0d0)
dSdPhi_boson_site=(0d0,0d0)
dSdPhi_boson_link=(0d0,0d0)
dSdPhi_fermion=(0d0,0d0)
dSdA=(0d0,0d0)
dSdA_boson_face=(0d0,0d0)
dSdA_boson_link=(0d0,0d0)
dSdA_fermion=(0d0,0d0)
!! force for Phi from boson
call Make_bosonic_force_Phi_mass(dSdPhi_boson_mass,PhiMat)
call Make_bosonic_force_Phi_site(dSdPhi_boson_site,PhiMat)
call Make_bosonic_force_Phi_link(dSdPhi_boson_link,UMAT,PhiMat)
!! force from fermion
call Make_fermionic_force(dSdPhi_fermion,dSdA_fermion,UMAT,PhiMat,PF,info)

!! force for A from boson
call Make_bosonic_force_A_link(dSdA_boson_link,UMAT,PhiMat)
call Make_bosonic_force_A_face(dSdA_boson_face,UMAT)
!dSdA_boson_face=(0d0,0d0)

dSdPhi= dSdPhi_boson_mass &   ! mass part
+ dSdPhi_boson_site  & ! site part
+ dSdPhi_boson_link  &  ! link part
+ dSdPhi_fermion   ! link part

dSdA = dSdA_boson_link  &
+ dSdA_boson_face &
+ dSdA_fermion

end subroutine Make_force


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Forces from mass terms, Sites, Links, Faces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! dS/dPhi from mass action
subroutine Make_bosonic_force_Phi_mass(dSdPhi_boson_mass,PhiMat)
implicit none

complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(out) :: dSdPhi_boson_mass(1:NMAT,1:NMAT,1:num_sites)
integer :: s,i,j

dSdPhi_boson_mass=(0d0,0d0)
do s=1,num_sites
  do j=1,NMAT
    do i=1,NMAT
  dSdPhi_boson_mass(i,j,s) = dcmplx(mass_square_phi) &!*  dconjg(Phi(a,s)) &
    * dconjg(PhiMat(j,i,s)) * dcmplx(overall_factor)
    !* cmplx(one_ov_2g2N)
    enddo
  enddo
enddo
end subroutine Make_bosonic_force_Phi_mass

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! dS/Dphi from site action
!complex(kind(0d0)) function dSdPhi_boson_site(a,s)
subroutine Make_bosonic_force_Phi_site(dSdPhi_boson_site,PhiMat)
use matrix_functions, only : matrix_commutator
implicit none

complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(out) :: dSdPhi_boson_site(1:NMAT,1:NMAT,1:num_sites)
integer :: a,s
complex(kind(0d0)) :: comm1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: comm2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmp(1:dimG), tmp1,tmp2
integer :: i,j

dSdPhi_boson_site=(0d0,0d0)

do s=1,num_sites
  call matrix_commutator(comm1,phimat(:,:,s),phimat(:,:,s),'N','C')
  call matrix_commutator(comm2,phimat(:,:,s),comm1,'C','N')
  !do a=1,dimG
  do j=1,NMAT
    do i=1,NMAT
    dSdPhi_boson_site(i,j,s)=dcmplx(0.5d0*alpha_s(s))*comm2(i,j)& 
      * cmplx( overall_factor )
      !* cmplx( one_ov_2g2N )
    enddo
  enddo
enddo

end subroutine Make_bosonic_force_Phi_site

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! dS/Dphi from link action
!complex(kind(0d0)) function dSdPhi_boson_link(a,s)
subroutine Make_bosonic_force_Phi_link(dSdPhi_boson_link,UMAT,PhiMat)
use matrix_functions, only : matrix_3_product
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(out) :: dSdPhi_boson_link(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: dPhi(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT), tmpmat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Ud_dPhibar_U(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: dPhibar(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: trace
integer :: a,s
integer :: l, l_i
integer :: i,j,k

dSdPhi_boson_link=(0d0,0d0)
do l=1,num_links
! diff( \Phi_s ) for the link l
  call Make_diff_PhiMat(dPhi,l,UMAT,PhiMat)
! U_l^\dagger diff( \bar\Phi_s ) U_l for all links 
  call matrix_3_product(Ud_dPhibar_U(:,:,l),UMAT(:,:,l),dPhi,UMAT(:,:,l),'C','C','N')

  do j=1,NMAT
    do i=1,NMAT
      dPhibar(i,j,l)=conjg(dPhi(j,i))
    enddo
  enddo
enddo

dSdPhi_boson_link=(0d0,0d0)
do s=1,num_sites
  do i=1,linkorg_to_s(s)%num_
    l_i=linkorg_to_s(s)%Labels_(i)
    do k=1,NMAT
      do j=1,NMAT
        dSdPhi_boson_link(j,k,s)=dSdPhi_boson_link(j,k,s) &
          + alpha_l(l_i) *  Ud_dPhibar_U(j,k,l_i)
      enddo
    enddo
  enddo

  do i=1,linktip_from_s(s)%num_
    l_i=linktip_from_s(s)%Labels_(i)
    do k=1,NMAT
      do j=1,NMAT
        dSdPhi_boson_link(j,k,s)=dSdPhi_boson_link(j,k,s) - alpha_l(l_i) &
          * dPhibar(j,k,l_i)
      enddo
    enddo
  enddo
enddo

dSdPhi_boson_link = dSdPhi_boson_link  &
  * cmplx( overall_factor )

end subroutine Make_bosonic_force_Phi_link

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! dS/DA from link action
subroutine Make_bosonic_force_A_link(dSdA_boson_link,UMAT,PhiMat)
use SUN_generators, only : trace_MTa, make_traceless_matrix_from_modes
use matrix_functions, only : matrix_3_product, matrix_commutator
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
!double precision, intent(out) :: dSdA_boson_link_org(1:dimG,1:num_links)
complex(kind(0d0)), intent(out) :: dSdA_boson_link(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: dPhi(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Phi_tip(1:NMAT,1:NMAT)
complex(kind(0d0)) :: UPhiUinv(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Comm(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmp
integer :: a,l,i,j

dSdA_boson_link=(0d0,0d0)
do l=1,num_links
  call Make_diff_PhiMAT(dPhi, l,UMAT,PhiMat)

  Phi_tip=(0d0,1d0)*PhiMat(:,:,link_tip(l))

  ! U_l.Phi_tip
!  ! i*U_l.Phi_tip.U_l^\dagger
  call matrix_3_product(UPhiUinv,UMAT(:,:,l),Phi_tip,UMAT(:,:,l),'N','N','C')

  ! i[ UPhiUinv, dSPhi^\dagger ]
  call matrix_commutator(Comm,UPhiUinv,dPhi,'N','C')

  do j=1,NMAT
    do i=1,NMAT
      dSdA_boson_link(i,j,l) = alpha_l(l) * ( Comm(i,j) + conjg( Comm(j,i)) )
    enddo
  enddo
enddo

dSdA_boson_link=dSdA_boson_link * overall_factor !*one_ov_2g2N

end subroutine Make_bosonic_force_A_link


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! dS/DA from face action
subroutine Make_bosonic_force_A_face(dSdA_boson_face,UMAT)
use SUN_generators, only : trace_mta, make_traceless_matrix_from_modes
use matrix_functions, only : matrix_power,matrix_product
!use Dirac_operator, only : calc_Sinu_and_CosUinv
!use diferential_Dirac, only :  calc_dCosUinvdA_dSinUdA
implicit none
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
double precision dSdA_boson_face_org(1:dimG,1:num_links)
complex(kind(0d0)), intent(out) :: dSdA_boson_face(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT,1:num_faces), Ufm(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: Omega(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: diff_Omega_org(1:NMAT,1:NMAT,1:dimG)
!!!!!!!!!
complex(kind(0d0)) :: diff_Omega(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
complex(kind(0d0)) :: sinU(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: CosUinv(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: dCosUinvdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
complex(kind(0d0)) :: dSinUdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
!!!!!!!!!
!integer :: FaceSize
!integer, allocatable :: sites(:),link_labels(:),link_dirs(:)
!integer, allocatable :: faces_l(:)
complex(kind(0d0)) :: trace,tmp
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: tmp2(1:dimG)
!complex(kind(0d0)) :: comp(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:num_faces,1:num_links)
integer :: a,l,f!,ll
integer :: l_label
integer :: i,j,k,ii,jj

!Omega=(0d0,0d0)
sinU=(0d0,0d0)
cosUinv=(0d0,0d0)
do f=1,num_faces
  call Make_face_variable(Uf(:,:,f),f,UMAT) 
  if( m_omega == 0 ) then
    call Make_moment_map0(Omega(:,:,f),Uf(:,:,f))
  else
    call matrix_power(Ufm(:,:,f),Uf(:,:,f),m_omega)
    call Make_moment_map(Omega(:,:,f),Ufm(:,:,f))
    call calc_sinU_and_cosUinv(sinU(:,:,f),cosUinv(:,:,f),Ufm(:,:,f))
  endif
enddo

dSdA_boson_face=(0d0,0d0)
do l=1,num_links
  do k=1,face_in_l(l)%num_
    f=face_in_l(l)%label_(k)
    do l_label=1,links_in_f(f)%num_
      if ( l == links_in_f(f)%link_labels_(l_label) ) exit
    enddo

  if( m_omega == 0 ) then 
    dSinUdA=(0d0,0d0)
    call calc_dSinUdA(dSinUdA(:,:,:,:),Uf(:,:,f),UMAT,f,l_label)
    diff_Omega=(0d0,-1d0)*dSinUdA
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! make diff_Omega traceless
    tmpmat=(0d0,0d0)
    do ii=1,NMAT
      tmpmat=tmpmat+diff_Omega(:,:,ii,ii)
    enddo
    do ii=1,NMAT
      diff_Omega(:,:,ii,ii)=diff_Omega(:,:,ii,ii)-tmpmat/cmplx(dble( NMAT ))
    enddo
    do jj=1,NMAT
      do ii=1,NMAT
       !!!!!!!!!!!!!!!!!!!!!
       ! make traceless
       if( NMAT > 2 ) then
         trace=(0d0,0d0)
         do i=1,NMAT
           trace=trace+diff_Omega(i,i,ii,jj)
         enddo
         do i=1,NMAT
           diff_Omega(i,i,ii,jj)=diff_Omega(i,i,ii,jj) &
             - trace / cmplx(dble( NMAT ))
         enddo
       endif
       !!!!!!!!!!!!!!!!!!!!!
      enddo
    enddo
      
    do jj=1,NMAT
      do ii=1,NMAT
        do j=1,NMAT
          do i=1,NMAT
            dSdA_boson_face(ii,jj,l)=dSdA_boson_face(ii,jj,l) &
              + (0.5d0,0d0) &
                *cmplx( overall_factor * alpha_f(f)*beta_f(f)*beta_f(f) )&
                *diff_Omega(i,j,ii,jj) *  Omega(j,i,f) 
                    !+  Omega(i,j,f) *  diff_Omega(j,i,ii,jj) )
            !comp(i,j,ii,jj,f,l)=diff_Omega(i,j,ii,jj)
          enddo
        enddo
      enddo
    enddo
 
  else
    dCosUinvdA=(0d0,0d0)
    dSinUdA=(0d0,0d0)
    call calc_dCosUinvdA_dSinUdA(&
      dCosUinvdA(:,:,:,:),dSinUdA(:,:,:,:),&
      cosUinv(:,:,f),Uf(:,:,f),UMAT,f,l_label)

    diff_Omega=(0d0,0d0)
    tmpmat=(0d0,0d0)
    do jj=1,NMAT
      do ii=1,NMAT
         call matrix_product(&
           diff_Omega(:,:,ii,jj),dSinUdA(:,:,ii,jj),cosUinv(:,:,f))
         call zgemm('N','N',NMAT,NMAT,NMAT,(0d0,-1d0), &
           SinU(:,:,f), NMAT, &
           dCosUinvdA(:,:,ii,jj), NMAT, &
           (0d0,-1d0), diff_Omega(:,:,ii,jj), NMAT)
         call zgemm('N','N',NMAT,NMAT,NMAT,(0d0,-1d0), &
           dCosUinvdA(:,:,ii,jj), NMAT, &
           SinU(:,:,f), NMAT, &
           (1d0,0d0), diff_Omega(:,:,ii,jj), NMAT)
         call zgemm('N','N',NMAT,NMAT,NMAT,(0d0,-1d0), &
           CosUinv(:,:,f), NMAT, &
           dSinUdA(:,:,ii,jj), NMAT, &
           (1d0,0d0), diff_Omega(:,:,ii,jj), NMAT)
         if( ii==jj ) then
           tmpmat=tmpmat+diff_Omega(:,:,ii,ii)
         endif
       enddo
     enddo
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! make diff_Omega traceless
     do ii=1,NMAT
       diff_Omega(:,:,ii,ii)=diff_Omega(:,:,ii,ii)-tmpmat/cmplx(dble( NMAT ))
     enddo
     do jj=1,NMAT
       do ii=1,NMAT
        !!!!!!!!!!!!!!!!!!!!!
        ! make traceless
        if( NMAT > 2 ) then
          trace=(0d0,0d0)
          do i=1,NMAT
            trace=trace+diff_Omega(i,i,ii,jj)
          enddo
          do i=1,NMAT
            diff_Omega(i,i,ii,jj)=diff_Omega(i,i,ii,jj) &
              - trace / cmplx(dble( NMAT ))
          enddo
        endif
        !!!!!!!!!!!!!!!!!!!!!
      enddo
    enddo
      
      do jj=1,NMAT
        do ii=1,NMAT
          do j=1,NMAT
            do i=1,NMAT
              dSdA_boson_face(ii,jj,l)=dSdA_boson_face(ii,jj,l) &
                + (0.5d0,0d0) &
                  *cmplx( overall_factor * alpha_f(f)*beta_f(f)*beta_f(f) / dble(m_omega)) &
                  *diff_Omega(i,j,ii,jj) *  Omega(j,i,f) 
                      !+  Omega(i,j,f) *  diff_Omega(j,i,ii,jj) )
              !comp(i,j,ii,jj,f,l)=diff_Omega(i,j,ii,jj)
            enddo
          enddo
        enddo
      enddo
    endif
  enddo
enddo

end subroutine Make_bosonic_force_A_face

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! dS/dPhi and dS/dA from fermion part
subroutine Make_fermionic_force(dSdPhi_fermion,dSdA_fermion,UMAT,Phimat,PF,info)
use rational_algorithm
!use Dirac_operator
!use differential_Dirac
implicit none 

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(in) :: PF(1:sizeD)
complex(kind(0d0)), intent(out) :: dSdPhi_fermion(1:NMAT,1:NMAT,1:num_sites)
!complex(kind(0d0)) :: dSdPhi_fermion_org(1:NMAT,1:NMAT,1:num_sites)
!double precision, intent(out) :: dSdA_fermion_org(1:dimG,1:num_links)
complex(kind(0d0)), intent(out) :: dSdA_fermion(1:NMAT,1:NMAT,1:num_links)
integer, intent(inout) :: info

! vec_r = (D\dagger D + \beta_r)^{-1} F
complex(kind(0d0)) :: vec(1:sizeD,1:N_Remez4)
! Dvec_r = D.vec_r
complex(kind(0d0)) :: Dvec(1:sizeD,1:N_Remez4)
complex(kind(0d0)) :: dvec_eta(1:NMAT,1:NMAT,1:num_sites,1:N_Remez4)
complex(kind(0d0)) :: dvec_lambda(1:NMAT,1:NMAT,1:num_links,1:N_Remez4)
complex(kind(0d0)) :: dvec_chi(1:NMAT,1:NMAT,1:num_faces,1:N_Remez4)
! [\frac{dD^\dagger D}{dPhi_s^a}]_{ij} \vec_j
!complex(kind(0d0)) :: dDdagD_dPhi_vec(1:sizeD, 1:dimG,1:num_sites,1:N_Remez4)
! [\frac{dD^\dagger D}{dA_l^a}]_{ij} \vec_j
!complex(kind(0d0)) :: dDdagD_dA_vec(1:sizeD, 1:dimG,1:num_links,1:N_Remez4)
! 
complex(kind(0d0)) :: dDdPhi_vec(1:sizeD,1:dimG,1:num_sites,1:N_Remez4)
complex(kind(0d0)) :: dDdbPhi_vec(1:sizeD,1:dimG,1:num_sites,1:N_Remez4)
complex(kind(0d0)) :: dDdA_vec(1:sizeD,1:dimG,1:num_links,1:N_Remez4)
!!
complex(kind(0d0)) :: dDdPhi_eta(1:NMAT,1:NMAT,1:num_sites,1:NMAT,1:NMAT,1:num_sites,1:N_Remez4)
complex(kind(0d0)) :: dDdPhi_lambda(1:NMAT,1:NMAT,1:num_links,1:NMAT,1:NMAT,1:num_sites,1:N_Remez4)
complex(kind(0d0)) :: dDdPhi_chi(1:NMAT,1:NMAT,1:num_faces,1:NMAT,1:NMAT,1:num_sites,1:N_Remez4)
complex(kind(0d0)) :: dDdbPhi_lambda(1:NMAT,1:NMAT,1:num_links,1:NMAT,1:NMAT,1:num_sites,1:N_Remez4)
complex(kind(0d0)) :: dDdA_eta(1:NMAT,1:NMAT,1:num_sites,1:NMAT,1:NMAT,1:num_links,1:N_Remez4)
complex(kind(0d0)) :: dDdA_lambda(1:NMAT,1:NMAT,1:num_links,1:NMAT,1:NMAT,1:num_links,1:N_Remez4)
complex(kind(0d0)) :: dDdA_chi(1:NMAT,1:NMAT,1:num_faces,1:NMAT,1:NMAT,1:num_links,1:N_Remez4)

complex(kind(0d0)) :: eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: chi(1:NMAT,1:NMAT,1:num_faces)

complex(kind(0d0)) :: tmp
double precision :: rtmp
integer :: CGite
integer :: r,i,j,s,l,a,b,ss,f,ll,ii,jj

dSdPhi_fermion=(0d0,0d0)
dSdA_fermion=(0d0,0d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! comutation by CG
! compute vec_r(i,r) =  (D^\dagger D + \beta_r)^{-1} F(i)
!write(*,*) "   mmBiCG start"
!call check_Dirac(UMAT,Phi)
call mmBiCG( vec, PF, Remez_beta4, sizeD, N_Remez4, epsilon, &
             CG_max, info, CGite, UMAT, PhiMat, Prod_DdagD )
!write(*,*) "   mmBiCG end", cgite
do r=1,N_Remez4
  call vec_to_mat(eta,lambda,chi,vec(:,r))
  call prod_Dirac(Dvec(:,r),vec(:,r),sizeD,UMAT,PhiMat)
  
  !call prod_dDdPhi(dDdPhi_vec(:,:,:,r),vec(:,r),UMAT)
  call prod_dDdPhi(&
    dDdPhi_eta(:,:,:,:,:,:,r),&
    dDdPhi_chi(:,:,:,:,:,:,r), &
    eta,chi,UMAT)
  call prod_dDdbPhi(&
    dDdbPhi_lambda(:,:,:,:,:,:,r),&
    lambda,UMAT)
  !call prod_dDdA(dDdA_vec(:,:,:,r),eta,lambda,chi,UMAT,PhiMat)
  call prod_dDdA(&
    dDdA_eta(:,:,:,:,:,:,r),&
    dDdA_lambda(:,:,:,:,:,:,r), &
    dDdA_chi(:,:,:,:,:,:,r),&
    eta,lambda,chi,UMAT,PhiMat)
enddo

do r=1,N_Remez4
  call vec_to_mat( &
    Dvec_eta(:,:,:,r),&
    Dvec_lambda(:,:,:,r),&
    Dvec_chi(:,:,:,r), &
    Dvec(:,r) )
enddo

do s=1,num_sites
  do ii=1,NMAT
    do jj=1,NMAT
      do r=1,N_Remez4
        tmp=(0d0,0d0)
        do ss=1,num_sites
          do i=1,NMAT
            do j=1,NMAT
              tmp=tmp&
                +conjg( Dvec_eta(i,j,ss,r) )*dDdPhi_eta(i,j,ss,ii,jj,s,r) !&
                !+conjg( dDdbPhi_eta(i,j,ss,ii,jj,s,r) )*Dvec_eta(i,j,ss,r)
            enddo
          enddo
        enddo
        do l=1,num_links
          do i=1,NMAT
            do j=1,NMAT
              tmp=tmp&
                !+conjg( Dvec_lambda(i,j,l,r) )*dDdPhi_lambda(i,j,l,ii,jj,s,r) &
                +conjg( dDdbPhi_lambda(i,j,l,jj,ii,s,r) )*Dvec_lambda(i,j,l,r)
            enddo
          enddo
        enddo
        do f=1,num_faces
          do i=1,NMAT
            do j=1,NMAT
              tmp=tmp&
                +conjg( Dvec_chi(i,j,f,r) )*dDdPhi_chi(i,j,f,ii,jj,s,r)! &
                !+conjg( dDdbPhi_chi(i,j,f,ii,jj,s,r) )*Dvec_chi(i,j,f,r)
            enddo
          enddo
        enddo
        dSdPhi_fermion(ii,jj,s)=dSdPhi_fermion(ii,jj,s) &
          - dcmplx(Remez_alpha4(r))*tmp
      enddo
    enddo
  enddo
enddo



do ll=1,num_links
  do ii=1,NMAT
    do jj=1,NMAT
      do r=1,N_Remez4
        tmp=(0d0,0d0)
        do ss=1,num_sites
          do i=1,NMAT
            do j=1,NMAT
              tmp=tmp &
                + conjg( Dvec_eta(i,j,ss,r) ) * dDdA_eta(i,j,ss,ii,jj,ll,r)&
                + conjg( dDdA_eta(i,j,ss,jj,ii,ll,r) ) * Dvec_eta(i,j,ss,r) 
            enddo
          enddo
        enddo
        do l=1,num_links
          do i=1,NMAT
            do j=1,NMAT
              tmp=tmp &
                +conjg( Dvec_lambda(i,j,l,r) ) * dDdA_lambda(i,j,l,ii,jj,ll,r)&
                + conjg( dDdA_lambda(i,j,l,jj,ii,ll,r) ) * Dvec_lambda(i,j,l,r) 
            enddo
          enddo
        enddo
        do f=1,num_faces
          do i=1,NMAT
            do j=1,NMAT
              tmp=tmp &
                + conjg( Dvec_chi(i,j,f,r) ) * dDdA_chi(i,j,f,ii,jj,ll,r) &
                + conjg( dDdA_chi(i,j,f,jj,ii,ll,r) ) * Dvec_chi(i,j,f,r) 
            enddo
          enddo
        enddo
        dSdA_fermion(ii,jj,ll)=dSdA_fermion(ii,jj,ll) & 
          - Remez_alpha4(r) * tmp 
      enddo
    enddo
  enddo
enddo

!do l=1,num_links
  !do a=1,dimG
    !call trace_MTa(tmp,dSdA_fermion(:,:,l),a,NMAT)
    !dSdA_fermion_org(a,l)=dble(tmp)
  !enddo
!enddo
end subroutine Make_fermionic_force



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! End Force
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Begin Observables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ver.01: correct the normalization of the PCSC relation 
!! ver.04: bug fix for mass part of PCSC
!! ver.05: include WT id. in naive quench
!! ver.06: added compensator for SU(2) (triple cover version)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_bosonic_action(Sb,UMAT,PhiMat)
implicit none

double precision, intent(out) :: Sb
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)

double precision :: SB_S,SB_L,SB_F

SB_S=0d0
SB_L=0d0
SB_F=0d0
call bosonic_action_site(SB_S,PhiMat)
call bosonic_action_link(SB_L,UMAT,PhiMat)
call bosonic_action_face(SB_F,UMAT)

Sb=SB_S+SB_L+SB_F

end subroutine calc_bosonic_action


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_TrX2(TrX2, PhiMat)
implicit none

double precision, intent(out) :: TrX2
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)

integer :: s,i,j

TrX2=0d0
do s=1,num_sites
  do i=1,NMAT
    do j=1,NMAT
      TrX2=TrX2+ dble( PhiMat(i,j,s)*conjg( PhiMat(i,j,s) ) )
    enddo
  enddo
enddo

TrX2 = TrX2 / (2d0 * dble(num_sites) )

end subroutine calc_TrX2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_eigenvalues_Dirac(eigenvalues,UMAT,PhiMat)
!use Dirac_operator, only : make_Dirac
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
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

call make_Dirac(Dirac,UMAT,PhiMat)
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
subroutine calc_min_and_max_of_eigenvalues_Dirac(minimal,maximal,UMAT,PhiMat)
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(inout) :: minimal, maximal

complex(kind(0d0)) :: eigenvalues(1:sizeD)

call calc_eigenvalues_Dirac(eigenvalues,UMAT,PhiMat)

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End Observables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Begin OUTPUT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_header(seed,UMAT,PhiMat)
implicit none

integer, intent(in) :: seed 
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)

complex(kind(0d0)) :: min_eigen,max_eigen
integer :: output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compute the eigenvalues of D
call calc_smallset_and_largest_eigenvalues_of_D(min_eigen,max_eigen,UMAT,PhiMat)

!!!!!!!!!!!!!!!!!!!!!!!!!
!! for standard output
output=6 
call write_header_to(output,min_eigen,max_eigen)
call write_observable_list(output)

!!!!!!!!!!!!!!!!!!!!!!!!!
!! for output file
output=OUTPUT_FILE
call write_header_to(output,min_eigen,max_eigen)
!!!
if( read_alpha == 0 ) then
  write(output,*) "# alpha and beta for SC are set to 1.0"
else
  write(output,*) "# alpha and beta are set in", trim(ALPHA_BETA)
endif
write(output,*)
!!!
if( new_config == 0 ) then
  write(output,*) "# configs read from ", trim(Fconfigin)
  if( fix_seed == 0 ) then
    write(output,*) "# random seed is succeeded from the previous simulation"
  elseif( fix_seed == 1 ) then
    write(output,*) "# random seed is fixed to seed=",seed
  elseif( fix_seed == 2 ) then
    write(output,*) "# random seed is determined by the system time"
  endif
else
  write(output,*) "# new configs"
  if( fix_seed == 1 ) then
    write(output,*) "# random seed is fixed to seed=",seed
  else
    write(output,*) "# random seed is determined by the system time"
  endif
endif
!!!
call write_observable_list(output)

end subroutine write_header

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_header_to(output,min_eigen,max_eigen)
implicit none

complex(kind(0d0)), intent(in) :: min_eigen,max_eigen
integer, intent(in) :: output

write(output,'(a,a)') "# simplicial complex : ",trim(SC_FILE_NAME)
write(output,'(A,F6.4)') "# lattice spacing (\lambda=1)= ",LatticeSpacing
write(output,'(A,I5)') "# NMAT= ",NMAT
write(output,'(a)') "#"
write(output,'(A,I5)') "# Ntau= ",Ntau
write(output,'(A,F10.8)') "# Dtau= ",Dtau
write(output,'(A,I5)') "# iterations= ", num_ite
write(output,'(a)') "#"
write(output,'(A,F6.4)') "# factor of Dtau for Phi= ",R_phi
write(output,'(A,F6.4)') "# factor of Dtau for UMAT= ",R_A
write(output,'(A,I5)') "# m_omega= ",m_omega
write(output,'(A,E12.5)') "# mass_square_phi= ",mass_square_phi
write(output,'(A,E12.5)') "# phys_mass_square_phi= ",phys_mass_square_phi
write(output,'(A,E12.5)') "# mass_f= ",mass_f
write(output,'(A,E12.5)') "# Remez_factor4= ",Remez_factor4
write(output,'(A,E12.5)') "# Remez_factor8= ",Remez_factor8
write(output,'(A,E12.5)') "# minimal eigenvalue of DD^\dagger= ",dble(min_eigen*conjg(min_eigen))
write(output,'(A,E12.5)') "# maximal eigenvalue of DD^\dagger= ",dble(max_eigen*conjg(max_eigen))

end subroutine write_header_to


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_observable_list(output)
implicit none

integer, intent(in) :: output
integer :: i

!obs_name(1)="Sb"

write(output,'(a)',advance='no') "# 1) iteration, 2) delta Hamiltonian, 3) max||Uf-1||/max 4) CG ite, "
do i=1,num_obs
  write(output,'(I3,a1,a,a1)',advance='no') i+3,")",trim(obs_name(i) ),","
enddo
write(output,'(I3,a)') num_obs+4,") acceptance rate"

end subroutine write_observable_list

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_observables(PhiMat,UMAT,ite,accept,delta_Ham,total_ite,CGite,ratio)
use SUN_generators, only : trace_MTa
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: Phi(1:dimG,1:num_sites)
double precision, intent(in) :: delta_Ham, ratio
integer, intent(in) :: ite, accept, total_ite, CGite
integer :: output,i,s,a

call calc_bosonic_action(OBS(1),UMAT,PhiMat)
call calc_TrX2(OBS(2),PhiMat)

!! for standard output
call write_observables_to(6,ite,total_ite,accept,delta_Ham,CGite,ratio)
!! for output file
call write_observables_to(OUTPUT_FILE,ite,total_ite,accept,delta_Ham,CGite,ratio)

end subroutine write_observables

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_observables_to(output,ite,total_ite,accept,delta_Ham,CGite,ratio)
implicit none

integer, intent(in) :: output,ite,accept,total_ite
double precision, intent(in) :: delta_Ham, ratio
integer, intent(in) :: CGite
integer i

write(output,'(I6,2X)',advance='no') ite
write(output,'(f12.5,2X)',advance='no') delta_Ham
write(output,'(f6.2,2X)',advance='no') ratio
write(output,'(I6,2X)',advance='no') CGite
do i=1,num_obs
  write(output,'(E12.5,2X)',advance='no') OBS(i)
enddo
write(output,"(F6.4)") dble(accept)/dble(ite-total_ite)

end subroutine write_observables_to

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_basic_info_to_medfile(MED_CONF_FILE)
implicit none

integer, intent(in) :: MED_CONF_FILE

write(MED_CONF_FILE) NMAT
write(MED_CONF_FILE) LatticeSpacing
write(MED_CONF_FILE) SC_FILE_NAME
write(MED_CONF_FILE) ALPHA_BETA
write(MED_CONF_FILE) test_mode
write(MED_CONF_FILE) new_config
write(MED_CONF_FILE) fix_seed
write(MED_CONF_FILE) read_alpha
write(MED_CONF_FILE) save_med_step
write(MED_CONF_FILE) obs_step
write(MED_CONF_FILE) m_omega
write(MED_CONF_FILE) mass_square_phi
write(MED_CONF_FILE) mass_f
write(MED_CONF_FILE) Remez_factor4
write(MED_CONF_FILE) Remez_factor8
write(MED_CONF_FILE) epsilon
write(MED_CONF_FILE) CG_max
write(MED_CONF_FILE) num_ite
write(MED_CONF_FILE) Ntau
write(MED_CONF_FILE) Dtau
write(MED_CONF_FILE) R_phi
write(MED_CONF_FILE) R_A
write(MED_CONF_FILE) Fconfigin
write(MED_CONF_FILE) Foutput
write(MED_CONF_FILE) Fconfigout
write(MED_CONF_FILE) Fmedconf

end subroutine write_basic_info_to_medfile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_smallset_and_largest_eigenvalues_of_D(min_eigen,max_eigen,UMAT,PhiMat)
implicit none

complex(kind(0d0)), intent(out) :: min_eigen, max_eigen
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)

complex(kind(0d0)) :: eigenvalues(1:sizeD)

call calc_eigenvalues_Dirac(eigenvalues,UMAT,PhiMat)

min_eigen=eigenvalues(1)
max_eigen=eigenvalues(sizeD)

end subroutine calc_smallset_and_largest_eigenvalues_of_D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End OUTPUT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Begin Dirac operators
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! make Dirac matrix
subroutine make_Dirac(Dirac,UMAT,PhiMat)
use SUN_generators, only : trace_MTa
implicit none

complex(kind(0d0)), intent(inout) :: Dirac(1:sizeD,1:sizeD)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
!complex(kind(0d0)) :: Phi(1:dimG,1:num_sites)

complex(kind(0d0)) :: unitvec(1:sizeD),dirac_vec(1:sizeD)
integer i,j
integer s,a

!do s=1,num_sites
  !do a=1,dimG
    !call trace_MTa(Phi(a,s),PhiMat(:,:,s),a,NMAT)
  !enddo
!enddo
do i=1,sizeD
  unitvec=(0d0,0d0)
  unitvec(i)=(1d0,0d0)
  call prod_Dirac(Dirac(:,i),unitvec,sizeD,UMAT,PhiMat)
enddo

end subroutine make_Dirac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! make D^\dagger D matrix
subroutine make_DdagD(DdagD,UMAT,PhiMat)
implicit none

complex(kind(0d0)), intent(inout) :: DdagD(1:sizeD,1:sizeD)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)


complex(kind(0d0)) :: unitvec(1:sizeD),dirac_vec(1:sizeD)
integer i,j

do i=1,sizeD
  unitvec=(0d0,0d0)
  unitvec(i)=(1d0,0d0)
  call prod_DdagD(DdagD(:,i),unitvec,sizeD,UMAT,Phimat)
enddo

end subroutine make_DdagD


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compute D^dagger.D.vec
!!   1) D_{ij}=-D_{ji}
!!   2) D^\dag)_{ij} = (D^*)_{ji} = - D^*_{ij}
!! ==>
!!   [D^\dagger D v]_i 
!!     = conjg(D_{ji}) D_{jk} v_k 
!!     = -conjg(D_{ij}) D_{jk} v_k                 
!! w=Dv => D^\dagger w = -( D.conjg(w) )^\dagger 
!!                     = -conjg[ D.conjg( D v ) ]
subroutine Prod_DdagD(DdagD_vec, vec, Vsize,UMAT,Phimat)
implicit none

integer, intent(in) :: Vsize !! vecsize must be sizeD
complex(kind(0d0)), intent(in) :: vec(1:Vsize)
complex(kind(0d0)), intent(inout) :: DdagD_vec(1:Vsize)
complex(kind(0d0)), intent(in) :: UMAT(:,:,:),PhiMat(:,:,:)

complex(kind(0d0)) :: tmpvec(1:Vsize),tmpvec2(1:Vsize)

integer :: i



call Prod_Dirac(tmpvec,vec,Vsize,UMAT,PhiMat)
do i=1,Vsize
  tmpvec2(i)=-dconjg(tmpvec(i))
enddo

call Prod_Dirac(tmpvec,tmpvec2,Vsize,UMAT,PhiMat)

do i=1,Vsize
  DdagD_vec(i)=dconjg(tmpvec(i))!*0.25d0 !! HERE!!
enddo


end subroutine Prod_DdagD



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compute D.vec
!!  S_f = 1/2 \Psi^T D \Psi
subroutine Prod_Dirac(D_vec, vec, Vsize,UMAT,PhiMat)
!use matrix_functions, only : matrix_3_product
use matrix_functions
implicit none

integer, intent(in) :: Vsize !! vecsize must be sizeD
complex(kind(0d0)), intent(in) :: vec(1:Vsize)
complex(kind(0d0)), intent(inout) :: D_vec(1:Vsize)
complex(kind(0d0)), intent(in) :: UMAT(:,:,:)
complex(kind(0d0)), intent(in) :: PhiMat(:,:,:)

!complex(kind(0d0)), parameter :: mass_f=1d0

complex(kind(0d0)) :: eta_mat(1:NMAT,1:NMAT,1:num_sites)
!complex(kind(0d0)) :: eta_ele(1:dimG,1:num_sites)
complex(kind(0d0)) :: lambda_mat(1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)) :: lambda_ele(1:dimG,1:num_links)
complex(kind(0d0)) :: chi_mat(1:NMAT,1:NMAT,1:num_faces)
!complex(kind(0d0)) :: chi_ele(1:dimG,1:num_faces)

complex(kind(0d0)) :: DF_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: DF_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: DF_chi(1:NMAT,1:NMAT,1:num_faces)



!complex(kind(0d0)) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: bPhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: r_site(1:dimG,1:num_sites)
complex(kind(0d0)) :: r_link(1:dimG,1:num_links)
complex(kind(0d0)) :: r_face(1:dimG,1:num_faces)


complex(kind(0d0)) :: Uinv(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT,1:num_faces),Ufm(1:NMAT,1:NMAT,1:num_faces)
!complex(kind(0d0)) :: diff_Omega(1:NMAT,1:NMAT,1:dimG)

complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat3(1:NMAT,1:NMAT)
complex(kind(0d0)) :: comm(1:NMAT,1:NMAT)
complex(kind(0d0)) :: acomm(1:NMAT,1:NMAT)
complex(kind(0d0)) :: line1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: line2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: trace,tmp

complex(kind(0d0)) :: AMAT(1:NMAT,1:NMAT)
complex(kind(0d0)) :: BMAT(1:NMAT,1:NMAT)
complex(kind(0d0)) :: sinU(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: cosUinv(1:NMAT,1:NMAT,1:num_faces)


integer :: s,l,f
integer :: i,j,k
integer :: a,b,c
integer :: r
integer :: label

!! for test
!complex(kind(0d0)) :: tmp_diff_Omega(1:NMAT,1:NMAT,1:dimG)
!complex(kind(0d0)) :: tmp_diff_Omega2(1:NMAT,1:NMAT,1:dimG)
integer :: ii,jj

!! preparation
D_vec=(0d0,0d0)
!r_site=(0d0,0d0)
!r_link=(0d0,0d0)
!r_face=(0d0,0d0)

call vec_to_mat(eta_mat,lambda_mat,chi_mat,vec)

do s=1,num_sites
  do i=1,NMAT
    do j=1,NMAT
  !call make_traceless_matrix_from_modes(PhiMat(:,:,s),NMAT,Phi(:,s))
  bPhiMat(i,j,s)=conjg(PhiMat(j,i,s))
  !call make_traceless_matrix_from_modes(bPhiMat(:,:,s),NMAT,dconjg(Phi(:,s)))
    enddo
  enddo
enddo

DF_eta=(0d0,0d0)
DF_lambda=(0d0,0d0)
DF_chi=(0d0,0d0)

if ( p1 == 0 ) then 
    !write(*,*) p1
!! (1) site action 
do s=1,num_sites
  call matrix_commutator(tmpmat1,PhiMat(:,:,s),eta_mat(:,:,s))
  DF_eta(:,:,s)=DF_eta(:,:,s) &
      +dcmplx(alpha_s(s))*(-0.5d0,0d0) *overall_factor*tmpmat1
enddo

endif



if (p2==0) then
    !write(*,*) p2
!! (2) link action 1
do s=1,num_sites
  tmpmat2=(0d0,0d0)
  do k=1,linkorg_to_s(s)%num_
    ! tmpmat2 = -i \sum_l( \alpha_l U_l^{-1}.\lambda_l.U_l )
    l=linkorg_to_s(s)%labels_(k)
    call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
      lambda_mat(:,:,l), NMAT, &
      UMAT(:,:,l), NMAT, &
      (0d0,0d0), tmpmat1, NMAT)
    call ZGEMM('C','N',NMAT,NMAT,NMAT,(0d0,-1d0)*dcmplx(alpha_l(l)), &
      UMAT(:,:,l), NMAT, &
      tmpmat1, NMAT, &
      (1d0,0d0), tmpmat2, NMAT)  
  enddo
  DF_eta(:,:,s)=DF_eta(:,:,s) &
    + overall_factor * tmpmat2
  !!!
  do k=1,linktip_from_s(s)%num_
    l=linktip_from_s(s)%labels_(k)
    DF_eta(:,:,s)=DF_eta(:,:,s) &
      + (0d0,1d0) * dcmplx(alpha_l(l)) * overall_factor * lambda_mat(:,:,l)
  enddo
enddo

do l=1,num_links
  !tmpmat2=i alpha_l U_l \lambda_l U_l^{-1}
  s=link_tip(l)
  call ZGEMM('N','C',NMAT,NMAT,NMAT,(1d0,0d0), &
    eta_mat(:,:,s), NMAT, &
    UMAT(:,:,l), NMAT, &
    (0d0,0d0), tmpmat1, NMAT)
  call ZGEMM('N','N',NMAT,NMAT,NMAT,(0d0,1d0)*dcmplx(alpha_l(l)), &
    UMAT(:,:,l), NMAT, &
    tmpmat1, NMAT, &
    (0d0,0d0), tmpmat2, NMAT)
  DF_lambda(:,:,l) = DF_lambda(:,:,l) &
    + overall_factor * tmpmat2
  s=link_org(l)
  DF_lambda(:,:,l) = DF_lambda(:,:,l) &
      - (0d0,1d0) * dcmplx(alpha_l(l)) * overall_factor * eta_mat(:,:,s)
enddo
endif

    
if(p3==0) then
!! (3) link action 2 
do l=1,num_links
  ! compute Ul.bPhi_tip(l).Ul^\dagger + bPhi_org(l)
  tmpmat2=bPhiMat(:,:,link_org(l))
  s=link_tip(l) 
  call ZGEMM('N','C',NMAT,NMAT,NMAT,(1d0,0d0), &
    bPhiMat(:,:,s), NMAT, &
    UMAT(:,:,l), NMAT, &
    (0d0,0d0), tmpmat1, NMAT)
  call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
    UMAT(:,:,l), NMAT, &
    tmpmat1, NMAT, &
    (1d0,0d0), tmpmat2, NMAT)
  !!!
  ! compute the commutator 
  call Matrix_Commutator(comm,tmpmat2,lambda_mat(:,:,l))
  !!
  DF_lambda(:,:,l)=DF_lambda(:,:,l) &
      + dcmplx(alpha_l(l)) * overall_factor * comm
enddo
endif

if(p4==0) then
!! (4) face action 1
do f=1,num_faces
  s=sites_in_f(f)%label_(1)
  call matrix_commutator(comm,PhiMat(:,:,s),chi_mat(:,:,f))
  DF_chi(:,:,f)=DF_chi(:,:,f) &
      + (-2d0,0d0) * dcmplx(alpha_f(f)) * overall_factor * comm
enddo
endif

if( p5==0 ) then
!! (5) face action 2
do f=1,num_faces
    call Make_face_variable(Uf(:,:,f),f,UMAT) 
  if( m_omega .ne. 0) then 
    call matrix_power(Ufm(:,:,f),Uf(:,:,f),m_omega)
    call calc_sinU_and_cosUinv(sinU(:,:,f),cosUinv(:,:,f),Ufm(:,:,f))
  endif
enddo


do f=1,num_faces
  do i=1,links_in_f(f)%num_
    l=links_in_f(f)%link_labels_(i)

    if( m_omega == 0 ) then
      call calc_Amat(Amat,f,i,1,Uf(:,:,f),UMAT)
      call calc_Bmat(Bmat,f,i,1,Uf(:,:,f),UMAT)

      call matrix_3_product(tmpmat1,Amat,lambda_mat(:,:,l),Bmat)
      call matrix_3_product(tmpmat2,Bmat,lambda_mat(:,:,l),Amat,'C','N','C')

      DF_chi(:,:,f)=DF_chi(:,:,f)&
        + cmplx(dble(links_in_f(f)%link_dirs_(i)))*(0d0,1d0)&
          * cmplx(overall_factor) * cmplx(alpha_f(f)*beta_f(f)) & 
          * (tmpmat1 + tmpmat2)
    else
      line1=(0d0,0d0)
      line2=(0d0,0d0)
      do k=1,m_omega
        call calc_Amat(Amat,f,i,k,Uf(:,:,f),UMAT)
        call calc_Bmat(Bmat,f,i,k,Uf(:,:,f),UMAT)
  
        ! tmpmat2=A.lambda.B
        call matrix_3_product(tmpmat2,Amat,lambda_mat(:,:,l),Bmat)
        ! tmpmat3=B^dag.lambda.A^dag 
        call matrix_3_product(tmpmat3,Bmat,lambda_mat(:,:,l),Amat,'C','N','C')
  
        ! line1 = A.lambda.B+Bdag.lambda.Adag
        line1=line1+tmpmat2+tmpmat3
        ! line2 = cosUinv.(A.lambda.B-Bdag.lambda.Adag).cosUinv
        call matrix_3_product(tmpmat1,cosUinv(:,:,f),tmpmat2-tmpmat3,cosUinv(:,:,f))
        line2=line2+tmpmat1
      enddo
  
      ! line1 = {A.lambda.B+B^dag.lambda.A^dag , (Um+Uminv)^{-1}
      tmpmat1=line1
      call matrix_AntiCommutator(line1,tmpmat1,cosUinv(:,:,f))
      ! line2 = {cosUinv.(A.lambda.B-Bdag.lambda.Adag).cosUinv, sinU}
      tmpmat1=line2
      call matrix_AntiCommutator(line2,tmpmat1,sinU(:,:,f))
      
      DF_chi(:,:,f)=DF_chi(:,:,f)&
        + cmplx(dble(links_in_f(f)%link_dirs_(i)))*(0d0,1d0)&
          * cmplx(overall_factor) / cmplx(dble(m_omega))&
          * cmplx(alpha_f(f)*beta_f(f)) * (line1-line2)
    endif
  enddo
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! make traceless for (i,j)
  trace=(0d0,0d0)
  do j=1,NMAT
    trace=trace+DF_chi(j,j,f)
  enddo
  do j=1,NMAT
    DF_chi(j,j,f)=DF_chi(j,j,f)-trace/cmplx(dble(NMAT))
  enddo
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
enddo

do l=1,num_links
  do i=1,face_in_l(l)%num_       
    f=face_in_l(l)%label_(i)
    ! j: position of the link l in the face l
    do j=1,links_in_f(f)%num_
      if ( l == links_in_f(f)%link_labels_(j) ) exit
    enddo

    if( m_omega == 0 ) then 
      call calc_Amat(Amat,f,j,1,Uf(:,:,f),UMAT)
      call calc_Bmat(Bmat,f,j,1,Uf(:,:,f),UMAT)

      call matrix_3_product(tmpmat1,Bmat,chi_mat(:,:,f),Amat)
      call matrix_3_product(tmpmat2,Amat,chi_mat(:,:,f),Bmat,'C','N','C')

      DF_lambda(:,:,l)=DF_lambda(:,:,l) &
        + cmplx(dble(links_in_f(f)%link_dirs_(j)))*(0d0,-1d0)&
          * cmplx(overall_factor) * cmplx(alpha_f(f)*beta_f(f)) &
          * (tmpmat1+tmpmat2)
    else
      ! tmpmat1= { (Uf+Ufinv)^{-1} , \chi_f }
      call matrix_anticommutator(tmpmat1,cosUinv(:,:,f),chi_mat(:,:,f))
  
      ! tmpmat2= (Uf+Ufinv)^{-1}.{ Uf-Ufinv , \chi_f }.(Uf+Ufinv)^{-1} 
      call matrix_anticommutator(tmpmat3,sinU(:,:,f),chi_mat(:,:,f))
      call matrix_3_product(tmpmat2,cosUinv(:,:,f),tmpmat3,cosUinv(:,:,f))
  
      line1=tmpmat1-tmpmat2
      line2=tmpmat1+tmpmat2
  
      acomm=(0d0,0d0)
      do k=1,m_omega
        call calc_Amat(Amat,f,j,k,Uf(:,:,f),UMAT)
        call calc_Bmat(Bmat,f,j,k,Uf(:,:,f),UMAT)
        ! tmpmat1= B.(tmpmat1-tmpmat2).A
        call matrix_3_product(tmpmat1,BMAT,line1,AMAT)
        ! tmpmat2 = Adag.(tmpmat1+tmpmat2).Bdag
        call matrix_3_product(tmpmat2, AMAT,line2,BMAT,'C','N','C')
  
        acomm=acomm+tmpmat1+tmpmat2
      enddo
  
      DF_lambda(:,:,l)=DF_lambda(:,:,l) &
        + cmplx(dble(links_in_f(f)%link_dirs_(j)))*(0d0,1d0)&
          * cmplx(-overall_factor) / cmplx(dble(m_omega))&
          * cmplx(alpha_f(f)*beta_f(f)) * acomm
    endif
  enddo
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! make traceless for (i,j)
   trace=(0d0,0d0)
   do j=1,NMAT
     trace=trace+DF_lambda(j,j,l)
   enddo
   do j=1,NMAT
     DF_lambda(j,j,l)=DF_lambda(j,j,l)-trace/cmplx(dble(NMAT))
   enddo
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
enddo
endif

call mat_to_vec(D_vec,DF_eta,DF_lambda,DF_chi)

do i=1,sizeD,2
  D_vec(i)=D_vec(i) + dcmplx(mass_f)*vec(i+1)
  D_vec(i+1)=D_vec(i+1) - dcmplx(mass_f)*vec(i)
enddo
!D_vec=D_vec*overall_factor
!write(*,*) D_vec
end subroutine Prod_Dirac


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! check anti-symmetricity of D
subroutine check_Dirac(UMAT,PhiMat)
implicit none 

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)

complex(kind(0d0)) :: PF(1:sizeD)
complex(kind(0d0)) :: Dirac(1:sizeD,1:sizeD), Det
double precision :: rtmp,dtmp
integer :: seed
integer :: i,j
integer :: s,l,a,b


!seed=12345
!call genrand_init( put=seed )
! produce pseudo-fermion
!call make_pseudo_fermion(PF,UMAT,Phi)
!! calculate Hamiltonian 

call make_Dirac(Dirac,UMAT,PhiMat)

rtmp=0d0
dtmp=0d0
do i=1,sizeD-1
  rtmp=rtmp+dble(Dirac(i,i)*dconjg(Dirac(i,i)))
  do j=i+1,sizeD
    rtmp=rtmp+dble( (Dirac(i,j)+Dirac(j,i)) * dconjg(Dirac(i,j)+Dirac(j,i)) ) 
  enddo
enddo
write(*,*) "# D's diagonal elements?", dtmp+dble(Dirac(sizeD,sizeD)*dconjg(Dirac(sizeD,sizeD)))
write(*,*) "# Is D anti-symmetric?", rtmp

call make_DdagD(Dirac,UMAT,PhiMat)
rtmp=0d0
do i=1,sizeD
  do j=i,sizeD
    rtmp=rtmp&
        +dble ( &
          ( Dirac(i,j) - dconjg( Dirac(j,i) ) ) &
          *dconjg( Dirac(i,j) - dconjg( Dirac(j,i) ) ) )
  enddo
enddo
write(*,*) "# Is D^\dag D hermitian?", rtmp

end subroutine check_Dirac



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to return product of d/d\Phi D
subroutine prod_dDdPhi(dDdPhi_eta,dDdPhi_chi,eta_mat,chi_mat,UMAT)
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links) 
complex(kind(0d0)), intent(in) :: eta_mat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(in) :: chi_mat(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)), intent(out) :: dDdPhi_eta(1:NMAT,1:NMAT,1:num_sites,1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(out) :: dDdPhi_chi(1:NMAT,1:NMAT,1:num_faces,1:NMAT,1:NMAT,1:num_sites)


complex(kind(0d0)) :: tmp
integer :: s,f,i,j,ii,jj
integer :: a,b,c,r

!! site: (s,a) <--> r=a+dimG*(s-1)
!! link: (l,a) <--> r=a+dimG*(num_sites + l -1)
!! face: (f,a) <--> r=a+dimG*(num_sites + num_links + f -1 )
!! preparation
!dDdPhi_vec=(0d0,0d0)
dDdPhi_eta=(0d0,0d0)
dDdPhi_chi=(0d0,0d0)

if( p1 == 0 ) then
    !write(*,*) p1
!! (1) Dirac from site
do s=1,num_sites
  do jj=1,NMAT
    i=jj
    do ii=1,NMAT
      do j=1,NMAT 
        dDdPhi_eta(i,j,s,ii,jj,s)= dDdPhi_eta(i,j,s,ii,jj,s) &
          + (-0.5d0,0d0)*cmplx( alpha_s(s)*overall_factor ) &
            * eta_mat(ii,j,s)
      enddo
    enddo
    !!
    do ii=1,NMAT
      j=ii
      do i=1,NMAT
        dDdPhi_eta(i,j,s,ii,jj,s)= dDdPhi_eta(i,j,s,ii,jj,s) &
          - (-0.5d0,0d0)*cmplx( alpha_s(s)*overall_factor ) &
            * eta_mat(i,jj,s)
      enddo
    enddo
  enddo
enddo
endif

!! (2) Dirac from link 1
! there is no contribution
    
!! (3) Dirac from link 2
! there is no contribution

if ( p4 == 0 ) then
!! (4) Dirac from face 1
do f=1,num_faces
  s=sites_in_f(f)%label_(1)
  do jj=1,NMAT
    i=jj
    do ii=1,NMAT
      do j=1,NMAT 
        dDdPhi_chi(i,j,f,ii,jj,s)= dDdPhi_chi(i,j,f,ii,jj,s) &
          + (-2d0,0d0)*cmplx( alpha_f(f)*overall_factor ) &
            * chi_mat(ii,j,f)
      enddo
    enddo
    !!
    do ii=1,NMAT
      j=ii
      do i=1,NMAT
        dDdPhi_chi(i,j,f,ii,jj,s)= dDdPhi_chi(i,j,f,ii,jj,s) &
          - (-2d0,0d0)*cmplx( alpha_f(f)*overall_factor ) &
            * chi_mat(i,jj,f)
      enddo
    enddo
  enddo
enddo
endif


end subroutine prod_dDdPhi


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to return product of d/d\Phi D^\dagger
subroutine prod_dDdbPhi(dDdbPhi_lambda,lambda_mat,UMAT)
use matrix_functions
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links) 
complex(kind(0d0)), intent(in) :: lambda_mat(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(out) :: dDdbPhi_lambda(1:NMAT,1:NMAT,1:num_links,1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: dDdbPhi_vec(1:sizeD,1:dimG,1:num_sites)


complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmp

integer :: s,l
integer :: i,j,ii,jj
integer :: a,b

!! preparation
dDdbPhi_vec=(0d0,0d0)
dDdbPhi_lambda=(0d0,0d0)


! (1) Dirac from site
! no contribution

!! (2) Dirac from link 1
! no contribution

if ( p3 == 0 ) then
!! (3) Dirac from link 2
do l=1,num_links
  call matrix_product(tmpmat1,Umat(:,:,l),lambda_mat(:,:,l),'C','N')
  call matrix_product(tmpmat2,lambda_mat(:,:,l),Umat(:,:,l))
  do jj=1,NMAT
    do ii=1,NMAT
      do j=1,NMAT
        do i=1,NMAT
          dDdbPhi_lambda(i,j,l,ii,jj,link_tip(l)) &
           = dDdbPhi_lambda(i,j,l,ii,jj,link_tip(l)) &
             + cmplx( alpha_l(l) * overall_factor ) &
               * ( UMAT(i,jj,l) * tmpmat1(ii,j) &
                   - tmpmat2(i,jj) * conjg( UMAT(j,ii,l) ))
          if ( i == jj ) then
            dDdbPhi_lambda(i,j,l,ii,jj,link_org(l)) &
             = dDdbPhi_lambda(i,j,l,ii,jj,link_org(l)) &
               + cmplx( alpha_l(l) * overall_factor ) &
                 * lambda_mat(ii,j,l)
          endif
          if ( ii == j ) then
            dDdbPhi_lambda(i,j,l,ii,jj,link_org(l)) &
             = dDdbPhi_lambda(i,j,l,ii,jj,link_org(l)) &
               - cmplx( alpha_l(l) * overall_factor ) &
                 * lambda_mat(i,jj,l)
          endif
        enddo
      enddo
    enddo
  enddo
enddo

endif

!! (4) Dirac from face 1
!   no contribution

!dDdbPhi_chi=dDdbPhi_chi*overall_factor
end subroutine prod_dDdbPhi


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to return product of d/dA D . vec
subroutine prod_dDdA(dDdA_eta,dDdA_lambda,dDdA_chi,&
    eta_mat,lambda_mat,chi_mat,UMAT,PhiMat)
use matrix_functions
implicit none

! d/dA_{ll,ii,jj} (D\Psi)_{s,i,j}=dDdA_eta(i,j,s,ii,jj,ll)
complex(kind(0d0)), intent(inout) :: dDdA_eta(1:NMAT,1:NMAT,1:num_sites,1:NMAT,1:NMAT,1:num_links)
! d/dA_{ll,ii,jj} (D\Psi)_{l,i,j}=dDdA_lambda(i,j,l,ii,jj,ll)
complex(kind(0d0)), intent(inout) :: dDdA_lambda(1:NMAT,1:NMAT,1:num_links,1:NMAT,1:NMAT,1:num_links)
! d/dA_{ll,ii,jj} (D\Psi)_{f,i,j}=dDdA_chi(i,j,f,ii,jj,ll)
complex(kind(0d0)),intent(inout) :: dDdA_chi(1:NMAT,1:NMAT,1:num_faces,1:NMAT,1:NMAT,1:num_links)

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links) 
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(in) :: eta_mat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(in) :: lambda_mat(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: chi_mat(1:NMAT,1:NMAT,1:num_faces)

complex(kind(0d0)) :: trace,tmp
complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT),tmpmat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat3(1:NMAT,1:NMAT),tmpmat4(1:NMAT,1:NMAT)


integer :: s,l,f,ll
integer :: i,j,k,nl,r,ii,jj
integer :: a,b,c,d,e

!type diffdiff_by_linkvals_in_face
  !complex(kind(0d0)), allocatable :: val(:,:,:,:,:,:)
!end type diffdiff_by_linkvals_in_face

!type(diffdiff_by_linkvals_in_face) :: diffdiff_Omega2(1:num_faces)

!do f=1,num_faces
  !allocate( diffdiff_Omega2(f)%val(1:NMAT,1:NMAT,1:dimG,1:dimG,&
    !1:links_in_f(f)%num_, 1:links_in_f(f)%num_) )
!enddo

!! preparation
!dDdA_vec=(0d0,0d0)
dDdA_eta=(0d0,0d0)
dDdA_lambda=(0d0,0d0)
dDdA_chi=(0d0,0d0)

! for test
!call make_SUN_generators(T,NMAT)

!! (1) Dirac from site
!   no contribution

if ( p2==0 ) then
!! (2) Dirac from link 1
do s=1,num_sites
  do k=1,linkorg_to_s(s)%num_
    l=linkorg_to_s(s)%labels_(k)
    ! tmpmat1 = lambda_l.U_l
    call matrix_product(tmpmat1,lambda_mat(:,:,l),UMAT(:,:,l))
    ! tmpmat2 = Ul^{-1}.lambda_l
    call matrix_product(tmpmat2,UMAT(:,:,l),lambda_mat(:,:,l),'C','N')
    do jj=1,NMAT
      do ii=1,NMAT
        do j=1,NMAT
          do i=1,NMAT
            dDdA_eta(i,j,s,ii,jj,l)= dDdA_eta(i,j,s,ii,jj,l) &
              - cmplx(alpha_l(l))*conjg(UMAT(jj,i,l))*tmpmat1(ii,j) &
              + cmplx(alpha_l(l))*tmpmat2(i,jj)*UMAT(ii,j,l)
          enddo
        enddo
      enddo
    enddo
  enddo
enddo    

do l=1,num_links
  s=link_tip(l)
  ! tmpmat2 = Ul.eta_l.Ul^{-1}
  call matrix_product(tmpmat1,UMAT(:,:,l),eta_mat(:,:,s))
  call matrix_product(tmpmat2,tmpmat1,UMAT(:,:,l),'N','C')
  do ii=1,NMAT
    do i=1,NMAT
      jj=ii
      j=i
      do k=1,NMAT
        dDdA_lambda(k,j,l,ii,k,l) = dDdA_lambda(k,j,l,ii,k,l) & 
          - cmplx(alpha_l(l)) * tmpmat2(ii,j)
        dDdA_lambda(i,k,l,k,jj,l) = dDdA_lambda(i,k,l,k,jj,l) & 
          + cmplx(alpha_l(l)) * tmpmat2(i,jj)
      enddo
    enddo
  enddo
enddo
endif 

if ( p3==0 ) then
!! (3) Dirac from link 2
do l=1,num_links
  s=link_tip(l)

  ! tmpmat1 = Ul.\bar\Phi.Ul^{-1}
  ! tmpmat2 = Ul.\bar\Phi.Ul^{-1}.lambda_l
  ! tmpmat3 = lambda_l.Ul.\bar\Phi.Ul^{-1}
  call matrix_product(tmpmat2,UMAT(:,:,l),PhiMat(:,:,s),'N','C')
  call matrix_product(tmpmat1,tmpmat2,UMAT(:,:,l),'N','C')
  call matrix_product(tmpmat2,tmpmat1,lambda_mat(:,:,l))
  call matrix_product(tmpmat3,lambda_mat(:,:,l),tmpmat1)
  do jj=1,NMAT
    do ii=1,NMAT
      do j=1,NMAT
        do i=1,NMAT
          if ( i==jj ) then
            dDdA_lambda(i,j,l,ii,jj,l) = dDdA_lambda(i,j,l,ii,jj,l) & 
              +(0d0,1d0)*cmplx(alpha_l(l))*tmpmat2(ii,j)
          endif
          if ( j==ii ) then
            dDdA_lambda(i,j,l,ii,jj,l) = dDdA_lambda(i,j,l,ii,jj,l) & 
              +(0d0,1d0)*cmplx(alpha_l(l))*tmpmat3(i,jj)
          endif
          dDdA_lambda(i,j,l,ii,jj,l) = dDdA_lambda(i,j,l,ii,jj,l) & 
            -(0d0,1d0)*cmplx(alpha_l(l))*lambda_mat(i,jj,l)*tmpmat1(ii,j) &
            -(0d0,1d0)*cmplx(alpha_l(l))*lambda_mat(ii,j,l)*tmpmat1(i,jj)
        enddo
      enddo
    enddo
  enddo
enddo
endif

!! (4) Dirac from face 1
!   no contribution

if( p5 == 0 ) then 
!! (5) Dirac from face 2
! preparation
  call calc_fermion_force_from_omega&
      (dDdA_lambda,dDdA_chi,UMAT,lambda_mat,chi_mat)
endif 
dDdA_eta=dDdA_eta*cmplx(overall_factor)
dDdA_lambda=dDdA_lambda*cmplx(overall_factor)
dDdA_chi=dDdA_chi*cmplx(overall_factor)

end subroutine prod_dDdA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_fermion_force_from_omega&
    (dDdA_lambda,dDdA_chi,UMAT,lambda_mat,chi_mat)
!use Dirac_operator, only : calc_sinu_and_CosUinv, calc_Amat, Calc_Bmat
use matrix_functions
implicit none 
complex(kind(0d0)), intent(inout) :: dDdA_lambda&
  (1:NMAT,1:NMAT,1:num_links, 1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(inout) :: dDdA_chi&
  (1:NMAT,1:NMAT,1:num_faces, 1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: lambda_mat(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: chi_mat(1:NMAT,1:NMAT,1:num_faces)

complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: Ufm(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: sinU(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: cosUinv(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: Amat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Bmat(1:NMAT,1:NMAT)

complex(kind(0d0)) :: dCosUinvdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
complex(kind(0d0)) :: dSinUdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
complex(kind(0d0)) :: dAmatdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
complex(kind(0d0)) :: dBmatdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)

complex(kind(0d0)) :: prodmat1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: prodmat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: globalmat1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: globalmat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: globalmat3(1:NMAT,1:NMAT)
complex(kind(0d0)) :: globalmat4(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat3(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat4(1:NMAT,1:NMAT)
complex(kind(0d0)) :: line(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
integer :: l_label,ll_label,k,l,ll
integer :: i,j,ii,jj,f
complex(kind(0d0)) :: trace

! label of face
do f=1,num_faces
  call Make_face_variable(Uf(:,:,f),f,UMAT) 
  if (m_omega .ne. 0) then 
    call matrix_power(Ufm(:,:,f),Uf(:,:,f),m_omega)
    call calc_sinU_and_cosUinv(sinU(:,:,f),cosUinv(:,:,f),Ufm(:,:,f))
  endif
enddo

! label of face
do f=1,num_faces
! label of link fermion \lambda_l
do l_label=1,links_in_f(f)%num_
  l=links_in_f(f)%link_labels_(l_label)

  ! label of link to differentiate
  do ll_label=1,links_in_f(f)%num_
    ll=links_in_f(f)%link_labels_(ll_label)
!subroutine calc_dCosUinvdA_dSinUdA(dCosUinvdA,dSinUdA,&
!    cosUinv,sinU,Uf,UMAT,f,ll_label)

    line=(0d0,0d0)
    if( m_omega == 0 ) then
      call calc_dABmatdA(dAmatdA,dBmatdA,&
        Uf(:,:,f),UMAT,ll_label,f,l_label,1)
      call calc_Amat(Amat,f,l_label,1,Uf(:,:,f),UMAT)
      call calc_Bmat(Bmat,f,l_label,1,Uf(:,:,f),UMAT)

      
      do jj=1,NMAT
      do ii=1,NMAT
        call matrix_3_product(tmpmat1, &
          dAmatdA(:,:,ii,jj),lambda_mat(:,:,l),BMAT)
          !dAmatdA(:,:,ii,jj),lambda_mat(:,:,l),BMAT)
        call matrix_3_product(tmpmat2, &
          Amat,lambda_mat(:,:,l),dBmatdA(:,:,ii,jj))
        call matrix_3_product(tmpmat3, &
          dBmatdA(:,:,jj,ii),lambda_mat(:,:,l),Amat,'C','N','C')
        call matrix_3_product(tmpmat4, &
          Bmat,lambda_mat(:,:,l),dAmatdA(:,:,jj,ii),'C','N','C')
        line(:,:,ii,jj)=tmpmat1+tmpmat2+tmpmat3+tmpmat4
      enddo
      enddo
        !call matrix_product(tmpmat1,Amat,lambda_mat(:,:,l))
        !call zgemm('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
        !  tmpmat1, NMAT, &
        !  dBmatdA(:,:,ii,jj), NMAT, &
        !  (1d0,0d0), line(:,:,ii,jj), NMAT)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! make traceless for (ii,jj)
      tmpmat1=(0d0,0d0)
      do ii=1,NMAT
        tmpmat1=tmpmat1+line(:,:,ii,ii)
      enddo
      do ii=1,NMAT
        line(:,:,ii,ii)=line(:,:,ii,ii)-tmpmat1/cmplx(dble(NMAT))
      enddo
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do jj=1,NMAT
        do ii=1,NMAT
          dDdA_chi(:,:,f,ii,jj,ll)=dDdA_chi(:,:,f,ii,jj,ll) &
            + cmplx(dble(links_in_f(f)%link_dirs_(l_label) ))*(0d0,1d0) &
             * cmplx( alpha_f(f)*beta_f(f) ) &
             * line(:,:,ii,jj) 
        enddo
      enddo
    else
      call calc_dCosUinvdA_dSinUdA(dCosUinvdA,dSinUdA,&
        cosUinv(:,:,f),Uf(:,:,f),UMAT,f,ll_label)
      do k=1,m_omega
        call calc_dABmatdA(dAmatdA,dBmatdA,&
          Uf(:,:,f),UMAT,ll_label,f,l_label,k)
        call calc_Amat(Amat,f,l_label,k,Uf(:,:,f),UMAT)
        call calc_Bmat(Bmat,f,l_label,k,Uf(:,:,f),UMAT)
  
        ! globalmat1 = A lambda B
        call matrix_3_product(globalmat1,Amat,lambda_mat(:,:,l),Bmat)
  
        ! globalmat2 = B^dag lambda A^dag
        call matrix_3_product(globalmat2,Bmat,lambda_mat(:,:,l),Amat,'C','N','C')
  
        do jj=1,NMAT
        do ii=1,NMAT
          ! prodat1 = ¥delta_A lambda B + A lambda ¥delta_B
          call matrix_3_product(prodmat1, &
            dAmatdA(:,:,ii,jj),lambda_mat(:,:,l),BMAT)
          call matrix_product(tmpmat1,Amat,lambda_mat(:,:,l))
          call zgemm('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat1, NMAT, &
            dBmatdA(:,:,ii,jj), NMAT, &
            (1d0,0d0), prodmat1, NMAT)
  
          ! prodmat2 = ¥delta_B^dag lambda A^dag + B^dag lambda ¥delta_A^dag
          call matrix_3_product(prodmat2,&
            dBmatdA(:,:,jj,ii),lambda_mat(:,:,l),Amat,'C','N','C') 
          call matrix_product(tmpmat1,Bmat,lambda_mat(:,:,l),'C','N')
          call zgemm('N','C',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat1, NMAT, &
            dAmatdA(:,:,jj,ii), NMAT, & 
            (1d0,0d0), prodmat2, NMAT)
  
          ! line1
          call matrix_anticommutator(tmpmat1, prodmat1+prodmat2, CosUinv(:,:,f))
          line(:,:,ii,jj)=line(:,:,ii,jj)+tmpmat1
          ! line2
          call matrix_anticommutator(tmpmat1,globalmat1+globalmat2, dCosUinvdA(:,:,ii,jj))
          line(:,:,ii,jj) = line(:,:,ii,jj) + tmpmat1
          ! line3
          call matrix_product(tmpmat1,CosUinv(:,:,f),globalmat1-globalmat2)
          call matrix_product(tmpmat2,tmpmat1,CosUinv(:,:,f)) 
          call matrix_anticommutator(tmpmat3,dSinUdA(:,:,ii,jj),tmpmat2)
          line(:,:,ii,jj) = line(:,:,ii,jj) - tmpmat3
          ! line4
           ! 4-1 collect in tmpmat2
          call matrix_product(tmpmat3,dCosUinvdA(:,:,ii,jj),globalmat1-globalmat2)
          call matrix_product(tmpmat2,tmpmat3,CosUinv(:,:,f))
           ! 4-2 tmpmat1=CosUinv.(globalmat1-globalmat2) here
          call zgemm('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat1, NMAT, &
            dCosUinvdA(:,:,ii,jj), NMAT, &
            (1d0,0d0), tmpmat2, NMAT)
           ! 4-3
          call matrix_product(tmpmat3,CosUinv(:,:,f),prodmat1-prodmat2)
          call zgemm('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat3, NMAT, &
            CosUinv(:,:,f), NMAT, &
            (1d0,0d0), tmpmat2, NMAT)
           ! take anti-commutator
          call matrix_anticommutator(tmpmat1,tmpmat2,SinU(:,:,f))
          line(:,:,ii,jj) = line(:,:,ii,jj) - tmpmat1
  
        enddo ! ii
        enddo ! jj
      enddo ! k
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! make traceless for (ii,jj)
      tmpmat1=(0d0,0d0)
      do ii=1,NMAT
        tmpmat1=tmpmat1+line(:,:,ii,ii)
      enddo
      do ii=1,NMAT
        line(:,:,ii,ii)=line(:,:,ii,ii)-tmpmat1/cmplx(dble(NMAT))
      enddo
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do jj=1,NMAT
        do ii=1,NMAT
          dDdA_chi(:,:,f,ii,jj,ll)=dDdA_chi(:,:,f,ii,jj,ll) &
            + cmplx(dble(links_in_f(f)%link_dirs_(l_label) ))*(0d0,1d0) &
             * cmplx( alpha_f(f)*beta_f(f)/dble(m_omega) ) &
             * line(:,:,ii,jj)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! make traceless for (i,j)
          trace=(0d0,0d0)
          do j=1,NMAT
            trace=trace+dDdA_chi(j,j,f,ii,jj,ll)
          enddo
          do j=1,NMAT
            dDdA_chi(j,j,f,ii,jj,ll)=dDdA_chi(j,j,f,ii,jj,ll)-trace/cmplx(dble(NMAT))
          enddo
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        enddo
      enddo
    endif
  enddo ! ll
enddo ! l
enddo ! f

do l=1,num_links
do i=1,face_in_l(l)%num_       
  f=face_in_l(l)%label_(i)
  ! j: position of the link l in the face l
  do l_label=1,links_in_f(f)%num_
    if ( l == links_in_f(f)%link_labels_(l_label) ) exit
  enddo

  do ll_label=1,links_in_f(f)%num_
    ll=links_in_f(f)%link_labels_(ll_label)

    if (m_omega == 0) then 
      call calc_dABmatdA(dAmatdA,dBmatdA,&
        Uf(:,:,f),UMAT,ll_label,f,l_label,1)
      call calc_Amat(Amat,f,l_label,1,Uf(:,:,f),UMAT)
      call calc_Bmat(Bmat,f,l_label,1,Uf(:,:,f),UMAT)

      do jj=1,NMAT
        do ii=1,NMAT
          call matrix_3_product(tmpmat1,dBmatdA(:,:,ii,jj),chi_mat(:,:,f),Amat)
          call matrix_3_product(tmpmat2,Bmat,chi_mat(:,:,f),dAmatdA(:,:,ii,jj))
          call matrix_3_product(tmpmat3,dAmatdA(:,:,jj,ii),chi_mat(:,:,f),Bmat,'C','N','C')
          call matrix_3_product(tmpmat4,Amat,chi_mat(:,:,f),dBmatdA(:,:,jj,ii), 'C','N','C')
          line(:,:,ii,jj)=tmpmat1+tmpmat2+tmpmat3+tmpmat4
        enddo
      enddo
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! make traceless for (ii,jj)
      tmpmat1=(0d0,0d0)
      do ii=1,NMAT
        tmpmat1=tmpmat1+line(:,:,ii,ii)
      enddo
      do ii=1,NMAT
        line(:,:,ii,ii)=line(:,:,ii,ii)-tmpmat1/cmplx(dble(NMAT))
      enddo
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
      do jj=1,NMAT
        do ii=1,NMAT
           dDdA_lambda(:,:,l,ii,jj,ll)=dDdA_lambda(:,:,l,ii,jj,ll) &
              - cmplx(dble(links_in_f(f)%link_dirs_(l_label) ))*(0d0,1d0) &
               * cmplx( alpha_f(f)*beta_f(f) ) &
               * line(:,:,ii,jj)
        enddo
      enddo
    else
      call calc_dCosUinvdA_dSinUdA(dCosUinvdA,dSinUdA,&
        cosUinv(:,:,f),Uf(:,:,f),UMAT,f,ll_label)
  
      globalmat3=(0d0,0d0)
      globalmat4=(0d0,0d0)
      ! globalmat3 = { CosUinv, chi_f }
      call matrix_anticommutator(globalmat3,CosUinv(:,:,f),chi_mat(:,:,f))
      ! globalmat4 = CosUinv.{ SinU, chi_f }.CosUinv
      call matrix_anticommutator(tmpmat1,SinU(:,:,f),chi_mat(:,:,f))
      call matrix_3_product(globalmat4,CosUinv(:,:,f),tmpmat1,CosUinv(:,:,f))
  
      line=(0d0,0d0)
      do k=1,m_omega
        call calc_dABmatdA(dAmatdA,dBmatdA,&
          Uf(:,:,f),UMAT,ll_label,f,l_label,k)
        call calc_Amat(Amat,f,l_label,k,Uf(:,:,f),UMAT)
        call calc_Bmat(Bmat,f,l_label,k,Uf(:,:,f),UMAT)
  
        do jj=1,NMAT
        do ii=1,NMAT
          ! prodmat1
          call matrix_anticommutator(prodmat1,dCosUinvdA(:,:,ii,jj),chi_mat(:,:,f))
          ! prodmat2
          ! 1st term
          call matrix_anticommutator(tmpmat1,dSinUdA(:,:,ii,jj),chi_mat(:,:,f))
          call matrix_3_product(prodmat2,CosUinv(:,:,f),tmpmat1,CosUinv(:,:,f))
          ! 2nd term
          call matrix_anticommutator(tmpmat1,SinU(:,:,f),chi_mat(:,:,f))
          call matrix_product(tmpmat2,dCosUinvdA(:,:,ii,jj),tmpmat1)
          call zgemm('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat2, NMAT, &
            CosUinv(:,:,f),NMAT, &
            (1d0,0d0), prodmat2, NMAT)
          ! 3rd term tmpmat1={ SinU, chi_mat } now
          call matrix_product(tmpmat2,CosUinv(:,:,f),tmpmat1)
          call zgemm('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat2, NMAT, &
            dCosUinvdA(:,:,ii,jj),NMAT, &
            (1d0,0d0), prodmat2, NMAT)
  
          ! line1
          call matrix_3_product(tmpmat1,&
            dBmatdA(:,:,ii,jj),globalmat3-globalmat4,Amat)
          line(:,:,ii,jj)=line(:,:,ii,jj)+tmpmat1
          ! line2
          call matrix_product(tmpmat1,Bmat,globalmat3-globalmat4)
          call zgemm('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat1, NMAT, &
            dAmatdA(:,:,ii,jj),NMAT, &
            (1d0,0d0), line(:,:,ii,jj), NMAT)
          ! line3
          call matrix_product(tmpmat1,dAmatdA(:,:,jj,ii),globalmat3+globalmat4,'C','N') 
          call zgemm('N','C',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat1, NMAT, &
            Bmat,NMAT, &
            (1d0,0d0), line(:,:,ii,jj), NMAT)
          ! line4
          call matrix_product(tmpmat1,Amat,globalmat3+globalmat4,'C','N')
          call zgemm('N','C',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat1, NMAT, &
            dBmatdA(:,:,jj,ii),NMAT, & 
            (1d0,0d0), line(:,:,ii,jj), NMAT)
          ! line5
          call matrix_product(tmpmat1,Bmat,prodmat1-prodmat2)
          call zgemm('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat1, NMAT, &
            Amat,NMAT, &
            (1d0,0d0), line(:,:,ii,jj), NMAT)
          ! line6
          call matrix_product(tmpmat1,Amat,prodmat1+prodmat2,'C','N')
          call zgemm('N','C',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat1, NMAT, &
            Bmat,NMAT, &
            (1d0,0d0), line(:,:,ii,jj), NMAT)
        enddo ! ii
        enddo ! jj
      enddo ! k
  
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! make traceless for (ii,jj)
      tmpmat1=(0d0,0d0)
      do ii=1,NMAT
        tmpmat1=tmpmat1+line(:,:,ii,ii)
      enddo
      do ii=1,NMAT
        line(:,:,ii,ii)=line(:,:,ii,ii)-tmpmat1/cmplx(dble(NMAT))
      enddo
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
      do jj=1,NMAT
        do ii=1,NMAT
           dDdA_lambda(:,:,l,ii,jj,ll)=dDdA_lambda(:,:,l,ii,jj,ll) &
              - cmplx(dble(links_in_f(f)%link_dirs_(l_label) ))*(0d0,1d0) &
               * cmplx( alpha_f(f)*beta_f(f)/dble(m_omega) ) &
               * line(:,:,ii,jj)
        enddo
      enddo
    endif 
    do jj=1,NMAT
      do ii=1,NMAT
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! make traceless for (i,j)
        trace=(0d0,0d0)
        do j=1,NMAT
          trace=trace+dDdA_lambda(j,j,l,ii,jj,ll)
        enddo
        do j=1,NMAT
          dDdA_lambda(j,j,l,ii,jj,ll)=dDdA_lambda(j,j,l,ii,jj,ll)-trace/cmplx(dble(NMAT))
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      enddo
    enddo
  enddo ! ll
enddo ! f
enddo ! l

end subroutine calc_fermion_force_from_omega

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End Dirac operators
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Begin other subroutines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Gaussian random number
!!  BoxMuller(gauss,N)
!!
!! generage ensemble exp(-1/2(x^2) and exp(-1/2(y^2))
!! 
!! output 2N gaussian randum numbers
!! gauss is an allocatable array.
!! It will be reallocated after calling this routine.
subroutine BoxMuller(gauss,N)
use mt95
implicit none

double precision, parameter :: PI=dacos(-1d0)
integer, intent(in) :: N
double precision, allocatable :: gauss(:)
double precision rand(1:2*N)
integer i

if( allocated(gauss) ) then 
    deallocate( gauss )
endif

allocate( gauss(1:2*N) )

call genrand_real3(rand)
!write(*,*) rand

do i=1,N
  gauss(2*i-1) = dsqrt(-2d0*dlog(rand(2*i-1)))*dsin(2d0*Pi*rand(2*i))
  gauss(2*i) = dsqrt(-2d0*dlog(rand(2*i-1)))*dcos(2d0*Pi*rand(2*i))
enddo
end subroutine BoxMuller


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Gaussian random number
!!  BoxMuller(gauss,N)
!!
!! output 2N gaussian randum numbers
subroutine BoxMuller2(gauss,N)
use mt95
implicit none

double precision, parameter :: PI=dacos(-1d0)
integer, intent(in) :: N
double precision :: gauss(1:2*N)
double precision :: rand(1:2*N)
integer i

call genrand_real3(rand)

do i=1,N
  gauss(2*i-1) = dsqrt(-2d0*dlog(rand(2*i-1)))*dsin(2d0*Pi*rand(2*i))
  gauss(2*i) = dsqrt(-2d0*dlog(rand(2*i-1)))*dcos(2d0*Pi*rand(2*i))
enddo
end subroutine BoxMuller2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute forward covariant difference of Phi
!!  DMAT = U_l Phi_tip(l) U_l^\dagger - Phi_org(l)
subroutine Make_diff_Phi(Dphi, l,UMAT,Phi)
use SUN_generators, only : Make_traceless_matrix_from_modes
implicit none

complex(kind(0d0)), intent(out) :: Dphi(1:NMAT,1:NMAT)
integer, intent(in) :: l
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
integer :: i,j
complex(kind(0d0)) :: Phi_tip(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: Phi_org(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)

call Make_traceless_matrix_from_modes(Phi_tip,NMAT,Phi(:,link_tip(l)))
!call Make_traceless_matrix_from_modes(Phi_org,NMAT,Phi(:,link_org(l)))
call Make_traceless_matrix_from_modes(Dphi,NMAT,Phi(:,link_org(l)))

! U_l.Phi_tip
call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
    UMAT(:,:,l), NMAT, &
    Phi_tip, NMAT, &
    (0d0,0d0), tmpmat1, NMAT)
! U_l.Phi_tip.U_l^\dagger
call ZGEMM('N','C',NMAT,NMAT,NMAT,(1d0,0d0), &
    tmpmat1, NMAT, &
    UMAT(:,:,l), NMAT, &
    (-1d0,0d0), Dphi, NMAT)
! U_l.Phi_tip.U_l^\dagger - Phi_org
!Dphi = Dphi - Phi_org

end subroutine Make_diff_Phi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute forward covariant difference of Phi
!!  DMAT = U_l Phi_tip(l) U_l^\dagger - Phi_org(l)
subroutine Make_diff_PhiMat(Dphi, l,UMAT,PhiMat)
use SUN_generators, only : Make_traceless_matrix_from_modes
implicit none

complex(kind(0d0)), intent(out) :: Dphi(1:NMAT,1:NMAT)
integer, intent(in) :: l
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
integer :: i,j
!complex(kind(0d0)) :: Phi_tip(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: Phi_org(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)

!call Make_traceless_matrix_from_modes(Phi_tip,NMAT,Phi(:,link_tip(l)))
!call Make_traceless_matrix_from_modes(Phi_org,NMAT,Phi(:,link_org(l)))
!call Make_traceless_matrix_from_modes(Dphi,NMAT,Phi(:,link_org(l)))

Dphi=PhiMat(:,:,link_org(l))

! U_l.Phi_tip
call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
    UMAT(:,:,l), NMAT, &
    PhiMat(:,:,link_tip(l)), NMAT, &
    (0d0,0d0), tmpmat1, NMAT)
! U_l.Phi_tip.U_l^\dagger
call ZGEMM('N','C',NMAT,NMAT,NMAT,(1d0,0d0), &
    tmpmat1, NMAT, &
    UMAT(:,:,l), NMAT, &
    (-1d0,0d0), Dphi, NMAT)
! U_l.Phi_tip.U_l^\dagger - Phi_org
!Dphi = Dphi - Phi_org

end subroutine Make_diff_PhiMat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to make moment map \Omega(Uf)
subroutine Make_moment_map0(Omega,Uf)
use matrix_functions, only : matrix_inverse
implicit none

complex(kind(0d0)), intent(out) :: Omega(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT)
integer :: i,j

do i=1,NMAT
  do j=1,NMAT
    Omega(i,j)=-im_unit*( Uf(i,j) - dconjg( Uf(j,i) ) )
  enddo
enddo

end subroutine Make_moment_map0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to make moment map \Omega(Uf)
subroutine Make_moment_map(Omega,Ufm)
use matrix_functions, only : matrix_inverse
implicit none

complex(kind(0d0)), intent(out) :: Omega(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: Ufm(1:NMAT,1:NMAT)
complex(kind(0d0)) :: SMAT(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Cinv(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
integer :: i,j,k
complex(kind(0d0)) :: trace

!call matrix_power(NMAT,Uf,m_omega,Ufm)

do i=1,NMAT
  do j=1,NMAT
    SMAT(i,j)=-im_unit*( Ufm(i,j) - dconjg( Ufm(j,i) ) )
    Cinv(i,j)= Ufm(i,j) + dconjg( Ufm(j,i) ) 
  enddo
enddo
! CMAT --> CMAT^{-1}
call matrix_inverse(Cinv)

!call ZGEMM('N','N',NMAT,NMAT,NMAT,dcmplx(1d0/dble(m_omega)), &
call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0),&
  SMAT,NMAT,&
  Cinv,NMAT,&
  (0d0,0d0),Omega,NMAT)
call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0),&
  Cinv,NMAT,&
  SMAT,NMAT,&
  (1d0,0d0),Omega,NMAT)

!Omega=(0d0,0d0)
!do i=1,NMAT
!  do j=1,NMAT
!    Omega(i,j)=tmpmat(i,j)+dconjg(tmpmat(j,i))
!  enddo
!enddo

Omega = Omega / dcmplx( dble(m_omega) )

if (NMAT > 2) then
  trace=(0d0,0d0)
  do i=1,NMAT
    trace=trace+Omega(i,i)
  enddo
  do i=1,NMAT
    Omega(i,i)=Omega(i,i)-trace / cmplx(dble( NMAT ))
  enddo
endif

end subroutine Make_moment_map


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to make plaquette variable of face f
subroutine Make_face_variable(Uf,f,UMAT)
!use simplicial_complex, only : get_links_in_face_sc
implicit none

complex(kind(0d0)), intent(out) :: Uf(1:NMAT,1:NMAT)
integer, intent(in) :: f
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
!integer :: FaceSize
!integer, allocatable :: sites(:),link_labels(:),link_dirs(:)
character(1) :: char1
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
integer :: l,i,j


if(links_in_f(f)%link_dirs_(1)==1) then
  Uf=UMAT(:,:,links_in_f(f)%link_labels_(1))
elseif(links_in_f(f)%link_dirs_(1)==-1) then
  do i=1,NMAT
  do j=1,NMAT
    UF(i,j)=dconjg(UMAT(j,i,links_in_f(f)%link_labels_(1)))
  enddo
  enddo
else
  write(*,*) "There is a blanck link in the face",f
  stop 1
endif
!
do l=2,links_in_f(f)%num_
  tmpmat=Uf
  if(links_in_f(f)%link_dirs_(l)==1) then 
    char1='N'
  elseif(links_in_f(f)%link_dirs_(l)==-1) then
    char1='C'
  else
    write(*,*) "There is a blanck link in the face",f
    stop
  endif
  !! Uf=tmpmat*Uf or tmpmat*(Uf)^\dagger
  call ZGEMM('N',char1,NMAT,NMAT,NMAT,(1d0,0d0), &
    tmpmat,NMAT, &
    UMAT(:,:,links_in_f(f)%link_labels_(l)),NMAT, &
    (0d0,0d0), Uf, NMAT)
enddo

end subroutine Make_face_variable


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to compute \frac{d}{d A_l^a} \Omega_m
!!  Uf and Uf^m must be iput matrices.
subroutine calc_diff_omega(diff_Omega, Uf, Ufm, f, l,UMAT)
use matrix_functions, only : matrix_inverse
implicit none

complex(kind(0d0)), intent(out) :: diff_Omega(1:NMAT,1:NMAT,1:dimG)
integer, intent(in) :: f, l
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT), Ufm(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: diff_S(1:NMAT,1:NMAT,1:dimG), diff_C(1:NMAT,1:NMAT,1:dimG)
complex(kind(0d0)) :: SMAT(1:NMAT,1:NMAT), Cinv(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: SCinv(1:NMAT,1:NMAT)
integer :: i,j,a


do i=1,NMAT
  do j=1,NMAT
    SMAT(i,j)=-im_unit*( Ufm(i,j) - dconjg( Ufm(j,i) ) )
    Cinv(i,j)=( Ufm(i,j) + dconjg( Ufm(j,i) ) )
  enddo
enddo
call matrix_inverse(Cinv)

call calc_diff_SandC(diff_S,diff_C,Uf,f,l,UMAT)
! S.C^{-1}
  call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
    SMAT,NMAT, &
    Cinv,NMAT, &
    (0d0,0d0), SCinv, NMAT)

do a=1,dimG
  ! tmpmat3 = diff_S . C^{-1}
  call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
    diff_S(:,:,a),NMAT, &
    Cinv,NMAT, &
    (0d0,0d0), tmpmat, NMAT)
  
  ! tmpmat2 = S.C^{-1}.diff_C
  call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
    SCinv,NMAT, &
    diff_C(:,:,a),NMAT, &
    (0d0,0d0), tmpmat2, NMAT)
  ! tmpmat3 = diff_S.C^{-1} - S.C^{-1}.diff_C.C^{-1}
  call ZGEMM('N','N',NMAT,NMAT,NMAT,(-1d0,0d0), &
    tmpmat2,NMAT, &
    Cinv,NMAT, &
    (1d0,0d0), tmpmat, NMAT)

  ! diff_Omega = 1/m_omega * ( tmpmat3 + (tmpmat3)^\dagger ) 
  do i=1,NMAT
    do j=1,NMAT
      diff_Omega(i,j,a)= &
        ( tmpmat(i,j) + dconjg(tmpmat(j,i) ) ) ! / dcmplx(dble(m_omega))
    enddo
  enddo
enddo
diff_Omega=diff_Omega / dcmplx(dble(m_omega)) 

end subroutine calc_diff_omega


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to compute \frac{d}{d A_l^a} S and \frac{d}{d A_l^a} C
subroutine calc_diff_SandC(diff_S,diff_C, Uf, f, l,UMAT)

complex(kind(0d0)), intent(out) :: diff_S(1:NMAT,1:NMAT,1:dimG)
complex(kind(0d0)), intent(out) :: diff_C(1:NMAT,1:NMAT,1:dimG)
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
integer, intent(in) :: f, l

complex(kind(0d0)) :: diff_Ufm(1:NMAT,1:NMAT,1:dimG)
integer a,i,j

call calc_diff_Ufm(diff_Ufm, Uf,f,l,UMAT)
  
diff_S=(0d0,0d0)
diff_C=(0d0,0d0)
do a=1,dimG
  do j=1,NMAT
    do i=1,NMAT
      diff_S(i,j,a)=(0d0,-1d0)*( diff_Ufm(i,j,a) - dconjg( diff_Ufm(j,i,a) ) )
      !diff_S(i,j,a) = diff_Ufm(i,j,a) + dconjg( diff_Ufm(j,i,a) ) 
      diff_C(i,j,a) = diff_Ufm(i,j,a) + dconjg( diff_Ufm(j,i,a) )
    enddo
  enddo
enddo

end subroutine calc_diff_SandC



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to compute d/dA(a,l) Ufm
subroutine calc_diff_Ufm(diff_Ufm, Uf,f, l,UMAT)
use matrix_functions, only : matrix_power
implicit none

complex(kind(0d0)), intent(out) :: diff_Ufm(1:NMAT,1:NMAT,1:dimG)
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: XMAT(1:NMAT,1:NMAT), YMAT(1:NMAT,1:NMAT)
integer, intent(in) :: l,f

integer :: nl
complex(kind(0d0)) :: Ufla(1:NMAT,1:NMAT,1:dimG)
complex(kind(0d0)) :: Uf_power1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Uf_power2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
character(1) :: char1

integer :: i,j,k,ll,a

!HERE
call calc_diff_Uf(Ufla, Uf,f, l,UMAT)

if ( m_omega == 1 ) then 
  diff_Ufm=Ufla
  return
endif

!!! In the following, we can assume m_omega >= 2.
do a=1,dimG
! k=0
  !call matrix_power(NMAT,Uf,m_omega-1,Uf_power1)
  call matrix_power(Uf_power1,Uf,m_omega-1)
  call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
    Ufla(:,:,a),NMAT, &
    Uf_power1,NMAT, &
    (0d0,0d0), diff_Ufm(:,:,a), NMAT)
enddo
! k=1,...,m_omega-2
if ( m_omega >= 3 ) then
  do k=1, m_omega-2
    call matrix_power(Uf_power1,Uf,k)
    call matrix_power(Uf_power2,Uf,m_omega-k-1)
    do a=1,dimG
      call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
        Uf_power1,NMAT, &
        Ufla(:,:,a),NMAT, &
        (0d0,0d0), tmpmat, NMAT)
      call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
        tmpmat,NMAT, &
        Uf_power2,NMAT, &
        (1d0,0d0), diff_Ufm(:,:,a), NMAT)
    enddo
  enddo
endif
! k=m_omega-1
call matrix_power(Uf_power1,Uf,m_omega-1)
do a=1,dimG
  call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
    Uf_power1,NMAT, &
    Ufla(:,:,a),NMAT, &
    (1d0,0d0), diff_Ufm(:,:,a), NMAT)
enddo

end subroutine calc_diff_Ufm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to compute d/dA(a,l) Uf
subroutine calc_diff_Uf(diff_Uf, Uf,f, l,UMAT)
implicit none

complex(kind(0d0)), intent(out) :: diff_Uf(1:NMAT,1:NMAT,1:dimG)
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
integer, intent(in) :: l,f

complex(kind(0d0)) :: XMAT(1:NMAT,1:NMAT)
complex(kind(0d0)) :: YMAT(1:NMAT,1:NMAT)
complex(kind(0d0)) :: diff_UlAl(1:NMAT,1:NMAT,1:dimG)

integer :: a,link_place
complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)


! find the place of l1 and l2
do link_place=1,links_in_f(f)%num_
  if( links_in_f(f)%link_labels_(link_place) == l ) then 
    exit
  endif
enddo

! d/dA(a1,l1) d/dA(a2,l2) Ul
call calc_diff_Ul_in_face(diff_UlAl,f,link_place,UMAT)
! X=U(1)...U(l-1)
call calc_prodUl_from_n1_to_n2_in_Uf(XMAT,f,1,link_place-1,UMAT)
! Z=U(l+1)...U(f_size)
call calc_prodUl_from_n1_to_n2_in_Uf(YMAT,f,link_place+1,links_in_f(f)%num_,UMAT)
! X.diffdiff.Y
do a=1,dimG
  call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
     XMAT,NMAT, &
     diff_UlAl(:,:,a),NMAT, &
     (0d0,0d0), tmpmat1, NMAT)
  call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
     tmpmat1,NMAT, &
     YMAT,NMAT, &
     (0d0,0d0), diff_Uf(:,:,a), NMAT)
enddo

end subroutine calc_diff_Uf




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to compute 
!!   d/dA(a1,l1) d/dA(a2,l2) \Omega_m
!!  Uf and Uf^m must be iput matrices.
!!  This is aranged in order to reduce the amount of computation.
subroutine calc_diffdiff_omega(diffdiff_Omega, Uf, Ufm, f, l1,l2, UMAT)
use matrix_functions, only : hermitianmatrix_inverse
implicit none

complex(kind(0d0)), intent(out) :: diffdiff_Omega(1:NMAT,1:NMAT,1:dimG,1:dimG)
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT), Ufm(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_faces)
integer, intent(in) :: f, l1,l2

complex(kind(0d0)) :: SMAT(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Cinv(1:NMAT,1:NMAT)
complex(kind(0d0)) :: SCinv(1:NMAT,1:NMAT)
complex(kind(0d0)) :: diff_S_l1(1:NMAT,1:NMAT,1:dimG)
complex(kind(0d0)) :: diff_C_l1(1:NMAT,1:NMAT,1:dimG)
complex(kind(0d0)) :: diff_S_l2(1:NMAT,1:NMAT,1:dimG)
complex(kind(0d0)) :: diff_C_l2(1:NMAT,1:NMAT,1:dimG)
complex(kind(0d0)) :: diffdiff_S(1:NMAT,1:NMAT,1:dimG,1:dimG)
complex(kind(0d0)) :: diffdiff_C(1:NMAT,1:NMAT,1:dimG,1:dimG)
integer :: i,j,a1,a2

complex(kind(0d0)) :: tmp_dd_Omega(1:NMAT,1:NMAT,1:dimG,1:dimG)
complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat3(1:NMAT,1:NMAT)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! preparation
! SMAT and Cinv
do i=1,NMAT
  do j=1,NMAT
    SMAT(i,j)=-im_unit*( Ufm(i,j) - dconjg( Ufm(j,i) ) )
    Cinv(i,j)=( Ufm(i,j) + dconjg( Ufm(j,i) ) )
  enddo
enddo
!call matrix_inverse(NMAT,Cinv)
call hermitianmatrix_inverse(Cinv)

! diff_S_l1, diff_C_l1, diff_S_l2, diff_C_l2
call calc_diff_SandC(diff_S_l1,diff_C_l1,Uf,f,l1,UMAT)
call calc_diff_SandC(diff_S_l2,diff_C_l2,Uf,f,l2,UMAT)

! diffdiff_S and diffdiff_C
call calc_diffdiff_SandC(diffdiff_S, diffdiff_C, Uf, f, l1, l2, UMAT)

! S.Cinv
call ZHEMM('L','U',NMAT,NMAT,(1d0,0d0), &
  SMAT, NMAT, &
  Cinv, NMAT, &
  (0d0,0d0), SCinv, NMAT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! d/dA(a1,l1) d/dA(a2,l2) S . Cinv
do a1=1,dimG
  do a2=1,dimG
    call ZHEMM('R','U',NMAT,NMAT,(1d0,0d0), &
      Cinv, NMAT, &
      diffdiff_S(:,:,a1,a2), NMAT, &
      (0d0,0d0), tmp_dd_Omega(:,:,a1,a2), NMAT)
    !call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
    !  diffdiff_S(:,:,a1,a2), NMAT, &
    !  Cinv, NMAT, &
    !  (0d0,0d0), tmp_dd_Omega(:,:,a1,a2), NMAT)
  enddo
enddo

! add -S.Cinv. (d/dA(a1,l1) d/dA(a2,l2) C). Cinv
do a1=1,dimG
  do a2=1,dimG
    call ZHEMM('R','U',NMAT,NMAT,(1d0,0d0), &
      Cinv, NMAT, &
      diffdiff_C(:,:,a1,a2), NMAT, &
      (0d0,0d0), tmpmat2, NMAT)
    !call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
    !  diffdiff_C(:,:,a1,a2), NMAT, &
    !  Cinv, NMAT, &
    !  (0d0,0d0), tmpmat2, NMAT)
    !!!
    call ZGEMM('N','N',NMAT,NMAT,NMAT,(-1d0,0d0), &
      SCinv, NMAT, &
      tmpmat2, NMAT, &
      (1d0,0d0), tmp_dd_Omega(:,:,a1,a2), NMAT)
  enddo
enddo

do a1=1,dimG
  do a2=1,dimG
    ! d/dA(a1,l1)SMAT . Cinv
    call ZHEMM('R','U',NMAT,NMAT,(1d0,0d0), &
      Cinv, NMAT, &
      diff_S_l1(:,:,a1), NMAT, &
      (0d0,0d0), tmpmat1, NMAT)
    !call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
    !  diff_S_l1(:,:,a1), NMAT, &
    !  Cinv, NMAT, &
    !  (0d0,0d0), tmpmat1, NMAT)
    ! d/dA(a2,l2)CMAT . Cinv
    call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
      diff_C_l2(:,:,a2), NMAT, &
      Cinv, NMAT, &
      (0d0,0d0), tmpmat2, NMAT)
    ! 
    call ZGEMM('N','N',NMAT,NMAT,NMAT,(-1d0,0d0), &
      tmpmat1, NMAT, &
      tmpmat2, NMAT, &
      (1d0,0d0), tmp_dd_Omega(:,:,a1,a2), NMAT)
    !!!!!!!!!!!!!!!!
    ! d/dA(a2,l2)SMAT . Cinv
    call ZHEMM('R','U',NMAT,NMAT,(1d0,0d0), &
      Cinv, NMAT, &
      diff_S_l2(:,:,a2), NMAT, &
      (0d0,0d0), tmpmat1, NMAT)
    ! d/dA(a1,l1)CMAT . Cinv
    call ZHEMM('R','U',NMAT,NMAT,(1d0,0d0), &
      Cinv, NMAT, &
      diff_C_l1(:,:,a1), NMAT, &
      (0d0,0d0), tmpmat2, NMAT)
    ! 
    call ZGEMM('N','N',NMAT,NMAT,NMAT,(-1d0,0d0), &
      tmpmat1, NMAT, &
      tmpmat2, NMAT, &
      (1d0,0d0), tmp_dd_Omega(:,:,a1,a2), NMAT)
  enddo
enddo

do a1=1,dimG
  do a2=1,dimG
    ! tmpmat1 = S.Cinv.dC/dA(a1,l1).Cinv
    call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
      SCinv, NMAT, &
      diff_C_l1(:,:,a1), NMAT, &
      (0d0,0d0), tmpmat2, NMAT)
    call ZHEMM('R','U',NMAT,NMAT,(1d0,0d0), &
      Cinv, NMAT, &
      tmpmat2, NMAT, &
      (0d0,0d0), tmpmat1, NMAT)
    ! tmpmat2 = dC/dA(a2,l2).Cinv
    call ZHEMM('R','U',NMAT,NMAT,(1d0,0d0), &
      Cinv, NMAT, &
      diff_C_l2(:,:,a2), NMAT, &
      (0d0,0d0), tmpmat2, NMAT)
    ! add S.Cinv.dC/dA(a1,l1).Cinv.dC/dA(a2,l2).Cinv
    call ZGEMM('N','N',NMAT,NMAT,NMAT,(+1d0,0d0), &
      tmpmat1, NMAT, &
      tmpmat2, NMAT, &
      (1d0,0d0), tmp_dd_Omega(:,:,a1,a2), NMAT)
    !!!!!!!!!!!!!
    ! tmpmat1 = S.Cinv.dC/dA(a2,l2).Cinv
    call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
      SCinv, NMAT, &
      diff_C_l2(:,:,a2), NMAT, &
      (0d0,0d0), tmpmat2, NMAT)
    call ZHEMM('R','U',NMAT,NMAT,(1d0,0d0), &
      Cinv, NMAT, &
      tmpmat2, NMAT, &
      (0d0,0d0), tmpmat1, NMAT)
    ! tmpmat2 = dC/dA(a1,l1).Cinv
    call ZHEMM('R','U',NMAT,NMAT,(1d0,0d0), &
      Cinv, NMAT, &
      diff_C_l1(:,:,a1), NMAT, &
      (0d0,0d0), tmpmat2, NMAT)
    ! add S.Cinv.dC/dA(a2,l2).Cinv.dC/dA(a1,l1).Cinv
    call ZGEMM('N','N',NMAT,NMAT,NMAT,(+1d0,0d0), &
      tmpmat1, NMAT, &
      tmpmat2, NMAT, &
      (1d0,0d0), tmp_dd_Omega(:,:,a1,a2), NMAT)
  enddo
enddo


do a2=1,dimG
  do a1=1,dimG
    do j=1,NMAT
      do i=1,NMAT
        diffdiff_Omega(i,j,a1,a2)=&
          ( tmp_dd_Omega(i,j,a1,a2) + dconjg( tmp_dd_Omega(j,i,a1,a2) ) ) &
           / dcmplx(dble(m_omega))
      enddo
    enddo
  enddo
enddo

end subroutine calc_diffdiff_omega

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to compute \frac{d}{d A_l^a} S and \frac{d}{d A_l^a} C
subroutine calc_diffdiff_SandC(diffdiff_S,diffdiff_C, Uf,f,l1,l2,UMAT)

complex(kind(0d0)), intent(out) :: diffdiff_S(1:NMAT,1:NMAT,1:dimG,1:dimG)
complex(kind(0d0)), intent(out) :: diffdiff_C(1:NMAT,1:NMAT,1:dimG,1:dimG)
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
integer, intent(in) :: f,l1,l2

complex(kind(0d0)) :: diffdiff_Ufm(1:NMAT,1:NMAT,1:dimG,1:dimG)
integer :: i,j,a1,a2

call calc_diffdiff_Ufm(diffdiff_Ufm,Uf,f,l1,l2,UMAT)

do a2=1,dimG
  do a1=1,dimG
    do j=1,NMAT
      do i=1,NMAT
        diffdiff_S(i,j,a1,a2) = (0d0,-1d0)  &
          * ( diffdiff_Ufm(i,j,a1,a2) - dconjg( diffdiff_Ufm(j,i,a1,a2) ) )
        !diffdiff_S(i,j,a1,a2) = &
        !  diffdiff_Ufm(i,j,a1,a2) + dconjg( diffdiff_Ufm(j,i,a1,a2) ) 
        diffdiff_C(i,j,a1,a2) = &
          diffdiff_Ufm(i,j,a1,a2) + dconjg( diffdiff_Ufm(j,i,a1,a2) ) 
      enddo
    enddo
  enddo
enddo


!do i=1,NMAT
!  do j=1,NMAT
!    do a1=1,dimG
!      do a2=1,dimG
!        write(*,*) diffdiff_S(i,j,a1,a2) - dconjg( diffdiff_S(j,i,a1,a2) ) 
!    enddo
!enddo
!enddo
!enddo

end subroutine calc_diffdiff_SandC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to compute d/dA(a2,l2) d/dA(a1,l1) Ufm
!! <notation>
!!   diffdiff_Ufm(:,;,a1,a2) = d/dA(a1,l1) d/dA(a2,l2) Ufm
subroutine calc_diffdiff_Ufm(diffdiff_Ufm, Uf,f,l1,l2,UMAT)
use matrix_functions, only : matrix_power
implicit none

complex(kind(0d0)), intent(out) :: diffdiff_Ufm(1:NMAT,1:NMAT,1:dimG,1:dimG)
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
integer, intent(in) :: f,l1,l2

complex(kind(0d0)) :: diffdiff_Uf(1:NMAT,1:NMAT,1:dimG,1:dimG)
complex(kind(0d0)) :: diff_Uf_l1(1:NMAT,1:NMAT,1:dimG)
complex(kind(0d0)) :: diff_Uf_l2(1:NMAT,1:NMAT,1:dimG)
complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT),tmpmat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat3(1:NMAT,1:NMAT)
integer :: j,k
integer :: a1,a2

!! diffdiff_Uf(:,:,a1,a2)=d/dA(a1,l1) d/dA(a2,l2) Uf
call calc_diffdiff_Uf(diffdiff_Uf,Uf,f,l1,l2,UMAT)
!! diff_Uf_l1(:,:,a1)=d/dA(a1,l1) Uf
call calc_diff_Uf(diff_Uf_l1, Uf, f, l1, UMAT)
!! diff_Uf_l2(:,:,a2)=d/dA(a2,l2) Uf
call calc_diff_Uf(diff_Uf_l2, Uf, f, l2, UMAT)

diffdiff_Ufm=(0d0,0d0)

if (m_omega == 1) then 
  diffdiff_Ufm=diffdiff_Uf
  return
endif

if ( m_omega >= 2 ) then
  do k=0, m_omega-1
    if ( k .ne. 0 ) then 
      call matrix_power(tmpmat1,Uf,k)
    endif
    if( k .ne. m_omega-1 ) then 
      call matrix_power(tmpmat2,Uf,m_omega-k-1)
    endif
    do a1=1,dimG
      do a2=1,dimG
        if( k.ne.0 ) then 
          call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat1, NMAT, &
            diffdiff_Uf(:,:,a1,a2), NMAT, &
            (0d0,0d0), tmpmat3, NMAT)
        else
          tmpmat3=diffdiff_Uf(:,:,a1,a2)
        endif
        !!!
        if( k .ne. m_omega-1 ) then
          call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat3, NMAT, &
            tmpmat2, NMAT, &
            (1d0,0d0), diffdiff_Ufm(:,:,a1,a2), NMAT)
        else
          diffdiff_Ufm(:,:,a1,a2)=diffdiff_Ufm(:,:,a1,a2)+tmpmat3
        endif
      enddo
    enddo
  enddo

  do k=0,m_omega-2
    do j=0,m_omega-k-2
      do a1=1,dimG
        do a2=1,dimG
          ! Uf^j
          call matrix_power(tmpmat1,Uf,j)
          ! Uf^j.diff_Uf_(a1,l1)
          call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat1, NMAT, &
            diff_Uf_l1(:,:,a1), NMAT, &
            (0d0,0d0), tmpmat2, NMAT)
          ! Uf^k
          call matrix_power(tmpmat1,Uf,k)
          ! Uf^j.diff_Uf_(a1,l1).Uf^k
          call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat2, NMAT, &
            tmpmat1, NMAT, &
            (0d0,0d0), tmpmat3, NMAT)
          ! Uf^j.diff_Uf_l1(a1).Uf^k.diff_Uf_l2(a2)
          call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat3, NMAT, &
            diff_Uf_l2(:,:,a2), NMAT, &
            (0d0,0d0), tmpmat2, NMAT)
          ! Uf^{m-j-k-2}
            call matrix_power(tmpmat1,Uf,m_omega-j-k-2)
          ! Uf^j.diff_Uf_(a1).Uf^k.diff_Uf_l2(a2).Uf^{m-j-k-2}
          call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat2, NMAT, &
            tmpmat1, NMAT, &
            (1d0,0d0), diffdiff_Ufm(:,:,a1,a2), NMAT)
        enddo
      enddo
      !!!
      do a1=1,dimG
        do a2=1,dimG
          ! Uf^j
          call matrix_power(tmpmat1,Uf,j)
          ! Uf^j.diff_Uf_(a2,l2)
          call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat1, NMAT, &
            diff_Uf_l2(:,:,a2), NMAT, &
            (0d0,0d0), tmpmat2, NMAT)
          ! Uf^k
          call matrix_power(tmpmat1,Uf,k)
          ! Uf^j.diff_Uf_(a2,l2).Uf^k
          call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat2, NMAT, &
            tmpmat1, NMAT, &
            (0d0,0d0), tmpmat3, NMAT)
          ! Uf^j.diff_Uf_l2(a2).Uf^k.diff_Uf_l1(a1)
          call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat3, NMAT, &
            diff_Uf_l1(:,:,a1), NMAT, &
            (0d0,0d0), tmpmat2, NMAT)
          ! Uf^{m-j-k-2}
            call matrix_power(tmpmat1,Uf,m_omega-j-k-2)
          ! Uf^j.diff_Uf_(a1).Uf^k.diff_Uf_l2(a2).Uf^{m-j-k-2}
          call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat2, NMAT, &
            tmpmat1, NMAT, &
            (1d0,0d0), diffdiff_Ufm(:,:,a1,a2), NMAT)
        enddo
      enddo
    enddo
  enddo
!do a1=1,dimG
!  do a2=1,dimG
!    write(*,*)  diffdiff_Ufm(:,:,a1,a2)-diffdiff_Ufm(:,:,a2,a1)
!  enddo
!enddo
endif



end subroutine calc_diffdiff_Ufm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to compute d/dA(a1,l1) d/dA(a2,l2) Uf
!! <notation>
!!   diffdiff_Uf(:,;,a1,a2) = d/dA(a1,l1) d/dA(a2,l2) Uf
subroutine calc_diffdiff_Uf(diffdiff_Uf, Uf, f, l1, l2, UMAT)

complex(kind(0d0)), intent(out) :: diffdiff_Uf(1:NMAT,1:NMAT,1:dimG,1:dimG)
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
integer, intent(in) :: l1,l2,f

integer :: place_l1, place_l2
complex(kind(0d0)) :: XMAT(1:NMAT,1:NMAT)
complex(kind(0d0)) :: YMAT(1:NMAT,1:NMAT)
complex(kind(0d0)) :: ZMAT(1:NMAT,1:NMAT)
complex(kind(0d0)) :: diff_UlAl1(1:NMAT,1:NMAT,1:dimG)
complex(kind(0d0)) :: diff_UlAl2(1:NMAT,1:NMAT,1:dimG)
complex(kind(0d0)) :: diffdiff_UlAl(1:NMAT,1:NMAT,1:dimG,1:dimG)
complex(kind(0d0)) :: tmp_diffdiff_Uf(1:NMAT,1:NMAT,1:dimG,1:dimG)
complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT), tmpmat2(1:NMAT,1:NMAT)
character(1) :: char1
integer :: link_place, link_place_l1, link_place_l2
integer :: nl
integer :: a1,a2,info

!! when l1=l2
if ( l1 == l2 ) then 
  ! find the place of l1 and l2
  do link_place_l1=1,links_in_f(f)%num_
    if( links_in_f(f)%link_labels_(link_place_l1) == l1 ) then 
      exit
    endif
  enddo
  ! d/dA(a1,l1) d/dA(a2,l2) Ul
  call calc_diffdiff_Ul_in_face(diffdiff_UlAl,f,link_place_l1,UMAT)
  ! X=U(1)...U(l1-1)
  call calc_prodUl_from_n1_to_n2_in_Uf(XMAT,f,1,link_place_l1-1,UMAT)
  ! Z=U(l1+1)...U(f_size)
  call calc_prodUl_from_n1_to_n2_in_Uf(ZMAT,f,link_place_l1+1,links_in_f(f)%num_,UMAT)
  ! X.diffdiff.Z
  do a1=1,dimG
    do a2=1,dimG
      call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
         XMAT,NMAT, &
         diffdiff_UlAl(:,:,a1,a2),NMAT, &
         (0d0,0d0), tmpmat1, NMAT)
      call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
         tmpmat1,NMAT, &
         ZMAT,NMAT, &
         (0d0,0d0), diffdiff_Uf(:,:,a1,a2), NMAT)
    enddo
  enddo
!do a1=1,dimG
!  do a2=1,dimG
!    write(*,*) diffdiff_Uf(:,:,a1,a2)-diffdiff_Uf(:,:,a2,a1)
!  enddo
!enddo
  return
!!!!!!!!!!!!!!
! when l_1 != l_2
else
! find the places of l1 and l2
  info=0
  do nl=1,links_in_f(f)%num_
    if( links_in_f(f)%link_labels_(nl) == l1 ) then 
      link_place_l1=nl
      info=info+1
    elseif( links_in_f(f)%link_labels_(nl) == l2 ) then 
      link_place_l2=nl
      info=info+1
    endif
    if(info==2) exit
  enddo
  if(info*(info-1)==0) then 
    write(*,*) "l_1 and/or l_2 are not included in the face f"
    stop
  endif

  call calc_diff_Ul_in_face(diff_UlAl1,f,link_place_l1,UMAT)
  call calc_diff_Ul_in_face(diff_UlAl2,f,link_place_l2,UMAT)

  !!!!!!!!!!!!!!!!!!
  if(link_place_l1 < link_place_l2) then 
    call calc_prodUl_from_n1_to_n2_in_Uf(XMAT,f,1,link_place_l1-1,UMAT)
    call calc_prodUl_from_n1_to_n2_in_Uf(YMAT,f,link_place_l1+1,link_place_l2-1,UMAT)
    call calc_prodUl_from_n1_to_n2_in_Uf(ZMAT,f,link_place_l2+1,links_in_f(f)%num_,UMAT)
  do a1=1,dimG
    do a2=1,dimG
      call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
         XMAT,NMAT, &
         diff_UlAl1(:,:,a1),NMAT, &
         (0d0,0d0), tmpmat1, NMAT)
      call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
         tmpmat1,NMAT, &
         YMAT,NMAT, &
         (0d0,0d0), tmpmat2, NMAT)
      call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
         tmpmat2,NMAT, &
         diff_UlAl2(:,:,a2),NMAT, &
         (0d0,0d0), tmpmat1, NMAT)
      call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
         tmpmat1,NMAT, &
         ZMAT,NMAT, &
         (0d0,0d0),  diffdiff_Uf(:,:,a1,a2), NMAT)
    enddo
  enddo
  !!!!!!!!!!!!!!!!!!1
  else
    call calc_prodUl_from_n1_to_n2_in_Uf(XMAT,f,1,link_place_l2-1,UMAT)
    call calc_prodUl_from_n1_to_n2_in_Uf(YMAT,f,link_place_l2+1,link_place_l1-1,UMAT)
    call calc_prodUl_from_n1_to_n2_in_Uf(ZMAT,f,link_place_l1+1,links_in_f(f)%num_,UMAT)
  do a1=1,dimG
    do a2=1,dimG
      call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
         XMAT,NMAT, &
         diff_UlAl2(:,:,a2),NMAT, &
         (0d0,0d0), tmpmat1, NMAT)
      call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
         tmpmat1,NMAT, &
         YMAT,NMAT, &
         (0d0,0d0), tmpmat2, NMAT)
      call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
         tmpmat2,NMAT, &
         diff_UlAl1(:,:,a1),NMAT, &
         (0d0,0d0), tmpmat1, NMAT)
      call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
         tmpmat1,NMAT, &
         ZMAT,NMAT, &
         (0d0,0d0),  diffdiff_Uf(:,:,a1,a2), NMAT)
    enddo
  enddo
  endif 
endif

return

end subroutine calc_diffdiff_Uf


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! caluculate product of Ul's from l1 to l2 in the face f.
subroutine calc_prodUl_from_n1_to_n2_in_Uf(ProdU,f,n1,n2,UMAT)
use matrix_functions, only : make_unit_matrix
implicit none

complex(kind(0d0)), intent(out) :: ProdU(1:NMAT,1:NMAT)
integer, intent(in) :: f,n1,n2
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)

complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
character(1) :: char1
integer :: link_place,i,j

!! if there is no links in this period, return the unit matrix
if ( n2 <= n1-1 ) then 
  call make_unit_matrix(ProdU)
  return
endif

if(links_in_f(f)%link_dirs_(n1)==1) then 
  ProdU=UMAT(:,:,links_in_f(f)%link_labels_(n1))
else
  do i=1,NMAT
    do j=1,NMAT
      ProdU(i,j)=dconjg(UMAT(j,i,links_in_f(f)%link_labels_(n1)))
    enddo
  enddo
endif
if ( n2 >= n1+1 ) then 
  do link_place=n1+1,n2
    tmpmat=ProdU
    if(links_in_f(f)%link_dirs_(link_place)==1) then 
      char1='N'
    elseif(links_in_f(f)%link_dirs_(link_place)==-1) then
      char1='C'
    endif
    !! Uf=tmpmat*Uf or tmpmat*(Uf)^\dagger
    call ZGEMM('N',char1,NMAT,NMAT,NMAT,(1d0,0d0), &
      tmpmat,NMAT, &
      UMAT(:,:,links_in_f(f)%link_labels_(link_place)),NMAT, &
      (0d0,0d0), ProdU, NMAT)
  enddo
endif

end subroutine calc_prodUl_from_n1_to_n2_in_Uf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate d/dA(a,l) U_l or U_l^{-1} in Uf for all a
!!  where l is the n-th link in f
subroutine calc_diff_Ul_in_face(diff_UlAl,f,link_place,UMAT)
use SUN_generators, only : MtimesT, TtimesM
use SUN_generators, only : make_SUN_generators ! for test
implicit none

complex(kind(0d0)), intent(out) :: diff_UlAl(1:NMAT,1:NMAT,1:dimG)
integer, intent(in) :: f,link_place
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)

complex(kind(0d0)) :: Uinv_times_minus_i(1:NMAT,1:NMAT)
integer :: l,a,i,j


l=links_in_f(f)%link_labels_(link_place)

if ( links_in_f(f)%link_dirs_(link_place) == 1 ) then 
  do a=1,dimG
    call TtimesM(diff_UlAl(:,:,a), (0d0,1d0)*UMAT(:,:,l), a, NMAT) ! NO BUG
  enddo
else
  do i=1,NMAT
    do j=1,NMAT
      Uinv_times_minus_i(i,j)=(0d0,-1d0)*dconjg(UMAT(j,i,l))
    enddo
  enddo
  do a=1,dimG
    call MtimesT(diff_UlAl(:,:,a), Uinv_times_minus_i, a, NMAT) ! NO BUG
  enddo
endif

end subroutine calc_diff_Ul_in_face


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate d/dA(a1,l) d/dA(a2,l) U_l or U_l^{-1} in Uf
!!  for all a1 and a2,  where l is the n-th link in f
subroutine calc_diffdiff_Ul_in_face(diffdiff_UlAl,f,link_place,UMAT)
use SUN_generators, only : TtimesM, MtimesT
implicit none

complex(kind(0d0)), intent(out) :: diffdiff_UlAl(1:NMAT,1:NMAT,1:dimG,1:dimG)
integer, intent(in) :: f,link_place
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)

complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Uinv_times_minus_1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmp_diffdiff_UlAl(1:NMAT,1:NMAT,1:dimG,1:dimG)
integer :: l,a1,a2,i,j

l=links_in_f(f)%link_labels_(link_place)

if ( links_in_f(f)%link_dirs_(link_place) == 1 ) then 
  do a2=1,dimG
      call TtimesM(tmpmat, UMAT(:,:,l), a2, NMAT)
    do a1=1,dimG
      call TtimesM(diffdiff_UlAl(:,:,a2,a1),-tmpmat,a1,NMAT)
    enddo
  enddo
else
  do i=1,NMAT
    do j=1,NMAT
      Uinv_times_minus_1(i,j)=-dconjg(UMAT(j,i,l))
    enddo
  enddo
  do a2=1,dimG
      call MtimesT(tmpmat, Uinv_times_minus_1(:,:), a2, NMAT)
    do a1=1,dimG
      call MtimesT(diffdiff_UlAl(:,:,a2,a1),tmpmat,a1,NMAT)
    enddo
  enddo
endif

end subroutine calc_diffdiff_Ul_in_face


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!
integer function site_index(a,s) 
implicit none

integer, intent(in) :: a,s

site_index=a+dimG*(s-1)

end function site_index

!!!!!!!!!!!!!
integer function link_index(a,l) 
implicit none

integer, intent(in) :: a,l

link_index=a+dimG*(num_sites + l - 1)

end function link_index

!!!!!!!!!!!!!
integer function face_index(a,f)
implicit none

integer, intent(in) :: a,f

face_index=a+dimG*(num_sites + num_links + f - 1)

end function face_index


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to compute d/dA(a,l) Uf
subroutine tmp_calc_diff_Uf(Ufla, Uf,f, l,UMAT)
use SUN_generators, only : TtimesM, MtimesT, MTN
implicit none

complex(kind(0d0)), intent(out) :: Ufla(1:NMAT,1:NMAT,1:dimG)
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: XMAT(1:NMAT,1:NMAT), YMAT(1:NMAT,1:NMAT)
integer, intent(in) :: l,f

integer :: nl
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
character(1) :: char1
integer :: i,j,k,ll,a

do nl=1,links_in_f(f)%num_
  if( links_in_f(f)%link_labels_(nl) == l ) then 
    exit
  endif
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute \frac{d}{d A_l^a} U_f
!! when Ul is the first link
if( nl == 1 ) then 
!!!
! X=1, Y=U_l^{-1} U_f ( dir=1 )
  if( links_in_f(f)%link_dirs_(nl) == 1 ) then
    do a=1,dimG
      call TtimesM(Ufla(:,:,a), Uf, a, NMAT)
    enddo
    Ufla = (0d0,1d0)*Ufla
! X=1, Y=U_l U_f      ( dir=-1 )
  elseif( links_in_f(f)%link_dirs_(nl) == -1 ) then
    do i=1,NMAT
      do j=1,NMAT
        ! X \equiv U_l^{-1}
        XMAT(i,j)=dconjg( UMAT(j,i,l) )
      enddo
    enddo
    ! Y=U_l.U_f
    call ZGEMM('N','N',NMAT,NMAT,NMAT,(0d0,-1d0), &
      UMAT(:,:,l),NMAT, &
      Uf, NMAT, &
      (0d0,0d0), YMAT, NMAT)
    ! Uf_l^a = -i Ul^{-1}.T_a.Y
    do a=1,dimG
      call MTN(Ufla(:,:,a),XMAT,YMAT,a,NMAT)
    enddo
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! when Ul is the last link
elseif( nl == links_in_f(f)%num_  ) then
  if( links_in_f(f)%link_dirs_(nl) == 1) then
!  X=U_f U_l^{-1}, Y=1 (dir=1)
    call ZGEMM('N','C',NMAT,NMAT,NMAT,(0d0,1d0), &
      Uf,NMAT, &
      UMAT(:,:,l), NMAT, &
      (0d0,0d0), XMAT, NMAT)
    do a=1,dimG
      call MTN(Ufla(:,:,a), XMAT, UMAT(:,:,l),a, NMAT)
    enddo
  elseif( links_in_f(f)%link_dirs_(nl) == -1 ) then
!  X=U_f U_l,      Y=1 (dir=-1)
    do a=1,dimG
      call MtimesT(Ufla(:,:,a),Uf,a,NMAT)
    enddo
    Ufla=(0d0,-1d0)*Ufla
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Other cases
else
! X=U_1 ... U_{nl-1}
  if(links_in_f(f)%link_dirs_(1)==1) then
    XMAT=UMAT(:,:,links_in_f(f)%link_labels_(1))
  elseif(links_in_f(f)%link_dirs_(1)==-1) then
    do i=1,NMAT
    do j=1,NMAT
      XMAT(i,j)=dconjg(UMAT(j,i,links_in_f(f)%link_labels_(1)))
    enddo
    enddo
  endif
  !
  if(nl >= 3) then
  do ll=2,nl-1
    tmpmat=XMAT
    if(links_in_f(f)%link_dirs_(ll)==1) then 
      char1='N'
    elseif(links_in_f(f)%link_dirs_(ll)==-1) then
      char1='C'
    endif
    !! Uf=tmpmat*Uf or tmpmat*(Uf)^\dagger
    call ZGEMM('N',char1,NMAT,NMAT,NMAT,(1d0,0d0), &
      tmpmat,NMAT, &
      UMAT(:,:,links_in_f(f)%link_labels_(ll)),NMAT, &
      (0d0,0d0), XMAT, NMAT)
  enddo
  endif
! Y=U_{nl+1}...U_{links_in_f(f)%num_ }
  if(links_in_f(f)%link_dirs_(nl+1)==1) then
    YMAT=UMAT(:,:,links_in_f(f)%link_labels_(nl+1))
  elseif(links_in_f(f)%link_dirs_(nl+1)==-1) then
    do i=1,NMAT
    do j=1,NMAT
      YMAT(i,j)=dconjg(UMAT(j,i,links_in_f(f)%link_labels_(nl+1)))
    enddo
    enddo
  endif
  !
  if(nl <= links_in_f(f)%num_ -2) then
  do ll=nl+2,links_in_f(f)%num_ 
    tmpmat=YMAT
    if(links_in_f(f)%link_dirs_(ll)==1) then 
      char1='N'
    elseif(links_in_f(f)%link_dirs_(ll)==-1) then
      char1='C'
    endif
    !! Uf=tmpmat*Uf or tmpmat*(Uf)^\dagger
    call ZGEMM('N',char1,NMAT,NMAT,NMAT,(1d0,0d0), &
      tmpmat,NMAT, &
      UMAT(:,:,links_in_f(f)%link_labels_(ll)),NMAT, &
      (0d0,0d0), YMAT, NMAT)
  enddo
  endif

  if( links_in_f(f)%link_dirs_(nl) == 1 ) then
    ! tmpmat = i U_l Y
    call ZGEMM('N','N',NMAT,NMAT,NMAT,(0d0,1d0), &
      UMAT(:,:,l),NMAT, &
      YMAT,NMAT, &
      (0d0,0d0), tmpmat, NMAT)
    do a=1,dimG
      call MTN(Ufla(:,:,a),XMAT,tmpmat,a,NMAT)
    enddo
  elseif( links_in_f(f)%link_dirs_(nl) == -1 ) then
    ! tmpmat = -i X U_l^{-1}
    call ZGEMM('N','C',NMAT,NMAT,NMAT,(0d0,-1d0), &
      XMAT,NMAT, &
      UMAT(:,:,l),NMAT, &
      (0d0,0d0), tmpmat, NMAT)
    do a=1,dimG
      call MTN(Ufla(:,:,a),tmpmat,YMAT,a,NMAT)
    enddo
  endif
endif

end subroutine tmp_calc_diff_Uf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to compute d/dA(b,l_2) d/dA(a,l_1) Uf
!! <notation>
!!   diffdiff_Uf(:,;,b,a) = d/dA(b,l_2) d/dA(a,l_1) Uf
subroutine tmp_calc_diffdiff_Uf(diffdiff_Uf, Uf, f, l_2, l_1, UMAT)
use matrix_functions, only : make_unit_matrix
use SUN_generators, only : TtimesM, MtimesT, MTN
implicit none

! 
complex(kind(0d0)), intent(out) :: diffdiff_Uf(1:NMAT,1:NMAT,1:dimG,1:dimG)
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
integer, intent(in) :: l_1,l_2,f
complex(kind(0d0)) :: XMAT(1:NMAT,1:NMAT), YMAT(1:NMAT,1:NMAT), ZMAT(1:NMAT,1:NMAT)

integer :: nl
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT), tmpmat2(1:NMAT,1:NMAT)
character(1) :: char1
integer :: i,j,k,ll,l,a,b,a1,a2
integer :: nl1, nl2, l1, l2, info

!! when l_1 = l_2 
if( l_1 == l_2 ) then 
  l=l_1

  do nl=1,links_in_f(f)%num_
    if( links_in_f(f)%link_labels_(nl) == l ) then 
      exit
    endif
  enddo
  !! when Ul is the first link
  if( nl == 1 ) then 
  !!!
  ! X=1, Y=U_l^{-1} U_f ( dir=1 )
    if( links_in_f(f)%link_dirs_(nl) == 1 ) then
      do a=1,dimG
        call TtimesM(tmpmat, Uf, a, NMAT)
        do b=1,dimG
          call TtimesM(diffdiff_Uf(:,:,b,a), tmpmat, b, NMAT)
        enddo
      enddo
      diffdiff_Uf = (-1d0,0d0)*diffdiff_Uf
  ! X=1, Y=U_l U_f      ( dir=-1 )
    elseif( links_in_f(f)%link_dirs_(nl) == -1 ) then
      do i=1,NMAT
        do j=1,NMAT
          ! X \equiv U_l^{-1}
          XMAT(i,j)=dconjg( UMAT(j,i,l) )
        enddo
      enddo
      !Y \equiv -U_l.U_f 
      call ZGEMM('N','N',NMAT,NMAT,NMAT,(-1d0,0d0), &
        UMAT(:,:,l),NMAT, &
        Uf, NMAT, &
        (0d0,0d0), YMAT, NMAT)
      ! 
      do a=1,dimG
        call TtimesM(tmpmat,YMAT,b,NMAT)
        do b=1,dimG
          call MTN(diffdiff_Uf(:,:,b,a),XMAT,tmpmat,a,NMAT)
        enddo
      enddo
    endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! when Ul is the last link
  elseif( nl == links_in_f(f)%num_  ) then
    if( links_in_f(f)%link_dirs_(nl) == 1) then
  !  X=-U_f U_l^{-1}, Y=1 (dir=1)
      call ZGEMM('N','C',NMAT,NMAT,NMAT,(-1d0,0d0), &
        Uf,NMAT, &
        UMAT(:,:,l), NMAT, &
        (0d0,0d0), XMAT, NMAT)
      do a=1,dimG
        ! tmpmat = Ta . UMAT
        call TtimesM(tmpmat, UMAT(:,:,l),a, NMAT)
        do b=1,dimG
          call MTN(diffdiff_Uf(:,:,b,a),XMAT,tmpmat,b,NMAT)
        enddo
      enddo
    elseif( links_in_f(f)%link_dirs_(nl) == -1 ) then
  !  X=U_f U_l,      Y=1 (dir=-1)
      do a=1,dimG
        call MtimesT(tmpmat,Uf,a,NMAT)
        do b=1,dimG
          call MtimesT(diffdiff_Uf(:,:,b,a),tmpmat,b,NMAT)
        enddo
      enddo
      diffdiff_Uf=(-1d0,0d0)*diffdiff_Uf
    endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Other cases
  else
  ! X=U_1 ... U_{nl-1}
    if(links_in_f(f)%link_dirs_(1)==1) then
      XMAT=UMAT(:,:,links_in_f(f)%link_labels_(1))
    elseif(links_in_f(f)%link_dirs_(1)==-1) then
      do i=1,NMAT
      do j=1,NMAT
        XMAT(i,j)=dconjg(UMAT(j,i,links_in_f(f)%link_labels_(1)))
      enddo
      enddo
    endif
    !
    if(nl >= 3) then
    do ll=2,nl-1
      tmpmat=XMAT
      if(links_in_f(f)%link_dirs_(ll)==1) then 
        char1='N'
      elseif(links_in_f(f)%link_dirs_(ll)==-1) then
        char1='C'
      endif
      !! Uf=tmpmat*Uf or tmpmat*(Uf)^\dagger
      call ZGEMM('N',char1,NMAT,NMAT,NMAT,(1d0,0d0), &
        tmpmat,NMAT, &
        UMAT(:,:,links_in_f(f)%link_labels_(ll)),NMAT, &
        (0d0,0d0), XMAT, NMAT)
    enddo
    endif
  ! Y=U_{nl+1}...U_{links_in_f(f)%num_ }
    if(links_in_f(f)%link_dirs_(nl+1)==1) then
      YMAT=UMAT(:,:,links_in_f(f)%link_labels_(nl+1))
    elseif(links_in_f(f)%link_dirs_(nl+1)==-1) then
      do i=1,NMAT
      do j=1,NMAT
        YMAT(i,j)=dconjg(UMAT(j,i,links_in_f(f)%link_labels_(nl+1)))
      enddo
      enddo
    endif
    !
    if(nl <= links_in_f(f)%num_ -2) then
    do ll=nl+2,links_in_f(f)%num_ 
      tmpmat=YMAT
      if(links_in_f(f)%link_dirs_(ll)==1) then 
        char1='N'
      elseif(links_in_f(f)%link_dirs_(ll)==-1) then
        char1='C'
      endif
      !! Uf=tmpmat*Uf or tmpmat*(Uf)^\dagger
      call ZGEMM('N',char1,NMAT,NMAT,NMAT,(1d0,0d0), &
        tmpmat,NMAT, &
        UMAT(:,:,links_in_f(f)%link_labels_(ll)),NMAT, &
        (0d0,0d0), YMAT, NMAT)
    enddo
    endif
  
    if( links_in_f(f)%link_dirs_(nl) == 1 ) then
      ! tmpmat = - U_l Y
      call ZGEMM('N','N',NMAT,NMAT,NMAT,(-1d0,0d0), &
        UMAT(:,:,l),NMAT, &
        YMAT,NMAT, &
        (0d0,0d0), tmpmat, NMAT)
      ! tmpmat2 = -Ta U_l Y
      do a=1,dimG
        call TtimesM(tmpmat2,tmpmat,a,NMAT)
        do b=1,dimG
          call MTN(diffdiff_Uf(:,:,b,a),XMAT,tmpmat2,b,NMAT)
        enddo
      enddo
    elseif( links_in_f(f)%link_dirs_(nl) == -1 ) then
      ! tmpmat = - X U_l^{-1}
      call ZGEMM('N','C',NMAT,NMAT,NMAT,(-1d0,0d0), &
        XMAT,NMAT, &
        UMAT(:,:,l),NMAT, &
        (0d0,0d0), tmpmat, NMAT)
      ! tmpmat2 = -X U_l^{-1} Ta
      do a=1,dimG
      call MtimesT(tmpmat2,tmpmat,a,NMAT)
        do b=1,dimG
          call MTN(diffdiff_Uf(:,:,b,a),tmpmat2,YMAT,b,NMAT)
        enddo
      enddo
    endif
  endif
  
!!!!!!!!!!!!!!
! when l_1 != l_2
else
! set l1 and l2 in the order of l_1 and l_2 in the face f
  info=0
  do nl=1,links_in_f(f)%num_
    if( links_in_f(f)%link_labels_(nl) == l_1 ) then 
      if(info==0) then
        nl1=nl
        l1=l_1
        info=info+1
      elseif(info==1) then
        nl2=nl
        l2=l_1
        info=info+1
      endif
    elseif( links_in_f(f)%link_labels_(nl) == l_2 ) then 
      if(info==0) then
        nl1=nl
        l1=l_2
        info=info+1
      elseif(info==1) then
        nl2=nl
        l2=l_2
        info=info+1
      endif
    endif
    if(info==2) exit
  enddo
  if(info*(info-1)==0) then 
    write(*,*) "l_1 and/or l_2 are not included in the face f"
    stop
  endif
! now 
!   l1: former link in f 
!   nl1: position of l1 in the face
!
!   l2: later link in f
!   nl2: position of l2 in the face

!! X=U_l1...U_{nl1-1}
if( nl1==1 ) then 
  call make_unit_matrix(XMAT)
else
  if(links_in_f(f)%link_dirs_(1)==1) then
    XMAT=UMAT(:,:,links_in_f(f)%link_labels_(1))
  elseif(links_in_f(f)%link_dirs_(1)==-1) then
    do i=1,NMAT
    do j=1,NMAT
      XMAT(i,j)=dconjg(UMAT(j,i,links_in_f(f)%link_labels_(1)))
    enddo
    enddo
  endif
endif
!
if(nl1 >= 3) then
  do ll=2,nl1-1
    tmpmat=XMAT
    if(links_in_f(f)%link_dirs_(ll)==1) then 
      char1='N'
    elseif(links_in_f(f)%link_dirs_(ll)==-1) then
      char1='C'
    endif
    !! Uf=tmpmat*Uf or tmpmat*(Uf)^\dagger
    call ZGEMM('N',char1,NMAT,NMAT,NMAT,(1d0,0d0), &
      tmpmat,NMAT, &
      UMAT(:,:,links_in_f(f)%link_labels_(ll)),NMAT, &
      (0d0,0d0), XMAT, NMAT)
  enddo
endif
!! dir(U_{l1}) = -1 --> X -> -i X.U_{l1}^{-1}
if( links_in_f(f)%link_dirs_(nl1)==-1) then
  tmpmat=XMAT
  call ZGEMM('N','C',NMAT,NMAT,NMAT,(0d0,-1d0), &
    tmpmat,NMAT, &
    UMAT(:,:,l1),NMAT, &
    (0d0,0d0), XMAT, NMAT)
endif

!! Y=U_{nl1+1}...U_{nl2-1}
if( nl2 == nl1+1 ) then 
  call make_unit_matrix(YMAT)
else
  if(links_in_f(f)%link_dirs_(nl1+1)==1) then
    YMAT=UMAT(:,:,links_in_f(f)%link_labels_(nl1+1))
  elseif(links_in_f(f)%link_dirs_(nl1+1)==-1) then
    do i=1,NMAT
    do j=1,NMAT
      YMAT(i,j)=dconjg(UMAT(j,i,links_in_f(f)%link_labels_(nl1+1)))
    enddo
    enddo
  endif
endif
!
if(nl2-1 >= nl1+2) then
do ll=nl1+2,nl2-1
  tmpmat=YMAT
  if(links_in_f(f)%link_dirs_(ll)==1) then 
    char1='N'
  elseif(links_in_f(f)%link_dirs_(ll)==-1) then
    char1='C'
  endif
  !! Uf=tmpmat*Uf or tmpmat*(Uf)^\dagger
  call ZGEMM('N',char1,NMAT,NMAT,NMAT,(1d0,0d0), &
    tmpmat,NMAT, &
    UMAT(:,:,links_in_f(f)%link_labels_(ll)),NMAT, &
    (0d0,0d0), YMAT, NMAT)
  enddo
endif
!! dir(U_{l1}) = 1 --> Y -> i U_{l1}.Y
if( links_in_f(f)%link_dirs_(nl1)==1) then
  tmpmat=YMAT
  call ZGEMM('N','N',NMAT,NMAT,NMAT,(0d0,1d0), &
    UMAT(:,:,l1),NMAT, &
    tmpmat,NMAT, &
    (0d0,0d0), YMAT, NMAT)
endif
!! dir(U_{l2}) = -1 --> Y -> -i Y.U_{l2}^{-1}
if( links_in_f(f)%link_dirs_(nl2)==-1) then
  tmpmat=YMAT
  call ZGEMM('N','C',NMAT,NMAT,NMAT,(0d0,-1d0), &
    tmpmat,NMAT, &
    UMAT(:,:,l2),NMAT, &
    (0d0,0d0), YMAT, NMAT)
endif



!! Z=U_{nl2+1}...U_{num_links(f))
if( nl2 == links_in_f(f)%num_) then 
  call make_unit_matrix(ZMAT)
else 
  if(links_in_f(f)%link_dirs_(nl2+1)==1) then
    ZMAT=UMAT(:,:,links_in_f(f)%link_labels_(nl2+1))
  elseif(links_in_f(f)%link_dirs_(nl2+1)==-1) then
    do i=1,NMAT
    do j=1,NMAT
      ZMAT(i,j)=dconjg(UMAT(j,i,links_in_f(f)%link_labels_(nl2+1)))
    enddo
    enddo
  endif
endif
!
if(links_in_f(f)%num_ >= nl2+2) then
  do ll=nl2+2, links_in_f(f)%num_
    tmpmat=ZMAT
    if(links_in_f(f)%link_dirs_(ll)==1) then 
      char1='N'
    elseif(links_in_f(f)%link_dirs_(ll)==-1) then
      char1='C'
    endif
    !! Uf=tmpmat*Uf or tmpmat*(Uf)^\dagger
    call ZGEMM('N',char1,NMAT,NMAT,NMAT,(1d0,0d0), &
      tmpmat,NMAT, &
      UMAT(:,:,links_in_f(f)%link_labels_(ll)),NMAT, &
      (0d0,0d0), ZMAT, NMAT)
  enddo
endif
!! dir(U_{l2}) = 1 --> Z -> i U_{l1}.Z
if( links_in_f(f)%link_dirs_(nl2)==1) then
  tmpmat=ZMAT
  call ZGEMM('N','N',NMAT,NMAT,NMAT,(0d0,1d0), &
    tmpmat,NMAT, &
    UMAT(:,:,l2),NMAT, &
    (0d0,0d0), ZMAT, NMAT)
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   diffdiff_Uf(:,;,b,a) = d/dA(b,l_2) d/dA(a,l_1) Uf
! if l_1 is the first link, Ta is between X and Y 
  if( l1 == l_1 ) then 
    do a=1,dimG
      call MTN(tmpmat,XMAT,YMAT,a,NMAT)
      do b=1,dimG
        call MTN(diffdiff_Uf(:,:,b,a),tmpmat,ZMAT,b,NMAT)
      enddo
    enddo
! if l_2 is the first link, Ta is between Y and Z
  else 
    do a=1,dimG
      call MTN(tmpmat,YMAT,ZMAT,a,NMAT)
      do b=1,dimG
        call MTN(diffdiff_Uf(:,:,b,a),XMAT,tmpmat,b,NMAT)
      enddo
    enddo
  endif
endif

end subroutine tmp_calc_diffdiff_Uf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate Amat used in fermionic action and force
subroutine calc_Amat(Amat,f,l_label,k,Uf,UMAT)
use matrix_functions, only : matrix_power,matrix_product
implicit none

complex(kind(0d0)), intent(out) :: Amat(1:NMAT,1:NMAT)
integer, intent(in) :: f,l_label,k
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)

complex(kind(0d0)):: tmpmat1(1:NMAT,1:NMAT)
complex(kind(0d0)):: tmpmat2(1:NMAT,1:NMAT)
!character :: C1
integer :: l,i,label

l=links_in_f(f)%link_labels_(l_label)

! tmpmat=U_f^{k-1}
if (k.ne.1) then
  call matrix_power(tmpmat1,Uf,k-1)
endif

! tmpmat2=U_f^{k-1} U_{l1}^{e1}...U_{l-1}^{el-1}
if ( links_in_f(f)%link_dirs_(l_label) == 1 ) then
  label = l_label - 1
else
  label = l_label
endif 
call calc_prodUl_from_n1_to_n2_in_Uf(tmpmat2,f,1,label,UMAT)

if (k.ne.1) then
  call matrix_product(Amat,tmpmat1,tmpmat2)
else
  Amat=tmpmat2
endif


end subroutine calc_Amat


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate Bmat used in fermionic action and force
subroutine calc_Bmat(Bmat,f,l_label,k,Uf,UMAT)
use matrix_functions, only : matrix_power,matrix_product
implicit none

complex(kind(0d0)), intent(out) :: Bmat(1:NMAT,1:NMAT)
integer, intent(in) :: f,l_label,k
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)

complex(kind(0d0)):: tmpmat1(1:NMAT,1:NMAT)
complex(kind(0d0)):: tmpmat2(1:NMAT,1:NMAT)
!character :: C1
integer :: l,i,label

l=links_in_f(f)%link_labels_(l_label)

! tmpmat1=U_f^{m-k-1}
if( m_omega-k >= 0 ) then 
  call matrix_power(tmpmat1,Uf,m_omega-k)
endif

! tmpmat=U_f^{k-1} U_{l1}^{e1}...U_{l-1}^{el-1}
if ( links_in_f(f)%link_dirs_(l_label) == 1 ) then
  label = l_label
else
  label = l_label+1
endif 
call calc_prodUl_from_n1_to_n2_in_Uf(tmpmat2,f,label,links_in_f(f)%num_ ,UMAT)

if( m_omega-k >= 0 ) then 
  call matrix_product(Bmat, tmpmat2, tmpmat1)
else
  Bmat=tmpmat2
endif 

end subroutine calc_Bmat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calclate U^m - U^{-m}
!!
!! CAUTION
!! Now sinU = Uf^m - Uf^{-1}
!! THERE IS NO -i FACTOR
subroutine calc_sinU_and_cosUinv(sinU,cosUinv,Ufm)
use matrix_functions, only : matrix_inverse,matrix_power,matrix_dagger_power
implicit none

complex(kind(0d0)), intent(out) :: sinU(NMAT,NMAT)
complex(kind(0d0)), intent(out) :: cosUinv(NMAT,NMAT)
complex(kind(0d0)), intent(in) :: Ufm(NMAT,NMAT)
complex(kind(0d0)):: tmpmat(NMAT,NMAT)
integer :: i,j

do i=1,NMAT
  do j=1,i
    if (i==j) then
      sinU(i,i)=Ufm(i,i) - conjg(Ufm(i,i))
      cosUinv(i,i)=Ufm(i,i) + conjg(Ufm(i,i))
    else
      sinU(i,j)=Ufm(i,j) - conjg(Ufm(j,i))
      sinU(j,i)=Ufm(j,i) - conjg(Ufm(i,j))
      cosUinv(i,j)=Ufm(i,j) + conjg(Ufm(j,i))
      cosUinv(j,i)=Ufm(j,i) + conjg(Ufm(i,j))
    endif 
  enddo
enddo
call matrix_inverse(cosUinv)

end subroutine calc_sinU_and_cosUinv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_dSinUdA(dSinUdA, Uf,UMAT,f,ll_label)
use matrix_functions, only : matrix_product
!use Dirac_operator, only : calc_Amat, calc_Bmat
implicit none

! for given f and ll
! d cosUinv(i,j,f) / dA_{ii,jj,ll)
complex(kind(0d0)), intent(out) :: dSinUdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
!complex(kind(0d0)), intent(in) :: sinU(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: UMat(1:NMAT,1:NMAT,1:num_links)
integer, intent(in) :: f,ll_label

complex(kind(0d0)) :: Amat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Bmat(1:NMAT,1:NMAT)

integer :: i,j,ii,jj

dSinUdA=(0d0,0d0)
call calc_Amat(Amat,f,ll_label,1,Uf,Umat)
call calc_Bmat(Bmat,f,ll_label,1,Uf,Umat)

do jj=1,NMAT
  do ii=1,NMAT
    do j=1,NMAT
      do i=1,NMAT
        dSinUdA(i,j,ii,jj)=dSinUdA(i,j,ii,jj) &
          + Amat(i,jj)*Bmat(ii,j) &
          + conjg(Bmat(jj,i))*conjg(Amat(j,ii))
      enddo
    enddo
  enddo
enddo

dSinUdA=dSinUdA * cmplx(dble(links_in_f(f)%link_dirs_(ll_label)))*(0d0,1d0)

end subroutine calc_dSinUdA



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to return product of d/dA D . vec
subroutine calc_dCosUinvdA_dSinUdA(dCosUinvdA,dSinUdA,&
    cosUinv,Uf,UMAT,f,ll_label)
use matrix_functions, only : matrix_product
!use Dirac_operator, only : calc_Amat, calc_Bmat
implicit none

! for given f and ll
! d cosUinv(i,j,f) / dA_{ii,jj,ll)
complex(kind(0d0)), intent(out) :: dCosUinvdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
complex(kind(0d0)), intent(out) :: dSinUdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
!complex(kind(0d0)), intent(in) :: sinU(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: cosUinv(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: UMat(1:NMAT,1:NMAT,1:num_links)
integer, intent(in) :: f,ll_label

complex(kind(0d0)) :: Amat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Bmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: UinvA(1:NMAT,1:NMAT), BUinv(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: AdagUinv(1:NMAT,1:NMAT), UinvBdag(1:NMAT,1:NMAT)


integer :: k,i,j,ii,jj
complex(kind(0d0)) :: factor

dCosUinvdA=(0d0,0d0)
dSinUdA=(0d0,0d0)
do k=1,m_omega
  call calc_Amat(Amat,f,ll_label,k,Uf,Umat)
  call calc_Bmat(Bmat,f,ll_label,k,Uf,Umat)

  call matrix_product(UinvA,cosUinv,Amat)
  call matrix_product(BUinv,Bmat,cosUinv)
  !call matrix_product(AdagUinv,Amat,cosUinv,'C','N')
  !call matrix_product(UinvBdag,cosUinv,Bmat,'N','C')

  do jj=1,NMAT
    do ii=1,NMAT
      do j=1,NMAT
        do i=1,NMAT
          dCosUinvdA(i,j,ii,jj)=dCosUinvdA(i,j,ii,jj) &
            + UinvA(i,jj)*BUinv(ii,j) &
            !- UinvBdag(i,jj)*AdagUinv(ii,j)
            - conjg( BUinv(jj,i) ) * conjg( UinvA(j,ii) )
          dSinUdA(i,j,ii,jj)=dSinUdA(i,j,ii,jj) &
            + Amat(i,jj)*Bmat(ii,j) &
            + conjg(Bmat(jj,i))*conjg(Amat(j,ii))
        enddo
      enddo
    enddo
  enddo
enddo

factor=cmplx(dble(links_in_f(f)%link_dirs_(ll_label)))*(0d0,1d0)
dCosUinvdA=dCosUinvdA*(-factor)
dSinUdA=dSinUdA * factor

end subroutine calc_dCosUinvdA_dSinUdA


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to return product of d/dA D . vec
subroutine calc_dABmatdA(dAmatdA,dBmatdA,Uf,UMAT,ll_label,f,l_label,k)
use matrix_functions, only : matrix_product, matrix_3_product, matrix_power,make_unit_matrix
implicit none

! for given f,l,k and ll
! d Amat(i,j,f,l,k) / dA_{ii,jj,ll)
complex(kind(0d0)), intent(out) :: dAmatdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
complex(kind(0d0)), intent(out) :: dBmatdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: UMat(1:NMAT,1:NMAT,1:num_links)
integer, intent(in) :: f,ll_label,l_label,k

complex(kind(0d0)) :: UU_initial_to_ll(1:NMAT,1:NMAT) ! 1..ll
complex(kind(0d0)) :: UU_initial_to_l(1:NMAT,1:NMAT) ! 1..l
complex(kind(0d0)) :: UU_ll_to_n(1:NMAT,1:NMAT) ! ll..n
complex(kind(0d0)) :: UU_l_to_n(1:NMAT,1:NMAT) ! l..n
complex(kind(0d0)) :: UU_ll_to_l(1:NMAT,1:NMAT) ! ll..l
complex(kind(0d0)) :: UU_l_to_ll(1:NMAT,1:NMAT) ! l..ll

complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat3(1:NMAT,1:NMAT)
complex(kind(0d0)) :: epsilon_r
character :: C1
integer :: n,label1,label2
integer :: i,j,ii,jj,kk

dAmatdA=(0d0,0d0)
dBmatdA=(0d0,0d0)
epsilon_r=cmplx(dble( links_in_f(f)%link_dirs_(ll_label) ))*(0d0,1d0)

n=links_in_f(f)%num_
!!!!!!!!!!!!!!
if ( links_in_f(f)%link_dirs_(ll_label) == 1 ) then
  label2 = ll_label - 1
else
  label2 = ll_label
endif 
call calc_prodUl_from_n1_to_n2_in_Uf(UU_initial_to_ll,f,1,label2,UMAT)
!!!!!!!!!!!!!!
if ( links_in_f(f)%link_dirs_(l_label) == 1 ) then
  label2 = l_label-1
else
  label2 = l_label
endif 
call calc_prodUl_from_n1_to_n2_in_Uf(UU_initial_to_l,f,1,label2,UMAT)
!!!!!!!!!!!!!!
if ( links_in_f(f)%link_dirs_(ll_label) == 1 ) then
  label1 = ll_label 
else
  label1 = ll_label+1
endif 
call calc_prodUl_from_n1_to_n2_in_Uf(UU_ll_to_n,f,label1,n,UMAT)
!!!!!!!!!!!!!!
if ( links_in_f(f)%link_dirs_(l_label) == 1 ) then
  label1 = l_label
else
  label1 = l_label+1
endif 
call calc_prodUl_from_n1_to_n2_in_Uf(UU_l_to_n,f,label1,n,UMAT)
!!!!!!!!!!!!!!
if ( links_in_f(f)%link_dirs_(ll_label) == 1 ) then
  label1 = ll_label 
else
  label1 = ll_label+1
endif 
if ( links_in_f(f)%link_dirs_(l_label) == 1 ) then
  label2 = l_label-1
else
  label2 = l_label
endif 
call calc_prodUl_from_n1_to_n2_in_Uf(UU_ll_to_l,f,label1,label2,UMAT)
!!!!!!!!!!!!!!
if ( links_in_f(f)%link_dirs_(l_label) == 1 ) then
  label1 = l_label 
else
  label1 = l_label+1
endif 
if ( links_in_f(f)%link_dirs_(ll_label) == 1 ) then
  label2 = ll_label-1
else
  label2 = ll_label
endif 
call calc_prodUl_from_n1_to_n2_in_Uf(UU_l_to_ll,f,label1,label2,UMAT)

        
!!!!!!!!!!!!!! dAmatdA !!!!!!!!!!!!!!!!!
if ( m_omega > 1 ) then
  do kk=1,k-1
    call matrix_power(tmpmat2,Uf,kk-1)
    call matrix_product(tmpmat1,tmpmat2,UU_initial_to_ll)

    call matrix_power(tmpmat3,Uf,k-kk-1)
    call matrix_3_product(tmpmat2,UU_ll_to_n,tmpmat3,UU_initial_to_l)

    do jj=1,NMAT
      do ii=1,NMAT
        do j=1,NMAT
          do i=1,NMAT
            dAmatdA(i,j,ii,jj)=dAmatdA(i,j,ii,jj) &
              + tmpmat1(i,jj) * tmpmat2(ii,j) * epsilon_r
          enddo
        enddo
      enddo
    enddo

  enddo
endif

if ( ll_label < l_label ) then
  call make_unit_matrix(tmpmat2)
  if( k-1 >= 0) call matrix_power(tmpmat2,Uf,k-1)
  call matrix_product(tmpmat1,tmpmat2,UU_initial_to_ll)

  do jj=1,NMAT
    do ii=1,NMAT
      do j=1,NMAT
        do i=1,NMAT
          dAmatdA(i,j,ii,jj)=dAmatdA(i,j,ii,jj) &
            + tmpmat1(i,jj) * UU_ll_to_l(ii,j) * epsilon_r
        enddo
      enddo
    enddo
  enddo
endif

if (ll_label == l_label) then
  if ( links_in_f(f)%link_dirs_(ll_label) == -1 ) then
    call make_unit_matrix(tmpmat2)
    if( k-1 >= 0 ) call matrix_power(tmpmat2,Uf,k-1)
    call matrix_product(tmpmat1, tmpmat2,UU_initial_to_l)
    do jj=1,NMAT
      do ii=1,NMAT
        j=ii
        do i=1,NMAT
          dAmatdA(i,j,ii,jj)=dAmatdA(i,j,ii,jj) &
            - (0d0,1d0)* tmpmat1(i,jj)
        enddo
      enddo
    enddo
  endif
endif

!!!!!!!!!!! dBmatdA !!!!!!!!!!!!!!!
if (m_omega > 1) then
  do kk=1,m_omega-k

    call make_unit_matrix(tmpmat2)
    if( kk-1 > 0) call matrix_power(tmpmat2, Uf, kk-1)
    call matrix_3_product(tmpmat1,UU_l_to_n,tmpmat2,UU_initial_to_ll)

    call make_unit_matrix(tmpmat3)
    if( m_omega-k-kk > 0) call matrix_power(tmpmat3,Uf,m_omega-k-kk)
    call matrix_product(tmpmat2,UU_ll_to_n,tmpmat3)
!
    do jj=1,NMAT
      do ii=1,NMAT
        do j=1,NMAT
          do i=1,NMAT
            dBmatdA(i,j,ii,jj)=dBmatdA(i,j,ii,jj) &
              + tmpmat1(i,jj) * tmpmat2(ii,j) * epsilon_r
          enddo
        enddo
      enddo
    enddo

  enddo
endif

if ( ll_label > l_label) then
  call make_unit_matrix(tmpmat1)
  if(m_omega-k > 0) call matrix_power(tmpmat1,Uf,m_omega-k)
  call matrix_product(tmpmat2,UU_ll_to_n,tmpmat1)

  do jj=1,NMAT
    do ii=1,NMAT
      do j=1,NMAT
        do i=1,NMAT
          dBmatdA(i,j,ii,jj)=dBmatdA(i,j,ii,jj) &
            + UU_l_to_ll(i,jj) * tmpmat2(ii,j) * epsilon_r
        enddo
      enddo
    enddo
  enddo
endif

if ( ll_label == l_label) then
  if ( links_in_f(f)%link_dirs_(ll_label) == 1) then
    call make_unit_matrix(tmpmat2)
    if( m_omega-k >= 0) call matrix_power(tmpmat2,Uf,m_omega-k)
    call matrix_product(tmpmat1,UU_ll_to_n,tmpmat2)
    do jj=1,NMAT
      do ii=1,NMAT
        i=jj
        do j=1,NMAT
          dBmatdA(i,j,ii,jj)=dBmatdA(i,j,ii,jj) &
            + (0d0,1d0) * tmpmat1(ii,j) 
        enddo
      enddo
    enddo
  endif
endif

end subroutine calc_dABmatdA



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
subroutine vec_to_mat(eta,lambda,chi,vec)
use SUN_generators, only : make_traceless_matrix_from_modes
implicit none

complex(kind(0d0)), intent(out) :: eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(out) :: lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(out) :: chi(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)), intent(in) :: vec(1:sizeD)

complex(kind(0d0)) :: ele(1:dimG)
integer :: s,l,f,a

do s=1,num_sites
  do a=1,dimG
    ele(a)=vec(site_index(a,s))
  enddo
  call make_traceless_matrix_from_modes(eta(:,:,s),NMAT,ele)
enddo

do l=1,num_links
  do a=1,dimG
    ele(a)=vec(link_index(a,l))
  enddo
  call make_traceless_matrix_from_modes(lambda(:,:,l),NMAT,ele)
enddo

do f=1,num_faces
  do a=1,dimG
    ele(a)=vec(face_index(a,f))
  enddo
  call make_traceless_matrix_from_modes(chi(:,:,f),NMAT,ele)
enddo

end subroutine vec_to_mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
subroutine mat_to_vec(vec,eta,lambda,chi)
use SUN_generators, only : trace_MTa
implicit none

complex(kind(0d0)), intent(in) :: eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(in) :: lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: chi(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)), intent(out) :: vec(1:sizeD)

complex(kind(0d0)) :: trace
integer :: s,l,f,a

do s=1,num_sites
  do a=1,dimG
    call trace_MTa(trace,eta(:,:,s),a,NMAT)
    vec(site_index(a,s))=vec(site_index(a,s))+trace
  enddo
enddo
do l=1,num_links
  do a=1,dimG
    call trace_MTa(trace,lambda(:,:,l),a,NMAT)
    vec(link_index(a,l))=vec(link_index(a,l))+trace
  enddo
enddo
do f=1,num_faces
  do a=1,dimG
    call trace_MTa(trace,chi(:,:,f),a,NMAT)
    vec(face_index(a,f))=vec(face_index(a,f))+trace
  enddo
enddo
end subroutine mat_to_vec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! check distance from 1 of Uf
subroutine  check_distance(info,ratio,UMAT)
use global_parameters
use matrix_functions, only : matrix_norm, make_unit_matrix
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
integer, intent(out) :: info
double precision, intent(out) :: ratio
double precision :: distance,tmp
integer l,i,j,f
!double precision norm
complex(kind(0d0)) UNITM(1:NMAT,1:NMAT)
!complex(kind(0d0)) tmp(1:NMAT,1:NMAT)
complex(kind(0d0)) Uf(1:NMAT,1:NMAT)
!double precision dist(1:NMAT-1)

!write(*,*) "===== check distance from 1 ==========="
!write(*,*) "theoretical dist. to the nearest center=",dsin(PI/dble(NMAT))*2d0

info=0
ratio=0d0
call make_unit_matrix(UNITM)
do f=1,num_faces
  call Make_face_variable(Uf,f,UMAT)
  call matrix_norm(distance,UNITM-Uf)
  tmp=distance/maximal_dist
  if( tmp > ratio ) then 
    ratio=tmp
  endif
  if ( distance > maximal_dist ) then
    info=info+1
  endif
enddo

end subroutine check_distance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! function to obtain the argument of a complex number z
real(8) function arg(z)
implicit none

complex(kind(0d0)), intent(in) :: z
complex(kind(0d0)) :: phase

phase = z / cmplx( abs(z) )
arg=atan2( dble((0d0,-1d0)*phase), dble(phase) )

end function arg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End other subroutines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module simulation


