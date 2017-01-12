module output
use global_parameters
implicit none

integer, parameter :: num_obs=2
character(10) :: obs_name(1:num_obs)
data obs_name/ "Sb","TrX2" /

double precision :: OBS(1:num_obs) ! 1) bosonic action SB

contains
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
use observables
use SUN_generators, only : trace_MTa
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: Phi(1:dimG,1:num_sites)
double precision, intent(in) :: delta_Ham, ratio
integer, intent(in) :: ite, accept, total_ite, CGite
integer :: output,i,s,a

!complex(kind(0d0)) :: min_eigen,max_eigen
!call calc_smallset_and_largest_eigenvalues_of_D(min_eigen,max_eigen,UMAT,Phi)

!do s=1,num_sites
  !do a=1,dimG
    !call trace_MTa(Phi(a,s),PhiMat(:,:,s),a,NMAT)
  !enddo
!enddo

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
write(MED_CONF_FILE) config_step
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
use observables, only : calc_eigenvalues_Dirac
implicit none

complex(kind(0d0)), intent(out) :: min_eigen, max_eigen
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)

complex(kind(0d0)) :: eigenvalues(1:sizeD)

call calc_eigenvalues_Dirac(eigenvalues,UMAT,PhiMat)

min_eigen=eigenvalues(1)
max_eigen=eigenvalues(sizeD)

end subroutine calc_smallset_and_largest_eigenvalues_of_D
end module output
