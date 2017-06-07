!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! module for molecular dynamics
module simulation
use global_parameters
use global_subroutines
!use hamiltonian
!use forces
use mt95
#ifdef PARALLEL
use parallel
#endif
implicit none

! for output
integer, parameter :: num_obs=2
character(10) :: obs_name(1:num_obs)
data obs_name/ "Sb","TrX2" /

double precision :: OBS(1:num_obs) ! 1) bosonic action SB


contains
#include "MonteCarloSteps.f90"
#include "hamiltonian.f90"
#include "forces.f90"
#include "observables.f90"
#include "output.f90"

end module simulation
