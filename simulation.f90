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
#include "MonteCarloSteps.f90"

end module simulation
