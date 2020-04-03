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
!integer, parameter :: num_obs=4
integer, parameter :: num_obs=2
character(10) :: obs_name(1:num_obs)
data obs_name/ "Sb","TrX2" /
!data obs_name/ "Sb","TrX2","Re(WT_Tr)","Im(WT_Tr)" /

double precision :: OBS(1:num_obs) ! 1) bosonic action SB


contains
#include "MonteCarloSteps.f90"
#include "hamiltonian.f90"
#include "forces.f90"
#include "output.f90"
#include "check_QS.f90"

#include "Observables/div_rot.f90"
#include "Observables/trphi2.f90"
#include "Observables/trf2.f90"
#include "Observables/compensators.f90"
#include "Observables/bosonic_action.f90"
#include "Observables/U1V_current.f90"
#include "Observables/fermion_action_site.f90"
#include "Observables/fermion_action_link1.f90"
#include "Observables/fermion_action_link2.f90"
#include "Observables/fermion_action_face1.f90"
#include "Observables/fermion_action_face2.f90"
#include "Observables/fermionic_face_lagrangian.f90"
#include "Observables/phichi.f90"
#include "Observables/checkFF.f90"
#include "Observables/fermionic_operators.f90"
#include "Observables/trivialWT.f90"
#include "Observables/eigenvalues_of_Dirac.f90"
#include "Observables/exact_U1R.f90"
#include "Observables/mass_reweighting.f90"
#include "Observables/WT_mass_contribution_site.f90"
#include "Observables/WT_mass_contribution_link.f90"
#include "Observables/WT_mass_contribution_face.f90"
#include "Observables/make_Xi.f90"
#include "Observables/Qfermion.f90"
#include "Observables/QS_site.f90"
#include "Observables/QS_link.f90"
#include "Observables/QS_face.f90"
#include "Observables/Qexact_fermion_action_link1.f90"
#include "Observables/Qexact_fermion_action_link2.f90"
#include "Observables/Qexact_fermion_action_face2.f90"
#include "Observables/Qexact_fermionic_face_lagrangian.f90"
#include "Observables/Xisite_Dinv.f90"
#include "Observables/Xilink_Dinv.f90"
#include "Observables/Xiface_Dinv.f90"
#include "Observables/QS_3fermion_link.f90"
!#include "Observables/Xilink_QS.f90"
!#include "Observables/Xiface_QS.f90"
!#include "Observables/test_fermion_action_site.f90"
!#include "Observables/test_fermion_action_link1.f90"




end module simulation
