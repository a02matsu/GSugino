!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ver.01: correct the normalization of the PCSC relation 
!! ver.04: bug fix for mass part of PCSC
!! ver.05: include WT id. in naive quench
!! ver.06: added compensator for SU(2) (triple cover version)
!module observables
!use global_parameters
!use global_subroutines
!implicit none
!
!contains

#include "Observables/div_rot.f90"
#include "Observables/trphi2.f90"
#include "Observables/trf2.f90"
#include "Observables/compensators.f90"
#include "Observables/bosonic_action.f90"
#include "Observables/U1V_current.f90"
#include "Observables/fermion_action_link2.f90"
#include "Observables/phichi.f90"
#include "Observables/checkFF.f90"
#include "Observables/fermionic_operators.f90"
#include "Observables/trivialWT.f90"
#include "Observables/eigenvalues_of_Dirac.f90"
#include "Observables/exact_U1R.f90"
#include "Observables/mass_reweighting.f90"
#include "Observables/fermionic_face_lagrangian.f90"
#include "Observables/siteWT.f90"
#include "Observables/make_Xi.f90"


