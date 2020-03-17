!!!!!!!!!!!!!!!!!!!!!!
!! Qeta must be calculated throught make_Qfermion
subroutine calc_QS_site(QS_eta,PhiMat,Qeta)
use matrix_functions, only : Hermitian_Conjugate
use Dirac_operator, only : prod_Dirac_site
implicit none

complex(kind(0d0)), intent(out) :: QS_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(in) :: Qeta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

integer :: s
complex(kind(0d0)) :: Bforce(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: DQeta(1:NMAT,1:NMAT,1:num_sites)

Bforce=(0d0,0d0)
call Make_bosonic_force_Phi_site(Bforce,PhiMat)
DQeta=(0d0,0d0)
call prod_Dirac_site(DQeta,PhiMat,Qeta)

!QS_eta=(0d0,0d0)
do s=1,num_sites
  !! boson
  call Hermitian_Conjugate(QS_eta(:,:,s),Bforce(:,:,s))
  QS_eta(:,:,s)=QS_eta(:,:,s)*U1Rfactor_site(s)
  !! fermion
  QS_eta(:,:,s)=QS_eta(:,:,s)-DQeta(:,:,s)
enddo

end subroutine calc_QS_site







