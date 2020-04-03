!!!!!!!!!!!!!!!!!!!!!!
!! Qeta must be calculated throught make_Qfermion
subroutine calc_QS_link(QSL_eta,QSL_lambda,PhiMat,Umat,Qeta,Qlambda)
use matrix_functions, only : Hermitian_Conjugate
use Dirac_operator, only : prod_Dirac_link1, prod_Dirac_link2
implicit none

complex(kind(0d0)), intent(out) :: QSL_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(out) :: QSL_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: Qeta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: Qlambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: UMat(1:NMAT,1:NMAT,1:num_necessary_links)

integer :: l,s
complex(kind(0d0)) :: Bforce_phi(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: Bforce_A(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: DQeta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: DQlambda(1:NMAT,1:NMAT,1:num_links)

QSL_eta=(0d0,0d0)
QSL_lambda=(0d0,0d0)

Bforce_phi=(0d0,0d0)
Bforce_A=(0d0,0d0)
call Make_bosonic_force_Phi_link(Bforce_phi,Umat,PhiMat)
call Make_bosonic_force_A_link(Bforce_A,Umat,PhiMat)

DQeta=(0d0,0d0)
DQlambda=(0d0,0d0)
call prod_Dirac_link1(DQeta,DQlambda,PhiMat,Umat,Qeta,Qlambda)
call prod_Dirac_link2(DQlambda,PhiMat,Umat,Qlambda)

QSL_eta=(0d0,0d0)
QSL_lambda=(0d0,0d0)
do s=1,num_sites
  !! boson
  call Hermitian_Conjugate(QSL_eta(:,:,s),Bforce_phi(:,:,s))
  QSL_eta(:,:,s)=QSL_eta(:,:,s)*U1Rfactor_site(s)
  !! fermion
  QSL_eta(:,:,s)=QSL_eta(:,:,s)-DQeta(:,:,s)
enddo
!!
do l=1,num_links
  !! boson
  QSL_lambda(:,:,l)=Bforce_A(:,:,l)*U1Rfactor_site(link_org(l))
  !! fermion
  QSL_lambda(:,:,l)=QSL_lambda(:,:,l)-DQlambda(:,:,l)
enddo

end subroutine calc_QS_link







