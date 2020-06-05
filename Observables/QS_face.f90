!!!!!!!!!!!!!!!!!!!!!!
!! Qeta must be calculated throught make_Qfermion
subroutine calc_QS_face(QSF_lambda,QSF_chi,PhiMat,Umat,Qlambda,Qchi)
use Dirac_operator, only : prod_Dirac_face1, Dirac_Omega_m0, Dirac_Omega_adm, Dirac_Omega
implicit none

complex(kind(0d0)), intent(out) :: QSF_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(out) :: QSF_chi(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: Qlambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: Qchi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: UMat(1:NMAT,1:NMAT,1:num_necessary_links)

integer :: l,f
complex(kind(0d0)) :: Bforce_A(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: DQlambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: DQchi(1:NMAT,1:NMAT,1:num_faces)


QSF_lambda=(0d0,0d0)
QSF_chi=(0d0,0d0)

call Make_bosonic_force_A_face(Bforce_A,Umat)

DQlambda=(0d0,0d0)
DQchi=(0d0,0d0)
!! from Sf_face1
call prod_Dirac_face1(DQchi,PhiMat,Qchi)
!! from Sf_face2
if( m_omega == 0 ) then
  call Dirac_Omega_m0(DQchi,DQlambda,Umat,Qlambda,Qchi)
elseif( m_omega == -1 ) then
  call Dirac_Omega_adm(DQchi,DQlambda,Umat,Qlambda,Qchi)
else
  call Dirac_Omega(DQchi,DQlambda,Umat,Qlambda,Qchi)
endif

QSF_lambda=(0d0,0d0)
do l=1,num_links
  QSF_lambda(:,:,l)=&
    Bforce_A(:,:,l)*U1Rfactor_site(link_org(l)) & !! fermion
    -DQlambda(:,:,l) !! boson
enddo
!!
QSF_chi=(0d0,0d0)
do f=1,num_faces
  !! boson
  ! no terms
  !! fermion
  QSF_chi(:,:,f)=-DQchi(:,:,f)
enddo
!!

end subroutine calc_QS_face
