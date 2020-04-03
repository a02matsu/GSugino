!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to compute QS in the form that
!!   QS = Tr(eta. QS_eta) + Tr(lambda. QS_lambda) + Tr(chi . QS_chi)
subroutine calc_QS(QS_eta, QS_lambda, QS_chi, &
    PhiMat, Umat)
implicit none

complex(kind(0d0)), intent(out) :: QS_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(out) :: QS_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(out) :: QS_chi(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) :: QSS_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: QSL_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: QSL_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: QSF_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: QSF_chi(1:NMAT,1:NMAT,1:num_faces)


call make_Qfermion(Qeta,Qlambda,Qchi,Omega,Umat,PhiMat)
 !! QS = \eta.QS_eta + \lambda.QS_lambda + \chi.QS_chi
call calc_QS_site(QSS_eta,PhiMat,Qeta)
call calc_QS_link(QSL_eta,QSL_lambda,PhiMat,Umat,Qeta,Qlambda)
call calc_QS_face(QSF_lambda,QSF_chi,PhiMat,Umat,Qlambda,Qchi)
QS_eta=QSS_eta+QSL_eta
QS_lambda=QSL_lambda+QSF_lambda
QS_chi=QSF_chi





