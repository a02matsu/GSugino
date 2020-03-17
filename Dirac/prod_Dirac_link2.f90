subroutine prod_Dirac_link2(DF_lambda,PhiMat,Umat,lambda_mat)
implicit none

complex(kind(0d0)), intent(out) :: DF_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: lambda_mat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
integer :: l


do l=1,num_links
  call Hermitian_conjugate(tmpmat1,PhiMat(:,:,link_org(l)))
  !tmpmat1=bPhiMat(:,:,link_org(l))

  call matrix_3_product(tmpmat1,Umat(:,:,l),PhiMat(:,:,link_tip(l)),Umat(:,:,l),&
    'N','C','C',dconjg(U1Rfactor_link(l)**2d0 ),'ADD')

  call matrix_commutator(tmpmat2,tmpmat1,lambda_mat(:,:,l))
  DF_lambda(:,:,l)=DF_lambda(:,:,l) &
      + dcmplx( alpha_l(l) * overall_factor ) * tmpmat2
enddo


end subroutine prod_Dirac_link2
