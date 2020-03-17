subroutine prod_Dirac_link1(DF_eta,DF_lambda,PhiMat,Umat,eta_mat,lambda_mat)
implicit none

complex(kind(0d0)), intent(out) :: DF_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(out) :: DF_lambda(1:NMAT,1:NMAT,1:num_links)

complex(kind(0d0)), intent(in) :: eta_mat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: lambda_mat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
integer :: s,l,k



do s=1,num_sites
  tmpmat2=(0d0,0d0)
  do k=1,linkorg_to_s(s)%num_
    l=linkorg_to_s(s)%labels_(k)
    call matrix_3_product(tmpmat2,Umat(:,:,l),lambda_mat(:,:,l),Umat(:,:,l),&
      'C','N','N',&
      (0d0,-1d0)*dcmplx(alpha_l(l))*dconjg(U1Rfactor_link(l)),'ADD')
      !(0d0,-1d0)*dcmplx(alpha_l(l))*dconjg(U1Rfactor_link(l)*U1R_ratio(l)),'ADD')
  enddo
  do k=1,linktip_from_s(s)%num_
    l=linktip_from_s(s)%labels_(k)
    tmpmat2=tmpmat2 + (0d0,1d0)*dcmplx(alpha_l(l)) * lambda_mat(:,:,l)
  enddo
  DF_eta(:,:,s)=DF_eta(:,:,s) &
    + dcmplx(overall_factor) * tmpmat2
enddo

do l=1,num_links
  tmpmat2=(0d0,-1d0)*eta_mat(:,:,link_org(l))
  call matrix_3_product(tmpmat2,Umat(:,:,l),eta_mat(:,:,link_tip(l)),Umat(:,:,l),&
    'N','N','C',(0d0,1d0)*dconjg(U1Rfactor_link(l)),'ADD')
    !'N','N','C',(0d0,1d0)*dconjg(U1Rfactor_link(l)*U1R_ratio(l)),'ADD')
  DF_lambda(:,:,l) = DF_lambda(:,:,l) &
      + dcmplx(alpha_l(l) * overall_factor) * tmpmat2 
enddo

end subroutine prod_Dirac_link1
