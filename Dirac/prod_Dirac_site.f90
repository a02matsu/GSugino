subroutine prod_Dirac_site(DF_eta,PhiMat,eta_mat)
implicit none

complex(kind(0d0)), intent(out) :: DF_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(in) :: eta_mat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
integer :: s

!! (1) site action 
do s=1,num_sites
  call matrix_commutator(tmpmat1,PhiMat(:,:,s),eta_mat(:,:,s))
  DF_eta(:,:,s)=DF_eta(:,:,s) &
      +dcmplx(alpha_s(s))*(-0.5d0,0d0) *overall_factor*tmpmat1
enddo


end subroutine prod_Dirac_site
