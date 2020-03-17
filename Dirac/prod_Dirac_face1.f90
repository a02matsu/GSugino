subroutine prod_Dirac_face1(DF_chi,PhiMat,chi_mat)
implicit none

complex(kind(0d0)), intent(out) :: DF_chi(1:NMAT,1:NMAT,1:num_faces)

complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)),intent(in) :: chi_mat(1:NMAT,1:NMAT,1:num_necessary_faces)

complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)

integer :: f,s

do f=1,num_faces
  s=sites_in_f(f)%label_(1)
  call matrix_commutator(tmpmat1,PhiMat(:,:,s),chi_mat(:,:,f))
  DF_chi(:,:,f)=DF_chi(:,:,f) &
      + (-2d0,0d0) * dcmplx(alpha_f(f) * overall_factor) * tmpmat1
enddo


end subroutine prod_Dirac_face1
