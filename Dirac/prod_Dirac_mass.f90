subroutine prod_Dirac_mass(DF_chi,DF_lambda,DF_eta,eta_mat,lambda_mat,chi_mat)
implicit none

complex(kind(0d0)), intent(out) :: DF_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(out) :: DF_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(out) :: DF_chi(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)), intent(in) :: eta_mat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: lambda_mat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: chi_mat(1:NMAT,1:NMAT,1:num_necessary_faces)

!! for fermion mass term
complex(kind(0d0)) :: tmp_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: tmp_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: tmp_chi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), allocatable :: tmp_vec(:),tmp_Dvec(:)
integer :: sizeD
integer :: i

sizeD = dimG*(global_num_sites+global_num_links+global_num_faces)
allocate( tmp_vec(1:sizeD), tmp_Dvec(1:sizeD) )

#ifdef PARALLEL
  call localmat_to_globalvec(tmp_vec,eta_mat,lambda_mat,chi_mat)
  if(MYRANK==0) then 
    tmp_Dvec=(0d0,0d0)
    do i=1,sizeD/2
      tmp_Dvec(2*i-1)=tmp_Dvec(2*i-1) + overall_factor * dcmplx(mass_f)*tmp_vec(2*i) 
      tmp_Dvec(2*i)=tmp_Dvec(2*i) - overall_factor * dcmplx(mass_f)*tmp_vec(2*i-1) 
    enddo
  endif
  call MPI_BCAST(tmp_Dvec,sizeD,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
  call globalvec_to_localmat(tmp_eta,tmp_lambda,tmp_chi,tmp_Dvec)
#else
  call mat_to_vec(tmp_vec,eta_mat,lambda_mat,chi_mat)
  tmp_Dvec=(0d0,0d0)
  do i=1,sizeD/2
    tmp_Dvec(2*i-1)=tmp_Dvec(2*i-1) + overall_factor * dcmplx(mass_f)*tmp_vec(2*i) 
    tmp_Dvec(2*i)=tmp_Dvec(2*i) - overall_factor * dcmplx(mass_f)*tmp_vec(2*i-1) 
  enddo
  call vec_to_mat(tmp_eta,tmp_lambda,tmp_chi,tmp_Dvec)
#endif

  DF_eta=DF_eta+tmp_eta
  DF_lambda=DF_lambda+tmp_lambda
  DF_chi=DF_chi+tmp_chi

end subroutine prod_Dirac_mass
