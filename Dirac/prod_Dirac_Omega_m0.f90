!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
subroutine Dirac_Omega_m0(DF_chi,DF_lambda,Umat,lambda_mat,chi_mat)
use matrix_functions, only : make_unit_matrix
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: lambda_mat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: chi_mat(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), intent(out) :: DF_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(out) :: DF_chi(1:NMAT,1:NMAT,1:num_faces)


complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT,1:num_faces)
!complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Xmat(1:NMAT,1:NMAT),Ymat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
integer :: i,j,k,l,f,l_place
complex(kind(0d0)) :: dir_factor
complex(kind(0d0)) :: U1Rfactor_fl

!! m_omega=0 case
do f=1,num_necessary_faces
  do l_place=1,links_in_f(f)%num_
    l=links_in_f(f)%link_labels_(l_place)

    if( f <= num_faces .or. l <= num_links ) then 
      call calc_XYmat(Xmat,Ymat,f,l_place,UMAT)
      !call make_unit_matrix(Ymat)
       
      !call calc_U1Rfactor_fl(U1Rfactor_fl,f,l)
      call calc_U1Rfactor_fl_by_route(U1Rfactor_fl,f,l_place)

      !dir_factor=dcmplx(dble(links_in_f(f)%link_dirs_(l_place)) &
        !* alpha_f(f)*beta_f(f)*overall_factor )*(0d0,1d0)
      dir_factor=&
        (0d0,-1d0)*overall_factor &
        *dcmplx(links_in_f(f)%link_dirs_(l_place))&
        *dcmplx(alpha_f(f)*beta_f(f)) &
        *U1Rfactor_fl

      if( f <= num_faces ) then 
        call matrix_3_product(tmpmat1,Xmat,lambda_mat(:,:,l),Ymat)
        call matrix_3_product(tmpmat1,Ymat,lambda_mat(:,:,l),Xmat,'C','N','C',(1d0,0d0),'ADD')
        DF_chi(:,:,f)=DF_chi(:,:,f) + dir_factor * tmpmat1
      endif
  
      if( l <= num_links ) then 
        call matrix_3_product(tmpmat1,Ymat,chi_mat(:,:,f),Xmat)
        call matrix_3_product(tmpmat1,Xmat,chi_mat(:,:,f),Ymat,'C','N','C',(1d0,0d0),'ADD')
        DF_lambda(:,:,l)=DF_lambda(:,:,l) - dir_factor * tmpmat1
      endif
    endif
  enddo
enddo
end subroutine Dirac_Omega_m0

