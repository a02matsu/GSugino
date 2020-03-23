!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
subroutine Dirac_Omega(DF_chi,DF_lambda,Umat,lambda_mat,chi_mat)
use matrix_functions, only : make_matrix_traceless
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: lambda_mat(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: chi_mat(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)), intent(out) :: DF_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(out) :: DF_chi(1:NMAT,1:NMAT,1:num_faces)

complex(kind(0d0)) :: Ufm(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Uf0tom(1:NMAT,1:NMAT,0:m_omega-1)
complex(kind(0d0)) :: Xmat(1:NMAT,1:NMAT),Ymat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Cosinv(1:NMAT,1:NMAT)
complex(kind(0d0)) :: UK_Cosinv(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
integer :: i,j,k,l,f,a, l_place
complex(kind(0d0)) :: dir_factor
complex(kind(0d0)) :: im_over_2

integer :: ll,last_place
complex(kind(0d0)) :: U1Rfactor_fl


do f=1,num_necessary_faces
!! preparation( Cos^{-1} and Omega )
  call Make_face_variable(Uf(:,:),f,UMAT) 
  if( m_omega == 1) then 
    Ufm=Uf
  else
    call matrix_power(Ufm,Uf(:,:),m_omega)
  endif
  !! Sinmat and Cos^{-1}
  do i=1,NMAT
    do j=1,NMAT
      Cosinv(i,j) = Ufm(i,j) + dconjg(Ufm(j,i))
    enddo
  enddo
  call Matrix_inverse(Cosinv)

  do l_place=1,links_in_f(f)%num_
    l=links_in_f(f)%link_labels_(l_place)
    !! U1Rfacor
    !call calc_U1Rfactor_fl(U1Rfactor_fl,f,l)
    call calc_U1Rfactor_fl_by_route(U1Rfactor_fl,f,l_place)

    dir_factor=&
      (0d0,4d0)/dcmplx(m_omega)*overall_factor & !sign had been flipped (200318)
      *dcmplx(links_in_f(f)%link_dirs_(l_place))&
      *dcmplx(alpha_f(f)*beta_f(f)) &
      *U1Rfactor_fl


    if( f<=num_faces .or. l<=num_links ) then 

      call calc_XYmat(Xmat,Ymat,f,l_place,UMAT)
      !!
      do k=0,m_omega-1
        call matrix_power(Uf0tom(:,:,k),Uf(:,:),k)
      enddo

      do k=0,m_omega-1
        call matrix_product(UK_Cosinv,Uf0tom(:,:,k),Cosinv)
        !!!!  DF_chi
        if( f<=num_faces ) then 
          call matrix_product(tmpmat2,Ymat,UK_Cosinv,'N','N')
          call matrix_3_product(tmpmat1,tmpmat2,lambda_mat(:,:,l),tmpmat2,'C','N','N')
          !!
          call matrix_product(tmpmat2,UK_Cosinv,Xmat)
          call matrix_3_product(tmpmat1,tmpmat2,lambda_mat(:,:,l),tmpmat2,&
            'N','N','C', (1d0,0d0),'ADD')
          !!
          DF_chi(:,:,f)=DF_chi(:,:,f) +  dir_factor * tmpmat1 
        endif
  
        !!!!  DF_lambda
        if( l<=num_links ) then 
          call matrix_3_product(tmpmat2,UK_Cosinv,chi_mat(:,:,f),UK_Cosinv,'N','N','C')
          call matrix_3_product(tmpmat1,Ymat,tmpmat2,Ymat,'N','N','C')
          !!
          call matrix_3_product(tmpmat2,UK_Cosinv,chi_mat(:,:,f),UK_Cosinv,'C','N','N')
          call matrix_3_product(tmpmat1,Xmat,tmpmat2,Xmat,'C','N','N',&
            (1d0,0d0),'ADD')
  
          DF_lambda(:,:,l)=DF_lambda(:,:,l) - dir_factor * tmpmat1 
        endif
      enddo
    endif
  enddo
enddo


end subroutine Dirac_Omega


