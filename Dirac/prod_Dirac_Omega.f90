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
complex(kind(0d0)) :: Sinmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Omega(1:NMAT,1:NMAT)
complex(kind(0d0)) :: UXmat(1:NMAT,1:NMAT),YUmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat3(1:NMAT,1:NMAT)
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
      Sinmat(i,j) = Ufm(i,j) - dconjg(Ufm(j,i))
    enddo
  enddo
  call Matrix_inverse(Cosinv)

  !! Omega = Cosinv . Sinmat
  call matrix_product(Omega,Cosinv,Sinmat)
  !call make_matrix_traceless(Omega)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! for development
  !call make_unit_matrix(Cosinv)
  !call make_unit_matrix(Sinmat)
  !call make_unit_matrix(Omega)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do l_place=1,links_in_f(f)%num_
    l=links_in_f(f)%link_labels_(l_place)
    !! U1Rfacor
    !call calc_U1Rfactor_fl(U1Rfactor_fl,f,l)
    call calc_U1Rfactor_fl_by_route(U1Rfactor_fl,f,l_place)

    dir_factor=&
      (0d0,-2d0)/dcmplx(m_omega)*overall_factor &
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
        !! UXmat
        call matrix_product(UXmat,Uf0tom(:,:,k),Xmat)
        !! YUmat
        call matrix_product(YUmat,Ymat,Uf0tom(:,:,m_omega-k-1))
  
        !!!!  DF_chi
        if( f<=num_faces ) then 
          tmpmat1=(0d0,0d0)
          tmpmat2=(0d0,0d0)
          ! term 1
          ! tmpmat1 = UX.lambda.YU 
          call matrix_3_product(tmpmat1,UXmat,lambda_mat(:,:,l),YUmat)
          ! term 2
          ! tmpmat2 = YU^dag.lambda.UX^dag
          call matrix_3_product(tmpmat2,YUmat,lambda_mat(:,:,l),UXmat,'C','N','C')
  
          tmpmat3=tmpmat1+tmpmat2
          ! term 3 and term 4
          call matrix_product(tmpmat3,&
            tmpmat1-tmpmat2,Omega, &
            'N','N',(-1d0,0d0),'ADD')
  
          call matrix_product(tmpmat1,Cosinv,tmpmat3)

          DF_chi(:,:,f)=DF_chi(:,:,f) +  dir_factor * tmpmat1 
        endif
  
        !!!!  DF_lambda
        if( l<=num_links ) then 
          ! tmpmat1 = chi.Cosinv
          call matrix_product(tmpmat1,chi_mat(:,:,f),Cosinv)

          ! term 5 
          ! tmpmat3 = -YU.(chi.Cinv).UX - UX^dag.(chi.Cinv).(YU)^dag
          call matrix_3_product(tmpmat3,YUmat,-tmpmat1,UXmat)
          ! term 6
          call matrix_3_product(tmpmat3,UXmat,-tmpmat1,YUmat,'C','N','C',(1d0,0d0),'ADD')
  
          ! tmpmat2 = Sin.Cosinv.chi.Cosinv
          call matrix_product(tmpmat2,Omega,tmpmat1)
          
          call matrix_3_product(tmpmat3,YUmat,tmpmat2,UXmat,'N','N','N',&
            (1d0,0d0),'ADD')
          call matrix_3_product(tmpmat3,UXmat,tmpmat2,YUmat,'C','N','C',&
            (-1d0,0d0),'ADD')
  
          DF_lambda(:,:,l)=DF_lambda(:,:,l) + dir_factor * tmpmat3 
        endif
      enddo
    endif
  enddo
enddo


end subroutine Dirac_Omega


