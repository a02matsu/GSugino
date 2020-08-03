!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
subroutine Dirac_Omega_adm(DF_chi,DF_lambda,Umat,lambda_mat,chi_mat)
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: lambda_mat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: chi_mat(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), intent(inout) :: DF_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(inout) :: DF_chi(1:NMAT,1:NMAT,1:num_faces)

complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)) :: SinU(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Xmat(1:NMAT,1:NMAT),Ymat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
integer :: i,j,l,f,l_place
complex(kind(0d0)) :: dir_factor
complex(kind(0d0)) :: U1Rfactor_fl
complex(kind(0d0)) :: Bval
complex(kind(0d0)) :: trace

do f=1,num_necessary_faces
  !!!!!!!!!!!!!!!
  !! Uf and SinU
  call make_face_variable(Uf(:,:),f,Umat)
  do i=1,NMAT
    do j=1,NMAT
      sinU(i,j) = Uf(i,j) - dconjg( Uf(j,i) )
    enddo
  enddo
  !!!!!!!!!!!!!!!
  !! Bval
  Bval=(1d0,0d0)
  do i=1,NMAT
    Bval=Bval - ((2d0,0d0)-Uf(i,i)-dconjg(Uf(i,i)))/(e_max*e_max) 
  enddo

  do l_place=1,links_in_f(f)%num_
    l=links_in_f(f)%link_labels_(l_place)

    !! U1Rfactor
    call calc_U1Rfactor_fl_by_route(U1Rfactor_fl,f,l_place)

    if( f <= num_faces .or. l <= num_links ) then 
      dir_factor&
        =(0d0,1d0)*dcmplx(&  ! sign is flipped (2020/06/05)
          dble(links_in_f(f)%link_dirs_(l_place)) &
          * alpha_f(f) * beta_f(f) * overall_factor) &
          * U1Rfactor_fl

      !!!!!!!!!!!!!!!!!!!!!!!!
      !! X_(f,l) and Y_(f,l)
      call calc_XYmat(Xmat,Ymat,f,l_place,UMAT)
    endif

    if( f<= num_faces ) then
      !!!!!!!!!!!!!!!!!!!!!!!!
      !! tmpmat1 = X.lambda.Y  
      !! tmpmat2 = Y^\dag .lambda. X^\dag
      call matrix_3_product(tmpmat1,Xmat,lambda_mat(:,:,l),Ymat)
      call matrix_3_product(tmpmat2,Ymat,lambda_mat(:,:,l),Xmat, 'C','N','C')
      trace=(0d0,0d0)
      do i=1,NMAT
        trace=trace+tmpmat1(i,i)-tmpmat2(i,i)
      enddo

      Df_chi(:,:,f)=Df_chi(:,:,f) &
        + dir_factor/Bval * (tmpmat1+tmpmat2) &
        - dir_factor/(Bval*Bval*e_max*e_max) * trace * SinU
    endif

    if( l<= num_links ) then
      !!!!!!!!!!!!!!!!!!!!!!!!
      !! tmpmat1 = Y.chi.X + X^\dag .chi. Y^\dag
      call matrix_3_product(tmpmat1,Ymat,chi_mat(:,:,f),Xmat)
      call matrix_3_product(tmpmat1,Xmat,chi_mat(:,:,f),Ymat,&
        'C','N','C',(1d0,0d0),'ADD')
      !!!!!!!!!!!!!!!!!!!!!!!!
      !! tmpmat2 = Y.X - X^dag.Y^dag
      call matrix_product(tmpmat2,Ymat,Xmat)
      call matrix_product(tmpmat2,Xmat,Ymat,'C','C',(-1d0,0d0),'ADD')

      trace=(0d0,0d0)
      do i=1,NMAT
        do j=1,NMAT
          trace=trace+chi_mat(i,j,f)*SinU(j,i)
        enddo
      enddo

      Df_lambda(:,:,l)=Df_lambda(:,:,l) - dir_factor/Bval * tmpmat1
      Df_lambda(:,:,l)=Df_lambda(:,:,l) &
        + dir_factor/(Bval*Bval*e_max*e_max) * trace * tmpmat2
    endif
      
  enddo
enddo

end subroutine Dirac_Omega_adm


