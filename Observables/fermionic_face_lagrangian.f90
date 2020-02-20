!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SF^f = 1/2g^2 \sum_f \sum_{l\in Lf} \alpha_f\beta_f\epsilon_{l,f) Lff(f,l)
!! is essentially the rotation of Lff. 
!! In the continuum limit, this part is nothing but
!!  SFc^f = 1/2g^2 \int dx\sqrt{g} Tr( +2i \chi rot(\lambda) )
!! Namely,
!!   Lff ~ +2i \chi \lambda
subroutine fermionic_face_lagrangian_adm(Lff,lf,l_place,Glambda_chi,Umat)
use global_parameters
use matrix_functions, only : matrix_product
implicit none

complex(kind(0d0)), intent(out) :: Lff
integer, intent(in) :: lf      ! local face
integer, intent(in) :: l_place ! place of the link in lf
complex(kind(0d0)), intent(in) :: Glambda_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_faces) 
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)

complex(kind(0d0)) :: Xmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Ymat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Usin(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Tmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Bval
complex(kind(0d0)) :: tmp
integer :: gl
integer :: i,j,k,l

gl = global_link_of_local( links_in_f(lf)%link_labels_(l_place) )
call calc_XYmat(Xmat,Ymat,lf,l_place,UMAT)
call make_face_variable(Uf(:,:),lf,Umat)
call calc_Bval(Bval,Uf)

Lff=(0d0,0d0)
do l=1,NMAT
  do k=1,NMAT
    do j=1,NMAT
      do i=1,NMAT
        Lff=Lff - (1d0,0d0)/Bval &
          * Glambda_chi(i,j,k,l,gl,lf) &
          * (Ymat(j,k)*Xmat(l,i) + dconjg(Xmat(k,j)*Ymat(i,l)))
      enddo
    enddo
  enddo
enddo

!! Usin
do i=1,NMAT
  do j=1,NMAT
    Usin(i,j) = Uf(i,j) - dconjg(Uf(j,i))
  enddo
enddo
!! Tmat
Tmat=(0d0,0d0)
call matrix_product(Tmat,Ymat,Xmat)
call matrix_product(Tmat,Xmat,Ymat,'C','C',(-1d0,0d0),'ADD')
do i=1,NMAT
  do j=1,NMAT
    tmp=(0d0,0d0)
    do k=1,NMAT
      do l=1,NMAT
        tmp=tmp+Glambda_chi(i,j,k,l,gl,lf)*Usin(l,k)
      enddo
    enddo
    Lff = Lff + (1d0,0d0) / dcmplx(e_max*e_max*Bval*Bval) &
      * tmp * Tmat(j,i)
  enddo
enddo

Lff = Lff * (0d0,1d0)

end subroutine fermionic_face_lagrangian_adm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SF^f = 1/2g^2 \sum_f \sum_{l\in Lf} \alpha_f\beta_f\epsilon_{l,f) Lff(f,l)
subroutine fermionic_face_lagrangian(Lff,lf,l_place,Glambda_chi,Umat)
use global_parameters
use matrix_functions, only : matrix_product,Hermitian_conjugate,matrix_power,matrix_inverse
implicit none

complex(kind(0d0)), intent(out) :: Lff
integer, intent(in) :: lf      ! local face
integer, intent(in) :: l_place ! place of the link in lf
complex(kind(0d0)), intent(in) :: Glambda_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_faces) 
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)

complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT), Ufm(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Uf0tom(1:NMAT,1:NMAT,0:m_omega-1)
complex(kind(0d0)) :: Sinmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Cosinv(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Omega(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Xmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Ymat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: UXmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: YUmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Mat1(1:NMAT,1:NMAT), Mat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmp
integer :: gl
integer :: i,j,r

complex(kind(0d0)) :: U1Rfactor_fl

!! 
gl = global_link_of_local( links_in_f(lf)%link_labels_(l_place) )

!! U1Rfacor
call calc_U1Rfactor_fl_by_route(U1Rfactor_fl,lf,l_place)
!call calc_U1Rfactor_fl(U1Rfactor_fl,lf,links_in_f(lf)%link_labels_(l_place) )

!! Ufm
call make_face_variable(Uf(:,:),lf,Umat)
call matrix_power(Ufm,Uf,m_omega)
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

!!
call calc_XYmat(Xmat,Ymat,lf,l_place,UMAT)
do r=0,m_omega-1
  call matrix_power(Uf0tom(:,:,r),Uf(:,:),r)
enddo

Lff=(0d0,0d0)
do r=0,m_omega-1
  !! UXmat
  call matrix_product(UXmat,Uf0tom(:,:,r),Xmat)
  !! YUmat
  call matrix_product(YUmat,Ymat,Uf0tom(:,:,m_omega-r-1))

  !! term 1
  Mat1=-YUmat
  call matrix_product(Mat2,Cosinv,UXmat)
  call update_Lff(Lff, Glambda_chi, gl, lf, Mat1, Mat2)

  !! term 2
  call Hermitian_conjugate(Mat1, -UXmat)
  call matrix_product(Mat2,Cosinv,YUmat,'N','C')
  call update_Lff(Lff, Glambda_chi, gl, lf, Mat1, Mat2)

  !! term 3
  call matrix_product(Mat1,YUmat,Omega)
  call matrix_product(Mat2,Cosinv,UXmat)
  call update_Lff(Lff, Glambda_chi, gl, lf, Mat1, Mat2)

  !! term 4
  call matrix_product(Mat1,-UXmat,Omega,'C','N')
  call matrix_product(Mat2,Cosinv,YUmat,'N','C')
  call update_Lff(Lff, Glambda_chi, gl, lf, Mat1, Mat2)
enddo

Lff = Lff * (0d0,2d0)/dcmplx(m_omega) * U1Rfactor_fl


end subroutine fermionic_face_lagrangian

!!!!!!!!!!!!
subroutine update_Lff(Lff, Glambda_chi, gl, lf, Mat1, Mat2)
use global_parameters
implicit none

complex(kind(0d0)), intent(inout) :: Lff
complex(kind(0d0)), intent(in) :: Glambda_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_faces) 
integer, intent(in) :: gl,lf
complex(kind(0d0)), intent(in) :: Mat1(1:NMAT,1:NMAT), Mat2(1:NMAT,1:NMAT)

integer :: i,j,k,l

do l=1,NMAT
  do k=1,NMAT
    do j=1,NMAT
      do i=1,NMAT
        Lff=Lff + Glambda_chi(i,j,k,l,gl,lf) * Mat1(j,k) * Mat2(l,i)
      enddo
    enddo
  enddo
enddo

end subroutine update_Lff








