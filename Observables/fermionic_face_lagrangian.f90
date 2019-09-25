!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SF^f = 1/2g^2 \sum_f \sum_{l\in Lf} \alpha_f\beta_f\epsilon_{l,f) Lff(f,l)
!! is essentially the rotation of Lff. 
!! In the continuum limit, this part is nothing but
!!  SFc^f = 1/2g^2 \int dx\sqrt{g} Tr( +2i \chi rot(\lambda) )
!! Namely,
!!   Lff ~ +2i \chi \lambda
subroutine fermionic_face_lagrangian(Lff,lf,l_place,Glambda_chi,Umat)
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
    Lff = Lff + (1d0,0d0) / (e_max*e_max*Bval*Bval) &
      * tmp * Tmat(j,i)
  enddo
enddo

Lff = Lff * (0d0,1d0)

end subroutine fermionic_face_lagrangian

