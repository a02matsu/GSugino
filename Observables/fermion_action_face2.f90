!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SF^f = 1/2g^2 \sum_f \sum_{l\in Lf} \alpha_f\beta_f\epsilon_{l,f) Lff(f,l)
!! is essentially the rotation of Lff. 
!! In the continuum limit, this part is nothing but
!!  SFc^f = 1/2g^2 \int dx\sqrt{g} Tr( +2i \chi rot(\lambda) )
!! Namely,
!!   Lff ~ +2i \chi \lambda
subroutine calc_Sf_face2(Sf_face,Glambda_chi,Umat)
use global_parameters
use matrix_functions, only : matrix_product
implicit none

complex(kind(0d0)), intent(out) :: Sf_face
complex(kind(0d0)), intent(in) :: Glambda_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_faces) 
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)

complex(kind(0d0)) :: Lff, tmp_Sf_face
integer :: lf      ! local face
integer :: l_place ! place of the link in lf
integer :: dir
!complex(kind(0d0)) :: Xmat(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: Ymat(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: Usin(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: Tmat(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: Bval
!complex(kind(0d0)) :: tmp
!integer :: gl,dir
!integer :: i,j,k,l

Sf_face=(0d0,0d0)
do lf=1,num_faces
  do l_place=1,links_in_f(lf)%num_
    call fermionic_face_lagrangian(Lff,lf,l_place,Glambda_chi,Umat)
    dir = links_in_f(lf)%link_dirs_(l_place) 
    tmp_Sf_face = tmp_Sf_face &
      + (0d0,1d0)*dcmplx( dble(dir) * alpha_f(lf)*beta_f(lf) ) * Lff
  enddo
enddo

call MPI_REDUCE(tmp_Sf_face,Sf_face,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)

Sf_face=Sf_face * dcmplx( overall_factor )

end subroutine calc_Sf_face2
