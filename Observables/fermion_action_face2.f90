!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SF^f = 1/2g^2 \sum_f \sum_{l\in Lf} \alpha_f\beta_f\epsilon_{l,f) Lff(f,l)
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

Sf_face=(0d0,0d0)
tmp_Sf_face=(0d0,0d0)
do lf=1,num_faces
  do l_place=1,links_in_f(lf)%num_
    call fermionic_face_lagrangian(Lff,lf,l_place,Glambda_chi,Umat)
    write(*,*) Lff
    dir = links_in_f(lf)%link_dirs_(l_place) 
    tmp_Sf_face = tmp_Sf_face &
      + dcmplx( dble(dir) * alpha_f(lf)*beta_f(lf) ) * Lff
  enddo
enddo

call MPI_REDUCE(tmp_Sf_face,Sf_face,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)

Sf_face=Sf_face * dcmplx( overall_factor )

end subroutine calc_Sf_face2

