subroutine Qexact_bosonic_action_link(SB_L,UMAT,PhiMat)
use matrix_functions, only : matrix_3_product
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
double precision, intent(out) :: SB_L
integer :: i,j
integer :: l
double precision :: tmp
complex(kind(0d0)) :: dPhi(1:NMAT,1:NMAT)


SB_L=0d0
do l=1,num_links
  call matrix_3_product(dPhi,Umat(:,:,l),PhiMat(:,:,link_tip(l)),Umat(:,:,l),'N','N','C')
  dPhi=dPhi * U1Rfactor_link(l)**2d0 * U1R_ratio(l)**2d0
  dPhi=dPhi-PhiMat(:,:,link_org(l))

  tmp=0d0
  do i=1,NMAT
  do j=1,NMAT
    tmp=tmp+dble(dPhi(j,i)*dconjg(dPhi(j,i)))
  enddo
  enddo

  SB_L=SB_L+alpha_l(l)*tmp
enddo

SB_L=SB_L*overall_factor 

end subroutine Qexact_bosonic_action_link



