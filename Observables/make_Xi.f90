!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to make S_eta, S_lambda, S_chi of
!!  \Zia = Tr(eta S_eta) + Tr(lambda S_lambda) + Tr(chi S_chi)
subroutine make_XiVec(Xi_eta,Xi_lambda,Xi_chi,Umat,PhiMat)
use matrix_functions, only : matrix_commutator,matrix_3_product
use parallel
implicit none

complex(kind(0d0)), intent(out) :: Xi_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(out) :: Xi_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(out) :: Xi_chi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT)
integer :: s,l,f
integer :: i,j

do s=1,num_sites
  call matrix_commutator(tmpmat,PhiMat(:,:,s),Phimat(:,:,s),'N','C')
  Xi_eta(:,:,s)=dcmplx(alpha_s(s)*0.25d0*overall_factor)*tmpmat
enddo

do l=1,num_links
  do j=1,NMAT
    do i=1,NMAT
      tmpmat(i,j)=-dconjg( PhiMat(j,i,link_org(l)) )
    enddo
  enddo
  call matrix_3_product(tmpmat,&
    Umat(:,:,l),PhiMat(:,:,link_tip(l)),Umat(:,:,l),&
    'N','C','C',(1d0,0d0),'ADD') 
  Xi_lambda(:,:,l)=(0d0,-1d0)*dcmplx(alpha_l(l)*overall_factor)*tmpmat
enddo

do f=1,num_faces
  call make_face_variable(Uf,f,Umat)
  if(m_omega==0) call make_moment_map0(tmpmat,Uf)
  if(m_omega==-1) call make_moment_map_adm(tmpmat,Uf)
  Xi_chi(:,:,f)=(0d0,-0.5d0)*dcmplx(alpha_f(f)*beta_f(f)*overall_factor)*tmpmat
enddo

  call syncronize_sites(Xi_eta)
  call syncronize_links(Xi_lambda)
  call syncronize_faces(Xi_chi)
end subroutine make_XiVec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to make S_eta
!!  \Zia = Tr(eta S_eta) 
subroutine make_XiVec_site(Xi_eta,PhiMat)
use matrix_functions, only : matrix_commutator
use parallel
implicit none

complex(kind(0d0)), intent(out) :: Xi_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
integer :: s

do s=1,num_sites
  call matrix_commutator(tmpmat,PhiMat(:,:,s),Phimat(:,:,s),'N','C')
  Xi_eta(:,:,s)=dcmplx(alpha_s(s)*0.25d0*overall_factor)*tmpmat
enddo
call syncronize_sites(Xi_eta)

end subroutine make_XiVec_site
