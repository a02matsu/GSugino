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
  Xi_eta(:,:,s)=Xi_eta(:,:,s)*dconjg( U1Rfactor_site(s) )
enddo

do l=1,num_links
  do j=1,NMAT
    do i=1,NMAT
      tmpmat(i,j)=-dconjg( PhiMat(j,i,link_org(l)) )
    enddo
  enddo
  call matrix_3_product(tmpmat,&
    Umat(:,:,l),PhiMat(:,:,link_tip(l)),Umat(:,:,l),&
    'N','C','C',dconjg( U1Rfactor_link(l)**2 ),'ADD') 
  Xi_lambda(:,:,l)=(0d0,-1d0)*dcmplx(alpha_l(l)*overall_factor)*tmpmat
  Xi_lambda(:,:,l)=Xi_lambda(:,:,l)*dconjg( U1Rfactor_site(link_org(l)) )
enddo

do f=1,num_faces
  call make_face_variable(Uf,f,Umat)
  if(m_omega==0) call make_moment_map0(tmpmat,Uf)
  if(m_omega==-1) call make_moment_map_adm(tmpmat,Uf)
  Xi_chi(:,:,f)=(0d0,-0.5d0)*dcmplx(alpha_f(f)*beta_f(f)*overall_factor)*tmpmat
  Xi_chi(:,:,f)=Xi_chi(:,:,f)*dconjg( U1Rfactor_site(sites_in_f(f)%label_(1)) )
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
  Xi_eta(:,:,s)=Xi_eta(:,:,s)*dconjg( U1Rfactor_site(s) )
enddo
call syncronize_sites(Xi_eta)

end subroutine make_XiVec_site

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to make S_lambda
!!  \Zia = Tr(lambda S_lambda)
subroutine make_XiVec_link(Xi_lambda,Umat,PhiMat)
use matrix_functions, only : matrix_3_product
use parallel
implicit none

complex(kind(0d0)), intent(out) :: Xi_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
integer :: l
integer :: i,j

Xi_lambda=(0d0,0d0)
do l=1,num_links
  do j=1,NMAT
    do i=1,NMAT
      tmpmat(i,j)=-dconjg( PhiMat(j,i,link_org(l)) )
    enddo
  enddo
  call matrix_3_product(tmpmat,&
    Umat(:,:,l),PhiMat(:,:,link_tip(l)),Umat(:,:,l),&
    'N','C','C',dconjg( U1Rfactor_link(l)**2 ),'ADD') 
    !'N','C','C',(1d0,0d0),'ADD') 
  Xi_lambda(:,:,l)=(0d0,-1d0)*dcmplx(alpha_l(l)*overall_factor)*tmpmat
  Xi_lambda(:,:,l)=Xi_lambda(:,:,l)*dconjg( U1Rfactor_site( link_org(l) ) )
  !tmpmat=-PhiMat(:,:,link_org(l))
  !call matrix_3_product(tmpmat,&
  !  Umat(:,:,l),PhiMat(:,:,link_tip(l)),Umat(:,:,l),&
  !  'N','N','C',(1d0,0d0),'ADD') 
  !do j=1,NMAT
  !  do i=1,NMAT
  !    Xi_lambda(i,j,l)=(0d0,-1d0)*dcmplx(alpha_l(l)*overall_factor)*dconjg(tmpmat(j,i))
  !  enddo
  !enddo
enddo

call syncronize_links(Xi_lambda)
end subroutine make_XiVec_link

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to make S_chi of
!!  \Zia = Tr(eta S_eta) + Tr(lambda S_lambda) + Tr(chi S_chi)
subroutine make_XiVec_face(Xi_chi,Umat)
use matrix_functions, only : matrix_commutator,matrix_3_product,matrix_power
use parallel
implicit none

complex(kind(0d0)), intent(out) :: Xi_chi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)

complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT),Ufm(1:NMAT,1:NMAT)
integer :: s,l,f
integer :: i,j

do f=1,num_faces
  call make_face_variable(Uf,f,Umat)
  if(m_omega==0) then 
    call make_moment_map0(tmpmat,Uf)
  elseif(m_omega==-1) then
    call make_moment_map_adm(tmpmat,Uf)
  else
    call matrix_power(Ufm,Uf,m_omega)
    call make_moment_map(tmpmat,Ufm)
  endif
  Xi_chi(:,:,f)=(0d0,-0.5d0)*dcmplx(alpha_f(f)*beta_f(f)*overall_factor)*tmpmat
  Xi_chi(:,:,f)=Xi_chi(:,:,f)*dconjg( U1Rfactor_site(sites_in_f(f)%label_(1)) )
enddo

  call syncronize_faces(Xi_chi)
end subroutine make_XiVec_face

