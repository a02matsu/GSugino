!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate D^\mu J_\mu for U(1)_V current
!!  vec1 ~ 1/2 Tr(\lambda(l) \eta(s))
!!  vec2 ~ Tr(\lambda(l) \chi(f))
!! and
!!  DJ1 = d.vec1, DJ2 = d.vec2
!! where \chi(f) is associated with the link variable Uf
!! as in the fermionic part of the face action
subroutine calc_DJ_U1V(DJ1,DJ2,Glambda_eta,Glambda_chi,UMAT)
use global_parameters
!use initialization_calcobs
use parallel
use global_subroutines, only : syncronize_linkval
use matrix_functions, only :matrix_product, hermitian_conjugate, make_unit_matrix, matrix_exp, matrix_trace
implicit none

complex(kind(0d0)), intent(out) :: DJ1(1:num_faces)
complex(kind(0d0)), intent(out) :: DJ2(1:num_faces)
complex(kind(0d0)), intent(in) :: Glambda_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_sites) 
complex(kind(0d0)), intent(in) :: Glambda_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_faces) 
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)


complex(kind(0d0)) :: trvec1(1:num_necessary_links)
!complex(kind(0d0)) :: vec1(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: U_fs(1:NMAT,1:NMAT)
complex(kind(0d0)) :: U_sf(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)) :: mat1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: mat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmp
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
integer :: ls,ll,lf,gs,gl,gf,kk,dir
integer :: i,j,k,l,a,b


!!  DJ1 ~ 1/2 rot Tr(\lambda(l) \eta(s))
call make_trV1(trvec1,Glambda_eta)
DJ1=(0d0,0d0)
do lf=1,num_faces
  !call make_unit_matrix(mat1)
  !tmp=(1d0,0d0)
  do kk=1,links_in_f(lf)%num_
    ll=links_in_f(lf)%link_labels_(kk)
    dir=links_in_f(lf)%link_dirs_(kk)
    DJ1(lf)=DJ1(lf) + dcmplx(dir)*trvec1(ll)
  enddo
  DJ1(lf)=DJ1(lf)*beta_f(lf)
enddo

!! DJ2 ~ div Tr(\lambda(l) \chi(f))
call make_trV2(trvec2,Gchi_lambda,UMAT)
DJ2=(0d0,0d0)
do lf=1,num_faces
  do kk=1,sites_in_f(lf)%num_
    ls=sites_in_f(lf)%label_(kk)
    !! contribution of [link FROM gs]
    do a=1,linktip_from_s(ls)%num_
      ll=linktip_from_s(ls)%label_(a)
      DJ2(lf) = DJ2(lf) &
        + trvec2(ll) * dcmplx(alpha_l(ll))/dcmplx( num_faces_in_s(ls) )
    enddo
    !! contribution of [link TO gs]
    do a=1,linkorg_to_s(ls)%num_
      ll=linkorg_to_s(ls)%label_(a)
      DJ2(lf) = DJ2(lf) & 
        - trvec2(ll) * dcmplx(alpha_l(ll))/dcmplx( num_faces_in_s(ls) )
    enddo
  enddo
enddo

DJ1=(DJ1)/dcmplx(LatticeSpacing**4)
DJ2=(DJ2)/dcmplx(LatticeSpacing**4)

end subroutine calc_DJ_U1V

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate 1/2 Tr(\lambda(l) \eta(s(l)))
subroutine make_trV1(trvec1,Glambda_eta)
use global_parameters
use parallel
use global_subroutines, only : syncronize_linkval
implicit none

complex(kind(0d0)), intent(out) :: trvec1(1:num_necessary_links) ! 1/2 Tr(\lambda(l) \eta(s))
complex(kind(0d0)), intent(in) :: Glambda_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_sites) 

integer :: ll,ls
integer :: gl
integer :: org_ll
integer :: i,j

trvec1=(0d0,0d0)
do ll=1,num_links
  gl=global_link_of_local(ll)
  ls=link_org(ll)
  do j=1,NMAT
    do i=1,NMAT
      trvec1(ll)=trvec1(ll) + (0.5d0,0d0)*Glambda_eta(i,j,j,i,gl,ls)
    enddo
  enddo
enddo
call syncronize_linkval(trvec1)
end subroutine make_trV1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate 
subroutine make_trV2(trvec2,Gchi_lambda,UMAT)
use global_parameters
use parallel
use global_subroutines, only : syncronize_linkval, calc_prodUl_from_n1_to_n2_in_Uf
implicit none

complex(kind(0d0)), intent(out) :: trvec2(1:num_necessary_links) ! Tr(\lambda(l) \chi(f))
complex(kind(0d0)), intent(in) :: Gchi_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_links) 
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)

complex(kind(0d0)) :: Ucarry(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmp(1:num_faces)
integer :: ll,lf
integer :: gl,gf
integer :: org_ll
integer :: i,j,k,l

trvec2=(0d0,0d0)
do ll=1,num_links
  call find_origin_of_dual_link(lf,org_ll,ll)
  gf=global_face_of_local(lf)
  !write(*,*) gf,global_site_of_local(link_org(ll)),global_site_of_local(link_tip(ll))
  !! Ucarry = U1 ... U_orgll
  call calc_prodUl_from_n1_to_n2_in_Uf(Ucarry,lf,1,org_ll,Umat)
  do l=1,NMAT
    do k=1,NMAT
      do j=1,NMAT
        do i=1,NMAT
          trvec2(ll)=trvec2(ll)&
            + Ucarry(j,k) * dconjg(Ucarry(i,l)) * Gchi_lambda(i,j,k,l,gf,ll)
        enddo
      enddo
    enddo
  enddo  
  !!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!
  !! special treatment of the present discretization
  if(gf==1) trvec2(ll)=-trvec2(ll)
  !!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!
  !write(*,*) global_link_of_local(ll), dble(vec2(ll)), dble((0d0,-1d0)*vec2(ll))
enddo
call syncronize_linkval(trvec2)
end subroutine make_trV2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate D^\mu J_\mu for U(1)_V current
!!  vec1 ~ 1/2 Tr(\lambda(l) \eta(s))
!!  vec2 ~ Tr(\lambda(l) \chi(f))
!! and
!!  DJ1 = d.vec1, DJ2 = d.vec2
!! where \chi(f) is associated with the link variable Uf
!! as in the fermionic part of the face action
subroutine calc_DJ_U1V_org3(DJ1,DJ2,Glambda_eta,Glambda_chi,UMAT)
use global_parameters
!use initialization_calcobs
use parallel
use global_subroutines, only : syncronize_linkval
use matrix_functions, only :matrix_product, hermitian_conjugate, make_unit_matrix, matrix_exp, matrix_trace
implicit none

complex(kind(0d0)), intent(out) :: DJ1(1:num_faces)
complex(kind(0d0)), intent(out) :: DJ2(1:num_faces)
complex(kind(0d0)), intent(in) :: Glambda_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_sites) 
complex(kind(0d0)), intent(in) :: Glambda_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_faces) 
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)


complex(kind(0d0)) :: trvec1(1:num_necessary_links)
!complex(kind(0d0)) :: vec1(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: U_fs(1:NMAT,1:NMAT)
complex(kind(0d0)) :: U_sf(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)) :: mat1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: mat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmp
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
integer :: ls,ll,lf,gs,gl,gf,kk,dir
integer :: i,j,k,l,a,b


!!  DJ1 ~ 1/2 rot Tr(\lambda(l) \eta(s))
call make_trV1(trvec1,Glambda_eta)
!call make_V1(vec1,Glambda_eta)
DJ1=(0d0,0d0)
do lf=1,num_faces
  !call make_unit_matrix(mat1)
  !tmp=(1d0,0d0)
  do kk=1,links_in_f(lf)%num_
    ll=links_in_f(lf)%link_labels_(kk)
    dir=links_in_f(lf)%link_dirs_(kk)
    DJ1(lf)=DJ1(lf) + dcmplx(dir)*trvec1(ll)
  enddo
  DJ1(lf)=DJ1(lf)*beta_f(lf)
enddo

!! DJ2 ~ div Tr(\lambda(l) \chi(f))
DJ2=(0d0,0d0)
do lf=1,num_faces
  do kk=1,sites_in_f(lf)%num_
    ls=sites_in_f(lf)%label_(kk)
    gs=global_site_of_local(ls)
    call calc_prodUl_from_n1_to_n2_in_Uf(U_fs,lf,1,kk-1,Umat) ! 4's argument must be kk-1 
    call hermitian_conjugate(U_sf,U_fs)
    !! contribution of [link FROM gs]
    do a=1,global_linktip_from_s(gs)%num_
      tmp=(0d0,0d0)
      gl=global_linktip_from_s(gs)%labels_(a)
      do l=1,NMAT
        do k=1,NMAT
          do j=1,NMAT
            do i=1,NMAT
              tmp=tmp + Glambda_chi(i,j,k,l,gl,lf)*U_sf(j,k)*U_fs(l,i)
            enddo
          enddo
        enddo
      enddo
      DJ2(lf) = DJ2(lf) + tmp * dcmplx(global_alpha_l(gl))/dcmplx( num_faces_in_s(ls) )
    enddo
    !! contribution of [link TO gs]
    do a=1,global_linkorg_to_s(gs)%num_
      tmp=(0d0,0d0)
      gl=global_linkorg_to_s(gs)%labels_(a)
      do ll=1,num_necessary_links
        if( global_link_of_local(ll) == gl ) exit
      enddo
      call matrix_product(mat1,Umat(:,:,ll),U_sf)
      call matrix_product(mat2,U_fs,Umat(:,:,ll),'N','C')
      do l=1,NMAT
        do k=1,NMAT
          do j=1,NMAT
            do i=1,NMAT
              tmp=tmp + Glambda_chi(i,j,k,l,gl,lf)*mat1(j,k)*mat2(l,i)*dcmplx(global_alpha_l(gl))
            enddo
          enddo
        enddo
      enddo
    enddo
    DJ2(lf) = DJ2(lf) - tmp / dcmplx( num_faces_in_s(ls) )
  enddo
enddo

DJ1=(DJ1)/dcmplx(LatticeSpacing**4)
DJ2=(DJ2)/dcmplx(LatticeSpacing**4)

end subroutine calc_DJ_U1V_org3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate D^\mu J_\mu for U(1)_V current
!!  vec1 ~ 1/2 Tr(\lambda(l) \eta(s))
!!  vec2 ~ Tr(\lambda(l) \chi(f))
!! where \chi(f) is associated with the link variable Uf
!! as in the fermionic part of the face action
subroutine calc_divJ_U1Vorg2(divJ1,divJ2,Glambda_eta,Glambda_chi,UMAT)
use global_parameters
!use initialization_calcobs
use parallel
use global_subroutines, only : syncronize_linkval
implicit none


complex(kind(0d0)), intent(in) :: Glambda_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_sites) 
complex(kind(0d0)), intent(in) :: Glambda_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_faces) 
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)

complex(kind(0d0)) :: vec1(1:NMAT,1:NMAT,1:num_necessary_links) ! 1/2 Tr(\lambda(l) \eta(s))
!complex(kind(0d0)) :: trvec1(1:num_necessary_links) ! 1/2 Tr(\lambda(l) \eta(s))
complex(kind(0d0)) :: trvec2(1:num_necessary_links) ! Tr(\lambda(l) \chi(f))

complex(kind(0d0)) :: divJ1(1:num_faces)
complex(kind(0d0)) :: divJ2(1:num_faces)

!!  vec1 ~ 1/2 Tr(\lambda(l) \eta(s))
call make_V1(vec1,Glambda_eta)
!call make_trV1(vec1,Glambda_eta)
!!  vec2 ~ Tr(\lambda(l) \chi(f))
call make_trV2_likeSf(trvec2,Glambda_chi,UMAT)

call calc_trrot2(divJ1,vec1,Umat)
!call calc_trrot(divJ1,vec1)
call calc_trdiv(divJ2,trvec2)

divJ1=(divJ1)/dcmplx(LatticeSpacing**4)
divJ2=(divJ2)/dcmplx(LatticeSpacing**4)

end subroutine calc_divJ_U1Vorg2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate D^\mu J_\mu for U(1)_V current
subroutine calc_divJ_U1Vorg(divJ,Glambda_eta,Gchi_lambda,UMAT)
use global_parameters
!use initialization_calcobs
use parallel
use global_subroutines, only : syncronize_linkval, calc_prodUl_from_n1_to_n2_in_Uf
implicit none


complex(kind(0d0)), intent(out) :: divJ(1:num_faces)
complex(kind(0d0)), intent(in) :: Glambda_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_sites) 
complex(kind(0d0)), intent(in) :: Gchi_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_links) 
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)

complex(kind(0d0)) :: vec1(1:num_necessary_links) ! 1/2 Tr(\lambda(l) \eta(s))
complex(kind(0d0)) :: vec2(1:num_necessary_links) ! Tr(\lambda(l) \chi(f))
complex(kind(0d0)) :: divJ1(1:num_faces)
complex(kind(0d0)) :: divJ2(1:num_faces)


call make_trV1(vec1,Glambda_eta)
call make_trV2(vec2,Gchi_lambda,UMAT)

call calc_trrot(divJ1,vec1)
call calc_trdiv(divJ2,vec2)

divJ=(divJ1+divJ2)/dcmplx(LatticeSpacing**4)

end subroutine calc_divJ_U1Vorg


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate D^\mu J_\mu for U(1)_V current
subroutine calc_divJ_U1V_2(divJ1,divJ2,Glambda_eta,Gchi_lambda,UMAT)
use global_parameters
!use initialization_calcobs
use parallel
use global_subroutines, only : syncronize_linkval, calc_prodUl_from_n1_to_n2_in_Uf
implicit none


complex(kind(0d0)), intent(in) :: Glambda_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_sites) 
complex(kind(0d0)), intent(in) :: Gchi_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_links) 
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)

complex(kind(0d0)) :: vec1(1:num_necessary_links) ! 1/2 Tr(\lambda(l) \eta(s))
complex(kind(0d0)) :: vec2(1:num_necessary_links) ! Tr(\lambda(l) \chi(f))

complex(kind(0d0)) :: divJ1(1:num_faces)
complex(kind(0d0)) :: divJ2(1:num_faces)

call make_trV1(vec1,Glambda_eta)
call make_trV2(vec2,Gchi_lambda,UMAT)

call calc_trrot(divJ1,vec1)
call calc_trdiv2(divJ2,vec2)

divJ1=(divJ1)/dcmplx(LatticeSpacing**4)
divJ2=(divJ2)/dcmplx(LatticeSpacing**4)

end subroutine calc_divJ_U1V_2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate D^\mu J_\mu for U(1)_V current
!!  rotation is taken as \sum_{l\in f} Tr( X V_l Y )
!!  divergence is taken as summation over site-divergence/#{facce on s}
subroutine calc_divJ_U1V_3(divJ1,divJ2,Glambda_eta,Gchi_lambda,UMAT)
use global_parameters
!use initialization_calcobs
use parallel
implicit none


complex(kind(0d0)), intent(in) :: Glambda_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_sites) 
complex(kind(0d0)), intent(in) :: Gchi_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_links) 
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)

complex(kind(0d0)) :: vec1(1:NMAT,1:NMAT,1:num_necessary_links) ! 1/2 (\lambda(l) \eta(s))_{ij}
complex(kind(0d0)) :: trvec2(1:num_necessary_links) ! Tr(\lambda(l) \chi(f))

complex(kind(0d0)) :: divJ1(1:num_faces)
complex(kind(0d0)) :: divJ2(1:num_faces)

call make_V1(vec1,Glambda_eta)
call make_trV2(trvec2,Gchi_lambda,UMAT)

call calc_trrot2(divJ1,vec1,Umat)
call calc_trdiv(divJ2,trvec2)

divJ1=(divJ1)/dcmplx(LatticeSpacing**4)
divJ2=(divJ2)/dcmplx(LatticeSpacing**4)

end subroutine calc_divJ_U1V_3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate lambda.lambda(U.\bar{phi}.U^-1 + \bar{phi})
subroutine make_V1(vec1,Glambda_eta)
use global_parameters
use parallel
use global_subroutines, only : syncronize_links
implicit none

complex(kind(0d0)), intent(out) :: vec1(1:NMAT,1:NMAT,1:num_necessary_links) ! 1/2 Tr(\lambda(l) \eta(s))
complex(kind(0d0)), intent(in) :: Glambda_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_sites) 

integer :: ll,ls
integer :: gl
integer :: org_ll
integer :: i,j,k


vec1=(0d0,0d0)
do ll=1,num_links
  gl=global_link_of_local(ll)
  ls=link_org(ll)
  do i=1,NMAT
    do j=1,NMAT
      do k=1,NMAT
        vec1(i,j,ll)=vec1(i,j,ll) + (0.5d0,0d0)*Glambda_eta(i,k,k,j,gl,ls)
      enddo
    enddo
  enddo
  !vec1(:,:,ll)=vec1(:,:,ll)*alpha_l(ll)
enddo
call syncronize_links(vec1)


end subroutine make_V1



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate lambda.lambda(U.\bar{phi}.U^-1 + \bar{phi})
subroutine make_V2_bak(vec2,Gchi_lambda,UMAT)
use global_parameters
use parallel
use global_subroutines, only : syncronize_linkval, calc_prodUl_from_n1_to_n2_in_Uf
implicit none

complex(kind(0d0)), intent(out) :: vec2(1:num_necessary_links) ! Tr(\lambda(l) \chi(f))
complex(kind(0d0)), intent(in) :: Gchi_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_links) 
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)

complex(kind(0d0)) :: Ucarry(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmp(1:num_faces)
integer :: ll,lf,ls
integer :: gl,gf
integer :: org_ll
integer :: i,j,k,l,ii

vec2=(0d0,0d0)
do ll=1,num_links
  do ii=1,face_in_l(ll)%num_
    lf=face_in_l(ll)%label_(ii)
    gf=global_face_of_local(lf)
    !! find the place of ll in the face sharing ll
    do org_ll=1,links_in_f(lf)%num_
      if( links_in_f(lf)%link_labels_(org_ll) == ll ) then 
        exit
      endif
    enddo
    if( links_in_f(lf)%link_dirs_(org_ll) == 1 ) then
      org_ll = org_ll - 1
    endif
    !! Ucarry = U1 ... U_orgll
    call calc_prodUl_from_n1_to_n2_in_Uf(Ucarry,lf,1,org_ll,Umat)
    do l=1,NMAT
      do k=1,NMAT
        do j=1,NMAT
          do i=1,NMAT
            vec2(ll)=vec2(ll)&
              + Ucarry(j,k) * dconjg(Ucarry(i,l)) &
                * (0.5d0,0d0) * Gchi_lambda(i,j,k,l,gf,ll)
          enddo
        enddo
      enddo
    enddo
  enddo
enddo
call syncronize_linkval(vec2)
end subroutine make_V2_bak



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! inspierd by
!!  SF^f = 1/g^2 \int ( -2i \chi D\times\lambda )
!!       = 1/g^2 sum_f rot( V_{f,l}
!! V_{f,l} = (i\alpha_\beta_f) (1/B(\chi X \lambda Y +...) )
!! \lambda \chi ~ -1/(-2i) V_{f,l} = -i/2 V_{f,l}
!! 
subroutine make_trV2_likeSf(vec2,Glambda_chi,UMAT)
use global_parameters
use global_subroutines, only : syncronize_linkval
use parallel
implicit none

complex(kind(0d0)), intent(out) :: vec2(1:num_necessary_links) ! Tr(\lambda(l) \chi(f))
complex(kind(0d0)), intent(in) :: Glambda_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_faces) 
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)

complex(kind(0d0)) :: Lff
integer :: ll,lf,l_place,dir
integer :: gl,gf,rank_send,rank_recv,tag
integer :: org_ll
integer :: ii

vec2=(0d0,0d0)
do gl=1,global_num_links
  !! find the origin face of the dual link
  !!  gf: global face label of the origin of the dual link of gl
  !!      defined so that dir(gl) is the same with the direction of the link in gf
  do ii=1,global_face_in_l(gl)%num_
    gf=global_face_in_l(gl)%label_(ii)
    do org_ll=1,global_links_in_f(gf)%num_
      if( global_links_in_f(gf)%link_labels_(org_ll) == gl ) then 
        dir=global_links_in_f(gf)%link_dirs_(org_ll)
        exit
      endif
    enddo
    !!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!
    !! special treatment of the present discretization
    if(gf==1) dir=-dir
    !!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!
    if( dir==1 ) exit
  enddo

  rank_send = local_face_of_global(gf)%rank_
  rank_recv = local_link_of_global(gl)%rank_
  tag=global_num_faces*(gf-1) + gl -1

  !! send phase
  if( MYRANK == rank_send ) then
    lf=local_face_of_global(gf)%label_
    !! find the place of gl in lf
    do l_place=1,links_in_f(lf)%num_
      ll=links_in_f(lf)%link_labels_(l_place)
      if( global_link_of_local(ll) == gl ) exit
    enddo
    call fermionic_face_lagrangian(Lff,lf,l_place,Glambda_chi,Umat)
    Lff=Lff*(0d0,1d0)

    if( MYRANK /= rank_recv ) then 
      call MPI_SEND(Lff,1,MPI_DOUBLE_COMPLEX,rank_recv,tag,MPI_COMM_WORLD,IERR)
    endif
  endif
  !! recv phase
  if( MYRANK == rank_recv ) then
    ll=local_link_of_global(gl)%label_
    if( MYRANK /= rank_send ) then
      call MPI_RECV(vec2(ll),1,MPI_DOUBLE_COMPLEX,rank_send,tag,MPI_COMM_WORLD,ISTATUS,IERR)
    else
      vec2(ll) = Lff
    endif
    !!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!
    !! special treatment of the present discretization
    if(gf==1) vec2(ll)=-vec2(ll)
    !!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!
  endif

  Lff = Lff * (0d0,-0.5d0)
enddo
call syncronize_linkval(vec2)
end subroutine make_trV2_likeSf



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! find the origin face of the dual link
!!  lf     : local face label of the origin of the dual link of ll
!!  org_ll : position of the origin of the link ll in the face lf
subroutine find_origin_of_dual_link(lf,org_ll,ll)
use global_parameters
use parallel
implicit none

integer, intent(out) :: lf, org_ll
integer, intent(in) :: ll

integer :: gf
integer :: ii
integer :: dir

do ii=1,face_in_l(ll)%num_
  lf=face_in_l(ll)%label_(ii)
  gf=global_face_of_local(lf)

  do org_ll=1,links_in_f(lf)%num_
    if( links_in_f(lf)%link_labels_(org_ll) == ll ) then 
      dir=links_in_f(lf)%link_dirs_(org_ll)
      exit
    endif
  enddo
  !!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!
  !! special treatment of the present discretization
  if(gf==1) dir=-dir
  !!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!
  if( dir==1 ) then 
    if( links_in_f(lf)%link_dirs_(org_ll) == 1 ) then
      org_ll = org_ll - 1
    endif
    exit
  endif 
enddo

end subroutine find_origin_of_dual_link


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate lambda.lambda(U.\bar{phi}.U^-1 + \bar{phi})
subroutine calc_SfL2(SfL2, PhiMat, Umat, Glambda_lambda)
use global_parameters
use parallel
use matrix_functions, only : matrix_3_product, trace_MM
use SUN_generators, only : trace_MTa, make_Mijkl_from_modes
implicit none

complex(kind(0d0)), intent(out) :: SfL2
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: Glambda_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_links)


complex(kind(0d0)) :: MMAT(1:NMAT,1:NMAT),tmpmat(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: DP(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
!complex(kind(0d0)) :: modes(1:dimG,1:dimG)
complex(kind(0d0)) :: trace, tmpSfL2, tmp
integer :: gl,ll
integer :: a,b,i,j,k

SfL2=(0d0,0d0)
tmpSfL2=(0d0,0d0)
do ll=1,num_links
  gl=global_link_of_local(ll)
  !!!
  do j=1, NMAT
    do i=1, NMAT
      MMat(i,j)=dconjg( PhiMat(j,i,link_org(ll)) )
    enddo
  enddo
  call matrix_3_product(MMat,Umat(:,:,ll),PhiMat(:,:,link_tip(ll)),Umat(:,:,ll),&
    'N','C','C',(1d0,0d0),'ADD')
  MMat= MMat * alpha_l(ll)
  !!!!
  tmpmat=(0d0,0d0)
  do j=1,NMAT
    do i=1,NMAT
      tmp=(0d0,0d0)
      do k=1,NMAT
        tmp = tmp + Glambda_lambda(i,k,k,j,gl,ll) !DP(i,k,k,j)
      enddo
      trace = trace + tmp*MMAT(j,i)
    enddo
  enddo
  tmpSfL2=tmpSfL2+trace
enddo

call MPI_REDUCE(tmpSfL2,SfL2,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)

if( MYRANK == 0 ) then
  SfL2 = SfL2 * overall_factor
endif

end subroutine calc_SfL2

