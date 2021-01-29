!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate D^\mu J_\mu for U(1)_V current
!!  vec1 ~ 1/2 Tr(\lambda(l) \eta(s))
!!  vec2 ~ Tr(\lambda(l) \chi(f))
!! and
!!  DJ1 = d.vec1, DJ2 = d.vec2
!! where \chi(f) is associated with the link variable Uf
!! as in the fermionic part of the face action
subroutine calc_DJ_U1V(DJ1,DJ2,Glambda_eta,Gchi_lambda,UMAT)
use global_parameters
!use initialization_calcobs
use parallel
use global_subroutines, only : syncronize_linkval
use matrix_functions, only :matrix_product, hermitian_conjugate, make_unit_matrix, matrix_exp, matrix_trace
implicit none

complex(kind(0d0)), intent(out) :: DJ1(1:num_faces)
complex(kind(0d0)), intent(out) :: DJ2(1:num_faces)
complex(kind(0d0)), intent(in) :: Glambda_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_sites) 
complex(kind(0d0)), intent(in) :: Gchi_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_links) 
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)


complex(kind(0d0)) :: trvec1(1:num_necessary_links)
complex(kind(0d0)) :: trvec2(1:num_necessary_links)
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
      ll=linktip_from_s(ls)%labels_(a)
      DJ2(lf) = DJ2(lf) &
        + trvec2(ll) * dcmplx(alpha_l(ll))/dcmplx( num_faces_in_s(ls) )
    enddo
    !! contribution of [link TO gs]
    do a=1,linkorg_to_s(ls)%num_
      ll=linkorg_to_s(ls)%labels_(a)
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
!! calculate Tr(\lambda(l) \chi(f) ) where f:sorce of dual link of l
subroutine make_trV2(trvec2,Gchi_lambda,UMAT)
use global_parameters
use parallel
use global_subroutines, only : syncronize_linkval, calc_prodUl_from_n1_to_n2_in_Uf
implicit none

complex(kind(0d0)), intent(out) :: trvec2(1:num_necessary_links) ! Tr(\lambda(l) \chi(f))
complex(kind(0d0)), intent(in) :: Gchi_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_necessary_links) 
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)

complex(kind(0d0)) :: Ucarry(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmp(1:num_faces)
integer :: ll,lf
integer :: gl,gf
integer :: org_ll
integer :: i,j,k,l
integer :: info

trvec2=(0d0,0d0)
do ll=1,num_links
  !call find_origin_of_dual_link(lf,org_ll,ll) ! org_ll:position of the origin of ll in the face fl
  !gf=global_face_of_local(lf)
  call find_global_origin_of_dual_local_link(gf,org_ll,ll) ! org_ll:position of the origin of ll in the face fl
  info=1
  do lf=1,num_necessary_faces
    if( global_face_of_local(lf) == gf ) then
      info=0
      exit
    endif
  enddo
  if( info==1 ) write(*,*) "There is no global face", gf, "in RANK", MYRANK
  call calc_prodUl_from_n1_to_n2_in_Uf(Ucarry,lf,1,org_ll-1,Umat)
  do l=1,NMAT
    do k=1,NMAT
      do j=1,NMAT
        do i=1,NMAT
          trvec2(ll)=trvec2(ll)&
            - Ucarry(j,k) * dconjg(Ucarry(i,l)) * Gchi_lambda(i,j,k,l,gf,ll)
        enddo
      enddo
    enddo
  enddo  
  !!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!
  !! special treatment of the present discretization
  !if(gf==1) trvec2(ll)=-trvec2(ll)
  !!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!
  !write(*,*) global_link_of_local(ll), dble(vec2(ll)), dble((0d0,-1d0)*vec2(ll))
enddo
call syncronize_linkval(trvec2)
end subroutine make_trV2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate 4-fermion part of <O(f_1) dJ(f_2)>
subroutine calc_OdJ_4F(\
  OdJ_4F, \
  Geta_eta,\
  Gchi_eta,\
  Geta_chi,\
  Geta_lambda,\
  Gchi_lambda,\
  Gchi_chi, \
  Phimat,Umat)
use global_parameters
!use initialization_calcobs
use parallel
use global_subroutines, only : syncronize_linkval
use matrix_functions, only :matrix_product, hermitian_conjugate, make_unit_matrix, matrix_exp, matrix_trace
implicit none

complex(kind(0d0)), intent(out) :: OdJ_4F(1:global_num_faces,1:global_num_faces)
complex(kind(0d0)), intent(in) :: Geta_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_necessary_sites) 
complex(kind(0d0)), intent(in) :: Gchi_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_necessary_sites) 
complex(kind(0d0)), intent(in) :: Geta_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_necessary_faces) 
complex(kind(0d0)), intent(in) :: Geta_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_necessary_links) 
complex(kind(0d0)), intent(in) :: Gchi_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_necessary_links) 
complex(kind(0d0)), intent(in) :: Gchi_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_necessary_faces) 
complex(kind(0d0)), intent(in) :: PHIMAT(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)

complex(kind(0d0)) :: GchiGchi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:global_num_faces) 
complex(kind(0d0)) :: GetaGchi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:global_num_faces) 
complex(kind(0d0)) :: tmpOdJ_4F(1:global_num_faces,1:num_faces)
complex(kind(0d0)) :: DD!(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
complex(kind(0d0)), allocatable :: phibar_p(:,:,:,:)
integer :: ratio,rank,tag
integer :: gs1,gl1,gf1,gs2,gl2,gf2,ls2,ll2,lf2,ls,gs,gfll2,fofll2
integer :: i,j,k,l,a,b,c,d,p,q,r

complex(kind(0d0)) :: Ucarry(1:NMAT,1:NMAT,1:num_necessary_links)
integer :: orgll, info

!! preparation
! ratio 
ratio = dimG*(global_num_sites-global_num_links+global_num_faces)/2
! phibar_p
allocate(phibar_p(1:NMAT,1:NMAT,0:ratio-1,1:global_num_sites))
do ls=1,num_sites
  gs=global_site_of_local(ls)
  call make_unit_matrix(phibar_p(:,:,0,gs))
  call hermitian_conjugate(phibar_p(:,:,1,gs),phimat(:,:,ls))
  do k=2,ratio-1
    call matrix_product(phibar_p(:,:,k,gs),phibar_p(:,:,k-1,gs),phimat(:,:,ls),'N','C')
  enddo
enddo
do gs=1,global_num_sites
  rank=local_site_of_global(gs)%rank_
  call MPI_BCAST(phibar_p(:,:,:,gs),NMAT*NMAT*ratio, MPI_DOUBLE_COMPLEX,rank,MPI_COMM_WORLD,IERR)
enddo
! dualU
call calc_dualU(Ucarry,Umat)
! GchiGchi
call construct_GFGF(GchiGchi,GetaGchi,Gchi_chi,Geta_chi)


OdJ_4F=(0d0,0d0)
tmpOdJ_4F=(0d0,0d0)
do gf1=1,global_num_faces
  gs1=global_sites_in_f(gf1)%label_(1)
  do lf2=1,num_faces
    do p=1,links_in_f(lf2)%num_
      ll2=links_in_f(lf2)%link_labels_(p)
      ls2=link_org(ll2)
      do l=1,NMAT
        do k=1,NMAT
          do j=1,NMAT
            do i=1,NMAT
              DD=(0d0,0d0)
              do b=1,NMAT
                do a=1,NMAT
                  DD=DD &
                    - Gchi_lambda(i,j,a,b,gf1,ll2)*Geta_eta(k,l,b,a,gs1,ls2) &
                    + Gchi_eta(i,j,b,a,gf1,ls2)*Geta_lambda(k,l,a,b,gs1,ll2)
                enddo
              enddo
              do q=0,ratio-1
                DD=DD*phibar_p(j,k,q,gs1)*phibar_p(l,i,ratio-q-1,gs1)
              enddo
              tmpOdJ_4F(gf1,lf2) = tmpOdJ_4F(gf1,lf2) &
                + DD * dcmplx( -beta_f(lf2)/(2d0*dble(NMAT)) )
            enddo
          enddo
        enddo
      enddo
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do p=1,sites_in_f(lf2)%num_
      ls2=sites_in_f(lf2)%label_(p)
      !! contribution of [link FROM ls2]
      do q=1,linktip_from_s(ls2)%num_
        ll2=linktip_from_s(ls2)%labels_(q)
        call find_global_origin_of_dual_global_link(gfll2,orgll,global_link_of_local(ll2))
        !! DD includes alpha_l(ll2)
        call calc_DD(DD,GchiGchi,Geta_lambda,Gchi_lambda,GetaGchi,&
          phibar_p,Ucarry(:,:,ll2),ratio,gf1,gs1,gfll2,ll2)
        tmpOdJ_4F(gf1,lf2) = tmpOdJ_4F(gf1,lf2) &
          + DD / dcmplx(NMAT *  num_faces_in_s(ls2) )
      enddo
      !! contribution of [link TO gs]
      do q=1,linkorg_to_s(ls2)%num_
        ll2=linkorg_to_s(ls2)%labels_(q)
        call find_global_origin_of_dual_global_link(gfll2,orgll,global_link_of_local(ll2))
        call calc_DD(DD,GchiGchi,Geta_lambda,Gchi_lambda,GetaGchi,phibar_p,Ucarry(:,:,ll2),ratio,gf1,gs1,gfll2,ll2)
        !! DD includes alpha_l(ll2)
        tmpOdJ_4F(gf1,lf2) = tmpOdJ_4F(gf1,lf2) &
          - DD / dcmplx(NMAT *  num_faces_in_s(ls2) )
      enddo
    enddo
  enddo
enddo

tmpOdJ_4F = tmpOdJ_4F / dcmplx( LatticeSpacing**(ratio+6) )

do gf2=1,global_num_faces
  lf2=local_face_of_global(gf2)%label_
  rank=local_face_of_global(gf2)%rank_
  tag=gf2
  if( MYRANK==0 ) then
    if( MYRANK==rank) then
      OdJ_4F(:,gf2)=tmpOdJ_4F(:,lf2)
    else
      call MPI_RECV(OdJ_4F(:,gf2),global_num_sites,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,ISTATUS,IERR)
    endif
  else
    if( MYRANK==rank ) then
      call MPI_SEND(tmpOdJ_4F(:,lf2),global_num_sites,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,IERR)
    endif
  endif
enddo

contains
  subroutine calc_DD(DD,GchiGchi,Geta_lambda,Gchi_lambda,GetaGchi,phibar_p,Ucarry,ratio,gf1,gs1,gfll2,ll2)
  use global_parameters
  implicit none
  complex(kind(0d0)), intent(out) :: DD
  integer, intent(in) :: ratio
  integer, intent(in) :: gf1,gs1,gfll2,ll2
  complex(kind(0d0)), intent(in) :: GchiGchi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:global_num_faces) 
  complex(kind(0d0)), intent(in) :: Geta_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_necessary_links) 
  complex(kind(0d0)), intent(in) :: Gchi_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_necessary_links) 
  complex(kind(0d0)), intent(in) :: GetaGchi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:global_num_faces) 

  complex(kind(0d0)), intent(in) :: Phibar_p(1:NMAT,1:NMAT,0:ratio-1,1:global_num_sites)
  complex(kind(0d0)), intent(in) :: Ucarry(1:NMAT,1:NMAT)

  integer :: i,j,k,l,a,b,c,d,p
  complex(kind(0d0)) :: tmp

  DD=(0d0,0d0)
  do i=1,NMAT
    do j=1,NMAT
      do k=1,NMAT
        do l=1,NMAT
          do a=1,NMAT
            do b=1,NMAT
              do c=1,NMAT
                do d=1,NMAT
                  tmp = - GchiGchi(i,j,a,b,gf1,gfll2)*Geta_lambda(k,l,c,d,gs1,ll2) &
                    + Gchi_lambda(i,j,c,d,gf1,ll2)*GetaGchi(k,l,a,b,gs1,gfll2) 
                  do p=0,ratio-1
                    DD = DD &
                      + tmp * phibar_p(j,k,p,gs1) * phibar_p(l,i,ratio-p-1,gs1) &
                        * Ucarry(b,c) * dconjg( Ucarry(a,d) )
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
  DD = DD * dcmplx( alpha_l(ll2) )

  end subroutine calc_DD

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Udual : Wilson line from f(l) to org(l) 
  !!          where f(l) is the origin of the dual link of l
  subroutine calc_dualU(dualU, Umat)
  use global_parameters
  implicit none
  
  complex(kind(0d0)), intent(out) :: dualU(1:NMAT,1:NMAT,1:num_necessary_links)
  complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
  
  integer :: ll, lf, gl, gf ,org_ll, info
  
  do ll=1,num_links
    call find_global_origin_of_dual_local_link(gf,org_ll,ll)
    info=1
    do lf=1,num_necessary_faces
      if( global_face_of_local(lf) == gf ) then
        info=0
        exit
      endif
    enddo
    if( info==1 ) then
      write(*,*) "[calc_dualU] There is no global face", gf, "in RANK", MYRANK
    endif
    call calc_prodUl_from_n1_to_n2_in_Uf(dualU(:,:,ll),lf,1,org_ll-1,Umat)
  enddo
  
  call syncronize_links(dualU)
  
  end subroutine calc_dualU

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! 
  subroutine construct_GFGF(GchiGchi,GetaGchi,Gchi_chi,Geta_chi)
  use global_parameters
  implicit none
  
  complex(kind(0d0)), intent(out) :: GchiGchi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:global_num_faces)
  complex(kind(0d0)), intent(out) :: GetaGchi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:global_num_faces)
  complex(kind(0d0)), intent(in) :: Gchi_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_necessary_faces)
  complex(kind(0d0)), intent(in) :: Geta_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_necessary_faces)
  
  integer :: gf,lf,rank
  
  do gf=1,global_num_faces
    rank = local_face_of_global(gf)%rank_
    lf = local_face_of_global(gf)%label_

    if(MYRANK==rank) then
      GchiGchi(:,:,:,:,:,gf)=Gchi_chi(:,:,:,:,:,lf)
      GetaGchi(:,:,:,:,:,gf)=Geta_chi(:,:,:,:,:,lf)
    endif
    call MPI_BCAST(GchiGchi(:,:,:,:,:,gf), NMAT**4*global_num_faces, MPI_DOUBLE_COMPLEX,rank,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(GetaGchi(:,:,:,:,:,gf), NMAT**4*global_num_sites, MPI_DOUBLE_COMPLEX,rank,MPI_COMM_WORLD,IERR)

  enddo

  end subroutine construct_GFGF


end subroutine calc_OdJ_4F


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! find the origin face of the dual link
!!  gf     : global face label of the origin of the dual link of ll
!!  org_ll : position of the origin of the link ll in the face lf
!! The origin of dual link is defined so that 
!! the face who includes that link in the POSITIVE direction.
subroutine find_global_origin_of_dual_local_link(gf,org_ll,ll)
use global_parameters
use parallel
implicit none

integer, intent(out) :: gf, org_ll
integer, intent(in) :: ll

integer :: gl
integer :: ii
integer :: dir
integer :: info

info=1
gl=global_link_of_local(ll)
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
  if( dir==1 ) then 
    if(gf==1) then
      org_ll = org_ll + 1
      if( org_ll > global_links_in_f(gf)%num_ ) org_ll = 1
    endif
    info=0
    return
  endif 
enddo

if(info==1) write(*,*) "ERROR in finding the dual link of", global_link_of_local(ll), "in RANK=", MYRANK

end subroutine find_global_origin_of_dual_local_link


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! find the origin face of the dual link
!!  gl     : global link
!!  gf     : global face label of the origin of the dual link of global gl
!!  org_gl : position of org(gl) in the face gf
!! The origin of dual link is defined so that 
!! the face who includes that link in the POSITIVE direction.
subroutine find_global_origin_of_dual_global_link(gf,org_gl,gl)
use global_parameters
use parallel
implicit none

integer, intent(out) :: gf, org_gl
integer, intent(in) :: gl

integer :: ii
integer :: dir
integer :: info

info=1
do ii=1,global_face_in_l(gl)%num_
  gf=global_face_in_l(gl)%label_(ii)

  do org_gl=1,global_links_in_f(gf)%num_
    if( global_links_in_f(gf)%link_labels_(org_gl) == gl ) then 
      dir=global_links_in_f(gf)%link_dirs_(org_gl)
      exit
    endif
  enddo
  !!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!
  !! special treatment of the present discretization
  if(gf==1) dir=-dir
  !!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!
  if( dir==1 ) then 
    if(gf==1) then
      org_gl = org_gl + 1
      if( org_gl > global_links_in_f(gf)%num_ ) org_gl = 1
    endif
    info=0
    return
  endif 
enddo
if(info==1) write(*,*) "ERROR in finding the dual link of", gl

end subroutine find_global_origin_of_dual_global_link


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate lambda.lambda(U.\bar{phi}.U^-1 + \bar{phi})
!subroutine calc_SfL2(SfL2, PhiMat, Umat, Glambda_lambda)
!use global_parameters
!use parallel
!use matrix_functions, only : matrix_3_product, trace_MM
!use SUN_generators, only : trace_MTa, make_Mijkl_from_modes
!implicit none
!
!complex(kind(0d0)), intent(out) :: SfL2
!complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
!complex(kind(0d0)), intent(in) :: Glambda_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_links)
!
!
!complex(kind(0d0)) :: MMAT(1:NMAT,1:NMAT),tmpmat(1:NMAT,1:NMAT)
!!complex(kind(0d0)) :: DP(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
!!complex(kind(0d0)) :: modes(1:dimG,1:dimG)
!complex(kind(0d0)) :: trace, tmpSfL2, tmp
!integer :: gl,ll
!integer :: a,b,i,j,k
!
!SfL2=(0d0,0d0)
!tmpSfL2=(0d0,0d0)
!do ll=1,num_links
!  gl=global_link_of_local(ll)
!  !!!
!  do j=1, NMAT
!    do i=1, NMAT
!      MMat(i,j)=dconjg( PhiMat(j,i,link_org(ll)) )
!    enddo
!  enddo
!  call matrix_3_product(MMat,Umat(:,:,ll),PhiMat(:,:,link_tip(ll)),Umat(:,:,ll),&
!    'N','C','C',(1d0,0d0),'ADD')
!  MMat= MMat * alpha_l(ll)
!  !!!!
!  tmpmat=(0d0,0d0)
!  do j=1,NMAT
!    do i=1,NMAT
!      tmp=(0d0,0d0)
!      do k=1,NMAT
!        tmp = tmp + Glambda_lambda(i,k,k,j,gl,ll) !DP(i,k,k,j)
!      enddo
!      trace = trace + tmp*MMAT(j,i)
!    enddo
!  enddo
!  tmpSfL2=tmpSfL2+trace
!enddo
!
!call MPI_REDUCE(tmpSfL2,SfL2,1,MPI_DOUBLE_COMPLEX, &
!  MPI_SUM,0,MPI_COMM_WORLD,IERR)
!
!if( MYRANK == 0 ) then
!  SfL2 = SfL2 * overall_factor
!endif
!
!end subroutine calc_SfL2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! calculate D^\mu J_\mu for U(1)_V current
!!!  vec1 ~ 1/2 Tr(\lambda(l) \eta(s))
!!!  vec2 ~ Tr(\lambda(l) \chi(f))
!!! and
!!!  DJ1 = d.vec1, DJ2 = d.vec2
!!! where \chi(f) is associated with the link variable Uf
!!! as in the fermionic part of the face action
!subroutine calc_DJ_U1V_org3(DJ1,DJ2,Glambda_eta,Glambda_chi,UMAT)
!use global_parameters
!!use initialization_calcobs
!use parallel
!use global_subroutines, only : syncronize_linkval
!use matrix_functions, only :matrix_product, hermitian_conjugate, make_unit_matrix, matrix_exp, matrix_trace
!implicit none
!
!complex(kind(0d0)), intent(out) :: DJ1(1:num_faces)
!complex(kind(0d0)), intent(out) :: DJ2(1:num_faces)
!complex(kind(0d0)), intent(in) :: Glambda_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_sites) 
!complex(kind(0d0)), intent(in) :: Glambda_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_faces) 
!complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
!
!
!complex(kind(0d0)) :: trvec1(1:num_necessary_links)
!!complex(kind(0d0)) :: vec1(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)) :: U_fs(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: U_sf(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: mat1(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: mat2(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: tmp
!complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
!integer :: ls,ll,lf,gs,gl,gf,kk,dir
!integer :: i,j,k,l,a,b
!
!
!!!  DJ1 ~ 1/2 rot Tr(\lambda(l) \eta(s))
!call make_trV1(trvec1,Glambda_eta)
!!call make_V1(vec1,Glambda_eta)
!DJ1=(0d0,0d0)
!do lf=1,num_faces
!  !call make_unit_matrix(mat1)
!  !tmp=(1d0,0d0)
!  do kk=1,links_in_f(lf)%num_
!    ll=links_in_f(lf)%link_labels_(kk)
!    dir=links_in_f(lf)%link_dirs_(kk)
!    DJ1(lf)=DJ1(lf) + dcmplx(dir)*trvec1(ll)
!  enddo
!  DJ1(lf)=DJ1(lf)*beta_f(lf)
!enddo
!
!!! DJ2 ~ div Tr(\lambda(l) \chi(f))
!DJ2=(0d0,0d0)
!do lf=1,num_faces
!  do kk=1,sites_in_f(lf)%num_
!    ls=sites_in_f(lf)%label_(kk)
!    gs=global_site_of_local(ls)
!    call calc_prodUl_from_n1_to_n2_in_Uf(U_fs,lf,1,kk-1,Umat) ! 4's argument must be kk-1 
!    call hermitian_conjugate(U_sf,U_fs)
!    !! contribution of [link FROM gs]
!    do a=1,global_linktip_from_s(gs)%num_
!      tmp=(0d0,0d0)
!      gl=global_linktip_from_s(gs)%labels_(a)
!      do l=1,NMAT
!        do k=1,NMAT
!          do j=1,NMAT
!            do i=1,NMAT
!              tmp=tmp + Glambda_chi(i,j,k,l,gl,lf)*U_sf(j,k)*U_fs(l,i)
!            enddo
!          enddo
!        enddo
!      enddo
!      DJ2(lf) = DJ2(lf) + tmp * dcmplx(global_alpha_l(gl))/dcmplx( num_faces_in_s(ls) )
!    enddo
!    !! contribution of [link TO gs]
!    do a=1,global_linkorg_to_s(gs)%num_
!      tmp=(0d0,0d0)
!      gl=global_linkorg_to_s(gs)%labels_(a)
!      do ll=1,num_necessary_links
!        if( global_link_of_local(ll) == gl ) exit
!      enddo
!      call matrix_product(mat1,Umat(:,:,ll),U_sf)
!      call matrix_product(mat2,U_fs,Umat(:,:,ll),'N','C')
!      do l=1,NMAT
!        do k=1,NMAT
!          do j=1,NMAT
!            do i=1,NMAT
!              tmp=tmp + Glambda_chi(i,j,k,l,gl,lf)*mat1(j,k)*mat2(l,i)*dcmplx(global_alpha_l(gl))
!            enddo
!          enddo
!        enddo
!      enddo
!    enddo
!    DJ2(lf) = DJ2(lf) - tmp / dcmplx( num_faces_in_s(ls) )
!  enddo
!enddo
!
!DJ1=(DJ1)/dcmplx(LatticeSpacing**4)
!DJ2=(DJ2)/dcmplx(LatticeSpacing**4)
!
!end subroutine calc_DJ_U1V_org3
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! calculate D^\mu J_\mu for U(1)_V current
!!!  vec1 ~ 1/2 Tr(\lambda(l) \eta(s))
!!!  vec2 ~ Tr(\lambda(l) \chi(f))
!!! where \chi(f) is associated with the link variable Uf
!!! as in the fermionic part of the face action
!subroutine calc_divJ_U1Vorg2(divJ1,divJ2,Glambda_eta,Glambda_chi,UMAT)
!use global_parameters
!!use initialization_calcobs
!use parallel
!use global_subroutines, only : syncronize_linkval
!implicit none
!
!
!complex(kind(0d0)), intent(in) :: Glambda_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_sites) 
!complex(kind(0d0)), intent(in) :: Glambda_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_faces) 
!complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
!
!complex(kind(0d0)) :: vec1(1:NMAT,1:NMAT,1:num_necessary_links) ! 1/2 Tr(\lambda(l) \eta(s))
!!complex(kind(0d0)) :: trvec1(1:num_necessary_links) ! 1/2 Tr(\lambda(l) \eta(s))
!complex(kind(0d0)) :: trvec2(1:num_necessary_links) ! Tr(\lambda(l) \chi(f))
!
!complex(kind(0d0)) :: divJ1(1:num_faces)
!complex(kind(0d0)) :: divJ2(1:num_faces)
!
!!!  vec1 ~ 1/2 Tr(\lambda(l) \eta(s))
!call make_V1(vec1,Glambda_eta)
!!call make_trV1(vec1,Glambda_eta)
!!!  vec2 ~ Tr(\lambda(l) \chi(f))
!call make_trV2_likeSf(trvec2,Glambda_chi,UMAT)
!
!call calc_trrot2(divJ1,vec1,Umat)
!!call calc_trrot(divJ1,vec1)
!call calc_trdiv(divJ2,trvec2)
!
!divJ1=(divJ1)/dcmplx(LatticeSpacing**4)
!divJ2=(divJ2)/dcmplx(LatticeSpacing**4)
!
!end subroutine calc_divJ_U1Vorg2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! calculate D^\mu J_\mu for U(1)_V current
!subroutine calc_divJ_U1Vorg(divJ,Glambda_eta,Gchi_lambda,UMAT)
!use global_parameters
!!use initialization_calcobs
!use parallel
!use global_subroutines, only : syncronize_linkval, calc_prodUl_from_n1_to_n2_in_Uf
!implicit none
!
!
!complex(kind(0d0)), intent(out) :: divJ(1:num_faces)
!complex(kind(0d0)), intent(in) :: Glambda_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_sites) 
!complex(kind(0d0)), intent(in) :: Gchi_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_links) 
!complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
!
!complex(kind(0d0)) :: vec1(1:num_necessary_links) ! 1/2 Tr(\lambda(l) \eta(s))
!complex(kind(0d0)) :: vec2(1:num_necessary_links) ! Tr(\lambda(l) \chi(f))
!complex(kind(0d0)) :: divJ1(1:num_faces)
!complex(kind(0d0)) :: divJ2(1:num_faces)
!
!
!call make_trV1(vec1,Glambda_eta)
!call make_trV2(vec2,Gchi_lambda,UMAT)
!
!call calc_trrot(divJ1,vec1)
!call calc_trdiv(divJ2,vec2)
!
!divJ=(divJ1+divJ2)/dcmplx(LatticeSpacing**4)
!
!end subroutine calc_divJ_U1Vorg
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! calculate D^\mu J_\mu for U(1)_V current
!subroutine calc_divJ_U1V_2(divJ1,divJ2,Glambda_eta,Gchi_lambda,UMAT)
!use global_parameters
!!use initialization_calcobs
!use parallel
!use global_subroutines, only : syncronize_linkval, calc_prodUl_from_n1_to_n2_in_Uf
!implicit none
!
!
!complex(kind(0d0)), intent(in) :: Glambda_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_sites) 
!complex(kind(0d0)), intent(in) :: Gchi_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_links) 
!complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
!
!complex(kind(0d0)) :: vec1(1:num_necessary_links) ! 1/2 Tr(\lambda(l) \eta(s))
!complex(kind(0d0)) :: vec2(1:num_necessary_links) ! Tr(\lambda(l) \chi(f))
!
!complex(kind(0d0)) :: divJ1(1:num_faces)
!complex(kind(0d0)) :: divJ2(1:num_faces)
!
!call make_trV1(vec1,Glambda_eta)
!call make_trV2(vec2,Gchi_lambda,UMAT)
!
!call calc_trrot(divJ1,vec1)
!call calc_trdiv2(divJ2,vec2)
!
!divJ1=(divJ1)/dcmplx(LatticeSpacing**4)
!divJ2=(divJ2)/dcmplx(LatticeSpacing**4)
!
!end subroutine calc_divJ_U1V_2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! calculate D^\mu J_\mu for U(1)_V current
!!!  rotation is taken as \sum_{l\in f} Tr( X V_l Y )
!!!  divergence is taken as summation over site-divergence/#{facce on s}
!subroutine calc_divJ_U1V_3(divJ1,divJ2,Glambda_eta,Gchi_lambda,UMAT)
!use global_parameters
!!use initialization_calcobs
!use parallel
!implicit none
!
!
!complex(kind(0d0)), intent(in) :: Glambda_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_sites) 
!complex(kind(0d0)), intent(in) :: Gchi_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_links) 
!complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
!
!complex(kind(0d0)) :: vec1(1:NMAT,1:NMAT,1:num_necessary_links) ! 1/2 (\lambda(l) \eta(s))_{ij}
!complex(kind(0d0)) :: trvec2(1:num_necessary_links) ! Tr(\lambda(l) \chi(f))
!
!complex(kind(0d0)) :: divJ1(1:num_faces)
!complex(kind(0d0)) :: divJ2(1:num_faces)
!
!call make_V1(vec1,Glambda_eta)
!call make_trV2(trvec2,Gchi_lambda,UMAT)
!
!call calc_trrot2(divJ1,vec1,Umat)
!call calc_trdiv(divJ2,trvec2)
!
!divJ1=(divJ1)/dcmplx(LatticeSpacing**4)
!divJ2=(divJ2)/dcmplx(LatticeSpacing**4)
!
!end subroutine calc_divJ_U1V_3
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! calculate lambda.lambda(U.\bar{phi}.U^-1 + \bar{phi})
!subroutine make_V1(vec1,Glambda_eta)
!use global_parameters
!use parallel
!use global_subroutines, only : syncronize_links
!implicit none
!
!complex(kind(0d0)), intent(out) :: vec1(1:NMAT,1:NMAT,1:num_necessary_links) ! 1/2 Tr(\lambda(l) \eta(s))
!complex(kind(0d0)), intent(in) :: Glambda_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_sites) 
!
!integer :: ll,ls
!integer :: gl
!integer :: org_ll
!integer :: i,j,k
!
!
!vec1=(0d0,0d0)
!do ll=1,num_links
!  gl=global_link_of_local(ll)
!  ls=link_org(ll)
!  do i=1,NMAT
!    do j=1,NMAT
!      do k=1,NMAT
!        vec1(i,j,ll)=vec1(i,j,ll) + (0.5d0,0d0)*Glambda_eta(i,k,k,j,gl,ls)
!      enddo
!    enddo
!  enddo
!  !vec1(:,:,ll)=vec1(:,:,ll)*alpha_l(ll)
!enddo
!call syncronize_links(vec1)
!
!
!end subroutine make_V1
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! calculate lambda.lambda(U.\bar{phi}.U^-1 + \bar{phi})
!subroutine make_V2_bak(vec2,Gchi_lambda,UMAT)
!use global_parameters
!use parallel
!use global_subroutines, only : syncronize_linkval, calc_prodUl_from_n1_to_n2_in_Uf
!implicit none
!
!complex(kind(0d0)), intent(out) :: vec2(1:num_necessary_links) ! Tr(\lambda(l) \chi(f))
!complex(kind(0d0)), intent(in) :: Gchi_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_links) 
!complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
!
!complex(kind(0d0)) :: Ucarry(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: tmp(1:num_faces)
!integer :: ll,lf,ls
!integer :: gl,gf
!integer :: org_ll
!integer :: i,j,k,l,ii
!
!vec2=(0d0,0d0)
!do ll=1,num_links
!  do ii=1,face_in_l(ll)%num_
!    lf=face_in_l(ll)%label_(ii)
!    gf=global_face_of_local(lf)
!    !! find the place of ll in the face sharing ll
!    do org_ll=1,links_in_f(lf)%num_
!      if( links_in_f(lf)%link_labels_(org_ll) == ll ) then 
!        exit
!      endif
!    enddo
!    if( links_in_f(lf)%link_dirs_(org_ll) == 1 ) then
!      org_ll = org_ll - 1
!    endif
!    !! Ucarry = U1 ... U_orgll
!    call calc_prodUl_from_n1_to_n2_in_Uf(Ucarry,lf,1,org_ll,Umat)
!    do l=1,NMAT
!      do k=1,NMAT
!        do j=1,NMAT
!          do i=1,NMAT
!            vec2(ll)=vec2(ll)&
!              + Ucarry(j,k) * dconjg(Ucarry(i,l)) &
!                * (0.5d0,0d0) * Gchi_lambda(i,j,k,l,gf,ll)
!          enddo
!        enddo
!      enddo
!    enddo
!  enddo
!enddo
!call syncronize_linkval(vec2)
!end subroutine make_V2_bak
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! inspierd by
!!!  SF^f = 1/g^2 \int ( -2i \chi D\times\lambda )
!!!       = 1/g^2 sum_f rot( V_{f,l}
!!! V_{f,l} = (i\alpha_\beta_f) (1/B(\chi X \lambda Y +...) )
!!! \lambda \chi ~ -1/(-2i) V_{f,l} = -i/2 V_{f,l}
!!! 
!subroutine make_trV2_likeSf(vec2,Glambda_chi,UMAT)
!use global_parameters
!use global_subroutines, only : syncronize_linkval
!use parallel
!implicit none
!
!complex(kind(0d0)), intent(out) :: vec2(1:num_necessary_links) ! Tr(\lambda(l) \chi(f))
!complex(kind(0d0)), intent(in) :: Glambda_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_faces) 
!complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
!
!complex(kind(0d0)) :: Lff
!integer :: ll,lf,l_place,dir
!integer :: gl,gf,rank_send,rank_recv,tag
!integer :: org_ll
!integer :: ii
!
!vec2=(0d0,0d0)
!do gl=1,global_num_links
!  !! find the origin face of the dual link
!  !!  gf: global face label of the origin of the dual link of gl
!  !!      defined so that dir(gl) is the same with the direction of the link in gf
!  do ii=1,global_face_in_l(gl)%num_
!    gf=global_face_in_l(gl)%label_(ii)
!    do org_ll=1,global_links_in_f(gf)%num_
!      if( global_links_in_f(gf)%link_labels_(org_ll) == gl ) then 
!        dir=global_links_in_f(gf)%link_dirs_(org_ll)
!        exit
!      endif
!    enddo
!    !!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!
!    !! special treatment of the present discretization
!    if(gf==1) dir=-dir
!    !!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!
!    if( dir==1 ) exit
!  enddo
!
!  rank_send = local_face_of_global(gf)%rank_
!  rank_recv = local_link_of_global(gl)%rank_
!  tag=global_num_faces*(gf-1) + gl -1
!
!  !! send phase
!  if( MYRANK == rank_send ) then
!    lf=local_face_of_global(gf)%label_
!    !! find the place of gl in lf
!    do l_place=1,links_in_f(lf)%num_
!      ll=links_in_f(lf)%link_labels_(l_place)
!      if( global_link_of_local(ll) == gl ) exit
!    enddo
!    call fermionic_face_lagrangian(Lff,lf,l_place,Glambda_chi,Umat)
!    Lff=Lff*(0d0,1d0)
!
!    if( MYRANK /= rank_recv ) then 
!      call MPI_SEND(Lff,1,MPI_DOUBLE_COMPLEX,rank_recv,tag,MPI_COMM_WORLD,IERR)
!    endif
!  endif
!  !! recv phase
!  if( MYRANK == rank_recv ) then
!    ll=local_link_of_global(gl)%label_
!    if( MYRANK /= rank_send ) then
!      call MPI_RECV(vec2(ll),1,MPI_DOUBLE_COMPLEX,rank_send,tag,MPI_COMM_WORLD,ISTATUS,IERR)
!    else
!      vec2(ll) = Lff
!    endif
!    !!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!
!    !! special treatment of the present discretization
!    if(gf==1) vec2(ll)=-vec2(ll)
!    !!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!
!  endif
!
!  Lff = Lff * (0d0,-0.5d0)
!enddo
!call syncronize_linkval(vec2)
!end subroutine make_trV2_likeSf
!
!
!
!
