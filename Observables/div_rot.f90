!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! divergence on a site
!!  div(V)(s) = \sum_{l\in <s,*>}  - \sum_{l\in <*,s>} Tr(V_l) 
subroutine calc_trdiv_in_site(trdiv,trvec)
use global_parameters
use global_subroutines, only : syncronize_siteval
use parallel
implicit none

complex(kind(0d0)), intent(out) :: trdiv(1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: trvec(1:num_necessary_links)

integer :: ls, ll
integer :: i

trdiv=(0d0,0d0)
do ls=1,num_sites
  do i=1,linktip_from_s(ls)%num_
    ll=linktip_from_s(ls)%labels_(i) 
    trdiv(ls)=trdiv(ls) + dcmplx(alpha_l(ll))*trvec(ll)
  enddo
  !!!
  do i=1,linkorg_to_s(ls)%num_
    ll=linkorg_to_s(ls)%labels_(i) 
    trdiv(ls)=trdiv(ls) - dcmplx(alpha_l(ll))*trvec(ll)
  enddo
  !write(*,*) global_site_of_local(ls), dble(trdiv(ls)), dble((0d0,-1d0)*trdiv(ls))
enddo

call syncronize_siteval(trdiv)

end subroutine calc_trdiv_in_site

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Tr( div(V) )
subroutine calc_trdiv(trdiv_f,trvec)
use global_parameters
use parallel
implicit none

complex(kind(0d0)), intent(out) :: trdiv_f(1:num_faces)
complex(kind(0d0)), intent(in) :: trvec(1:num_necessary_links)
complex(kind(0d0)) :: trdiv_s(1:num_necessary_sites)
integer :: lf

complex(kind(0d0)) :: tmp
integer :: ls, ll, gs, gl
integer :: i,j

trdiv_f=(0d0,0d0)
call calc_trdiv_in_site(trdiv_s,trvec)
do lf=1,num_faces
  do i=1,sites_in_f(lf)%num_
    ls=sites_in_f(lf)%label_(i)
    trdiv_f(lf)=trdiv_f(lf) + trdiv_s(ls)/dcmplx(dble(num_faces_in_s(ls)))
    !write(*,*) global_face_of_local(lf), dble(trdiv_s(ls)), dble((0d0,-1d0)*trdiv_s(ls))
  enddo
  trdiv_f(lf) = trdiv_f(lf) * beta_f(lf)
  !write(*,*) global_face_of_local(lf), dble(trdiv_f(lf)), dble((0d0,-1d0)*trdiv_f(lf))
enddo

end subroutine calc_trdiv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Tr( rot(V) )
subroutine calc_trrot(trrot,trvec)
use parallel
implicit none

complex(kind(0d0)), intent(out) :: trrot(1:num_faces)
complex(kind(0d0)), intent(in) :: trvec(1:num_necessary_links)

integer :: ls, ll, lf
integer :: i,j
integer :: dir

trrot=(0d0,0d0)
do lf=1,num_faces
  do i=1,links_in_f(lf)%num_
    ll=links_in_f(lf)%link_labels_(i)
    dir=links_in_f(lf)%link_dirs_(i)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if( global_face_of_local(lf)==1 ) dir=-dir
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    trrot(lf) = trrot(lf) + dcmplx(dble(dir))*trvec(ll)
  enddo
  trrot(lf) = trrot(lf) * beta_f(lf)
  !write(*,*) global_face_of_local(lf), dble(trrot(lf)), dble((0d0,-1d0)*trrot(lf))
enddo

end subroutine calc_trrot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Tr( rot(V) ) by using Xmat and Ymat
subroutine calc_trrot2(trrot,vec,Umat)
use global_parameters
use global_subroutines, only : calc_XYmat
use parallel
implicit none

complex(kind(0d0)), intent(out) :: trrot(1:num_faces)
complex(kind(0d0)), intent(in) :: vec(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)

complex(kind(0d0)) :: Xmat(1:NMAT,1:NMAT), Ymat(1:NMAT,1:NMAT)
integer :: ls, ll, lf
integer :: i,j,k
integer :: dir,l_place

trrot=(0d0,0d0)
do lf=1,num_faces
  do l_place=1,links_in_f(lf)%num_
    ll=links_in_f(lf)%link_labels_(l_place)
    dir=links_in_f(lf)%link_dirs_(l_place)
    call calc_XYmat(Xmat,Ymat,lf,l_place,Umat)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if( global_face_of_local(lf)==1 ) dir=-dir
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1,NMAT
      do j=1,NMAT
        do k=1,NMAT
          trrot(lf) = trrot(lf) &
            + (1.0d0,0d0)*dcmplx(dble(dir)) &
            *Xmat(i,j)*Ymat(k,i)*vec(j,k,ll)
            !+ (0.5d0,0d0)*dcmplx(dble(dir)) &
            !*(Xmat(i,j)*Ymat(k,i)+dconjg(Ymat(j,i)*Xmat(i,k)))*vec(j,k,ll)
        enddo
      enddo
    enddo
  enddo
  trrot(lf) = trrot(lf) * beta_f(lf)
enddo

end subroutine calc_trrot2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Tr( div(V) ) onother version
subroutine calc_trdiv2(trdiv,trvec)
use parallel
implicit none

complex(kind(0d0)), intent(out) :: trdiv(1:num_faces)
complex(kind(0d0)), intent(in) :: trvec(1:num_necessary_links)

complex(kind(0d0)) :: trdiv_s(1:num_necessary_sites)
integer :: ls, ll, lf
integer :: i,j
integer :: l1, l2
integer :: dir1, dir2

call calc_trdiv_in_site(trdiv_s,trvec)
trdiv=(0d0,0d0)
do lf=1,num_faces
  do i=1,sites_in_f(lf)%num_
    ls=sites_in_f(lf)%label_(i)
    trdiv(lf)=trdiv(lf) + trdiv_s(ls)
    !!!!!
    if( i==1 ) then
      l1=links_in_f(lf)%link_labels_(links_in_f(lf)%num_)
      dir1=links_in_f(lf)%link_dirs_(links_in_f(lf)%num_)
    else
      l1=links_in_f(lf)%link_labels_(i-1)
      dir1=links_in_f(lf)%link_dirs_(i-1)
    endif
    l2=links_in_f(lf)%link_labels_(i)
    dir2=links_in_f(lf)%link_dirs_(i)
    !!!!!
    trdiv(lf) = trdiv(lf) &
      + dcmplx(dble(dir1)*alpha_l(l1))*trvec(l1) &
      - dcmplx(dble(dir2)*alpha_l(l2))*trvec(l2)
  enddo
  trdiv(lf) = trdiv(lf) * beta_f(lf)
enddo
end subroutine calc_trdiv2
