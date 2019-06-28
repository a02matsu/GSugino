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
    trdiv(ls)=trdiv(ls) + alpha_l(ll)*trvec(ll)
  enddo
  !!!
  do i=1,linkorg_to_s(ls)%num_
    ll=linkorg_to_s(ls)%labels_(i) 
    trdiv(ls)=trdiv(ls) - alpha_l(ll)*trvec(ll)
  enddo
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
  enddo
  trdiv_f(lf) = trdiv_f(lf) * beta_f(lf)
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
    trrot(lf) = trrot(lf) + dcmplx(dble(dir))*trvec(ll)
  enddo
  trrot(lf) = trrot(lf) * beta_f(lf)
enddo

end subroutine calc_trrot



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Tr( rot(V) ) at face f
!! trvec is vectors defined in the face lf
subroutine calc_trrot_in_face(trrot,trvec,num,lf)
use parallel
implicit none

complex(kind(0d0)), intent(out) :: trrot
complex(kind(0d0)), intent(in) :: trvec(1:num)
integer, intent(in) :: num ! number of links in the face lf
integer, intent(in) :: lf ! label of local face

integer :: ls, ll
integer :: i,j
integer :: dir

trrot=(0d0,0d0)
do i=1,links_in_f(lf)%num_
  !ll=links_in_f(lf)%link_labels_(i)
  dir=links_in_f(lf)%link_dirs_(i)
  trrot = trrot + dcmplx(dble(dir))*trvec(i)
enddo
trrot = trrot * beta_f(lf)

end subroutine calc_trrot_in_face

