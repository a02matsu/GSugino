!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  compute div(K3)
!! *** eta, lambda, chi must be  ***
!! ***  1) initialized before    ***
!! ***  2) syncronized after     ***
subroutine make_divK3(eta,lambda,chi,F12,Dphi,CommPhi,Umat,PhiMat,fcurr)
use parallel
implicit none

complex(kind(0d0)), intent(inout) :: eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(inout) :: lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(inout) :: chi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), intent(in) :: F12(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), intent(in) :: DPhi(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: CommPhi(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: Phimat(1:NMAT,1:NMAT,1:num_necessary_sites)
integer, intent(in) :: fcurr ! global face

integer :: dir, Nfaces
integer :: s,l,f,i,j,k
integer :: rank
integer :: gs,gl,gf
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)

!! divergence part
do i=1, global_sites_in_f(fcurr)%num_
  gs=global_sites_in_f(fcurr)%label_(i)
  Nfaces=global_face_in_s(gs)%num_
  do j=1, global_linktip_from_s(gs)%num_
    gl=global_linktip_from_s(gs)%labels_(j)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    rank=local_link_of_global(gl)%rank_
    l=local_link_of_global(gl)%label_
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    if( MYRANK == rank ) then
      do k=1, face_in_l(l)%num_
        f=face_in_l(l)%label_(k)
        !write(*,*) MYRANK,"test1"
        call carry_MAT(tmpmat,F12(:,:,f),sites_in_f(f)%label_(1),link_org(l),f,Umat)
        lambda(:,:,l)=lambda(:,:,l) &
          + dcmplx(0.5d0 * alpha_l(l) / dble(Nfaces)) * tmpmat
      enddo
    endif
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do j=1, global_linkorg_to_s(gs)%num_
    gl=global_linkorg_to_s(gs)%labels_(j)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    rank=local_link_of_global(gl)%rank_
    l=local_link_of_global(gl)%label_
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    if( MYRANK == rank ) then
      do k=1, face_in_l(l)%num_
        f=face_in_l(l)%label_(k)
        !write(*,*) MYRANK,"test3"
        call carry_MAT(tmpmat,F12(:,:,f),sites_in_f(f)%label_(1),link_org(l),f,Umat)
        lambda(:,:,l)=lambda(:,:,l) &
          - dcmplx(0.5d0 * alpha_l(l) / dble(Nfaces)) * tmpmat
      enddo
    endif
  enddo
enddo

end subroutine make_divK3

