!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compute div(K4)
!! *** eta, lambda, chi must be  ***
!! ***  1) initialized before    ***
!! ***  2) syncronized after     ***
subroutine make_divK4(eta,lambda,chi,F12,Dphi,CommPhi,Umat,PhiMat,fcurr)
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
    do k=1, global_face_in_l(gl)%num_
      gf=global_face_in_l(gl)%label_(k)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      rank=local_face_of_global(gf)%rank_
      f=local_face_of_global(gf)%label_
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if( MYRANK == rank ) then
        do l=1,num_necessary_links
          if( global_link_of_local(l) == gl ) exit
        enddo
        call carry_MAT(tmpmat,DPhi(:,:,l),link_org(l),sites_in_f(f)%label_(1),f,Umat)
        chi(:,:,f)=chi(:,:,f) &
          - dcmplx( 0.5d0 * alpha_l(l) / dble(Nfaces) ) * tmpmat
      endif
    enddo
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do j=1, global_linkorg_to_s(gs)%num_
    gl=global_linkorg_to_s(gs)%labels_(j)
    do k=1, global_face_in_l(gl)%num_
      gf=global_face_in_l(gl)%label_(k)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      rank=local_face_of_global(gf)%rank_
      f=local_face_of_global(gf)%label_
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if( MYRANK == rank ) then
        do l=1,num_necessary_links
          if( global_link_of_local(l) == gl ) exit
        enddo
        !write(*,*) MYRANK,"test4"
        call carry_MAT(tmpmat,DPhi(:,:,l),link_org(l),sites_in_f(f)%label_(1),f,Umat)
        chi(:,:,f)=chi(:,:,f) &
          + dcmplx( 0.5d0 * alpha_l(l) / dble(Nfaces) ) * tmpmat
      endif
    enddo
  enddo
enddo

end subroutine make_divK4


