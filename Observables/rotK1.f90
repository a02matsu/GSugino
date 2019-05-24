!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  compute rotK1
!! *** eta, lambda, chi must be  ***
!! ***  1) initialized before    ***
!! ***  2) syncronized after     ***
subroutine make_rotK1(eta,lambda,chi,F12,Dphi,CommPhi,Umat,PhiMat,fcurr)
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

!! rotation part
do i=1, global_links_in_f(fcurr)%num_
  gl=global_links_in_f(fcurr)%link_labels_(i)
  dir=global_links_in_f(fcurr)%link_dirs_(i)
  if( dir==1 ) then 
    gs=global_link_org(gl)
  else
    gs=global_link_tip(gl)
  endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  rank=local_site_of_global(gs)%rank_
  s=local_site_of_global(gs)%label_
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if( MYRANK == rank ) then 
    do l=1,num_necessary_links
      if( global_link_of_local(l) == gl ) exit
    enddo
    do f=1,num_necessary_faces
      if( global_face_of_local(f) == fcurr ) exit
    enddo
    eta(:,:,s)=eta(:,:,s) &
      + dcmplx( 0.5d0 * beta_f(f) * dble(dir) ) * DPhi(:,:,l)
  endif
enddo



end subroutine make_rotK1

