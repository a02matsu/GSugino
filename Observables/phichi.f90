!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  compute rotK1
!! *** eta, lambda, chi must be  ***
!! ***  1) initialized before    ***
!! ***  2) syncronized after     ***
subroutine make_phichi(eta,lambda,chi,F12,Dphi,CommPhi,Umat,PhiMat,fcurr)
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

!! mass part
rank=local_face_of_global(fcurr)%rank_
f=local_face_of_global(fcurr)%label_
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if( MYRANK == rank ) then
  s=sites_in_f(f)%label_(1)
  chi(:,:,f)=chi(:,:,f) + dcmplx(mass_square_phi) * PhiMat(:,:,s)
endif

end subroutine make_phichi

