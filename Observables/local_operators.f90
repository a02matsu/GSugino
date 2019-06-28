!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_trphi2(trphi2,PhiMat)
use parallel
implicit none

double precision, intent(out) :: trphi2(1:num_sites)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) :: tmp
integer :: s,i,j

trphi2=(0d0,0d0)
do s=1,num_sites
  do j=1,NMAT
    do i=1,NMAT
      trphi2(s)=trphi2(s)+dble(PhiMat(i,j,s)*dconjg(PhiMat(i,j,s)))
    enddo
  enddo
enddo
trphi2=trphi2/(LatticeSpacing*LatticeSpacing*dble(NMAT))

end subroutine calc_trphi2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! trF2 = 1/N Tr( f^2 )*A_f
!! where
!!    f=1/2 E^{\mu\nu} F_{\mu\nu} = (Uf-Uf^\dagger) / 2ia^2 A_f
subroutine calc_trF2(trF2,UMat)
use parallel
use global_subroutines, only : make_face_variable
implicit none

double precision, intent(out) :: trF2(1:num_faces)
complex(kind(0d0)), intent(in) :: UMat(1:NMAT,1:NMAT,1:num_necessary_faces)

complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT), Fmat(1:NMAT,1:NMAT)
integer :: f,i,j

trF2=(0d0,0d0)
do f=1,num_faces
  call make_face_variable(Uf,f,UMAT)
  do j=1,NMAT
    do i=1,NMAT
      Fmat(i,j)= (0d0,-1d0)*(Uf(i,j)-dconjg(Uf(j,i)))
    enddo
  enddo
  !Fmat=Fmat * (0d0,-0.5d0) / dcmplx( LatticeSpacing*LatticeSpacing * alpha_f(f))

  do j=1,NMAT
    do i=1,NMAT
      trF2(f)=trF2(f) + dble(FMat(i,j)*dconjg(FMat(i,j)))
    enddo
  enddo
  trF2(f) = trF2(f) / ( alpha_f(f)**2 )
enddo
trF2 = trF2/(4d0 * LatticeSpacing**3 * dble(NMAT))

end subroutine calc_trF2


