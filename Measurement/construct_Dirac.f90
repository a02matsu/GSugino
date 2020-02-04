subroutine construct_Dirac(Dirac,Umat,PhiMat)
use global_parameters
use SUN_generators, only : make_SUN_generators
use global_subroutines, only : syncronize_sites, syncronize_links, syncronize_faces, localmat_to_globalvec
use Dirac_operator, only : Prod_Dirac
use parallel
implicit none


complex(kind(0d0)), intent(out) :: Dirac(:,:)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links) 
complex(kind(0d0)), intent(in) :: PHIMAT(1:NMAT,1:NMAT,1:num_necessary_sites) 

complex(kind(0d0)) :: PFeta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: PFlambda(1:NMAT,1:NMAT,1:num_necessary_links) 
complex(kind(0d0)) :: PFchi(1:NMAT,1:NMAT,1:num_necessary_faces) 
complex(kind(0d0)) :: DPFeta(1:NMAT,1:NMAT,1:num_sites) 
complex(kind(0d0)) :: DPFlambda(1:NMAT,1:NMAT,1:num_links) 
complex(kind(0d0)) :: DPFchi(1:NMAT,1:NMAT,1:num_faces) 
complex(kind(0d0)) :: T(1:NMAT,1:NMAT,1:dimG) 

integer :: gs,gl,gf
integer :: gs1,gl1,gf1
integer :: s,l,f
integer :: i,j,a,b
integer :: rank
integer :: numF

Dirac=(0d0,0d0)
numF= (global_num_sites+global_num_links+global_num_faces)*(NMAT*NMAT-1)

! SU(N) generators
call make_SUN_generators(T,NMAT)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! make Dirac matrix
J=0
do gs=1, global_num_sites
  s=local_site_of_global(gs)%label_
  rank=local_site_of_global(gs)%rank_

  PFlambda=(0d0,0d0)
  PFchi=(0d0,0d0)
  do a=1, dimG
    J=J+1
    PFeta=(0d0,0d0)
    if( MYRANK == rank ) then
      PFeta(:,:,s) = T(:,:,a)
    endif
    call syncronize_sites(PFeta)

    call Prod_Dirac(&
      DPFeta, DPFlambda, DPFchi, &
      PFeta,PFlambda,PFchi, UMAT,PhiMat)
    !!!
    call localmat_to_globalvec(Dirac(:,J),DPFeta,DPFlambda,DPFchi)
  enddo
enddo
!!!!!!!!!!!!!!
do gl=1, global_num_links
  l=local_link_of_global(gl)%label_
  rank=local_link_of_global(gl)%rank_

  PFeta=(0d0,0d0)
  PFchi=(0d0,0d0)
  do a=1, dimG
    J=J+1
    PFlambda=(0d0,0d0)
    if( MYRANK == rank ) then
      PFlambda(:,:,l) = T(:,:,a)
    endif
    call syncronize_links(PFlambda)

    call Prod_Dirac(&
      DPFeta, DPFlambda, DPFchi, &
      PFeta,PFlambda,PFchi, UMAT,PhiMat)
    !!!
    call localmat_to_globalvec(Dirac(:,J),DPFeta,DPFlambda,DPFchi)
  enddo
enddo
!!!!!!!!!!!!!!
do gf=1, global_num_faces
  f=local_face_of_global(gf)%label_
  rank=local_face_of_global(gf)%rank_

  PFeta=(0d0,0d0)
  PFlambda=(0d0,0d0)
  do a=1, dimG
    J=J+1
    PFchi=(0d0,0d0)
    if( MYRANK == rank ) then
      PFchi(:,:,f) = T(:,:,a)
    endif
    call syncronize_faces(PFchi)

    call Prod_Dirac(&
      DPFeta, DPFlambda, DPFchi, &
      PFeta,PFlambda,PFchi, UMAT,PhiMat)
    !!!
    call localmat_to_globalvec(Dirac(:,J),DPFeta,DPFlambda,DPFchi)
  enddo
enddo

end subroutine construct_Dirac
