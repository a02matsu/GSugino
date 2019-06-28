!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate D^\mu J_\mu for U(1)_V current
subroutine calc_divJ_U1V(divJ,Glambda_eta,Gchi_lambda,UMAT)
use global_parameters
!use initialization_calcobs
use parallel
use global_subroutines, only : syncronize_linkval, calc_prodUl_from_n1_to_n2_in_Uf
implicit none


complex(kind(0d0)), intent(out) :: divJ(1:num_faces)
complex(kind(0d0)), intent(in) :: Glambda_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_sites) 
!complex(kind(0d0)), intent(in) :: Glambda_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_faces) 
complex(kind(0d0)), intent(in) :: Gchi_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_links) 
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)

complex(kind(0d0)) :: vec1(1:num_necessary_links) ! 1/2 Tr(\lambda(l) \eta(s))
complex(kind(0d0)) :: vec2(1:num_necessary_links) ! Tr(\lambda(l) \chi(f))

!type VEC_IN_FACE
!  integer :: num_
!  complex(kind(0d0)), allocatable :: val(:)
!end type VEC_IN_FACE
!type(VEC_IN_FACE) :: vec2(1:num_faces) ! Tr(\lambda(l) \chi(f))

complex(kind(0d0)) :: divJ1(1:num_faces)
complex(kind(0d0)) :: divJ2(1:num_faces)
complex(kind(0d0)) :: Ucarry(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmp(1:num_faces)
integer :: ll,lf,ls
integer :: gl,gf
integer :: org_ll
integer :: i,j,k,l,ii

vec1=(0d0,0d0)
do ll=1,num_links
  gl=global_link_of_local(ll)
  ls=link_org(ll)
  do j=1,NMAT
    do i=1,NMAT
      vec1(ll)=vec1(ll) + (0.5d0,0d0)*Glambda_eta(i,j,j,i,gl,ls)
    enddo
  enddo
enddo
call syncronize_linkval(vec1)

vec2=(0d0,0d0)
do ll=1,num_links
  !gl=global_link_of_local(ll)
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
              + (0.5d0,0d0) * Ucarry(j,k) * dconjg(Ucarry(i,l)) &
                * Gchi_lambda(i,j,k,l,gf,ll)
                !* Glambda_chi(j,k,l,i,gl,lf)
          enddo
        enddo
      enddo
    enddo
  enddo
enddo
call syncronize_linkval(vec2)

!do lf=1,num_faces
!  !! set size of vec2
!  vec2(lf)%num_=links_in_f(lf)%num_
!  allocate( vec2(lf)%val(1:links_in_f(lf)%num_) )
!  vec2(lf)%val=(0d0,0d0)
!  !!
!  do ii=1,links_in_f(lf)%num_
!    ll=links_in_f(lf)%link_labels_(ii)
!    dir=links_in_f(lf)%link_dirs_(ii)
!    gl=global_link_of_local(ll)
!    if( dir == 1 ) then 
!      s_place=ii
!    else
!      s_place=ii+1
!      if(s_place == links_in_f(lf)%num_ + 1) s_place=1
!    endif
!    call calc_prodUl_from_n1_to_n2_in_Uf(Ucarry,lf,1,s_place-1,Umat)
!    do i=1,NMAT
!      do l=1,NMAT
!        do j=1,NMAT
!          do k=1,NMAT
!            vec2(lf)%val(ii)=vec2(lf)%val(ii) &
!              + Ucarry(i,j) * dconjg(Ucarry(l,k)) * Glambda_chi(j,k,l,i,gl,lf)
!          enddo
!        enddo
!      enddo
!    enddo
!  enddo
!enddo

divJ1=(0d0,0d0)
divJ2=(0d0,0d0)
call calc_trrot(divJ1,vec1)
call calc_trdiv(divJ2,vec2)

divJ=(divJ1+divJ2)/dcmplx(LatticeSpacing**3)

end subroutine calc_divJ_U1V

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

