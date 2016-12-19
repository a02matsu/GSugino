!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! module to construct sevral quantities used in simulation 
!!
!! DO NOT CHANGE THE VALUES OF VARIABLES IN GLOBASL_PARAMETERS
module global_subroutines
use global_parameters
!use simplicial_complex
implicit none
contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Gaussian random number
!!  BoxMuller(gauss,N)
!!
!! generage ensemble exp(-1/2(x^2) and exp(-1/2(y^2))
!! 
!! output 2N gaussian randum numbers
!! gauss is an allocatable array.
!! It will be reallocated after calling this routine.
subroutine BoxMuller(gauss,N)
use mt95
implicit none

double precision, parameter :: PI=dacos(-1d0)
integer, intent(in) :: N
double precision, allocatable :: gauss(:)
double precision rand(1:2*N)
integer i

if( allocated(gauss) ) then 
    deallocate( gauss )
endif

allocate( gauss(1:2*N) )

call genrand_real3(rand)
!write(*,*) rand

do i=1,N
  gauss(2*i-1) = dsqrt(-2d0*dlog(rand(2*i-1)))*dsin(2d0*Pi*rand(2*i))
  gauss(2*i) = dsqrt(-2d0*dlog(rand(2*i-1)))*dcos(2d0*Pi*rand(2*i))
enddo
end subroutine BoxMuller


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Gaussian random number
!!  BoxMuller(gauss,N)
!!
!! output 2N gaussian randum numbers
subroutine BoxMuller2(gauss,N)
use mt95
implicit none

double precision, parameter :: PI=dacos(-1d0)
integer, intent(in) :: N
double precision :: gauss(1:2*N)
double precision :: rand(1:2*N)
integer i

call genrand_real3(rand)

do i=1,N
  gauss(2*i-1) = dsqrt(-2d0*dlog(rand(2*i-1)))*dsin(2d0*Pi*rand(2*i))
  gauss(2*i) = dsqrt(-2d0*dlog(rand(2*i-1)))*dcos(2d0*Pi*rand(2*i))
enddo
end subroutine BoxMuller2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute forward covariant difference of Phi
!!  DMAT = U_l Phi_tip(l) U_l^\dagger - Phi_org(l)
subroutine Make_diff_Phi(Dphi, l,UMAT,Phi)
use SUN_generators, only : Make_traceless_matrix_from_modes
implicit none

complex(kind(0d0)), intent(out) :: Dphi(1:NMAT,1:NMAT)
integer, intent(in) :: l
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
integer :: i,j
complex(kind(0d0)) :: Phi_tip(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: Phi_org(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)

call Make_traceless_matrix_from_modes(Phi_tip,NMAT,Phi(:,link_tip(l)))
!call Make_traceless_matrix_from_modes(Phi_org,NMAT,Phi(:,link_org(l)))
call Make_traceless_matrix_from_modes(Dphi,NMAT,Phi(:,link_org(l)))

! U_l.Phi_tip
call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
    UMAT(:,:,l), NMAT, &
    Phi_tip, NMAT, &
    (0d0,0d0), tmpmat1, NMAT)
! U_l.Phi_tip.U_l^\dagger
call ZGEMM('N','C',NMAT,NMAT,NMAT,(1d0,0d0), &
    tmpmat1, NMAT, &
    UMAT(:,:,l), NMAT, &
    (-1d0,0d0), Dphi, NMAT)
! U_l.Phi_tip.U_l^\dagger - Phi_org
!Dphi = Dphi - Phi_org

end subroutine Make_diff_Phi


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to make moment map \Omega(Uf)
subroutine Make_moment_map(Omega,Ufm)
use matrix_functions, only : matrix_inverse
implicit none

complex(kind(0d0)), intent(out) :: Omega(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: Ufm(1:NMAT,1:NMAT)
complex(kind(0d0)) :: SMAT(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Cinv(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
integer :: i,j,k

!call matrix_power(NMAT,Uf,m_omega,Ufm)

do i=1,NMAT
  do j=1,NMAT
    SMAT(i,j)=-im_unit*( Ufm(i,j) - dconjg( Ufm(j,i) ) )
    Cinv(i,j)=( Ufm(i,j) + dconjg( Ufm(j,i) ) )
  enddo
enddo
! CMAT --> CMAT^{-1}
call matrix_inverse(Cinv)

!call ZGEMM('N','N',NMAT,NMAT,NMAT,dcmplx(1d0/dble(m_omega)), &
call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0),&
  SMAT,NMAT,&
  Cinv,NMAT,&
  (0d0,0d0),tmpmat,NMAT)

do i=1,NMAT
  do j=1,NMAT
    Omega(i,j)=tmpmat(i,j)+dconjg(tmpmat(j,i))
  enddo
enddo

Omega = Omega / dcmplx( dble(m_omega) )

end subroutine Make_moment_map


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to make plaquette variable of face f
subroutine Make_face_variable(Uf,f,UMAT)
!use simplicial_complex, only : get_links_in_face_sc
implicit none

complex(kind(0d0)), intent(out) :: Uf(1:NMAT,1:NMAT)
integer, intent(in) :: f
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
!integer :: FaceSize
!integer, allocatable :: sites(:),link_labels(:),link_dirs(:)
character(1) :: char1
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
integer :: l,i,j


if(links_in_f(f)%link_dirs_(1)==1) then
  Uf=UMAT(:,:,links_in_f(f)%link_labels_(1))
elseif(links_in_f(f)%link_dirs_(1)==-1) then
  do i=1,NMAT
  do j=1,NMAT
    UF(i,j)=dconjg(UMAT(j,i,links_in_f(f)%link_labels_(1)))
  enddo
  enddo
else
  write(*,*) "There is a blanck link in the face",f
  stop 1
endif
!
do l=2,links_in_f(f)%num_
  tmpmat=Uf
  if(links_in_f(f)%link_dirs_(l)==1) then 
    char1='N'
  elseif(links_in_f(f)%link_dirs_(l)==-1) then
    char1='C'
  else
    write(*,*) "There is a blanck link in the face",f
    stop
  endif
  !! Uf=tmpmat*Uf or tmpmat*(Uf)^\dagger
  call ZGEMM('N',char1,NMAT,NMAT,NMAT,(1d0,0d0), &
    tmpmat,NMAT, &
    UMAT(:,:,links_in_f(f)%link_labels_(l)),NMAT, &
    (0d0,0d0), Uf, NMAT)
enddo

end subroutine Make_face_variable


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to compute \frac{d}{d A_l^a} \Omega_m
!!  Uf and Uf^m must be iput matrices.
subroutine calc_diff_omega(diff_Omega, Uf, Ufm, f, l,UMAT)
use matrix_functions, only : matrix_inverse
implicit none

complex(kind(0d0)), intent(out) :: diff_Omega(1:NMAT,1:NMAT,1:dimG)
integer, intent(in) :: f, l
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT), Ufm(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: diff_S(1:NMAT,1:NMAT,1:dimG), diff_C(1:NMAT,1:NMAT,1:dimG)
complex(kind(0d0)) :: SMAT(1:NMAT,1:NMAT), Cinv(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: SCinv(1:NMAT,1:NMAT)
integer :: i,j,a


do i=1,NMAT
  do j=1,NMAT
    SMAT(i,j)=-im_unit*( Ufm(i,j) - dconjg( Ufm(j,i) ) )
    Cinv(i,j)=( Ufm(i,j) + dconjg( Ufm(j,i) ) )
  enddo
enddo
call matrix_inverse(Cinv)

call calc_diff_SandC(diff_S,diff_C,Uf,f,l,UMAT)
! S.C^{-1}
  call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
    SMAT,NMAT, &
    Cinv,NMAT, &
    (0d0,0d0), SCinv, NMAT)

do a=1,dimG
  ! tmpmat3 = diff_S . C^{-1}
  call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
    diff_S(:,:,a),NMAT, &
    Cinv,NMAT, &
    (0d0,0d0), tmpmat, NMAT)
  
  ! tmpmat2 = S.C^{-1}.diff_C
  call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
    SCinv,NMAT, &
    diff_C(:,:,a),NMAT, &
    (0d0,0d0), tmpmat2, NMAT)
  ! tmpmat3 = diff_S.C^{-1} - S.C^{-1}.diff_C.C^{-1}
  call ZGEMM('N','N',NMAT,NMAT,NMAT,(-1d0,0d0), &
    tmpmat2,NMAT, &
    Cinv,NMAT, &
    (1d0,0d0), tmpmat, NMAT)

  ! diff_Omega = 1/m_omega * ( tmpmat3 + (tmpmat3)^\dagger ) 
  do i=1,NMAT
    do j=1,NMAT
      diff_Omega(i,j,a)= &
        ( tmpmat(i,j) + dconjg(tmpmat(j,i) ) ) ! / dcmplx(dble(m_omega))
    enddo
  enddo
enddo
diff_Omega=diff_Omega / dcmplx(dble(m_omega)) 

end subroutine calc_diff_omega


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to compute \frac{d}{d A_l^a} S and \frac{d}{d A_l^a} C
subroutine calc_diff_SandC(diff_S,diff_C, Uf, f, l,UMAT)

complex(kind(0d0)), intent(out) :: diff_S(1:NMAT,1:NMAT,1:dimG)
complex(kind(0d0)), intent(out) :: diff_C(1:NMAT,1:NMAT,1:dimG)
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
integer, intent(in) :: f, l

complex(kind(0d0)) :: diff_Ufm(1:NMAT,1:NMAT,1:dimG)
integer a,i,j

call calc_diff_Ufm(diff_Ufm, Uf,f,l,UMAT)
  
diff_S=(0d0,0d0)
diff_C=(0d0,0d0)
do a=1,dimG
  do j=1,NMAT
    do i=1,NMAT
      diff_S(i,j,a)=(0d0,-1d0)*( diff_Ufm(i,j,a) - dconjg( diff_Ufm(j,i,a) ) )
      !diff_S(i,j,a) = diff_Ufm(i,j,a) + dconjg( diff_Ufm(j,i,a) ) 
      diff_C(i,j,a) = diff_Ufm(i,j,a) + dconjg( diff_Ufm(j,i,a) )
    enddo
  enddo
enddo

end subroutine calc_diff_SandC



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to compute d/dA(a,l) Ufm
subroutine calc_diff_Ufm(diff_Ufm, Uf,f, l,UMAT)
use matrix_functions, only : matrix_power
implicit none

complex(kind(0d0)), intent(out) :: diff_Ufm(1:NMAT,1:NMAT,1:dimG)
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: XMAT(1:NMAT,1:NMAT), YMAT(1:NMAT,1:NMAT)
integer, intent(in) :: l,f

integer :: nl
complex(kind(0d0)) :: Ufla(1:NMAT,1:NMAT,1:dimG)
complex(kind(0d0)) :: Uf_power1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Uf_power2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
character(1) :: char1

integer :: i,j,k,ll,a

!HERE
call calc_diff_Uf(Ufla, Uf,f, l,UMAT)

if ( m_omega == 1 ) then 
  diff_Ufm=Ufla
  return
endif

!!! In the following, we can assume m_omega >= 2.
do a=1,dimG
! k=0
  !call matrix_power(NMAT,Uf,m_omega-1,Uf_power1)
  call matrix_power(Uf_power1,Uf,m_omega-1)
  call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
    Ufla(:,:,a),NMAT, &
    Uf_power1,NMAT, &
    (0d0,0d0), diff_Ufm(:,:,a), NMAT)
enddo
! k=1,...,m_omega-2
if ( m_omega >= 3 ) then
  do k=1, m_omega-2
    call matrix_power(Uf_power1,Uf,k)
    call matrix_power(Uf_power2,Uf,m_omega-k-1)
    do a=1,dimG
      call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
        Uf_power1,NMAT, &
        Ufla(:,:,a),NMAT, &
        (0d0,0d0), tmpmat, NMAT)
      call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
        tmpmat,NMAT, &
        Uf_power2,NMAT, &
        (1d0,0d0), diff_Ufm(:,:,a), NMAT)
    enddo
  enddo
endif
! k=m_omega-1
call matrix_power(Uf_power1,Uf,m_omega-1)
do a=1,dimG
  call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
    Uf_power1,NMAT, &
    Ufla(:,:,a),NMAT, &
    (1d0,0d0), diff_Ufm(:,:,a), NMAT)
enddo

end subroutine calc_diff_Ufm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to compute d/dA(a,l) Uf
subroutine calc_diff_Uf(diff_Uf, Uf,f, l,UMAT)
implicit none

complex(kind(0d0)), intent(out) :: diff_Uf(1:NMAT,1:NMAT,1:dimG)
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
integer, intent(in) :: l,f

complex(kind(0d0)) :: XMAT(1:NMAT,1:NMAT)
complex(kind(0d0)) :: YMAT(1:NMAT,1:NMAT)
complex(kind(0d0)) :: diff_UlAl(1:NMAT,1:NMAT,1:dimG)

integer :: a,link_place
complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)


! find the place of l1 and l2
do link_place=1,links_in_f(f)%num_
  if( links_in_f(f)%link_labels_(link_place) == l ) then 
    exit
  endif
enddo

! d/dA(a1,l1) d/dA(a2,l2) Ul
call calc_diff_Ul_in_face(diff_UlAl,f,link_place,UMAT)
! X=U(1)...U(l-1)
call calc_prodUl_from_n1_to_n2_in_Uf(XMAT,f,1,link_place-1,UMAT)
! Z=U(l+1)...U(f_size)
call calc_prodUl_from_n1_to_n2_in_Uf(YMAT,f,link_place+1,links_in_f(f)%num_,UMAT)
! X.diffdiff.Y
do a=1,dimG
  call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
     XMAT,NMAT, &
     diff_UlAl(:,:,a),NMAT, &
     (0d0,0d0), tmpmat1, NMAT)
  call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
     tmpmat1,NMAT, &
     YMAT,NMAT, &
     (0d0,0d0), diff_Uf(:,:,a), NMAT)
enddo

end subroutine calc_diff_Uf




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to compute 
!!   d/dA(a1,l1) d/dA(a2,l2) \Omega_m
!!  Uf and Uf^m must be iput matrices.
!!  This is aranged in order to reduce the amount of computation.
subroutine calc_diffdiff_omega(diffdiff_Omega, Uf, Ufm, f, l1,l2, UMAT)
use matrix_functions, only : hermitianmatrix_inverse
implicit none

complex(kind(0d0)), intent(out) :: diffdiff_Omega(1:NMAT,1:NMAT,1:dimG,1:dimG)
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT), Ufm(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_faces)
integer, intent(in) :: f, l1,l2

complex(kind(0d0)) :: SMAT(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Cinv(1:NMAT,1:NMAT)
complex(kind(0d0)) :: SCinv(1:NMAT,1:NMAT)
complex(kind(0d0)) :: diff_S_l1(1:NMAT,1:NMAT,1:dimG)
complex(kind(0d0)) :: diff_C_l1(1:NMAT,1:NMAT,1:dimG)
complex(kind(0d0)) :: diff_S_l2(1:NMAT,1:NMAT,1:dimG)
complex(kind(0d0)) :: diff_C_l2(1:NMAT,1:NMAT,1:dimG)
complex(kind(0d0)) :: diffdiff_S(1:NMAT,1:NMAT,1:dimG,1:dimG)
complex(kind(0d0)) :: diffdiff_C(1:NMAT,1:NMAT,1:dimG,1:dimG)
integer :: i,j,a1,a2

complex(kind(0d0)) :: tmp_dd_Omega(1:NMAT,1:NMAT,1:dimG,1:dimG)
complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat3(1:NMAT,1:NMAT)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! preparation
! SMAT and Cinv
do i=1,NMAT
  do j=1,NMAT
    SMAT(i,j)=-im_unit*( Ufm(i,j) - dconjg( Ufm(j,i) ) )
    Cinv(i,j)=( Ufm(i,j) + dconjg( Ufm(j,i) ) )
  enddo
enddo
!call matrix_inverse(NMAT,Cinv)
call hermitianmatrix_inverse(Cinv)

! diff_S_l1, diff_C_l1, diff_S_l2, diff_C_l2
call calc_diff_SandC(diff_S_l1,diff_C_l1,Uf,f,l1,UMAT)
call calc_diff_SandC(diff_S_l2,diff_C_l2,Uf,f,l2,UMAT)

! diffdiff_S and diffdiff_C
call calc_diffdiff_SandC(diffdiff_S, diffdiff_C, Uf, f, l1, l2, UMAT)

! S.Cinv
call ZHEMM('L','U',NMAT,NMAT,(1d0,0d0), &
  SMAT, NMAT, &
  Cinv, NMAT, &
  (0d0,0d0), SCinv, NMAT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! d/dA(a1,l1) d/dA(a2,l2) S . Cinv
do a1=1,dimG
  do a2=1,dimG
    call ZHEMM('R','U',NMAT,NMAT,(1d0,0d0), &
      Cinv, NMAT, &
      diffdiff_S(:,:,a1,a2), NMAT, &
      (0d0,0d0), tmp_dd_Omega(:,:,a1,a2), NMAT)
    !call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
    !  diffdiff_S(:,:,a1,a2), NMAT, &
    !  Cinv, NMAT, &
    !  (0d0,0d0), tmp_dd_Omega(:,:,a1,a2), NMAT)
  enddo
enddo

! add -S.Cinv. (d/dA(a1,l1) d/dA(a2,l2) C). Cinv
do a1=1,dimG
  do a2=1,dimG
    call ZHEMM('R','U',NMAT,NMAT,(1d0,0d0), &
      Cinv, NMAT, &
      diffdiff_C(:,:,a1,a2), NMAT, &
      (0d0,0d0), tmpmat2, NMAT)
    !call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
    !  diffdiff_C(:,:,a1,a2), NMAT, &
    !  Cinv, NMAT, &
    !  (0d0,0d0), tmpmat2, NMAT)
    !!!
    call ZGEMM('N','N',NMAT,NMAT,NMAT,(-1d0,0d0), &
      SCinv, NMAT, &
      tmpmat2, NMAT, &
      (1d0,0d0), tmp_dd_Omega(:,:,a1,a2), NMAT)
  enddo
enddo

do a1=1,dimG
  do a2=1,dimG
    ! d/dA(a1,l1)SMAT . Cinv
    call ZHEMM('R','U',NMAT,NMAT,(1d0,0d0), &
      Cinv, NMAT, &
      diff_S_l1(:,:,a1), NMAT, &
      (0d0,0d0), tmpmat1, NMAT)
    !call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
    !  diff_S_l1(:,:,a1), NMAT, &
    !  Cinv, NMAT, &
    !  (0d0,0d0), tmpmat1, NMAT)
    ! d/dA(a2,l2)CMAT . Cinv
    call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
      diff_C_l2(:,:,a2), NMAT, &
      Cinv, NMAT, &
      (0d0,0d0), tmpmat2, NMAT)
    ! 
    call ZGEMM('N','N',NMAT,NMAT,NMAT,(-1d0,0d0), &
      tmpmat1, NMAT, &
      tmpmat2, NMAT, &
      (1d0,0d0), tmp_dd_Omega(:,:,a1,a2), NMAT)
    !!!!!!!!!!!!!!!!
    ! d/dA(a2,l2)SMAT . Cinv
    call ZHEMM('R','U',NMAT,NMAT,(1d0,0d0), &
      Cinv, NMAT, &
      diff_S_l2(:,:,a2), NMAT, &
      (0d0,0d0), tmpmat1, NMAT)
    ! d/dA(a1,l1)CMAT . Cinv
    call ZHEMM('R','U',NMAT,NMAT,(1d0,0d0), &
      Cinv, NMAT, &
      diff_C_l1(:,:,a1), NMAT, &
      (0d0,0d0), tmpmat2, NMAT)
    ! 
    call ZGEMM('N','N',NMAT,NMAT,NMAT,(-1d0,0d0), &
      tmpmat1, NMAT, &
      tmpmat2, NMAT, &
      (1d0,0d0), tmp_dd_Omega(:,:,a1,a2), NMAT)
  enddo
enddo

do a1=1,dimG
  do a2=1,dimG
    ! tmpmat1 = S.Cinv.dC/dA(a1,l1).Cinv
    call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
      SCinv, NMAT, &
      diff_C_l1(:,:,a1), NMAT, &
      (0d0,0d0), tmpmat2, NMAT)
    call ZHEMM('R','U',NMAT,NMAT,(1d0,0d0), &
      Cinv, NMAT, &
      tmpmat2, NMAT, &
      (0d0,0d0), tmpmat1, NMAT)
    ! tmpmat2 = dC/dA(a2,l2).Cinv
    call ZHEMM('R','U',NMAT,NMAT,(1d0,0d0), &
      Cinv, NMAT, &
      diff_C_l2(:,:,a2), NMAT, &
      (0d0,0d0), tmpmat2, NMAT)
    ! add S.Cinv.dC/dA(a1,l1).Cinv.dC/dA(a2,l2).Cinv
    call ZGEMM('N','N',NMAT,NMAT,NMAT,(+1d0,0d0), &
      tmpmat1, NMAT, &
      tmpmat2, NMAT, &
      (1d0,0d0), tmp_dd_Omega(:,:,a1,a2), NMAT)
    !!!!!!!!!!!!!
    ! tmpmat1 = S.Cinv.dC/dA(a2,l2).Cinv
    call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
      SCinv, NMAT, &
      diff_C_l2(:,:,a2), NMAT, &
      (0d0,0d0), tmpmat2, NMAT)
    call ZHEMM('R','U',NMAT,NMAT,(1d0,0d0), &
      Cinv, NMAT, &
      tmpmat2, NMAT, &
      (0d0,0d0), tmpmat1, NMAT)
    ! tmpmat2 = dC/dA(a1,l1).Cinv
    call ZHEMM('R','U',NMAT,NMAT,(1d0,0d0), &
      Cinv, NMAT, &
      diff_C_l1(:,:,a1), NMAT, &
      (0d0,0d0), tmpmat2, NMAT)
    ! add S.Cinv.dC/dA(a2,l2).Cinv.dC/dA(a1,l1).Cinv
    call ZGEMM('N','N',NMAT,NMAT,NMAT,(+1d0,0d0), &
      tmpmat1, NMAT, &
      tmpmat2, NMAT, &
      (1d0,0d0), tmp_dd_Omega(:,:,a1,a2), NMAT)
  enddo
enddo


do a2=1,dimG
  do a1=1,dimG
    do j=1,NMAT
      do i=1,NMAT
        diffdiff_Omega(i,j,a1,a2)=&
          ( tmp_dd_Omega(i,j,a1,a2) + dconjg( tmp_dd_Omega(j,i,a1,a2) ) ) &
           / dcmplx(dble(m_omega))
      enddo
    enddo
  enddo
enddo

end subroutine calc_diffdiff_omega

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to compute \frac{d}{d A_l^a} S and \frac{d}{d A_l^a} C
subroutine calc_diffdiff_SandC(diffdiff_S,diffdiff_C, Uf,f,l1,l2,UMAT)

complex(kind(0d0)), intent(out) :: diffdiff_S(1:NMAT,1:NMAT,1:dimG,1:dimG)
complex(kind(0d0)), intent(out) :: diffdiff_C(1:NMAT,1:NMAT,1:dimG,1:dimG)
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
integer, intent(in) :: f,l1,l2

complex(kind(0d0)) :: diffdiff_Ufm(1:NMAT,1:NMAT,1:dimG,1:dimG)
integer :: i,j,a1,a2

call calc_diffdiff_Ufm(diffdiff_Ufm,Uf,f,l1,l2,UMAT)

do a2=1,dimG
  do a1=1,dimG
    do j=1,NMAT
      do i=1,NMAT
        diffdiff_S(i,j,a1,a2) = (0d0,-1d0)  &
          * ( diffdiff_Ufm(i,j,a1,a2) - dconjg( diffdiff_Ufm(j,i,a1,a2) ) )
        !diffdiff_S(i,j,a1,a2) = &
        !  diffdiff_Ufm(i,j,a1,a2) + dconjg( diffdiff_Ufm(j,i,a1,a2) ) 
        diffdiff_C(i,j,a1,a2) = &
          diffdiff_Ufm(i,j,a1,a2) + dconjg( diffdiff_Ufm(j,i,a1,a2) ) 
      enddo
    enddo
  enddo
enddo


!do i=1,NMAT
!  do j=1,NMAT
!    do a1=1,dimG
!      do a2=1,dimG
!        write(*,*) diffdiff_S(i,j,a1,a2) - dconjg( diffdiff_S(j,i,a1,a2) ) 
!    enddo
!enddo
!enddo
!enddo

end subroutine calc_diffdiff_SandC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to compute d/dA(a2,l2) d/dA(a1,l1) Ufm
!! <notation>
!!   diffdiff_Ufm(:,;,a1,a2) = d/dA(a1,l1) d/dA(a2,l2) Ufm
subroutine calc_diffdiff_Ufm(diffdiff_Ufm, Uf,f,l1,l2,UMAT)
use matrix_functions, only : matrix_power
implicit none

complex(kind(0d0)), intent(out) :: diffdiff_Ufm(1:NMAT,1:NMAT,1:dimG,1:dimG)
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
integer, intent(in) :: f,l1,l2

complex(kind(0d0)) :: diffdiff_Uf(1:NMAT,1:NMAT,1:dimG,1:dimG)
complex(kind(0d0)) :: diff_Uf_l1(1:NMAT,1:NMAT,1:dimG)
complex(kind(0d0)) :: diff_Uf_l2(1:NMAT,1:NMAT,1:dimG)
complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT),tmpmat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat3(1:NMAT,1:NMAT)
integer :: j,k
integer :: a1,a2

!! diffdiff_Uf(:,:,a1,a2)=d/dA(a1,l1) d/dA(a2,l2) Uf
call calc_diffdiff_Uf(diffdiff_Uf,Uf,f,l1,l2,UMAT)
!! diff_Uf_l1(:,:,a1)=d/dA(a1,l1) Uf
call calc_diff_Uf(diff_Uf_l1, Uf, f, l1, UMAT)
!! diff_Uf_l2(:,:,a2)=d/dA(a2,l2) Uf
call calc_diff_Uf(diff_Uf_l2, Uf, f, l2, UMAT)

diffdiff_Ufm=(0d0,0d0)

if (m_omega == 1) then 
  diffdiff_Ufm=diffdiff_Uf
  return
endif

if ( m_omega >= 2 ) then
  do k=0, m_omega-1
    if ( k .ne. 0 ) then 
      call matrix_power(tmpmat1,Uf,k)
    endif
    if( k .ne. m_omega-1 ) then 
      call matrix_power(tmpmat2,Uf,m_omega-k-1)
    endif
    do a1=1,dimG
      do a2=1,dimG
        if( k.ne.0 ) then 
          call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat1, NMAT, &
            diffdiff_Uf(:,:,a1,a2), NMAT, &
            (0d0,0d0), tmpmat3, NMAT)
        else
          tmpmat3=diffdiff_Uf(:,:,a1,a2)
        endif
        !!!
        if( k .ne. m_omega-1 ) then
          call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat3, NMAT, &
            tmpmat2, NMAT, &
            (1d0,0d0), diffdiff_Ufm(:,:,a1,a2), NMAT)
        else
          diffdiff_Ufm(:,:,a1,a2)=diffdiff_Ufm(:,:,a1,a2)+tmpmat3
        endif
      enddo
    enddo
  enddo

  do k=0,m_omega-2
    do j=0,m_omega-k-2
      do a1=1,dimG
        do a2=1,dimG
          ! Uf^j
          call matrix_power(tmpmat1,Uf,j)
          ! Uf^j.diff_Uf_(a1,l1)
          call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat1, NMAT, &
            diff_Uf_l1(:,:,a1), NMAT, &
            (0d0,0d0), tmpmat2, NMAT)
          ! Uf^k
          call matrix_power(tmpmat1,Uf,k)
          ! Uf^j.diff_Uf_(a1,l1).Uf^k
          call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat2, NMAT, &
            tmpmat1, NMAT, &
            (0d0,0d0), tmpmat3, NMAT)
          ! Uf^j.diff_Uf_l1(a1).Uf^k.diff_Uf_l2(a2)
          call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat3, NMAT, &
            diff_Uf_l2(:,:,a2), NMAT, &
            (0d0,0d0), tmpmat2, NMAT)
          ! Uf^{m-j-k-2}
            call matrix_power(tmpmat1,Uf,m_omega-j-k-2)
          ! Uf^j.diff_Uf_(a1).Uf^k.diff_Uf_l2(a2).Uf^{m-j-k-2}
          call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat2, NMAT, &
            tmpmat1, NMAT, &
            (1d0,0d0), diffdiff_Ufm(:,:,a1,a2), NMAT)
        enddo
      enddo
      !!!
      do a1=1,dimG
        do a2=1,dimG
          ! Uf^j
          call matrix_power(tmpmat1,Uf,j)
          ! Uf^j.diff_Uf_(a2,l2)
          call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat1, NMAT, &
            diff_Uf_l2(:,:,a2), NMAT, &
            (0d0,0d0), tmpmat2, NMAT)
          ! Uf^k
          call matrix_power(tmpmat1,Uf,k)
          ! Uf^j.diff_Uf_(a2,l2).Uf^k
          call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat2, NMAT, &
            tmpmat1, NMAT, &
            (0d0,0d0), tmpmat3, NMAT)
          ! Uf^j.diff_Uf_l2(a2).Uf^k.diff_Uf_l1(a1)
          call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat3, NMAT, &
            diff_Uf_l1(:,:,a1), NMAT, &
            (0d0,0d0), tmpmat2, NMAT)
          ! Uf^{m-j-k-2}
            call matrix_power(tmpmat1,Uf,m_omega-j-k-2)
          ! Uf^j.diff_Uf_(a1).Uf^k.diff_Uf_l2(a2).Uf^{m-j-k-2}
          call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat2, NMAT, &
            tmpmat1, NMAT, &
            (1d0,0d0), diffdiff_Ufm(:,:,a1,a2), NMAT)
        enddo
      enddo
    enddo
  enddo
!do a1=1,dimG
!  do a2=1,dimG
!    write(*,*)  diffdiff_Ufm(:,:,a1,a2)-diffdiff_Ufm(:,:,a2,a1)
!  enddo
!enddo
endif



end subroutine calc_diffdiff_Ufm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to compute d/dA(a1,l1) d/dA(a2,l2) Uf
!! <notation>
!!   diffdiff_Uf(:,;,a1,a2) = d/dA(a1,l1) d/dA(a2,l2) Uf
subroutine calc_diffdiff_Uf(diffdiff_Uf, Uf, f, l1, l2, UMAT)

complex(kind(0d0)), intent(out) :: diffdiff_Uf(1:NMAT,1:NMAT,1:dimG,1:dimG)
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
integer, intent(in) :: l1,l2,f

integer :: place_l1, place_l2
complex(kind(0d0)) :: XMAT(1:NMAT,1:NMAT)
complex(kind(0d0)) :: YMAT(1:NMAT,1:NMAT)
complex(kind(0d0)) :: ZMAT(1:NMAT,1:NMAT)
complex(kind(0d0)) :: diff_UlAl1(1:NMAT,1:NMAT,1:dimG)
complex(kind(0d0)) :: diff_UlAl2(1:NMAT,1:NMAT,1:dimG)
complex(kind(0d0)) :: diffdiff_UlAl(1:NMAT,1:NMAT,1:dimG,1:dimG)
complex(kind(0d0)) :: tmp_diffdiff_Uf(1:NMAT,1:NMAT,1:dimG,1:dimG)
complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT), tmpmat2(1:NMAT,1:NMAT)
character(1) :: char1
integer :: link_place, link_place_l1, link_place_l2
integer :: nl
integer :: a1,a2,info

!! when l1=l2
if ( l1 == l2 ) then 
  ! find the place of l1 and l2
  do link_place_l1=1,links_in_f(f)%num_
    if( links_in_f(f)%link_labels_(link_place_l1) == l1 ) then 
      exit
    endif
  enddo
  ! d/dA(a1,l1) d/dA(a2,l2) Ul
  call calc_diffdiff_Ul_in_face(diffdiff_UlAl,f,link_place_l1,UMAT)
  ! X=U(1)...U(l1-1)
  call calc_prodUl_from_n1_to_n2_in_Uf(XMAT,f,1,link_place_l1-1,UMAT)
  ! Z=U(l1+1)...U(f_size)
  call calc_prodUl_from_n1_to_n2_in_Uf(ZMAT,f,link_place_l1+1,links_in_f(f)%num_,UMAT)
  ! X.diffdiff.Z
  do a1=1,dimG
    do a2=1,dimG
      call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
         XMAT,NMAT, &
         diffdiff_UlAl(:,:,a1,a2),NMAT, &
         (0d0,0d0), tmpmat1, NMAT)
      call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
         tmpmat1,NMAT, &
         ZMAT,NMAT, &
         (0d0,0d0), diffdiff_Uf(:,:,a1,a2), NMAT)
    enddo
  enddo
!do a1=1,dimG
!  do a2=1,dimG
!    write(*,*) diffdiff_Uf(:,:,a1,a2)-diffdiff_Uf(:,:,a2,a1)
!  enddo
!enddo
  return
!!!!!!!!!!!!!!
! when l_1 != l_2
else
! find the places of l1 and l2
  info=0
  do nl=1,links_in_f(f)%num_
    if( links_in_f(f)%link_labels_(nl) == l1 ) then 
      link_place_l1=nl
      info=info+1
    elseif( links_in_f(f)%link_labels_(nl) == l2 ) then 
      link_place_l2=nl
      info=info+1
    endif
    if(info==2) exit
  enddo
  if(info*(info-1)==0) then 
    write(*,*) "l_1 and/or l_2 are not included in the face f"
    stop
  endif

  call calc_diff_Ul_in_face(diff_UlAl1,f,link_place_l1,UMAT)
  call calc_diff_Ul_in_face(diff_UlAl2,f,link_place_l2,UMAT)

  !!!!!!!!!!!!!!!!!!
  if(link_place_l1 < link_place_l2) then 
    call calc_prodUl_from_n1_to_n2_in_Uf(XMAT,f,1,link_place_l1-1,UMAT)
    call calc_prodUl_from_n1_to_n2_in_Uf(YMAT,f,link_place_l1+1,link_place_l2-1,UMAT)
    call calc_prodUl_from_n1_to_n2_in_Uf(ZMAT,f,link_place_l2+1,links_in_f(f)%num_,UMAT)
  do a1=1,dimG
    do a2=1,dimG
      call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
         XMAT,NMAT, &
         diff_UlAl1(:,:,a1),NMAT, &
         (0d0,0d0), tmpmat1, NMAT)
      call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
         tmpmat1,NMAT, &
         YMAT,NMAT, &
         (0d0,0d0), tmpmat2, NMAT)
      call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
         tmpmat2,NMAT, &
         diff_UlAl2(:,:,a2),NMAT, &
         (0d0,0d0), tmpmat1, NMAT)
      call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
         tmpmat1,NMAT, &
         ZMAT,NMAT, &
         (0d0,0d0),  diffdiff_Uf(:,:,a1,a2), NMAT)
    enddo
  enddo
  !!!!!!!!!!!!!!!!!!1
  else
    call calc_prodUl_from_n1_to_n2_in_Uf(XMAT,f,1,link_place_l2-1,UMAT)
    call calc_prodUl_from_n1_to_n2_in_Uf(YMAT,f,link_place_l2+1,link_place_l1-1,UMAT)
    call calc_prodUl_from_n1_to_n2_in_Uf(ZMAT,f,link_place_l1+1,links_in_f(f)%num_,UMAT)
  do a1=1,dimG
    do a2=1,dimG
      call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
         XMAT,NMAT, &
         diff_UlAl2(:,:,a2),NMAT, &
         (0d0,0d0), tmpmat1, NMAT)
      call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
         tmpmat1,NMAT, &
         YMAT,NMAT, &
         (0d0,0d0), tmpmat2, NMAT)
      call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
         tmpmat2,NMAT, &
         diff_UlAl1(:,:,a1),NMAT, &
         (0d0,0d0), tmpmat1, NMAT)
      call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
         tmpmat1,NMAT, &
         ZMAT,NMAT, &
         (0d0,0d0),  diffdiff_Uf(:,:,a1,a2), NMAT)
    enddo
  enddo
  endif 
endif
!diffdiff_Uf=tmp_diffdiff_Uf
!! symmetrize
!do a1=1,dimG
!  do a2=a1,dimG
!    diffdiff_Uf(:,:,a1,a2)=&
!      (0.5d0,0d0)*(tmp_diffdiff_Uf(:,:,a1,a2)+tmp_diffdiff_Uf(:,:,a2,a1))
!    if ( a1 .ne. a2 ) then 
!      diffdiff_Uf(:,:,a2,a1)=diffdiff_Uf(:,:,a1,a2)
!    endif
!  enddo
!enddo

!do a1=1,dimG
!  do a2=1,dimG
!    write(*,*) diffdiff_Uf(:,:,a1,a2)-diffdiff_Uf(:,:,a2,a1)
!  enddo
!enddo


return

end subroutine calc_diffdiff_Uf


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! caluculate product of Ul's from l1 to l2 in the face f.
subroutine calc_prodUl_from_n1_to_n2_in_Uf(ProdU,f,n1,n2,UMAT)
use matrix_functions, only : make_unit_matrix
implicit none

complex(kind(0d0)), intent(out) :: ProdU(1:NMAT,1:NMAT)
integer, intent(in) :: f,n1,n2
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)

complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
character(1) :: char1
integer :: link_place,i,j

!! if there is no links in this period, return the unit matrix
if ( n2 == n1-1 ) then 
  call make_unit_matrix(ProdU)
  return
endif

if(links_in_f(f)%link_dirs_(n1)==1) then 
  ProdU=UMAT(:,:,links_in_f(f)%link_labels_(n1))
else
  do i=1,NMAT
    do j=1,NMAT
      ProdU(i,j)=dconjg(UMAT(j,i,links_in_f(f)%link_labels_(n1)))
    enddo
  enddo
endif
if ( n2 >= n1+1 ) then 
  do link_place=n1+1,n2
    tmpmat=ProdU
    if(links_in_f(f)%link_dirs_(link_place)==1) then 
      char1='N'
    elseif(links_in_f(f)%link_dirs_(link_place)==-1) then
      char1='C'
    endif
    !! Uf=tmpmat*Uf or tmpmat*(Uf)^\dagger
    call ZGEMM('N',char1,NMAT,NMAT,NMAT,(1d0,0d0), &
      tmpmat,NMAT, &
      UMAT(:,:,links_in_f(f)%link_labels_(link_place)),NMAT, &
      (0d0,0d0), ProdU, NMAT)
  enddo
endif

end subroutine calc_prodUl_from_n1_to_n2_in_Uf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate d/dA(a,l) U_l or U_l^{-1} in Uf for all a
!!  where l is the n-th link in f
subroutine calc_diff_Ul_in_face(diff_UlAl,f,link_place,UMAT)
use SUN_generators, only : MtimesT, TtimesM
use SUN_generators, only : make_SUN_generators ! for test
implicit none

complex(kind(0d0)), intent(out) :: diff_UlAl(1:NMAT,1:NMAT,1:dimG)
integer, intent(in) :: f,link_place
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)

complex(kind(0d0)) :: Uinv_times_minus_i(1:NMAT,1:NMAT)
integer :: l,a,i,j

!! for test
!complex(kind(0d0)) ::  T(1:NMAT,1:NMAT,1:NMAT**2-1) ! generators
!complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT),tmp
!call make_SUN_generators(T,NMAT)

l=links_in_f(f)%link_labels_(link_place)

if ( links_in_f(f)%link_dirs_(link_place) == 1 ) then 
  do a=1,dimG
    call TtimesM(diff_UlAl(:,:,a), (0d0,1d0)*UMAT(:,:,l), a, NMAT) ! NO BUG
    !tmpmat=diff_UlAl(:,:,a)-(0d0,1d0)*MATMUL(T(:,:,a),UMAT(:,:,l))
    !tmp=(0d0,0d0)
    !do i=1,NMAT
    !  do j=1,NMAT
    !    tmp=tmp+tmpmat(i,j)*dconjg(tmpmat(i,j))
    !  enddo
    !enddo
    !write(*,*) tmp
  enddo
else
  do i=1,NMAT
    do j=1,NMAT
      Uinv_times_minus_i(i,j)=(0d0,-1d0)*dconjg(UMAT(j,i,l))
    enddo
  enddo
  do a=1,dimG
    call MtimesT(diff_UlAl(:,:,a), Uinv_times_minus_i, a, NMAT) ! NO BUG
    !call ZGEMM('C','N',NMAT,NMAT,NMAT,(0d0,-1d0), &
    !  UMAT(:,:,l),NMAT, &
    !  T(:,:,a),NMAT, &
    !  (0d0,0d0), tmpmat, NMAT)
    !do i=1,NMAT
    !  do j=1,NMAT
    !    tmpmat(i,j)=diff_UlAl(i,j,a)-tmpmat(i,j)
    !  enddo
    !enddo
    !tmp=(0d0,0d0)
    !do i=1,NMAT
    !  do j=1,NMAT
    !    tmp=tmp+tmpmat(i,j)*dconjg(tmpmat(i,j))
    !  enddo
    !enddo
    !write(*,*) tmp
  enddo
endif

end subroutine calc_diff_Ul_in_face


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate d/dA(a1,l) d/dA(a2,l) U_l or U_l^{-1} in Uf
!!  for all a1 and a2,  where l is the n-th link in f
subroutine calc_diffdiff_Ul_in_face(diffdiff_UlAl,f,link_place,UMAT)
use SUN_generators, only : TtimesM, MtimesT
implicit none

complex(kind(0d0)), intent(out) :: diffdiff_UlAl(1:NMAT,1:NMAT,1:dimG,1:dimG)
integer, intent(in) :: f,link_place
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)

complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Uinv_times_minus_1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmp_diffdiff_UlAl(1:NMAT,1:NMAT,1:dimG,1:dimG)
integer :: l,a1,a2,i,j

l=links_in_f(f)%link_labels_(link_place)

if ( links_in_f(f)%link_dirs_(link_place) == 1 ) then 
  do a2=1,dimG
      call TtimesM(tmpmat, UMAT(:,:,l), a2, NMAT)
    do a1=1,dimG
      call TtimesM(diffdiff_UlAl(:,:,a2,a1),-tmpmat,a1,NMAT)
    enddo
  enddo
else
  do i=1,NMAT
    do j=1,NMAT
      Uinv_times_minus_1(i,j)=-dconjg(UMAT(j,i,l))
    enddo
  enddo
  do a2=1,dimG
      call MtimesT(tmpmat, Uinv_times_minus_1(:,:), a2, NMAT)
    do a1=1,dimG
      call MtimesT(diffdiff_UlAl(:,:,a2,a1),tmpmat,a1,NMAT)
    enddo
  enddo
endif
!diffdiff_UlAl=tmp_diffdiff_UlAl
!! symmetrize
!do a1=1,dimG
!  do a2=a1,dimG
!    diffdiff_UlAl(:,:,a1,a2)&
!      &=(tmp_diffdiff_UlAl(:,:,a1,a2)+tmp_diffdiff_UlAl(:,:,a2,a1))*(0.5d0,0d0)
!    if ( a1 .ne. a2 ) then
!      diffdiff_UlAl(:,:,a2,a1)=diffdiff_UlAl(:,:,a1,a2)
!    endif
!  enddo
!enddo


end subroutine calc_diffdiff_Ul_in_face


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!
integer function site_index(a,s) 
implicit none

integer, intent(in) :: a,s

site_index=a+dimG*(s-1)

end function site_index

!!!!!!!!!!!!!
integer function link_index(a,l) 
implicit none

integer, intent(in) :: a,l

link_index=a+dimG*(num_sites + l - 1)

end function link_index

!!!!!!!!!!!!!
integer function face_index(a,f)
implicit none

integer, intent(in) :: a,f

face_index=a+dimG*(num_sites + num_links + f - 1)

end function face_index


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to compute d/dA(a,l) Uf
subroutine tmp_calc_diff_Uf(Ufla, Uf,f, l,UMAT)
use SUN_generators, only : TtimesM, MtimesT, MTN
implicit none

complex(kind(0d0)), intent(out) :: Ufla(1:NMAT,1:NMAT,1:dimG)
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: XMAT(1:NMAT,1:NMAT), YMAT(1:NMAT,1:NMAT)
integer, intent(in) :: l,f

integer :: nl
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
character(1) :: char1
integer :: i,j,k,ll,a

do nl=1,links_in_f(f)%num_
  if( links_in_f(f)%link_labels_(nl) == l ) then 
    exit
  endif
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute \frac{d}{d A_l^a} U_f
!! when Ul is the first link
if( nl == 1 ) then 
!!!
! X=1, Y=U_l^{-1} U_f ( dir=1 )
  if( links_in_f(f)%link_dirs_(nl) == 1 ) then
    do a=1,dimG
      call TtimesM(Ufla(:,:,a), Uf, a, NMAT)
    enddo
    Ufla = (0d0,1d0)*Ufla
! X=1, Y=U_l U_f      ( dir=-1 )
  elseif( links_in_f(f)%link_dirs_(nl) == -1 ) then
    do i=1,NMAT
      do j=1,NMAT
        ! X \equiv U_l^{-1}
        XMAT(i,j)=dconjg( UMAT(j,i,l) )
      enddo
    enddo
    ! Y=U_l.U_f
    call ZGEMM('N','N',NMAT,NMAT,NMAT,(0d0,-1d0), &
      UMAT(:,:,l),NMAT, &
      Uf, NMAT, &
      (0d0,0d0), YMAT, NMAT)
    ! Uf_l^a = -i Ul^{-1}.T_a.Y
    do a=1,dimG
      call MTN(Ufla(:,:,a),XMAT,YMAT,a,NMAT)
    enddo
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! when Ul is the last link
elseif( nl == links_in_f(f)%num_  ) then
  if( links_in_f(f)%link_dirs_(nl) == 1) then
!  X=U_f U_l^{-1}, Y=1 (dir=1)
    call ZGEMM('N','C',NMAT,NMAT,NMAT,(0d0,1d0), &
      Uf,NMAT, &
      UMAT(:,:,l), NMAT, &
      (0d0,0d0), XMAT, NMAT)
    do a=1,dimG
      call MTN(Ufla(:,:,a), XMAT, UMAT(:,:,l),a, NMAT)
    enddo
  elseif( links_in_f(f)%link_dirs_(nl) == -1 ) then
!  X=U_f U_l,      Y=1 (dir=-1)
    do a=1,dimG
      call MtimesT(Ufla(:,:,a),Uf,a,NMAT)
    enddo
    Ufla=(0d0,-1d0)*Ufla
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Other cases
else
! X=U_1 ... U_{nl-1}
  if(links_in_f(f)%link_dirs_(1)==1) then
    XMAT=UMAT(:,:,links_in_f(f)%link_labels_(1))
  elseif(links_in_f(f)%link_dirs_(1)==-1) then
    do i=1,NMAT
    do j=1,NMAT
      XMAT(i,j)=dconjg(UMAT(j,i,links_in_f(f)%link_labels_(1)))
    enddo
    enddo
  endif
  !
  if(nl >= 3) then
  do ll=2,nl-1
    tmpmat=XMAT
    if(links_in_f(f)%link_dirs_(ll)==1) then 
      char1='N'
    elseif(links_in_f(f)%link_dirs_(ll)==-1) then
      char1='C'
    endif
    !! Uf=tmpmat*Uf or tmpmat*(Uf)^\dagger
    call ZGEMM('N',char1,NMAT,NMAT,NMAT,(1d0,0d0), &
      tmpmat,NMAT, &
      UMAT(:,:,links_in_f(f)%link_labels_(ll)),NMAT, &
      (0d0,0d0), XMAT, NMAT)
  enddo
  endif
! Y=U_{nl+1}...U_{links_in_f(f)%num_ }
  if(links_in_f(f)%link_dirs_(nl+1)==1) then
    YMAT=UMAT(:,:,links_in_f(f)%link_labels_(nl+1))
  elseif(links_in_f(f)%link_dirs_(nl+1)==-1) then
    do i=1,NMAT
    do j=1,NMAT
      YMAT(i,j)=dconjg(UMAT(j,i,links_in_f(f)%link_labels_(nl+1)))
    enddo
    enddo
  endif
  !
  if(nl <= links_in_f(f)%num_ -2) then
  do ll=nl+2,links_in_f(f)%num_ 
    tmpmat=YMAT
    if(links_in_f(f)%link_dirs_(ll)==1) then 
      char1='N'
    elseif(links_in_f(f)%link_dirs_(ll)==-1) then
      char1='C'
    endif
    !! Uf=tmpmat*Uf or tmpmat*(Uf)^\dagger
    call ZGEMM('N',char1,NMAT,NMAT,NMAT,(1d0,0d0), &
      tmpmat,NMAT, &
      UMAT(:,:,links_in_f(f)%link_labels_(ll)),NMAT, &
      (0d0,0d0), YMAT, NMAT)
  enddo
  endif

  if( links_in_f(f)%link_dirs_(nl) == 1 ) then
    ! tmpmat = i U_l Y
    call ZGEMM('N','N',NMAT,NMAT,NMAT,(0d0,1d0), &
      UMAT(:,:,l),NMAT, &
      YMAT,NMAT, &
      (0d0,0d0), tmpmat, NMAT)
    do a=1,dimG
      call MTN(Ufla(:,:,a),XMAT,tmpmat,a,NMAT)
    enddo
  elseif( links_in_f(f)%link_dirs_(nl) == -1 ) then
    ! tmpmat = -i X U_l^{-1}
    call ZGEMM('N','C',NMAT,NMAT,NMAT,(0d0,-1d0), &
      XMAT,NMAT, &
      UMAT(:,:,l),NMAT, &
      (0d0,0d0), tmpmat, NMAT)
    do a=1,dimG
      call MTN(Ufla(:,:,a),tmpmat,YMAT,a,NMAT)
    enddo
  endif
endif

end subroutine tmp_calc_diff_Uf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to compute d/dA(b,l_2) d/dA(a,l_1) Uf
!! <notation>
!!   diffdiff_Uf(:,;,b,a) = d/dA(b,l_2) d/dA(a,l_1) Uf
subroutine tmp_calc_diffdiff_Uf(diffdiff_Uf, Uf, f, l_2, l_1, UMAT)
use matrix_functions, only : make_unit_matrix
use SUN_generators, only : TtimesM, MtimesT, MTN
implicit none

! 
complex(kind(0d0)), intent(out) :: diffdiff_Uf(1:NMAT,1:NMAT,1:dimG,1:dimG)
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
integer, intent(in) :: l_1,l_2,f
complex(kind(0d0)) :: XMAT(1:NMAT,1:NMAT), YMAT(1:NMAT,1:NMAT), ZMAT(1:NMAT,1:NMAT)

integer :: nl
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT), tmpmat2(1:NMAT,1:NMAT)
character(1) :: char1
integer :: i,j,k,ll,l,a,b,a1,a2
integer :: nl1, nl2, l1, l2, info

!! when l_1 = l_2 
if( l_1 == l_2 ) then 
  l=l_1

  do nl=1,links_in_f(f)%num_
    if( links_in_f(f)%link_labels_(nl) == l ) then 
      exit
    endif
  enddo
  !! when Ul is the first link
  if( nl == 1 ) then 
  !!!
  ! X=1, Y=U_l^{-1} U_f ( dir=1 )
    if( links_in_f(f)%link_dirs_(nl) == 1 ) then
      do a=1,dimG
        call TtimesM(tmpmat, Uf, a, NMAT)
        do b=1,dimG
          call TtimesM(diffdiff_Uf(:,:,b,a), tmpmat, b, NMAT)
        enddo
      enddo
      diffdiff_Uf = (-1d0,0d0)*diffdiff_Uf
  ! X=1, Y=U_l U_f      ( dir=-1 )
    elseif( links_in_f(f)%link_dirs_(nl) == -1 ) then
      do i=1,NMAT
        do j=1,NMAT
          ! X \equiv U_l^{-1}
          XMAT(i,j)=dconjg( UMAT(j,i,l) )
        enddo
      enddo
      !Y \equiv -U_l.U_f 
      call ZGEMM('N','N',NMAT,NMAT,NMAT,(-1d0,0d0), &
        UMAT(:,:,l),NMAT, &
        Uf, NMAT, &
        (0d0,0d0), YMAT, NMAT)
      ! 
      do a=1,dimG
        call TtimesM(tmpmat,YMAT,b,NMAT)
        do b=1,dimG
          call MTN(diffdiff_Uf(:,:,b,a),XMAT,tmpmat,a,NMAT)
        enddo
      enddo
    endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! when Ul is the last link
  elseif( nl == links_in_f(f)%num_  ) then
    if( links_in_f(f)%link_dirs_(nl) == 1) then
  !  X=-U_f U_l^{-1}, Y=1 (dir=1)
      call ZGEMM('N','C',NMAT,NMAT,NMAT,(-1d0,0d0), &
        Uf,NMAT, &
        UMAT(:,:,l), NMAT, &
        (0d0,0d0), XMAT, NMAT)
      do a=1,dimG
        ! tmpmat = Ta . UMAT
        call TtimesM(tmpmat, UMAT(:,:,l),a, NMAT)
        do b=1,dimG
          call MTN(diffdiff_Uf(:,:,b,a),XMAT,tmpmat,b,NMAT)
        enddo
      enddo
    elseif( links_in_f(f)%link_dirs_(nl) == -1 ) then
  !  X=U_f U_l,      Y=1 (dir=-1)
      do a=1,dimG
        call MtimesT(tmpmat,Uf,a,NMAT)
        do b=1,dimG
          call MtimesT(diffdiff_Uf(:,:,b,a),tmpmat,b,NMAT)
        enddo
      enddo
      diffdiff_Uf=(-1d0,0d0)*diffdiff_Uf
    endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Other cases
  else
  ! X=U_1 ... U_{nl-1}
    if(links_in_f(f)%link_dirs_(1)==1) then
      XMAT=UMAT(:,:,links_in_f(f)%link_labels_(1))
    elseif(links_in_f(f)%link_dirs_(1)==-1) then
      do i=1,NMAT
      do j=1,NMAT
        XMAT(i,j)=dconjg(UMAT(j,i,links_in_f(f)%link_labels_(1)))
      enddo
      enddo
    endif
    !
    if(nl >= 3) then
    do ll=2,nl-1
      tmpmat=XMAT
      if(links_in_f(f)%link_dirs_(ll)==1) then 
        char1='N'
      elseif(links_in_f(f)%link_dirs_(ll)==-1) then
        char1='C'
      endif
      !! Uf=tmpmat*Uf or tmpmat*(Uf)^\dagger
      call ZGEMM('N',char1,NMAT,NMAT,NMAT,(1d0,0d0), &
        tmpmat,NMAT, &
        UMAT(:,:,links_in_f(f)%link_labels_(ll)),NMAT, &
        (0d0,0d0), XMAT, NMAT)
    enddo
    endif
  ! Y=U_{nl+1}...U_{links_in_f(f)%num_ }
    if(links_in_f(f)%link_dirs_(nl+1)==1) then
      YMAT=UMAT(:,:,links_in_f(f)%link_labels_(nl+1))
    elseif(links_in_f(f)%link_dirs_(nl+1)==-1) then
      do i=1,NMAT
      do j=1,NMAT
        YMAT(i,j)=dconjg(UMAT(j,i,links_in_f(f)%link_labels_(nl+1)))
      enddo
      enddo
    endif
    !
    if(nl <= links_in_f(f)%num_ -2) then
    do ll=nl+2,links_in_f(f)%num_ 
      tmpmat=YMAT
      if(links_in_f(f)%link_dirs_(ll)==1) then 
        char1='N'
      elseif(links_in_f(f)%link_dirs_(ll)==-1) then
        char1='C'
      endif
      !! Uf=tmpmat*Uf or tmpmat*(Uf)^\dagger
      call ZGEMM('N',char1,NMAT,NMAT,NMAT,(1d0,0d0), &
        tmpmat,NMAT, &
        UMAT(:,:,links_in_f(f)%link_labels_(ll)),NMAT, &
        (0d0,0d0), YMAT, NMAT)
    enddo
    endif
  
    if( links_in_f(f)%link_dirs_(nl) == 1 ) then
      ! tmpmat = - U_l Y
      call ZGEMM('N','N',NMAT,NMAT,NMAT,(-1d0,0d0), &
        UMAT(:,:,l),NMAT, &
        YMAT,NMAT, &
        (0d0,0d0), tmpmat, NMAT)
      ! tmpmat2 = -Ta U_l Y
      do a=1,dimG
        call TtimesM(tmpmat2,tmpmat,a,NMAT)
        do b=1,dimG
          call MTN(diffdiff_Uf(:,:,b,a),XMAT,tmpmat2,b,NMAT)
        enddo
      enddo
    elseif( links_in_f(f)%link_dirs_(nl) == -1 ) then
      ! tmpmat = - X U_l^{-1}
      call ZGEMM('N','C',NMAT,NMAT,NMAT,(-1d0,0d0), &
        XMAT,NMAT, &
        UMAT(:,:,l),NMAT, &
        (0d0,0d0), tmpmat, NMAT)
      ! tmpmat2 = -X U_l^{-1} Ta
      do a=1,dimG
      call MtimesT(tmpmat2,tmpmat,a,NMAT)
        do b=1,dimG
          call MTN(diffdiff_Uf(:,:,b,a),tmpmat2,YMAT,b,NMAT)
        enddo
      enddo
    endif
  endif
  
!!!!!!!!!!!!!!
! when l_1 != l_2
else
! set l1 and l2 in the order of l_1 and l_2 in the face f
  info=0
  do nl=1,links_in_f(f)%num_
    if( links_in_f(f)%link_labels_(nl) == l_1 ) then 
      if(info==0) then
        nl1=nl
        l1=l_1
        info=info+1
      elseif(info==1) then
        nl2=nl
        l2=l_1
        info=info+1
      endif
    elseif( links_in_f(f)%link_labels_(nl) == l_2 ) then 
      if(info==0) then
        nl1=nl
        l1=l_2
        info=info+1
      elseif(info==1) then
        nl2=nl
        l2=l_2
        info=info+1
      endif
    endif
    if(info==2) exit
  enddo
  if(info*(info-1)==0) then 
    write(*,*) "l_1 and/or l_2 are not included in the face f"
    stop
  endif
! now 
!   l1: former link in f 
!   nl1: position of l1 in the face
!
!   l2: later link in f
!   nl2: position of l2 in the face

!! X=U_l1...U_{nl1-1}
if( nl1==1 ) then 
  call make_unit_matrix(XMAT)
else
  if(links_in_f(f)%link_dirs_(1)==1) then
    XMAT=UMAT(:,:,links_in_f(f)%link_labels_(1))
  elseif(links_in_f(f)%link_dirs_(1)==-1) then
    do i=1,NMAT
    do j=1,NMAT
      XMAT(i,j)=dconjg(UMAT(j,i,links_in_f(f)%link_labels_(1)))
    enddo
    enddo
  endif
endif
!
if(nl1 >= 3) then
  do ll=2,nl1-1
    tmpmat=XMAT
    if(links_in_f(f)%link_dirs_(ll)==1) then 
      char1='N'
    elseif(links_in_f(f)%link_dirs_(ll)==-1) then
      char1='C'
    endif
    !! Uf=tmpmat*Uf or tmpmat*(Uf)^\dagger
    call ZGEMM('N',char1,NMAT,NMAT,NMAT,(1d0,0d0), &
      tmpmat,NMAT, &
      UMAT(:,:,links_in_f(f)%link_labels_(ll)),NMAT, &
      (0d0,0d0), XMAT, NMAT)
  enddo
endif
!! dir(U_{l1}) = -1 --> X -> -i X.U_{l1}^{-1}
if( links_in_f(f)%link_dirs_(nl1)==-1) then
  tmpmat=XMAT
  call ZGEMM('N','C',NMAT,NMAT,NMAT,(0d0,-1d0), &
    tmpmat,NMAT, &
    UMAT(:,:,l1),NMAT, &
    (0d0,0d0), XMAT, NMAT)
endif

!! Y=U_{nl1+1}...U_{nl2-1}
if( nl2 == nl1+1 ) then 
  call make_unit_matrix(YMAT)
else
  if(links_in_f(f)%link_dirs_(nl1+1)==1) then
    YMAT=UMAT(:,:,links_in_f(f)%link_labels_(nl1+1))
  elseif(links_in_f(f)%link_dirs_(nl1+1)==-1) then
    do i=1,NMAT
    do j=1,NMAT
      YMAT(i,j)=dconjg(UMAT(j,i,links_in_f(f)%link_labels_(nl1+1)))
    enddo
    enddo
  endif
endif
!
if(nl2-1 >= nl1+2) then
do ll=nl1+2,nl2-1
  tmpmat=YMAT
  if(links_in_f(f)%link_dirs_(ll)==1) then 
    char1='N'
  elseif(links_in_f(f)%link_dirs_(ll)==-1) then
    char1='C'
  endif
  !! Uf=tmpmat*Uf or tmpmat*(Uf)^\dagger
  call ZGEMM('N',char1,NMAT,NMAT,NMAT,(1d0,0d0), &
    tmpmat,NMAT, &
    UMAT(:,:,links_in_f(f)%link_labels_(ll)),NMAT, &
    (0d0,0d0), YMAT, NMAT)
  enddo
endif
!! dir(U_{l1}) = 1 --> Y -> i U_{l1}.Y
if( links_in_f(f)%link_dirs_(nl1)==1) then
  tmpmat=YMAT
  call ZGEMM('N','N',NMAT,NMAT,NMAT,(0d0,1d0), &
    UMAT(:,:,l1),NMAT, &
    tmpmat,NMAT, &
    (0d0,0d0), YMAT, NMAT)
endif
!! dir(U_{l2}) = -1 --> Y -> -i Y.U_{l2}^{-1}
if( links_in_f(f)%link_dirs_(nl2)==-1) then
  tmpmat=YMAT
  call ZGEMM('N','C',NMAT,NMAT,NMAT,(0d0,-1d0), &
    tmpmat,NMAT, &
    UMAT(:,:,l2),NMAT, &
    (0d0,0d0), YMAT, NMAT)
endif



!! Z=U_{nl2+1}...U_{num_links(f))
if( nl2 == links_in_f(f)%num_) then 
  call make_unit_matrix(ZMAT)
else 
  if(links_in_f(f)%link_dirs_(nl2+1)==1) then
    ZMAT=UMAT(:,:,links_in_f(f)%link_labels_(nl2+1))
  elseif(links_in_f(f)%link_dirs_(nl2+1)==-1) then
    do i=1,NMAT
    do j=1,NMAT
      ZMAT(i,j)=dconjg(UMAT(j,i,links_in_f(f)%link_labels_(nl2+1)))
    enddo
    enddo
  endif
endif
!
if(links_in_f(f)%num_ >= nl2+2) then
  do ll=nl2+2, links_in_f(f)%num_
    tmpmat=ZMAT
    if(links_in_f(f)%link_dirs_(ll)==1) then 
      char1='N'
    elseif(links_in_f(f)%link_dirs_(ll)==-1) then
      char1='C'
    endif
    !! Uf=tmpmat*Uf or tmpmat*(Uf)^\dagger
    call ZGEMM('N',char1,NMAT,NMAT,NMAT,(1d0,0d0), &
      tmpmat,NMAT, &
      UMAT(:,:,links_in_f(f)%link_labels_(ll)),NMAT, &
      (0d0,0d0), ZMAT, NMAT)
  enddo
endif
!! dir(U_{l2}) = 1 --> Z -> i U_{l1}.Z
if( links_in_f(f)%link_dirs_(nl2)==1) then
  tmpmat=ZMAT
  call ZGEMM('N','N',NMAT,NMAT,NMAT,(0d0,1d0), &
    tmpmat,NMAT, &
    UMAT(:,:,l2),NMAT, &
    (0d0,0d0), ZMAT, NMAT)
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   diffdiff_Uf(:,;,b,a) = d/dA(b,l_2) d/dA(a,l_1) Uf
! if l_1 is the first link, Ta is between X and Y 
  if( l1 == l_1 ) then 
    do a=1,dimG
      call MTN(tmpmat,XMAT,YMAT,a,NMAT)
      do b=1,dimG
        call MTN(diffdiff_Uf(:,:,b,a),tmpmat,ZMAT,b,NMAT)
      enddo
    enddo
! if l_2 is the first link, Ta is between Y and Z
  else 
    do a=1,dimG
      call MTN(tmpmat,YMAT,ZMAT,a,NMAT)
      do b=1,dimG
        call MTN(diffdiff_Uf(:,:,b,a),XMAT,tmpmat,b,NMAT)
      enddo
    enddo
  endif
endif

end subroutine tmp_calc_diffdiff_Uf


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! function to obtain the argument of a complex number z
real(8) function arg(z)
implicit none

complex(kind(0d0)), intent(in) :: z
complex(kind(0d0)) :: phase

phase = z / cmplx( abs(z) )
arg=atan2( dble((0d0,-1d0)*phase), dble(phase) )

end function arg

end module global_subroutines
