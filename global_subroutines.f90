!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! module to construct sevral quantities used in simulation 
!!
!! DO NOT CHANGE THE VALUES OF VARIABLES IN GLOBASL_PARAMETERS
module global_subroutines
use global_parameters
#ifdef PARALLEL
use parallel
#endif
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
!! Compute forward covariant difference of Phi
!!  DMAT = U_l Phi_tip(l) U_l^\dagger - Phi_org(l)
subroutine Make_diff_PhiMat(Dphi, l,UMAT,PhiMat)
use SUN_generators, only : Make_traceless_matrix_from_modes
implicit none

complex(kind(0d0)), intent(out) :: Dphi(1:NMAT,1:NMAT)
integer, intent(in) :: l
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
integer :: i,j
!complex(kind(0d0)) :: Phi_tip(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: Phi_org(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)

!call Make_traceless_matrix_from_modes(Phi_tip,NMAT,Phi(:,link_tip(l)))
!call Make_traceless_matrix_from_modes(Phi_org,NMAT,Phi(:,link_org(l)))
!call Make_traceless_matrix_from_modes(Dphi,NMAT,Phi(:,link_org(l)))

Dphi=PhiMat(:,:,link_org(l))

! U_l.Phi_tip
call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
    UMAT(:,:,l), NMAT, &
    PhiMat(:,:,link_tip(l)), NMAT, &
    (0d0,0d0), tmpmat1, NMAT)
! U_l.Phi_tip.U_l^\dagger
call ZGEMM('N','C',NMAT,NMAT,NMAT,(1d0,0d0), &
    tmpmat1, NMAT, &
    UMAT(:,:,l), NMAT, &
    (-1d0,0d0), Dphi, NMAT)
! U_l.Phi_tip.U_l^\dagger - Phi_org
!Dphi = Dphi - Phi_org

end subroutine Make_diff_PhiMat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to make moment map \Omega(Uf)
subroutine Make_moment_map0(Omega,Uf)
use matrix_functions, only : matrix_inverse
implicit none

complex(kind(0d0)), intent(out) :: Omega(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT)
integer :: i,j

do i=1,NMAT
  do j=1,NMAT
    Omega(i,j)=-im_unit*( Uf(i,j) - dconjg( Uf(j,i) ) )
  enddo
enddo

end subroutine Make_moment_map0


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
complex(kind(0d0)) :: trace

!call matrix_power(NMAT,Uf,m_omega,Ufm)

do i=1,NMAT
  do j=1,NMAT
    SMAT(i,j)=-im_unit*( Ufm(i,j) - dconjg( Ufm(j,i) ) )
    Cinv(i,j)= Ufm(i,j) + dconjg( Ufm(j,i) ) 
  enddo
enddo
! CMAT --> CMAT^{-1}
call matrix_inverse(Cinv)

!call ZGEMM('N','N',NMAT,NMAT,NMAT,dcmplx(1d0/dble(m_omega)), &
call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0),&
  SMAT,NMAT,&
  Cinv,NMAT,&
  (0d0,0d0),Omega,NMAT)
call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0),&
  Cinv,NMAT,&
  SMAT,NMAT,&
  (1d0,0d0),Omega,NMAT)

!Omega=(0d0,0d0)
!do i=1,NMAT
!  do j=1,NMAT
!    Omega(i,j)=tmpmat(i,j)+dconjg(tmpmat(j,i))
!  enddo
!enddo

Omega = Omega / dcmplx( dble(m_omega) )

if (NMAT > 2) then
  trace=(0d0,0d0)
  do i=1,NMAT
    trace=trace+Omega(i,i)
  enddo
  do i=1,NMAT
    Omega(i,i)=Omega(i,i)-trace / cmplx(dble( NMAT ))
  enddo
endif

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
if ( n2 <= n1-1 ) then 
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
!! calculate Amat used in fermionic action and force
subroutine calc_Amat(Amat,f,l_label,k,Uf,UMAT)
use matrix_functions, only : matrix_power,matrix_product
implicit none

complex(kind(0d0)), intent(out) :: Amat(1:NMAT,1:NMAT)
integer, intent(in) :: f,l_label,k
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)

complex(kind(0d0)):: tmpmat1(1:NMAT,1:NMAT)
complex(kind(0d0)):: tmpmat2(1:NMAT,1:NMAT)
!character :: C1
integer :: l,i,label

l=links_in_f(f)%link_labels_(l_label)

! tmpmat=U_f^{k-1}
if (k.ne.1) then
  call matrix_power(tmpmat1,Uf,k-1)
endif

! tmpmat2=U_f^{k-1} U_{l1}^{e1}...U_{l-1}^{el-1}
if ( links_in_f(f)%link_dirs_(l_label) == 1 ) then
  label = l_label - 1
else
  label = l_label
endif 
call calc_prodUl_from_n1_to_n2_in_Uf(tmpmat2,f,1,label,UMAT)

if (k.ne.1) then
  call matrix_product(Amat,tmpmat1,tmpmat2)
else
  Amat=tmpmat2
endif


end subroutine calc_Amat


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate Bmat used in fermionic action and force
subroutine calc_Bmat(Bmat,f,l_label,k,Uf,UMAT)
use matrix_functions, only : matrix_power,matrix_product
implicit none

complex(kind(0d0)), intent(out) :: Bmat(1:NMAT,1:NMAT)
integer, intent(in) :: f,l_label,k
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)

complex(kind(0d0)):: tmpmat1(1:NMAT,1:NMAT)
complex(kind(0d0)):: tmpmat2(1:NMAT,1:NMAT)
!character :: C1
integer :: l,i,label

l=links_in_f(f)%link_labels_(l_label)

! tmpmat1=U_f^{m-k-1}
if( m_omega-k >= 0 ) then 
  call matrix_power(tmpmat1,Uf,m_omega-k)
endif

! tmpmat=U_f^{k-1} U_{l1}^{e1}...U_{l-1}^{el-1}
if ( links_in_f(f)%link_dirs_(l_label) == 1 ) then
  label = l_label
else
  label = l_label+1
endif 
call calc_prodUl_from_n1_to_n2_in_Uf(tmpmat2,f,label,links_in_f(f)%num_ ,UMAT)

if( m_omega-k >= 0 ) then 
  call matrix_product(Bmat, tmpmat2, tmpmat1)
else
  Bmat=tmpmat2
endif 

end subroutine calc_Bmat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calclate U^m - U^{-m}
!!
!! CAUTION
!! Now sinU = Uf^m - Uf^{-1}
!! THERE IS NO -i FACTOR
subroutine calc_sinU_and_cosUinv(sinU,cosUinv,Ufm)
use matrix_functions, only : matrix_inverse,matrix_power,matrix_dagger_power
implicit none

complex(kind(0d0)), intent(out) :: sinU(NMAT,NMAT)
complex(kind(0d0)), intent(out) :: cosUinv(NMAT,NMAT)
complex(kind(0d0)), intent(in) :: Ufm(NMAT,NMAT)
complex(kind(0d0)):: tmpmat(NMAT,NMAT)
integer :: i,j

do i=1,NMAT
  do j=1,i
    if (i==j) then
      sinU(i,i)=Ufm(i,i) - conjg(Ufm(i,i))
      cosUinv(i,i)=Ufm(i,i) + conjg(Ufm(i,i))
    else
      sinU(i,j)=Ufm(i,j) - conjg(Ufm(j,i))
      sinU(j,i)=Ufm(j,i) - conjg(Ufm(i,j))
      cosUinv(i,j)=Ufm(i,j) + conjg(Ufm(j,i))
      cosUinv(j,i)=Ufm(j,i) + conjg(Ufm(i,j))
    endif 
  enddo
enddo
call matrix_inverse(cosUinv)

end subroutine calc_sinU_and_cosUinv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_dSinUdA(dSinUdA, Uf,UMAT,f,ll_label)
use matrix_functions, only : matrix_product
!use Dirac_operator, only : calc_Amat, calc_Bmat
implicit none

! for given f and ll
! d cosUinv(i,j,f) / dA_{ii,jj,ll)
complex(kind(0d0)), intent(out) :: dSinUdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
!complex(kind(0d0)), intent(in) :: sinU(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: UMat(1:NMAT,1:NMAT,1:num_links)
integer, intent(in) :: f,ll_label

complex(kind(0d0)) :: Amat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Bmat(1:NMAT,1:NMAT)

integer :: i,j,ii,jj

dSinUdA=(0d0,0d0)
call calc_Amat(Amat,f,ll_label,1,Uf,Umat)
call calc_Bmat(Bmat,f,ll_label,1,Uf,Umat)

do jj=1,NMAT
  do ii=1,NMAT
    do j=1,NMAT
      do i=1,NMAT
        dSinUdA(i,j,ii,jj)=dSinUdA(i,j,ii,jj) &
          + Amat(i,jj)*Bmat(ii,j) &
          + conjg(Bmat(jj,i))*conjg(Amat(j,ii))
      enddo
    enddo
  enddo
enddo

dSinUdA=dSinUdA * cmplx(dble(links_in_f(f)%link_dirs_(ll_label)))*(0d0,1d0)

end subroutine calc_dSinUdA



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to return product of d/dA D . vec
subroutine calc_dCosUinvdA_dSinUdA(dCosUinvdA,dSinUdA,&
    cosUinv,Uf,UMAT,f,ll_label)
use matrix_functions, only : matrix_product
!use Dirac_operator, only : calc_Amat, calc_Bmat
implicit none

! for given f and ll
! d cosUinv(i,j,f) / dA_{ii,jj,ll)
complex(kind(0d0)), intent(out) :: dCosUinvdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
complex(kind(0d0)), intent(out) :: dSinUdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
!complex(kind(0d0)), intent(in) :: sinU(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: cosUinv(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: UMat(1:NMAT,1:NMAT,1:num_links)
integer, intent(in) :: f,ll_label

complex(kind(0d0)) :: Amat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Bmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: UinvA(1:NMAT,1:NMAT), BUinv(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: AdagUinv(1:NMAT,1:NMAT), UinvBdag(1:NMAT,1:NMAT)


integer :: k,i,j,ii,jj
complex(kind(0d0)) :: factor

dCosUinvdA=(0d0,0d0)
dSinUdA=(0d0,0d0)
do k=1,m_omega
  call calc_Amat(Amat,f,ll_label,k,Uf,Umat)
  call calc_Bmat(Bmat,f,ll_label,k,Uf,Umat)

  call matrix_product(UinvA,cosUinv,Amat)
  call matrix_product(BUinv,Bmat,cosUinv)
  !call matrix_product(AdagUinv,Amat,cosUinv,'C','N')
  !call matrix_product(UinvBdag,cosUinv,Bmat,'N','C')

  do jj=1,NMAT
    do ii=1,NMAT
      do j=1,NMAT
        do i=1,NMAT
          dCosUinvdA(i,j,ii,jj)=dCosUinvdA(i,j,ii,jj) &
            + UinvA(i,jj)*BUinv(ii,j) &
            !- UinvBdag(i,jj)*AdagUinv(ii,j)
            - conjg( BUinv(jj,i) ) * conjg( UinvA(j,ii) )
          dSinUdA(i,j,ii,jj)=dSinUdA(i,j,ii,jj) &
            + Amat(i,jj)*Bmat(ii,j) &
            + conjg(Bmat(jj,i))*conjg(Amat(j,ii))
        enddo
      enddo
    enddo
  enddo
enddo

factor=cmplx(dble(links_in_f(f)%link_dirs_(ll_label)))*(0d0,1d0)
dCosUinvdA=dCosUinvdA*(-factor)
dSinUdA=dSinUdA * factor

end subroutine calc_dCosUinvdA_dSinUdA


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to return product of d/dA D . vec
subroutine calc_dABmatdA(dAmatdA,dBmatdA,Uf,UMAT,ll_label,f,l_label,k)
use matrix_functions, only : matrix_product, matrix_3_product, matrix_power,make_unit_matrix
implicit none

! for given f,l,k and ll
! d Amat(i,j,f,l,k) / dA_{ii,jj,ll)
complex(kind(0d0)), intent(out) :: dAmatdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
complex(kind(0d0)), intent(out) :: dBmatdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: UMat(1:NMAT,1:NMAT,1:num_links)
integer, intent(in) :: f,ll_label,l_label,k

complex(kind(0d0)) :: UU_initial_to_ll(1:NMAT,1:NMAT) ! 1..ll
complex(kind(0d0)) :: UU_initial_to_l(1:NMAT,1:NMAT) ! 1..l
complex(kind(0d0)) :: UU_ll_to_n(1:NMAT,1:NMAT) ! ll..n
complex(kind(0d0)) :: UU_l_to_n(1:NMAT,1:NMAT) ! l..n
complex(kind(0d0)) :: UU_ll_to_l(1:NMAT,1:NMAT) ! ll..l
complex(kind(0d0)) :: UU_l_to_ll(1:NMAT,1:NMAT) ! l..ll

complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat3(1:NMAT,1:NMAT)
complex(kind(0d0)) :: epsilon_r
character :: C1
integer :: n,label1,label2
integer :: i,j,ii,jj,kk

dAmatdA=(0d0,0d0)
dBmatdA=(0d0,0d0)
epsilon_r=cmplx(dble( links_in_f(f)%link_dirs_(ll_label) ))*(0d0,1d0)

n=links_in_f(f)%num_
!!!!!!!!!!!!!!
if ( links_in_f(f)%link_dirs_(ll_label) == 1 ) then
  label2 = ll_label - 1
else
  label2 = ll_label
endif 
call calc_prodUl_from_n1_to_n2_in_Uf(UU_initial_to_ll,f,1,label2,UMAT)
!!!!!!!!!!!!!!
if ( links_in_f(f)%link_dirs_(l_label) == 1 ) then
  label2 = l_label-1
else
  label2 = l_label
endif 
call calc_prodUl_from_n1_to_n2_in_Uf(UU_initial_to_l,f,1,label2,UMAT)
!!!!!!!!!!!!!!
if ( links_in_f(f)%link_dirs_(ll_label) == 1 ) then
  label1 = ll_label 
else
  label1 = ll_label+1
endif 
call calc_prodUl_from_n1_to_n2_in_Uf(UU_ll_to_n,f,label1,n,UMAT)
!!!!!!!!!!!!!!
if ( links_in_f(f)%link_dirs_(l_label) == 1 ) then
  label1 = l_label
else
  label1 = l_label+1
endif 
call calc_prodUl_from_n1_to_n2_in_Uf(UU_l_to_n,f,label1,n,UMAT)
!!!!!!!!!!!!!!
if ( links_in_f(f)%link_dirs_(ll_label) == 1 ) then
  label1 = ll_label 
else
  label1 = ll_label+1
endif 
if ( links_in_f(f)%link_dirs_(l_label) == 1 ) then
  label2 = l_label-1
else
  label2 = l_label
endif 
call calc_prodUl_from_n1_to_n2_in_Uf(UU_ll_to_l,f,label1,label2,UMAT)
!!!!!!!!!!!!!!
if ( links_in_f(f)%link_dirs_(l_label) == 1 ) then
  label1 = l_label 
else
  label1 = l_label+1
endif 
if ( links_in_f(f)%link_dirs_(ll_label) == 1 ) then
  label2 = ll_label-1
else
  label2 = ll_label
endif 
call calc_prodUl_from_n1_to_n2_in_Uf(UU_l_to_ll,f,label1,label2,UMAT)

        
!!!!!!!!!!!!!! dAmatdA !!!!!!!!!!!!!!!!!
if ( m_omega > 1 ) then
  do kk=1,k-1
    call matrix_power(tmpmat2,Uf,kk-1)
    call matrix_product(tmpmat1,tmpmat2,UU_initial_to_ll)

    call matrix_power(tmpmat3,Uf,k-kk-1)
    call matrix_3_product(tmpmat2,UU_ll_to_n,tmpmat3,UU_initial_to_l)

    do jj=1,NMAT
      do ii=1,NMAT
        do j=1,NMAT
          do i=1,NMAT
            dAmatdA(i,j,ii,jj)=dAmatdA(i,j,ii,jj) &
              + tmpmat1(i,jj) * tmpmat2(ii,j) * epsilon_r
          enddo
        enddo
      enddo
    enddo

  enddo
endif

if ( ll_label < l_label ) then
  call make_unit_matrix(tmpmat2)
  if( k-1 >= 0) call matrix_power(tmpmat2,Uf,k-1)
  call matrix_product(tmpmat1,tmpmat2,UU_initial_to_ll)

  do jj=1,NMAT
    do ii=1,NMAT
      do j=1,NMAT
        do i=1,NMAT
          dAmatdA(i,j,ii,jj)=dAmatdA(i,j,ii,jj) &
            + tmpmat1(i,jj) * UU_ll_to_l(ii,j) * epsilon_r
        enddo
      enddo
    enddo
  enddo
endif

if (ll_label == l_label) then
  if ( links_in_f(f)%link_dirs_(ll_label) == -1 ) then
    call make_unit_matrix(tmpmat2)
    if( k-1 >= 0 ) call matrix_power(tmpmat2,Uf,k-1)
    call matrix_product(tmpmat1, tmpmat2,UU_initial_to_l)
    do jj=1,NMAT
      do ii=1,NMAT
        j=ii
        do i=1,NMAT
          dAmatdA(i,j,ii,jj)=dAmatdA(i,j,ii,jj) &
            - (0d0,1d0)* tmpmat1(i,jj)
        enddo
      enddo
    enddo
  endif
endif

!!!!!!!!!!! dBmatdA !!!!!!!!!!!!!!!
if (m_omega > 1) then
  do kk=1,m_omega-k

    call make_unit_matrix(tmpmat2)
    if( kk-1 > 0) call matrix_power(tmpmat2, Uf, kk-1)
    call matrix_3_product(tmpmat1,UU_l_to_n,tmpmat2,UU_initial_to_ll)

    call make_unit_matrix(tmpmat3)
    if( m_omega-k-kk > 0) call matrix_power(tmpmat3,Uf,m_omega-k-kk)
    call matrix_product(tmpmat2,UU_ll_to_n,tmpmat3)
!
    do jj=1,NMAT
      do ii=1,NMAT
        do j=1,NMAT
          do i=1,NMAT
            dBmatdA(i,j,ii,jj)=dBmatdA(i,j,ii,jj) &
              + tmpmat1(i,jj) * tmpmat2(ii,j) * epsilon_r
          enddo
        enddo
      enddo
    enddo

  enddo
endif

if ( ll_label > l_label) then
  call make_unit_matrix(tmpmat1)
  if(m_omega-k > 0) call matrix_power(tmpmat1,Uf,m_omega-k)
  call matrix_product(tmpmat2,UU_ll_to_n,tmpmat1)

  do jj=1,NMAT
    do ii=1,NMAT
      do j=1,NMAT
        do i=1,NMAT
          dBmatdA(i,j,ii,jj)=dBmatdA(i,j,ii,jj) &
            + UU_l_to_ll(i,jj) * tmpmat2(ii,j) * epsilon_r
        enddo
      enddo
    enddo
  enddo
endif

if ( ll_label == l_label) then
  if ( links_in_f(f)%link_dirs_(ll_label) == 1) then
    call make_unit_matrix(tmpmat2)
    if( m_omega-k >= 0) call matrix_power(tmpmat2,Uf,m_omega-k)
    call matrix_product(tmpmat1,UU_ll_to_n,tmpmat2)
    do jj=1,NMAT
      do ii=1,NMAT
        i=jj
        do j=1,NMAT
          dBmatdA(i,j,ii,jj)=dBmatdA(i,j,ii,jj) &
            + (0d0,1d0) * tmpmat1(ii,j) 
        enddo
      enddo
    enddo
  endif
endif

end subroutine calc_dABmatdA



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
subroutine vec_to_mat(eta,lambda,chi,vec)
use SUN_generators, only : make_traceless_matrix_from_modes
implicit none

!complex(kind(0d0)), intent(out) :: eta(1:NMAT,1:NMAT,1:num_sites)
!complex(kind(0d0)), intent(out) :: lambda(1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)), intent(out) :: chi(1:NMAT,1:NMAT,1:num_faces)
!complex(kind(0d0)), intent(in) :: vec(1:sizeD)
complex(kind(0d0)), intent(out) :: eta(:,:,:)
complex(kind(0d0)), intent(out) :: lambda(:,:,:)
complex(kind(0d0)), intent(out) :: chi(:,:,:)
complex(kind(0d0)), intent(in) :: vec(:)

complex(kind(0d0)), allocatable :: ele(:)
integer :: s,l,f,a
integer :: dimG,sizeD,NMAT,num_sites,num_links,num_faces

NMAT=size(eta,1)
num_sites=size(eta,3)
num_links=size(lambda,3)
num_faces=size(chi,3)
sizeD=size(vec,1)
dimG=NMAT*NMAT-1

allocate( ele(1:dimG) )

do s=1,num_sites
  do a=1,dimG
    ele(a)=vec(site_index(a,s))
  enddo
  call make_traceless_matrix_from_modes(eta(:,:,s),NMAT,ele)
enddo

do l=1,num_links
  do a=1,dimG
    ele(a)=vec(link_index(a,l))
  enddo
  call make_traceless_matrix_from_modes(lambda(:,:,l),NMAT,ele)
enddo

do f=1,num_faces
  do a=1,dimG
    ele(a)=vec(face_index(a,f))
  enddo
  call make_traceless_matrix_from_modes(chi(:,:,f),NMAT,ele)
enddo

end subroutine vec_to_mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
subroutine mat_to_vec(vec,eta,lambda,chi)
use SUN_generators, only : trace_MTa
implicit none

!complex(kind(0d0)), intent(in) :: eta(1:NMAT,1:NMAT,1:num_sites)
!complex(kind(0d0)), intent(in) :: lambda(1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)), intent(in) :: chi(1:NMAT,1:NMAT,1:num_faces)
!complex(kind(0d0)), intent(out) :: vec(1:sizeD)
complex(kind(0d0)), intent(in) :: eta(:,:,:)
complex(kind(0d0)), intent(in) :: lambda(:,:,:)
complex(kind(0d0)), intent(in) :: chi(:,:,:)
complex(kind(0d0)), intent(out) :: vec(:)

complex(kind(0d0)) :: trace
integer :: s,l,f,a
integer :: dimG,sizeD,NMAT,num_sites,num_links,num_faces

NMAT=size(eta,1)
num_sites=size(eta,3)
num_links=size(lambda,3)
num_faces=size(chi,3)
sizeD=size(vec,1)
dimG=NMAT*NMAT-1

do s=1,num_sites
  do a=1,dimG
    call trace_MTa(trace,eta(:,:,s),a,NMAT)
    vec(site_index(a,s))=vec(site_index(a,s))+trace
  enddo
enddo
do l=1,num_links
  do a=1,dimG
    call trace_MTa(trace,lambda(:,:,l),a,NMAT)
    vec(link_index(a,l))=vec(link_index(a,l))+trace
  enddo
enddo
do f=1,num_faces
  do a=1,dimG
    call trace_MTa(trace,chi(:,:,f),a,NMAT)
    vec(face_index(a,f))=vec(face_index(a,f))+trace
  enddo
enddo
end subroutine mat_to_vec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! check distance from 1 of Uf
subroutine  check_distance(info,ratio,UMAT)
use global_parameters
use matrix_functions, only : matrix_norm, make_unit_matrix
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
integer, intent(out) :: info
double precision, intent(out) :: ratio
double precision :: distance,tmp
integer l,i,j,f
!double precision norm
complex(kind(0d0)) UNITM(1:NMAT,1:NMAT)
!complex(kind(0d0)) tmp(1:NMAT,1:NMAT)
complex(kind(0d0)) Uf(1:NMAT,1:NMAT)
!double precision dist(1:NMAT-1)

!write(*,*) "===== check distance from 1 ==========="
!write(*,*) "theoretical dist. to the nearest center=",dsin(PI/dble(NMAT))*2d0

info=0
ratio=0d0
call make_unit_matrix(UNITM)
do f=1,num_faces
  call Make_face_variable(Uf,f,UMAT)
  call matrix_norm(distance,UNITM-Uf)
  tmp=distance/maximal_dist
  if( tmp > ratio ) then 
    ratio=tmp
  endif
  if ( distance > maximal_dist ) then
    info=info+1
  endif
enddo

end subroutine check_distance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! function to obtain the argument of a complex number z
real(8) function arg(z)
implicit none

complex(kind(0d0)), intent(in) :: z
complex(kind(0d0)) :: phase

phase = z / cmplx( abs(z) )
arg=atan2( dble((0d0,-1d0)*phase), dble(phase) )

end function arg


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 必要な変数を通信するsubroutine
!subroutine variables_syncronizer(PhiMat, UMat)
!#ifdef PARALLEL
!use parallel
!
!integer s,l,f
!
!
!#endif
!end subroutine variables_syncronizer

end module global_subroutines
