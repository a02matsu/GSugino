!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! module for fermion matrix
module differential_Dirac
use global_parameters
use global_subroutines
use SUN_generators
use matrix_functions
#ifdef PARALLEL
use parallel
#endif
implicit none

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate force from dDdPhi
subroutine calc_force_from_dDdPhi(force,eta,lambda,chi,Deta,Dlambda,Dchi,Umat,ss)
implicit none

complex(kind(0d0)), intent(out) :: force(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: eta(1:NMAT,1:NMAT,1:num_necessary_sites,1:N_Remez4)
complex(kind(0d0)), intent(in) :: lambda(1:NMAT,1:NMAT,1:num_necessary_links,1:N_Remez4)
complex(kind(0d0)), intent(in) :: chi(1:NMAT,1:NMAT,1:num_necessary_faces,1:N_Remez4)
complex(kind(0d0)), intent(in) :: Deta(1:NMAT,1:NMAT,1:num_necessary_sites,1:N_Remez4)
complex(kind(0d0)), intent(in) :: Dlambda(1:NMAT,1:NMAT,1:num_necessary_links,1:N_Remez4)
complex(kind(0d0)), intent(in) :: Dchi(1:NMAT,1:NMAT,1:num_necessary_faces,1:N_Remez4)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links) 
integer, intent(in) :: ss

integer :: l,f,i,j,k,r
complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
integer :: ii,jj,info

force=(0d0,0d0)
!! (1) Dirac from site
if( p1 == 0 ) then
  do r=1,N_Remez4
    call matrix_Commutator(tmpmat1,eta(:,:,ss,r),Deta(:,:,ss,r),'N','C')
    force=force + (0.5d0,0d0)*dcmplx(Remez_alpha4(r)*alpha_s(ss))*tmpmat1
  enddo
endif

!! (2) Dirac from link 1
! there is no contribution
    
!! (3) Dirac from link 2
if( p3 == 0 ) then 
  do r=1,N_Remez4
    do k=1,linkorg_to_s(ss)%num_
      l=linkorg_to_s(ss)%labels_(k)
      call matrix_commutator(tmpmat1,Dlambda(:,:,l,r),lambda(:,:,l,r),'N','C')
      call matrix_3_product(tmpmat2,Umat(:,:,l),tmpmat1,Umat(:,:,l),'C','N','N')
      tmpmat2 = tmpmat2 * U1Rfactor(l)*U1Rfactor(l)
      force = force - dcmplx(Remez_alpha4(r)*alpha_l(l))*tmpmat2
    enddo
    !!
    do k=1,linktip_from_s(ss)%num_
      l=linktip_from_s(ss)%labels_(k)
      call matrix_commutator(tmpmat1,Dlambda(:,:,l,r),lambda(:,:,l,r),'N','C')
      force = force - dcmplx(Remez_alpha4(r)*alpha_l(l))*tmpmat1
    enddo
  enddo
endif

!! (4) Dirac from face 1
if ( p4 == 0 ) then
  do f=1,num_necessary_faces
    if( ss==sites_in_f(f)%label_(1) ) then
      do r=1,N_Remez4
        call matrix_commutator(tmpmat1,Dchi(:,:,f,r),chi(:,:,f,r),'C','N')
        force = force - dcmplx( 2d0*Remez_alpha4(r)*alpha_f(f) ) * tmpmat1
      enddo
    endif
  enddo
endif

  force = force * dcmplx(overall_factor)
end subroutine calc_force_from_dDdPhi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to return product of d/dA D . vec
subroutine calc_force_from_dDdA(&
    force,&
    eta,lambda,chi,Deta,Dlambda,Dchi,Umat,PhiMat,&
    Uf,Ufm,ll)
implicit none

complex(kind(0d0)), intent(out) :: force(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT,1:num_necessary_faces) 
complex(kind(0d0)), intent(in) :: Ufm(1:NMAT,1:NMAT,1:num_necessary_faces) 
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links) 
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: eta(1:NMAT,1:NMAT,1:num_necessary_sites,1:N_Remez4)
complex(kind(0d0)), intent(in) :: lambda(1:NMAT,1:NMAT,1:num_necessary_links,1:N_Remez4)
complex(kind(0d0)), intent(in) :: chi(1:NMAT,1:NMAT,1:num_necessary_faces,1:N_Remez4)
complex(kind(0d0)), intent(in) :: Deta(1:NMAT,1:NMAT,1:num_necessary_sites,1:N_Remez4)
complex(kind(0d0)), intent(in) :: Dlambda(1:NMAT,1:NMAT,1:num_necessary_links,1:N_Remez4)
complex(kind(0d0)), intent(in) :: Dchi(1:NMAT,1:NMAT,1:num_necessary_faces,1:N_Remez4)
integer, intent(in) :: ll

complex(kind(0d0)) :: trace,tmp,tmp_force(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT),tmpmat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat3(1:NMAT,1:NMAT),tmpmat4(1:NMAT,1:NMAT)


integer :: s,l,f
integer :: i,j,k,nl,r,ii,jj

!! preparation
!dDdA_vec=(0d0,0d0)
tmp_force=(0d0,0d0)

! for test
!call make_SUN_generators(T,NMAT)

!! (1) Dirac from site
!   no contribution

!! (2) Dirac from link 1
if ( p2==0 ) then
  tmp_force=(0d0,0d0)
  s=link_tip(ll)
  do r=1,N_Remez4
    call matrix_3_product(tmpmat2,Umat(:,:,ll),Deta(:,:,s,r),Umat(:,:,ll),'N','C','C')
    call matrix_commutator(tmpmat1,tmpmat2,lambda(:,:,ll,r))
    tmp_force= tmp_force &
    + dcmplx( (-Remez_alpha4(r))*alpha_l(ll) ) * tmpmat1 &
      * dconjg(U1Rfactor(ll))
    !!
    call matrix_3_product(tmpmat2,Umat(:,:,ll),eta(:,:,s,r),Umat(:,:,ll),'N','N','C')
    call matrix_commutator(tmpmat1,Dlambda(:,:,ll,r),tmpmat2,'C','N')
    tmp_force=tmp_force &
      + dcmplx( (-Remez_alpha4(r))*alpha_l(ll) ) * tmpmat1 &
      * dconjg(U1Rfactor(ll))
  enddo
endif


if ( p3==0 ) then
!! (3) Dirac from link 2
  s=link_tip(ll)
  call matrix_3_product(tmpmat2,Umat(:,:,ll),PhiMat(:,:,s),Umat(:,:,ll),'N','C','C')
  tmpmat2 = tmpmat2 * dconjg(U1Rfactor(ll)*U1Rfactor(ll))
  do r=1,N_Remez4
    call matrix_commutator(tmpmat3,lambda(:,:,ll,r),Dlambda(:,:,ll,r),'N','C')
    call matrix_commutator(tmpmat4,tmpmat2,tmpmat3)
    tmp_force=tmp_force-dcmplx(Remez_alpha4(r)*alpha_l(ll))*(0d0,1d0)*tmpmat4
  enddo
endif


!! (4) Dirac from face 1
!   no contribution

if( p5 == 0 ) then 
!! (5) Dirac from face 2
  !if ( p5_test == 0 ) then 
    if( m_omega == 0 ) then 
      call calc_fermion_force_from_omega_m0&
        (tmpmat1,lambda,chi,Dlambda,Dchi,Umat,ll)
    elseif( m_omega == -1 ) then 
      call calc_fermion_force_from_omega_adm&
        (tmpmat1,lambda,chi,Dlambda,Dchi,Umat,Uf,ll)
    else
      call calc_fermion_force_from_omega&
        (tmpmat1,lambda,chi,Dlambda,Dchi,Umat,Uf,ll)
    endif
  tmp_force=tmp_force+tmpmat1
endif 

!! fermion mass term
!if( p_mass == 0 ) then
!  s=link_tip(ll)
!  do r=1,N_Remez4
!    call matrix_3_product(tmpmat2,Umat(:,:,ll),Deta(:,:,s,r),Umat(:,:,ll),'N','C','C')
!    call matrix_commutator(tmpmat1,tmpmat2,lambda(:,:,ll,r))
!    tmp_force= tmp_force + dcmplx( (-Remez_alpha4(r))*alpha_l(ll)*mass_f ) * tmpmat1
!    !!
!    call matrix_3_product(tmpmat2,Umat(:,:,ll),eta(:,:,s,r),Umat(:,:,ll),'N','N','C')
!    call matrix_commutator(tmpmat1,Dlambda(:,:,ll,r),tmpmat2,'C','N')
!    tmp_force=tmp_force + dcmplx( (-Remez_alpha4(r))*alpha_l(ll)*mass_f ) * tmpmat1
!  enddo
!
! call calc_fermion_force_from_massf&
!    (tmpmat1,lambda,chi,Dlambda,Dchi,Umat,ll)
! tmp_force=tmp_force+tmpmat1
!endif

force=(0d0,0d0)
do ii=1,NMAT
  do jj=1,NMAT
    force(ii,jj)=force(ii,jj)&
      + (tmp_force(ii,jj)+dconjg(tmp_force(jj,ii)))
  enddo
enddo
force=force*dcmplx(overall_factor)

end subroutine calc_force_from_dDdA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_fermion_force_from_omega&
    (pre_force,lambda,chi,Dlambda,Dchi,Umat,Uf,ll)
implicit none

complex(kind(0d0)), intent(inout) :: pre_force(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: lambda(1:NMAT,1:NMAT,1:num_necessary_links,1:N_Remez4)
complex(kind(0d0)), intent(in) :: chi(1:NMAT,1:NMAT,1:num_necessary_faces,1:N_Remez4)
complex(kind(0d0)), intent(in) :: Dlambda(1:NMAT,1:NMAT,1:num_necessary_links,1:N_Remez4)
complex(kind(0d0)), intent(in) :: Dchi(1:NMAT,1:NMAT,1:num_necessary_faces,1:N_Remez4)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT,1:num_necessary_faces)
integer, intent(in) :: ll

complex(kind(0d0)) :: Uf0tom(1:NMAT,1:NMAT,0:m_omega)
complex(kind(0d0)) :: Ufm(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Xmat(1:NMAT,1:NMAT),Ymat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Cosinv(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Sinmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Omega(1:NMAT,1:NMAT)
complex(kind(0d0)) :: UXmat(1:NMAT,1:NMAT),YUmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat3(1:NMAT,1:NMAT)
complex(kind(0d0)) :: ini_F1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: F1_F2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: F2_fin(1:NMAT,1:NMAT)
integer :: X_last
integer :: l_place, ll_place

complex(kind(0d0)) :: dir_factor
complex(kind(0d0)) :: ll_dir_factor
!integer :: ll_dir

complex(kind(0d0)) :: dU_Mae(1:NMAT,1:NMAT),dU_Ushiro(1:NMAT,1:NMAT)
complex(kind(0d0)) :: dUX_Mae(1:NMAT,1:NMAT), dUX_Ushiro(1:NMAT,1:NMAT)
complex(kind(0d0)) :: dYU_Mae(1:NMAT,1:NMAT), dYU_Ushiro(1:NMAT,1:NMAT)
complex(kind(0d0)) :: dXY_Mae(1:NMAT,1:NMAT), dXY_Ushiro(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: Dsin_mae(1:NMAT,1:NMAT), Dsin_ushiro(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: Dcosinv_mae(1:NMAT,1:NMAT), Dcosinv_ushiro(1:NMAT,1:NMAT)

integer :: l,f
integer :: a,b,i,j,k,r,kk
integer :: info

complex(kind(0d0)) :: U1Rfactor_fl


pre_force=(0d0,0d0)
do a=1,face_in_l(ll)%num_
  f=face_in_l(ll)%label_(a)

  ! place of ll in f
  do ll_place=1,links_in_f(f)%num_
    if( links_in_f(f)%link_labels_(ll_place)==ll ) then
      ll_dir_factor = (0d0,1d0)*dcmplx(links_in_f(f)%link_dirs_(ll_place))
      exit
    endif
  enddo

  ! power seriese of Uf
  do k=0,m_omega-1
    call matrix_power(Uf0tom(:,:,k),Uf(:,:,f),k)
  enddo
  ! Uf^m
    call matrix_power(Ufm,Uf(:,:,f),m_omega)
  ! Sinmat and Cosinv
  do i=1,NMAT
    do j=1,NMAT
      Cosinv(i,j) = Ufm(i,j) + dconjg(Ufm(j,i))
      Sinmat(i,j) = Ufm(i,j) - dconjg(Ufm(j,i))
    enddo
  enddo
  call Matrix_inverse(Cosinv)
  !! Omega = (Ufm-Ufm^-1)/(Ufm+Ufm^-1)
  call matrix_product(Omega,Cosinv,Sinmat)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! for development
  !call make_unit_matrix(Cosinv)
  !call make_unit_matrix(Sinmat)
  !call make_unit_matrix(Omega)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !! d/dA U_f
  call calc_XYmat(dU_Mae,dU_Ushiro,f,ll_place,UMAT)

  do l_place=1,links_in_f(f)%num_
    l=links_in_f(f)%link_labels_(l_place)

    !! U1Rfactor
    call calc_U1Rfactor_fl(U1Rfactor_fl,f,l)


    dir_factor=&
      ll_dir_factor &
      * dcmplx(links_in_f(f)%link_dirs_(l_place)) &
      * (0d0,-2d0)/dcmplx(m_omega) & 
      * dcmplx(alpha_f(f) * beta_f(f)) &
      * U1Rfactor_fl


    !! Xmat and Ymat for l_place
    call calc_XYmat(Xmat,Ymat,f,l_place,UMAT)
    if( links_in_f(f)%link_dirs_(l_place) == 1 ) then
      X_last = l_place-1
    else
      X_last = l_place
    endif

    !! dX/dA and dY/dA
    if( ll_place <= l_place ) then
      call calc_dXdA(dXY_Mae,dXY_Ushiro,ll_dir_factor,f,l_place,ll_place,UMAT)
    else
      call calc_dYdA(dXY_Mae,dXY_Ushiro,ll_dir_factor,f,l_place,ll_place,UMAT)
    endif


    do k=0,m_omega-1
      !! UXmat
      call matrix_product(UXmat,Uf0tom(:,:,k),Xmat)
      !! YUmat
      call matrix_product(YUmat,Ymat,Uf0tom(:,:,m_omega-k-1))

      !! dUf/dA in UX part
      if( k>=1 ) then
        do kk=0,k-1
          !!!!!!!!!!!!!!!!!!!!!!!!!
          ! formar part
          call matrix_product( dUX_Mae, Uf0tom(:,:,kk), dU_Mae )
          ! latter part
          call matrix_3_product(dUX_Ushiro, dU_Ushiro, Uf0tom(:,:,k-kk-1), Xmat)

          ! term 1 and 5
          ini_F1=dUX_Ushiro
          F1_F2=YUmat
          call matrix_product(F2_fin,Cosinv,dUX_Mae)
          call update_preforce(pre_force,ini_F1,F1_F2,F2_fin,lambda,chi,Dlambda,Dchi,l,f,dir_factor,1)

          ! term 2 and 6
          call Hermitian_Conjugate(ini_F1,dUX_Mae)
          call matrix_product(F1_F2,Cosinv,YUmat,'N','C')
          call Hermitian_Conjugate(F2_fin,dUX_Ushiro)
          call update_preforce(pre_force,ini_F1,F1_F2,F2_fin,lambda,chi,Dlambda,Dchi,l,f,-dir_factor,2)

          ! term 3 and 7
          ini_F1=-dUX_Ushiro
          call matrix_product(F1_F2,YUmat,Omega)
          call matrix_product(F2_fin,Cosinv,dUX_Mae)
          call update_preforce(pre_force,ini_F1,F1_F2,F2_fin,lambda,chi,Dlambda,Dchi,l,f,dir_factor,1)

          ! term 4 and 8
          call matrix_product(ini_F1,dUX_Mae,Omega,'C','N')
          call matrix_product(F1_F2,Cosinv,YUmat,'N','C')
          call Hermitian_Conjugate(F2_fin,dUX_Ushiro)
          call update_preforce(pre_force,ini_F1,F1_F2,F2_fin,lambda,chi,Dlambda,Dchi,l,f,-dir_factor,2)

        enddo
      endif

      !! dX/dA in UX part
      if( ll_place <= X_last ) then 
        call calc_dXdA(dXY_Mae,dXY_Ushiro,ll_dir_factor,f,l_place,ll_place,UMAT)
        !!!!!!!!!!!!!!!!!!!!!!!!!
        ! formar part
        call matrix_product( dUX_Mae, Uf0tom(:,:,k), dXY_Mae )
        ! latter part
        dUX_Ushiro=dXY_Ushiro

        ! term 1 and 5
        ini_F1=dUX_Ushiro
        F1_F2=YUmat
        call matrix_product(F2_fin,Cosinv,dUX_Mae)
        call update_preforce(pre_force,ini_F1,F1_F2,F2_fin,lambda,chi,Dlambda,Dchi,l,f,dir_factor,1)

        ! term 2 and 6
        call Hermitian_Conjugate(ini_F1,dUX_Mae)
        call matrix_product(F1_F2,Cosinv,YUmat,'N','C')
        call Hermitian_Conjugate(F2_fin,dUX_Ushiro)
        call update_preforce(pre_force,ini_F1,F1_F2,F2_fin,lambda,chi,Dlambda,Dchi,l,f,-dir_factor,2)

        ! term 3 and 7
        ini_F1=-dUX_Ushiro
        call matrix_product(F1_F2,YUmat,Omega)
        call matrix_product(F2_fin,Cosinv,dUX_Mae)
        call update_preforce(pre_force,ini_F1,F1_F2,F2_fin,lambda,chi,Dlambda,Dchi,l,f,dir_factor,1)

        ! term 4 and 8
        call matrix_product(ini_F1,dUX_Mae,Omega,'C','N')
        call matrix_product(F1_F2,Cosinv,YUmat,'N','C')
        call Hermitian_Conjugate(F2_fin,dUX_Ushiro)
        call update_preforce(pre_force,ini_F1,F1_F2,F2_fin,lambda,chi,Dlambda,Dchi,l,f,-dir_factor,2)
      endif

      !! YU part
      if( m_omega-k-1>=1 ) then
        do kk=0,m_omega-k-2
          !!!!!!!!!!!!!!!!!!!!!!!!!
          ! formar part
          call matrix_3_product( dYU_Mae, Ymat, Uf0tom(:,:,kk), dU_Mae)
          ! latter part
          call matrix_product(dYU_Ushiro, dU_Ushiro, Uf0tom(:,:,m_omega-k-kk-2))

          ! term 1 and 5 
          ini_F1=dYU_Ushiro
          call matrix_product(F1_F2,Cosinv,UXmat)
          F2_fin=dYU_Mae
          call update_preforce(pre_force,ini_F1,F1_F2,F2_fin,lambda,chi,Dlambda,Dchi,l,f,dir_factor,2)

          ! term 2 and 6
          call Hermitian_Conjugate(ini_F1,dYU_Mae)
          call Hermitian_Conjugate(F1_F2,UXmat)
          call matrix_product(F2_fin,Cosinv,dYU_Ushiro,'N','C')
          call update_preforce(pre_force,ini_F1,F1_F2,F2_fin,lambda,chi,Dlambda,Dchi,l,f,-dir_factor,1)

          ! term 3 and 7
          call matrix_product(ini_F1,dYU_Ushiro,Omega)
          call matrix_product(F1_F2,Cosinv,UXmat)
          F2_fin=-dYU_Mae
          call update_preforce(pre_force,ini_F1,F1_F2,F2_fin,lambda,chi,Dlambda,Dchi,l,f,dir_factor,2)

          ! term 4 and 8
          call Hermitian_Conjugate(ini_F1,dYU_Mae)
          call matrix_product(F1_F2,UXmat,Omega,'C','N')
          call matrix_product(F2_fin,Cosinv,dYU_Ushiro,'N','C')
          call update_preforce(pre_force,ini_F1,F1_F2,F2_fin,lambda,chi,Dlambda,Dchi,l,f,-dir_factor,1)

        enddo
      endif

      !! dY/dA in UY part
      if( ll_place > X_last ) then 
        call calc_dYdA(dXY_Mae,dXY_Ushiro,ll_dir_factor,f,l_place,ll_place,UMAT)
        ! formar part
        dYU_Mae=dXY_Mae
        ! latter part
        call matrix_product(dYU_Ushiro, dXY_Ushiro, Uf0tom(:,:,m_omega-k-1))

        ! term 1 and 5 
        ini_F1=dYU_Ushiro
        call matrix_product(F1_F2,Cosinv,UXmat)
        F2_fin=dYU_Mae
        call update_preforce(pre_force,ini_F1,F1_F2,F2_fin,lambda,chi,Dlambda,Dchi,l,f,dir_factor,2)

        ! term 2 and 6
        call Hermitian_Conjugate(ini_F1,dYU_Mae)
        call Hermitian_Conjugate(F1_F2,UXmat)
        call matrix_product(F2_fin,Cosinv,dYU_Ushiro,'N','C')
        call update_preforce(pre_force,ini_F1,F1_F2,F2_fin,lambda,chi,Dlambda,Dchi,l,f,-dir_factor,1)

        ! term 3 and 7
        call matrix_product(ini_F1,dYU_Ushiro,Omega)
        call matrix_product(F1_F2,Cosinv,UXmat)
        F2_fin=-dYU_Mae
        call update_preforce(pre_force,ini_F1,F1_F2,F2_fin,lambda,chi,Dlambda,Dchi,l,f,dir_factor,2)

        ! term 4 and 8
        call Hermitian_Conjugate(ini_F1,dYU_Mae)
        call matrix_product(F1_F2,UXmat,Omega,'C','N')
        call matrix_product(F2_fin,Cosinv,dYU_Ushiro,'N','C')
        call update_preforce(pre_force,ini_F1,F1_F2,F2_fin,lambda,chi,Dlambda,Dchi,l,f,-dir_factor,1)

      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! cosinv part
      ! former part

      do kk=0,m_omega-1
        call matrix_product(tmpmat1, Uf0tom(:,:,kk),dU_Mae)
        call matrix_product(tmpmat2, dU_Ushiro,Uf0tom(:,:,m_omega-1-kk))

        ! F2_fin is common
        call matrix_product(F2_fin,Cosinv,-tmpmat1) ! take care for the sign

        ! term 1 and 5
        call matrix_3_product(ini_F1,tmpmat2,Cosinv,UXmat)
        F1_F2=YUmat
        call update_preforce(pre_force,ini_F1,F1_F2,F2_fin,lambda,chi,Dlambda,Dchi,l,f,dir_factor,1)
  
        ! term 2 and 6
        call matrix_3_product(ini_F1,tmpmat2,Cosinv,YUmat,'N','N','C')
        call Hermitian_conjugate(F1_F2,UXmat)
        call update_preforce(pre_force,ini_F1,F1_F2,F2_fin,lambda,chi,Dlambda,Dchi,l,f,dir_factor,1)

        ! term 3 and 7
        call matrix_3_product(ini_F1,tmpmat2,Cosinv,UXmat)
        call matrix_product(F1_F2,YUmat,Omega)
        call update_preforce(pre_force,ini_F1,F1_F2,F2_fin,lambda,chi,Dlambda,Dchi,l,f,-dir_factor,1)
  
        ! term 4 and 8
        call matrix_3_product(ini_F1,tmpmat2,Cosinv,YUmat,'N','N','C')
        call matrix_product(F1_F2,UXmat,Omega,'C','N')
        call update_preforce(pre_force,ini_F1,F1_F2,F2_fin,lambda,chi,Dlambda,Dchi,l,f,dir_factor,1)

        ! latter part
        ! F2_fin is common
        call matrix_product(F2_fin,Cosinv,tmpmat2,'N','C') 

        ! term 1 and 5
        call matrix_3_product(ini_F1,tmpmat1,Cosinv,UXmat,'C','N','N')
        F1_F2=YUmat
        call update_preforce(pre_force,ini_F1,F1_F2,F2_fin,lambda,chi,Dlambda,Dchi,l,f,dir_factor,1)
  
        ! term 2 and 6
        call matrix_3_product(ini_F1,tmpmat1,Cosinv,YUmat,'C','N','C')
        call Hermitian_conjugate(F1_F2,UXmat)
        call update_preforce(pre_force,ini_F1,F1_F2,F2_fin,lambda,chi,Dlambda,Dchi,l,f,dir_factor,1)

        ! term 3 and 7
        call matrix_3_product(ini_F1,tmpmat1,Cosinv,UXmat,'C','N','N')
        call matrix_product(F1_F2,YUmat,Omega)
        call update_preforce(pre_force,ini_F1,F1_F2,F2_fin,lambda,chi,Dlambda,Dchi,l,f,-dir_factor,1)
  
        ! term 4 and 8
        call matrix_3_product(ini_F1,tmpmat1,Cosinv,YUmat,'C','N','C')
        call matrix_product(F1_F2,UXmat,Omega,'C','N')
        call update_preforce(pre_force,ini_F1,F1_F2,F2_fin,lambda,chi,Dlambda,Dchi,l,f,dir_factor,1)

      enddo


      !! Omega part
      do kk=0, m_omega-1
        !!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! former part
        call matrix_3_product(tmpmat1,Cosinv,Uf0tom(:,:,kk),dU_Mae)

        ! ini_F1 is common
        call matrix_product(tmpmat2,dU_Ushiro, Uf0tom(:,:,m_omega-kk-1))
        ini_F1=tmpmat2
        call matrix_product(ini_F1,&
          tmpmat2, Omega, 'N','N',(-1d0,0d0),'ADD')

        ! term 3 and 7
        call matrix_product(F1_F2,-Cosinv,UXmat)
        call matrix_product(F2_fin,YUmat,tmpmat1)
        call update_preforce(pre_force,ini_F1,F1_F2,F2_fin,lambda,chi,Dlambda,Dchi,l,f,dir_factor,2)
  
        ! term 4 and 8
        call matrix_product(F1_F2,Cosinv,YUmat,'N','C')
        call matrix_product(F2_fin,UXmat,tmpmat1,'C','N')
        call update_preforce(pre_force,ini_F1,F1_F2,F2_fin,lambda,chi,Dlambda,Dchi,l,f,dir_factor,2)


        !!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! latter part
        call matrix_3_product(tmpmat1,&
          Cosinv,Uf0tom(:,:,m_omega-kk-1),dU_Ushiro,&
          'N','C','C')

        ! ini_F1 is common
        call matrix_product(tmpmat2, dU_Mae, Uf0tom(:,:,kk),'C','C')
        ini_F1=tmpmat2
        call matrix_product(ini_F1, &
          tmpmat2, Omega, 'N','N',(1d0,0d0),'ADD')

        ! term 3 and 7
        call matrix_product(F1_F2,-Cosinv,UXmat)
        call matrix_product(F2_fin,YUmat,tmpmat1)
        call update_preforce(pre_force,ini_F1,F1_F2,F2_fin,lambda,chi,Dlambda,Dchi,l,f,dir_factor,2)
  
        ! term 4 and 8
        call matrix_product(F1_F2,Cosinv,YUmat,'N','C')
        call matrix_product(F2_fin,UXmat,tmpmat1,'C','N')
        call update_preforce(pre_force,ini_F1,F1_F2,F2_fin,lambda,chi,Dlambda,Dchi,l,f,dir_factor,2)

      enddo

    enddo
  enddo
enddo

end subroutine calc_fermion_force_from_omega

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine update_preforce(pre_force,ini_F1,F1_F2,F2_fin,&
    lambda,chi,Dlambda,Dchi,l,f,dir_factor,order)
implicit none

complex(kind(0d0)), intent(inout) :: pre_force(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: ini_F1(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: F1_F2(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: F2_fin(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: lambda(1:NMAT,1:NMAT,1:num_necessary_links,1:N_Remez4)
complex(kind(0d0)), intent(in) :: chi(1:NMAT,1:NMAT,1:num_necessary_faces,1:N_Remez4)
complex(kind(0d0)), intent(in) :: Dlambda(1:NMAT,1:NMAT,1:num_necessary_links,1:N_Remez4)
complex(kind(0d0)), intent(in) :: Dchi(1:NMAT,1:NMAT,1:num_necessary_faces,1:N_Remez4)
integer, intent(in) :: l,f
complex(kind(0d0)), intent(in) :: dir_factor
integer, intent(in) :: order ! order=1: Dchi d(...)/dA lambda (...)
                             !      =2: Dchi (...) lambda d(...)/dA

complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)

integer :: r

if( order==1 ) then 
  do r=1,N_Remez4
    ! term in Dchi-lambda
    call matrix_3_product(tmpmat1,ini_F1,lambda(:,:,l,r),F1_F2)
    call matrix_3_product(tmpmat2,tmpmat1,Dchi(:,:,f,r),F2_fin,'N','C','N')
    pre_force=pre_force &
      -Remez_alpha4(r)*dir_factor*tmpmat2
    ! term in Dlambda-chi
    call matrix_3_product(tmpmat1,ini_F1,Dlambda(:,:,l,r),F1_F2,'N','C','N')
    call matrix_3_product(tmpmat2,tmpmat1,chi(:,:,f,r),F2_fin)
    pre_force=pre_force &
      +Remez_alpha4(r)*dir_factor*tmpmat2 
  enddo
elseif( order==2 ) then
  do r=1,N_Remez4
    ! term in Dchi-lambda
    call matrix_3_product(tmpmat1,ini_F1,Dchi(:,:,f,r),F1_F2,'N','C','N')
    call matrix_3_product(tmpmat2,tmpmat1,lambda(:,:,l,r),F2_fin)
    pre_force=pre_force &
      -Remez_alpha4(r)*dir_factor*tmpmat2
    ! term in Dlambda-chi
    call matrix_3_product(tmpmat1,ini_F1,chi(:,:,f,r),F1_F2)
    call matrix_3_product(tmpmat2,tmpmat1,Dlambda(:,:,l,r),F2_fin,'N','C','N')
    pre_force=pre_force &
      +Remez_alpha4(r)*dir_factor*tmpmat2 
  enddo
endif

end subroutine update_preforce

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_fermion_force_from_omega_adm&
    (pre_force,lambda,chi,Dlambda,Dchi,Umat,Uf,ll)
implicit none

complex(kind(0d0)), intent(inout) :: pre_force(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: lambda(1:NMAT,1:NMAT,1:num_necessary_links,1:N_Remez4)
complex(kind(0d0)), intent(in) :: chi(1:NMAT,1:NMAT,1:num_necessary_faces,1:N_Remez4)
complex(kind(0d0)), intent(in) :: Dlambda(1:NMAT,1:NMAT,1:num_necessary_links,1:N_Remez4)
complex(kind(0d0)), intent(in) :: Dchi(1:NMAT,1:NMAT,1:num_necessary_faces,1:N_Remez4)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT,1:num_necessary_faces)
intent(in) :: ll

complex(kind(0d0)) :: Xmat(1:NMAT,1:NMAT),Ymat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT),tmpmat2(1:NMAT,1:NMAT),tmpmat3(1:NMAT,1:NMAT)
complex(kind(0d0)) :: dBdA(1:NMAT,1:NMAT)
complex(kind(0d0)) :: sinU(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tr_DchiSin
complex(kind(0d0)) :: tr_chiSin
complex(kind(0d0)) :: tr_XlamY
complex(kind(0d0)) :: tr_XDlamY
!complex(kind(0d0)) :: dUfdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
complex(kind(0d0)) :: Mae(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Ushiro(1:NMAT,1:NMAT)
complex(kind(0d0)) :: factor
complex(kind(0d0)) :: XY_Mae(1:NMAT,1:NMAT)
complex(kind(0d0)) :: XY_Ushiro(1:NMAT,1:NMAT)
complex(kind(0d0)) :: XY_factor
!complex(kind(0d0)) :: dsinUdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
complex(kind(0d0)) :: dir_factor
complex(kind(0d0)) :: Bval,trace
complex(kind(0d0)) :: tmp,tmp1,tmp2
integer :: i,j,k,l,ll,f,a,b,l_place,ll_place,X_last,r!,kk
!complex(kind(0d0)) :: dUfdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)

pre_force=(0d0,0d0)
do k=1,face_in_l(ll)%num_
  f=face_in_l(ll)%label_(k)

  !!!!!!!!!!!!!!!!!!
  !! evaluate sin(U_f)
  do j=1,NMAT
    do i=1,NMAT
      sinU(i,j)=Uf(i,j,f)-dconjg( Uf(j,i,f) )
    enddo
  enddo
  !!!!!!!!!!!!!!!!!!
  !! evaluate B_f
  Bval=(1d0,0d0)
  do i=1,NMAT
    Bval=Bval - ((2d0,0d0)-Uf(i,i,f)-dconjg(Uf(i,i,f)))/(e_max*e_max) 
  enddo
  !!!!!!!!!!!!!!!!!!
  !! place of ll in f
  do ll_place=1,links_in_f(f)%num_
    if( links_in_f(f)%link_labels_(ll_place) == ll ) exit
  enddo
  
  !!!!!!!!!!!!!!!!!!
  !! Prepare dsin(U_f)/dA, dB/dA, dUf/DA
  call calc_dUfdA_dBdA(Mae,Ushiro,factor,dBdA,Umat,f,ll_place)

  do l_place=1,links_in_f(f)%num_
    l=links_in_f(f)%link_labels_(l_place)
    dir_factor=(0d0,-1d0)*dcmplx(&
      dble(links_in_f(f)%link_dirs_(l_place)) * alpha_f(f) * beta_f(f) )

    call calc_XYmat(Xmat,Ymat,f,l_place,UMAT)

    if( links_in_f(f)%link_dirs_(l_place) == 1 ) then
      X_last = l_place-1
    else
      X_last = l_place
    endif

    do r=1,N_Remez4
      !!!!!!!!!!!!!!!!!!
      !! tr( Dchi^\dag .(Uf-Ufdag) ), tr( chi.(Uf-Ufdag) )
      tr_DchiSin=(0d0,0d0)
      tr_chiSin=(0d0,0d0)
      do j=1,NMAT
        do i=1,NMAT
          tr_DchiSin=tr_DchiSin + dconjg(Dchi(i,j,f,r))*SinU(i,j)
          tr_chiSin=tr_chiSin +chi(j,i,f,r)*SinU(i,j)
        enddo
      enddo


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! dBdA contribution
      call matrix_3_product(tmpmat1,Xmat,lambda(:,:,l,r),Ymat)
      call matrix_3_product(tmpmat2,Ymat,lambda(:,:,l,r),Xmat,'C','N','C')
      !!
      trace=(0d0,0d0)
      tr_XlamY=(0d0,0d0)
      do j=1,NMAT
        tr_XlamY=tr_XlamY+tmpmat1(j,j)-tmpmat2(j,j)
        do i=1,NMAT
          trace=trace - dconjg(Dchi(i,j,f,r))&
            *(tmpmat1(i,j)+tmpmat2(i,j))/(Bval*Bval)
        enddo
      enddo
      !!
      call matrix_3_product(tmpmat1,Xmat,Dlambda(:,:,l,r),Ymat,'N','C','N')
      call matrix_3_product(tmpmat2,Ymat,Dlambda(:,:,l,r),Xmat,'C','C','C')
      tr_XDlamY=(0d0,0d0)
      do j=1,NMAT
        tr_XDlamY=tr_XDlamY+(tmpmat1(j,j) - tmpmat2(j,j))
        do i=1,NMAT
          trace=trace + chi(j,i,f,r) &
            *(tmpmat1(i,j)+tmpmat2(i,j))/(Bval*Bval)
        enddo
      enddo
      !!
      trace=trace+(2d0,0d0)*(tr_DchiSin*tr_XlamY-tr_chiSin*tr_XDlamY)&
        /(Bval**3 * e_max**2)
      !!
      !!!!!!!!!!!!!!!!!!
      pre_force=pre_force + trace*(-Remez_alpha4(r))*dir_factor*dBdA
      !!!!!!!!!!!!!!!!!!

      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! dUf/dA contribution
      tmp=(-Remez_alpha4(r))*dir_factor*(factor)/(Bval**2 * e_max**2)
      do i=1,NMAT
        do j=1,NMAT
          tmpmat1(i,j) = tmp * (-tr_XlamY * dconjg(Dchi(j,i,f,r)) &
            + tr_XDlamY * chi(i,j,f,r))
        enddo
      enddo
      !!!!!!!!!!!!!!!!!
      call matrix_3_product(pre_force,Ushiro,tmpmat1,Mae,'N','N','N',(1d0,0d0),'ADD')
      call matrix_3_product(pre_force,Mae,tmpmat1,Ushiro,'C','N','C',(1d0,0d0),'ADD')
      !!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! dX/dA contribution
      if( ll_place <= X_last ) then
        call calc_dXdA(XY_Mae,XY_Ushiro,XY_factor,f,l_place,ll_place,UMAT)

        !!!!!!!!!!!!!!!
        !! tmpmat1=lambda.Y
        !! tmpmat2=Dlambda^\dag . Y
        call matrix_product(tmpmat1,lambda(:,:,l,r),Ymat)
        call matrix_product(tmpmat2,Dlambda(:,:,l,r),Ymat,'C','N')

        call matrix_product(tmpmat3,tmpmat1,Dchi(:,:,f,r),'N','C',&
          XY_factor/Bval)               
        call matrix_product(tmpmat3,tmpmat2,chi(:,:,f,r),'N','N',&
          -XY_factor/Bval,'ADD')        
        tmpmat3=tmpmat3 + XY_factor/(Bval**2 * e_max**2) &
          *( -tr_DchiSin * tmpmat1 + tr_chiSin * tmpmat2 )
        !!!!!!!!!!!!!!!                 
        call matrix_3_product(pre_force,XY_Ushiro,tmpmat3,XY_Mae,'N','N','N',&
          -Remez_alpha4(r)*dir_factor,'ADD')
        !!!!!!!!!!!!!!!                 
                                        
                                        
        !!!!!!!!!!!!!!!                 
        !! tmpmat1=Ydag . lambda        
        !! tmpmat2=Ydag. Dlambda^dag    
        call matrix_product(tmpmat1,Ymat,lambda(:,:,l,r),'C','N')
        call matrix_product(tmpmat2,Ymat,Dlambda(:,:,l,r),'C','C')
                                        
        call matrix_product(tmpmat3,Dchi(:,:,f,r),tmpmat1,'C','N',&
          -XY_factor/Bval)              
        call matrix_product(tmpmat3,chi(:,:,f,r),tmpmat2,'N','N',&
          +XY_factor/Bval,'ADD')        
        tmpmat3=tmpmat3 + XY_factor/(Bval**2 * e_max**2) &
          *( -tr_DchiSin * tmpmat1 + tr_chiSin * tmpmat2 )
        !!!!!!!!!!!!!!!                 
        call matrix_3_product(pre_force,XY_Mae,tmpmat3,XY_Ushiro,'C','N','C',&
          -Remez_alpha4(r)*dir_factor,'ADD')
        !!!!!!!!!!!!!!!                 
                                        
                                        
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
      !! dY/dA contribution             
      else                              
        call calc_dYdA(XY_Mae,XY_Ushiro,XY_factor,f,l_place,ll_place,UMAT)
        !!!!!!!!!!!!!!!                 
        !! tmpmat1=X . lambda           
        !! tmpmat2=X. Dlambda^dag       
        call matrix_product(tmpmat1,Xmat,lambda(:,:,l,r))
        call matrix_product(tmpmat2,Xmat,Dlambda(:,:,l,r),'N','C')
                                        
        call matrix_product(tmpmat3,Dchi(:,:,f,r),tmpmat1,'C','N',&
          +XY_factor/Bval)              
        call matrix_product(tmpmat3,chi(:,:,f,r),tmpmat2,'N','N',&
          -XY_factor/Bval,'ADD')        
        tmpmat3 = tmpmat3 + XY_factor/(Bval**2 * e_max**2) &
          *( -tr_DchiSin * tmpmat1 + tr_chiSin * tmpmat2) 
        !!!!!!!!!!!!!!!                 
        call matrix_3_product(pre_force,XY_Ushiro,tmpmat3,XY_Mae,'N','N','N',&
          -Remez_alpha4(r)*dir_factor,'ADD')
        !!!!!!!!!!!!!!!                 
                                        
        !!!!!!!!!!!!!!!                 
        !! tmpmat1= lambda . Xdag       
        !! tmpmat2= Dlambda^dag . Xdag  
        call matrix_product(tmpmat1,lambda(:,:,l,r),Xmat,'N','C')
        call matrix_product(tmpmat2,Dlambda(:,:,l,r),Xmat,'C','C')
                                        
        call matrix_product(tmpmat3,tmpmat1,Dchi(:,:,f,r),'N','C',&
          -XY_factor/Bval)              
        call matrix_product(tmpmat3,tmpmat2,chi(:,:,f,r),'N','N',&
          +XY_factor/Bval,'ADD')        
        tmpmat3=tmpmat3 + XY_factor/(Bval**2 * e_max**2) &
          *( -tr_DchiSin * tmpmat1 + tr_chiSin * tmpmat2 )
        !!!!!!!!!!!!!!!                 
        call matrix_3_product(pre_force,XY_Mae,tmpmat3,XY_Ushiro,'C','N','C',&
          -Remez_alpha4(r)*dir_factor,'ADD')
        !!!!!!!!!!!!!!!                  
      endif                             
    enddo                                
  enddo                                   
enddo                                     
end subroutine calc_fermion_force_from_omega_adm
                                         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_fermion_force_from_omega_m0&
    (pre_force,lambda,chi,Dlambda,Dchi,Umat,ll)
use matrix_functions, only : make_unit_matrix
implicit none

complex(kind(0d0)), intent(out) :: pre_force(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: lambda(1:NMAT,1:NMAT,1:num_necessary_links,1:N_Remez4)
complex(kind(0d0)), intent(in) :: chi(1:NMAT,1:NMAT,1:num_necessary_faces,1:N_Remez4)
complex(kind(0d0)), intent(in) :: Dlambda(1:NMAT,1:NMAT,1:num_necessary_links,1:N_Remez4)
complex(kind(0d0)), intent(in) :: Dchi(1:NMAT,1:NMAT,1:num_necessary_faces,1:N_Remez4)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
integer, intent(in) :: ll

complex(kind(0d0)) :: Xmat(1:NMAT,1:NMAT),Ymat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: XY_Mae(1:NMAT,1:NMAT)
complex(kind(0d0)) :: XY_Ushiro(1:NMAT,1:NMAT)
complex(kind(0d0)) :: XY_factor
complex(kind(0d0)) :: dir_factor
integer :: k,l,f,l_place,ll_place,X_last,r!,kk

pre_force=(0d0,0d0)
do k=1,face_in_l(ll)%num_
  f=face_in_l(ll)%label_(k)

  !!!!!!!!!!!!!!!!!!
  !! place of ll in f
  do ll_place=1,links_in_f(f)%num_
    if( links_in_f(f)%link_labels_(ll_place) == ll ) exit
  enddo
  
  do l_place=1,links_in_f(f)%num_
    l=links_in_f(f)%link_labels_(l_place)
    dir_factor=(0d0,1d0)*dcmplx( &
      dble( links_in_f(f)%link_dirs_(l_place) ) * alpha_f(f) * beta_f(f) )

    !!!!!!!!!!!!!!!!!!
    !! Xmat and Ymat
    call calc_XYmat(Xmat,Ymat,f,l_place,UMAT)
    !call make_unit_matrix(Ymat)

    !!!!!!!!!!!!!!!!!!
    !! last entry of Xmat
    if( links_in_f(f)%link_dirs_(l_place) == 1 ) then
      X_last = l_place-1
    else
      X_last = l_place
    endif


    do r=1,N_Remez4
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! dX/dA contribution
      if( ll_place <= X_last ) then
        call calc_dXdA(XY_Mae,XY_Ushiro,XY_factor,f,l_place,ll_place,UMAT)

        !!!!!!!!!!!!!!!
        !! tmpmat1=lambda.Y.D\chi^\dagger + D\lambda^\dagger.Y.\chi
        call matrix_3_product(tmpmat1,lambda(:,:,l,r),Ymat,Dchi(:,:,f,r),&
          'N','N','C')!,XY_factor)
        call matrix_3_product(tmpmat1,Dlambda(:,:,l,r),Ymat,chi(:,:,f,r),&
          'C','N','N',(-1d0,0d0),'ADD')
        !!
        call matrix_3_product(pre_force,XY_Ushiro,tmpmat1,XY_Mae,'N','N','N',&
          -Remez_alpha4(r)*dir_factor*XY_factor,'ADD')
        !!!!!!!!!!!!!!
        !! tmpmat2=Dchi^\dag.Y^\dag.\lambda - \chi.Y^\dag.Dlambda^\dag
        call matrix_3_product(tmpmat1,Dchi(:,:,f,r),Ymat,lambda(:,:,l,r),&
          'C','C','N')!,-XY_factor)
        call matrix_3_product(tmpmat1,chi(:,:,f,r),Ymat,Dlambda(:,:,l,r),&
          'N','C','C',(-1d0,0d0),'ADD')
        !!
        call matrix_3_product(pre_force,XY_Mae,tmpmat1,XY_Ushiro,'C','N','C',&
          -Remez_alpha4(r)*dir_factor*(-XY_factor),'ADD')
        !!!!!!!!!!!!!!!                 
                                        
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
      !! dY/dA contribution             
      elseif( ll_place > X_last ) then
        call calc_dYdA(XY_Mae,XY_Ushiro,XY_factor,f,l_place,ll_place,UMAT)
        !!!!!!!!!!!!!!!                 
        !! tmpmat1=X . lambda           
        call matrix_3_product(tmpmat1,Dchi(:,:,f,r),Xmat,lambda(:,:,l,r),&
          'C','N','N')
        call matrix_3_product(tmpmat1,chi(:,:,f,r),Xmat,Dlambda(:,:,l,r),&
          'N','N','C',(-1d0,0d0),'ADD')
        !!
        call matrix_3_product(pre_force,XY_Ushiro,tmpmat1,XY_Mae,'N','N','N',&
          -Remez_alpha4(r)*dir_factor*XY_factor,'ADD')
        !!!!!!!!!!!!!!!                 
        !! tmpmat1= lambda . Xdag       
        call matrix_3_product(tmpmat1,lambda(:,:,l,r),Xmat,Dchi(:,:,f,r),&
          'N','C','C')
        call matrix_3_product(tmpmat1,Dlambda(:,:,l,r),Xmat,chi(:,:,f,r),&
          'C','C','N',(-1d0,0d0),'ADD')
        !!
        call matrix_3_product(pre_force,XY_Mae,tmpmat1,XY_Ushiro,'C','N','C',&
          -Remez_alpha4(r)*dir_factor*(-XY_factor),'ADD')
        !!!!!!!!!!!!!!!                  
      endif                             
    enddo                                
  enddo                                   
enddo                                     
end subroutine calc_fermion_force_from_omega_m0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_fermion_force_from_massf&
    (pre_force,lambda,chi,Dlambda,Dchi,Umat,ll)
implicit none

complex(kind(0d0)), intent(out) :: pre_force(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: lambda(1:NMAT,1:NMAT,1:num_necessary_links,1:N_Remez4)
complex(kind(0d0)), intent(in) :: chi(1:NMAT,1:NMAT,1:num_necessary_faces,1:N_Remez4)
complex(kind(0d0)), intent(in) :: Dlambda(1:NMAT,1:NMAT,1:num_necessary_links,1:N_Remez4)
complex(kind(0d0)), intent(in) :: Dchi(1:NMAT,1:NMAT,1:num_necessary_faces,1:N_Remez4)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
integer, intent(in) :: ll

complex(kind(0d0)) :: Xmat(1:NMAT,1:NMAT),Ymat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: XY_Mae(1:NMAT,1:NMAT)
complex(kind(0d0)) :: XY_Ushiro(1:NMAT,1:NMAT)
complex(kind(0d0)) :: XY_factor
complex(kind(0d0)) :: dir_factor
integer :: k,l,f,l_place,ll_place,X_last,r!,kk

pre_force=(0d0,0d0)
do k=1,face_in_l(ll)%num_
  f=face_in_l(ll)%label_(k)

  !!!!!!!!!!!!!!!!!!
  !! place of ll in f
  do ll_place=1,links_in_f(f)%num_
    if( links_in_f(f)%link_labels_(ll_place) == ll ) exit
  enddo
  
  do l_place=1,links_in_f(f)%num_
    l=links_in_f(f)%link_labels_(l_place)
    !dir_factor=(0d0,-1d0)*dcmplx(&
      !dble(links_in_f(f)%link_dirs_(l_place)) * alpha_f(f) * beta_f(f) )
    dir_factor=(0d0,-1d0)*dcmplx( alpha_f(f) * beta_f(f) * mass_f )

    !!!!!!!!!!!!!!!!!!
    !! Xmat and Ymat
    call calc_XYmat(Xmat,Ymat,f,l_place,UMAT)

    !!!!!!!!!!!!!!!!!!
    !! last entry of Xmat
    if( links_in_f(f)%link_dirs_(l_place) == 1 ) then
      X_last = l_place-1
    else
      X_last = l_place
    endif

    do r=1,N_Remez4
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! dX/dA contribution
      if( ll_place <= X_last ) then
        call calc_dXdA(XY_Mae,XY_Ushiro,XY_factor,f,l_place,ll_place,UMAT)

        !!!!!!!!!!!!!!!
        !! tmpmat1=lambda.Y
        !! tmpmat2=Dlambda^\dag . Y
        call matrix_3_product(tmpmat1,lambda(:,:,l,r),Ymat,Dchi(:,:,f,r),&
          'N','N','C',XY_factor)
        call matrix_3_product(tmpmat1,Dlambda(:,:,l,r),Ymat,chi(:,:,f,r),&
          'C','N','N',-XY_factor,'ADD')
        !!!!!!!!!!!!!!!                 
        call matrix_3_product(pre_force,XY_Ushiro,tmpmat1,XY_Mae,'N','N','N',&
          -Remez_alpha4(r)*dir_factor,'ADD')
        !!!!!!!!!!!!!!!                 
                                        
        call matrix_3_product(tmpmat1,Dchi(:,:,f,r),Ymat,lambda(:,:,l,r),&
          'C','C','N',-XY_factor)
        call matrix_3_product(tmpmat1,chi(:,:,f,r),Ymat,Dlambda(:,:,l,r),&
          'N','C','C',XY_factor,'ADD')
        !!!!!!!!!!!!!!!                 
        call matrix_3_product(pre_force,XY_Mae,tmpmat1,XY_Ushiro,'C','N','C',&
          -Remez_alpha4(r)*dir_factor,'ADD')
        !!!!!!!!!!!!!!!                 
                                        
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
      !! dY/dA contribution             
      else                              
        call calc_dYdA(XY_Mae,XY_Ushiro,XY_factor,f,l_place,ll_place,UMAT)
        !!!!!!!!!!!!!!!                 
        !! tmpmat1=X . lambda           
        !! tmpmat2=X. Dlambda^dag       
        call matrix_3_product(tmpmat1,Dchi(:,:,f,r),Xmat,lambda(:,:,l,r),&
          'C','N','N',XY_factor)
        call matrix_3_product(tmpmat1,chi(:,:,f,r),Xmat,Dlambda(:,:,l,r),&
          'N','N','C',-XY_factor,'ADD')
                                        
        !!!!!!!!!!!!!!!                 
        call matrix_3_product(pre_force,XY_Ushiro,tmpmat1,XY_Mae,'N','N','N',&
          -Remez_alpha4(r)*dir_factor,'ADD')
        !!!!!!!!!!!!!!!                 
                                        
        !!!!!!!!!!!!!!!                 
        !! tmpmat1= lambda . Xdag       
        !! tmpmat2= Dlambda^dag . Xdag  
        call matrix_3_product(tmpmat1,lambda(:,:,l,r),Xmat,Dchi(:,:,f,r),&
          'N','C','C',-XY_factor)
        call matrix_3_product(tmpmat1,Dlambda(:,:,l,r),Xmat,chi(:,:,f,r),&
          'C','C','N',XY_factor,'ADD')
                                        
        !!!!!!!!!!!!!!!                 
        call matrix_3_product(pre_force,XY_Mae,tmpmat1,XY_Ushiro,'C','N','C',&
          -Remez_alpha4(r)*dir_factor,'ADD')
        !!!!!!!!!!!!!!!                  
      endif                             
    enddo                                
  enddo                                   
enddo                                     
end subroutine calc_fermion_force_from_massf
                                         
                                         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to return product of d/d\Phi D
subroutine prod_dDdPhi(dDdPhi_eta,dDdPhi_chi,eta_mat,chi_mat,UMAT)
implicit none                            
                                         
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links) 
complex(kind(0d0)), intent(in) :: eta_mat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: chi_mat(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), intent(out) :: dDdPhi_eta(1:NMAT,1:NMAT,1:num_sites,1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(out) :: dDdPhi_chi(1:NMAT,1:NMAT,1:num_faces,1:NMAT,1:NMAT,1:num_necessary_sites)
                                         
                                         
complex(kind(0d0)) :: tmp                
integer :: s,f,i,j,ii,jj                 
integer :: a,b,c,r                       
                                         
!! site: (s,a) <--> r=a+dimG*(s-1)       
!! link: (l,a) <--> r=a+dimG*(num_sites + l -1)
!! face: (f,a) <--> r=a+dimG*(num_sites + num_links + f -1 )
!! preparation                           
!dDdPhi_vec=(0d0,0d0)                    
dDdPhi_eta=(0d0,0d0)                     
dDdPhi_chi=(0d0,0d0)                     
                                         
if( p1 == 0 ) then                       
!! (1) Dirac from site                   
do s=1,num_sites                         
  do jj=1,NMAT                           
    i=jj                                 
    do ii=1,NMAT                         
      do j=1,NMAT                        
        dDdPhi_eta(i,j,s,ii,jj,s)= dDdPhi_eta(i,j,s,ii,jj,s) &
          + (-0.5d0,0d0)*dcmplx( alpha_s(s)*overall_factor ) &
            * eta_mat(ii,j,s)            
      enddo                              
    enddo                                
    !!                                   
    do ii=1,NMAT                         
      j=ii                               
      do i=1,NMAT                        
        dDdPhi_eta(i,j,s,ii,jj,s)= dDdPhi_eta(i,j,s,ii,jj,s) &
          - (-0.5d0,0d0)*dcmplx( alpha_s(s)*overall_factor ) &
            * eta_mat(i,jj,s)            
      enddo                              
    enddo                                
  enddo                                  
enddo                                    
endif                                    
                                         
!! (2) Dirac from link 1                 
! there is no contribution               
                                         
!! (3) Dirac from link 2                 
! there is no contribution               
                                         
if ( p4 == 0 ) then                      
!! (4) Dirac from face 1                 
do f=1,num_faces                         
  s=sites_in_f(f)%label_(1)              
  do jj=1,NMAT                           
    i=jj                                 
    do ii=1,NMAT                         
      do j=1,NMAT                        
        dDdPhi_chi(i,j,f,ii,jj,s)= dDdPhi_chi(i,j,f,ii,jj,s) &
          + (-2d0,0d0)*dcmplx( alpha_f(f)*overall_factor ) &
            * chi_mat(ii,j,f)
      enddo
    enddo
    !!
    do ii=1,NMAT
      j=ii
      do i=1,NMAT
        dDdPhi_chi(i,j,f,ii,jj,s)= dDdPhi_chi(i,j,f,ii,jj,s) &
          - (-2d0,0d0)*dcmplx( alpha_f(f)*overall_factor ) &
            * chi_mat(i,jj,f)
      enddo
    enddo
  enddo
enddo
endif


end subroutine prod_dDdPhi


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to return product of d/d\Phi D^\dagger
subroutine prod_dDdbPhi(dDdbPhi_lambda,lambda_mat,UMAT)
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links) 
complex(kind(0d0)), intent(in) :: lambda_mat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(out) :: dDdbPhi_lambda(1:NMAT,1:NMAT,1:num_links,1:NMAT,1:NMAT,1:num_necessary_sites)
!complex(kind(0d0)) :: dDdbPhi_vec(1:sizeD,1:dimG,1:num_sites)


complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmp

integer :: s,l
integer :: i,j,ii,jj
integer :: a,b

!! preparation
!dDdbPhi_vec=(0d0,0d0)
dDdbPhi_lambda=(0d0,0d0)


! (1) Dirac from site
! no contribution

!! (2) Dirac from link 1
! no contribution

if ( p3 == 0 ) then
!! (3) Dirac from link 2
do l=1,num_links
  call matrix_product(tmpmat1,Umat(:,:,l),lambda_mat(:,:,l),'C','N')
  call matrix_product(tmpmat2,lambda_mat(:,:,l),Umat(:,:,l))
  do jj=1,NMAT
    do ii=1,NMAT
      do j=1,NMAT
        do i=1,NMAT
          dDdbPhi_lambda(i,j,l,ii,jj,link_tip(l)) &
           = dDdbPhi_lambda(i,j,l,ii,jj,link_tip(l)) &
             + dcmplx( alpha_l(l) * overall_factor ) &
               * ( UMAT(i,jj,l) * tmpmat1(ii,j) &
                   - tmpmat2(i,jj) * dconjg( UMAT(j,ii,l) ))
          !!
          if ( i == jj ) then
            dDdbPhi_lambda(i,j,l,ii,jj,link_org(l)) &
             = dDdbPhi_lambda(i,j,l,ii,jj,link_org(l)) &
               + dcmplx( alpha_l(l) * overall_factor ) &
                 * lambda_mat(ii,j,l)
          endif
          !!
          if ( ii == j ) then
            dDdbPhi_lambda(i,j,l,ii,jj,link_org(l)) &
             = dDdbPhi_lambda(i,j,l,ii,jj,link_org(l)) &
               - dcmplx( alpha_l(l) * overall_factor ) &
                 * lambda_mat(i,jj,l)
          endif
        enddo
      enddo
    enddo
  enddo
enddo

endif

!! (4) Dirac from face 1
!   no contribution

!dDdbPhi_chi=dDdbPhi_chi*overall_factor
end subroutine prod_dDdbPhi


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to return product of d/dA D . vec
subroutine prod_dDdA(dDdA_eta,dDdA_lambda,dDdA_chi,&
    eta_mat,lambda_mat,chi_mat,UMAT,PhiMat)
implicit none

! d/dA_{ll,ii,jj} (D\Psi)_{s,i,j}=dDdA_eta(i,j,s,ii,jj,ll)
complex(kind(0d0)), intent(inout) :: dDdA_eta(1:NMAT,1:NMAT,1:num_sites,1:NMAT,1:NMAT,1:num_necessary_links)
! d/dA_{ll,ii,jj} (D\Psi)_{l,i,j}=dDdA_lambda(i,j,l,ii,jj,ll)
complex(kind(0d0)), intent(inout) :: dDdA_lambda(1:NMAT,1:NMAT,1:num_links,1:NMAT,1:NMAT,1:num_necessary_links)
! d/dA_{ll,ii,jj} (D\Psi)_{f,i,j}=dDdA_chi(i,j,f,ii,jj,ll)
complex(kind(0d0)),intent(inout) :: dDdA_chi(1:NMAT,1:NMAT,1:num_faces,1:NMAT,1:NMAT,1:num_necessary_links)

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links) 
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: eta_mat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: lambda_mat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: chi_mat(1:NMAT,1:NMAT,1:num_necessary_faces)

complex(kind(0d0)) :: trace,tmp
complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT),tmpmat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat3(1:NMAT,1:NMAT),tmpmat4(1:NMAT,1:NMAT)

!for test
!complex(kind(0d0)) :: T(1:NMAT,1:NMAT,1:dimG)
!complex(kind(0d0)) :: tmp_diffdiff_Omega(1:NMAT,1:NMAT,1:dimG,1:dimG)
!complex(kind(0d0)) :: tmp_diffdiff_Omega2(1:NMAT,1:NMAT,1:dimG,1:dimG)
!complex(kind(0d0)) :: tmp_diff_Omega(1:NMAT,1:NMAT,1:dimG)
!complex(kind(0d0)) :: tmp_diff_Omega2(1:NMAT,1:NMAT,1:dimG)
!complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT,1:num_faces)
!complex(kind(0d0)) :: Ufm(1:NMAT,1:NMAT,1:num_faces)


integer :: s,l,f,ll
integer :: i,j,k,nl,r,ii,jj
integer :: a,b,c,d,e


!! preparation
!dDdA_vec=(0d0,0d0)
dDdA_eta=(0d0,0d0)
dDdA_lambda=(0d0,0d0)
dDdA_chi=(0d0,0d0)

! for test
!call make_SUN_generators(T,NMAT)

!! (1) Dirac from site
!   no contribution

if ( p2==0 ) then
!! (2) Dirac from link 1
do s=1,num_sites
  do k=1,linkorg_to_s(s)%num_
    l=linkorg_to_s(s)%labels_(k)
    !if( l<=num_links ) then 
      ! tmpmat1 = lambda_l.U_l
      call matrix_product(tmpmat1,lambda_mat(:,:,l),UMAT(:,:,l))
      ! tmpmat2 = Ul^{-1}.lambda_l
      call matrix_product(tmpmat2,UMAT(:,:,l),lambda_mat(:,:,l),'C','N')
      do jj=1,NMAT
        do ii=1,NMAT
          do j=1,NMAT
            do i=1,NMAT
              dDdA_eta(i,j,s,ii,jj,l)= dDdA_eta(i,j,s,ii,jj,l) &
                - dcmplx(alpha_l(l))*dconjg(UMAT(jj,i,l))*tmpmat1(ii,j) &
                + dcmplx(alpha_l(l))*tmpmat2(i,jj)*UMAT(ii,j,l)
            enddo
          enddo
        enddo
      enddo
    !endif
  enddo
enddo    

do l=1,num_links
  s=link_tip(l)
  ! tmpmat2 = Ul.eta_l.Ul^{-1}
  call matrix_product(tmpmat1,UMAT(:,:,l),eta_mat(:,:,s))
  call matrix_product(tmpmat2,tmpmat1,UMAT(:,:,l),'N','C')
  do ii=1,NMAT
    do i=1,NMAT
      jj=ii
      j=i
      do k=1,NMAT
        dDdA_lambda(k,j,l,ii,k,l) = dDdA_lambda(k,j,l,ii,k,l) & 
          - dcmplx(alpha_l(l)) * tmpmat2(ii,j)
        dDdA_lambda(i,k,l,k,jj,l) = dDdA_lambda(i,k,l,k,jj,l) & 
          + dcmplx(alpha_l(l)) * tmpmat2(i,jj)
      enddo
    enddo
  enddo
enddo
endif 

if ( p3==0 ) then
!! (3) Dirac from link 2
do l=1,num_links
  s=link_tip(l)

  ! tmpmat1 = Ul.\bar\Phi.Ul^{-1}
  ! tmpmat2 = Ul.\bar\Phi.Ul^{-1}.lambda_l
  ! tmpmat3 = lambda_l.Ul.\bar\Phi.Ul^{-1}
  call matrix_product(tmpmat2,UMAT(:,:,l),PhiMat(:,:,s),'N','C')
  call matrix_product(tmpmat1,tmpmat2,UMAT(:,:,l),'N','C')
  call matrix_product(tmpmat2,tmpmat1,lambda_mat(:,:,l))
  call matrix_product(tmpmat3,lambda_mat(:,:,l),tmpmat1)
  do jj=1,NMAT
    do ii=1,NMAT
      do j=1,NMAT
        do i=1,NMAT
          if ( i==jj ) then
            dDdA_lambda(i,j,l,ii,jj,l) = dDdA_lambda(i,j,l,ii,jj,l) & 
              +(0d0,1d0)*dcmplx(alpha_l(l))*tmpmat2(ii,j)
          endif
          if ( j==ii ) then
            dDdA_lambda(i,j,l,ii,jj,l) = dDdA_lambda(i,j,l,ii,jj,l) & 
              +(0d0,1d0)*dcmplx(alpha_l(l))*tmpmat3(i,jj)
          endif
          dDdA_lambda(i,j,l,ii,jj,l) = dDdA_lambda(i,j,l,ii,jj,l) & 
            -(0d0,1d0)*dcmplx(alpha_l(l))*lambda_mat(i,jj,l)*tmpmat1(ii,j) &
            -(0d0,1d0)*dcmplx(alpha_l(l))*lambda_mat(ii,j,l)*tmpmat1(i,jj)
        enddo
      enddo
    enddo
  enddo
enddo
endif

!! (4) Dirac from face 1
!   no contribution

if( p5 == 0 ) then 
!! (5) Dirac from face 2
! preparation
  !if ( p5_test == 0 ) then 
    !if( m_omega == 0 ) then 
      !call calc_fermion_force_from_omega_m0&
        !(dDdA_lambda,dDdA_chi,UMAT,lambda_mat,chi_mat)
    !elseif( m_omega == -1 ) then 
      !call calc_fermion_force_from_omega_adm&
        !(tmp_force,lambda,chi,Dlambda,Dchi,Umat,Uf,ll)
        !(dDdA_lambda,dDdA_chi,UMAT,lambda_mat,chi_mat)
    !else
      !call calc_fermion_force_from_omega&
        !(dDdA_lambda,dDdA_chi,UMAT,lambda_mat,chi_mat)
    !endif
  !else
  !  call calc_fermion_force_from_omega_test&
  !    (dDdA_lambda,dDdA_chi,UMAT,lambda_mat,chi_mat)
  !endif
endif 
dDdA_eta=dDdA_eta*dcmplx(overall_factor)
dDdA_lambda=dDdA_lambda*dcmplx(overall_factor)
dDdA_chi=dDdA_chi*dcmplx(overall_factor)

!! Make the force traceless
!trace=(0d0,0d0)
!do s=1,num_sites
!  do j=1,NMAT
!    do i=1,NMAT
!      do ll=1,num_links
!        do ii=1,NMAT
!          trace=trace+dDdA_eta(i,j,s,ii,ii,ll)
!        enddo
!        trace=trace/cmplx(dble(NMAT))
!        do ii=1,NMAT
!          dDdA_eta(i,j,s,ii,ii,ll)=dDdA_eta(i,j,s,ii,ii,ll)-trace
!        enddo
!      enddo
!    enddo
!  enddor
!enddo
!
!trace=(0d0,0d0)
!do l=1,num_links
!  do j=1,NMAT
!    do i=1,NMAT
!      do ll=1,num_links
!        do ii=1,NMAT
!          trace=trace+dDdA_lambda(i,j,l,ii,ii,ll)
!        enddo
!        trace=trace/cmplx(dble(NMAT))
!        do ii=1,NMAT
!          dDdA_lambda(i,j,l,ii,ii,ll)=dDdA_lambda(i,j,l,ii,ii,ll)-trace
!        enddo
!      enddo
!    enddo
!  enddo
!enddo
!
!trace=(0d0,0d0)
!do f=1,num_faces
!  do j=1,NMAT
!    do i=1,NMAT
!      do ll=1,num_links
!        do ii=1,NMAT
!          trace=trace+dDdA_chi(i,j,f,ii,ii,ll)
!        enddo
!        trace=trace/cmplx(dble(NMAT))
!        do ii=1,NMAT
!          dDdA_chi(i,j,f,ii,ii,ll)=dDdA_chi(i,j,f,ii,ii,ll)-trace
!        enddo
!      enddo
!    enddo
!  enddo
!enddo
    

end subroutine prod_dDdA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine calc_fermion_force_from_omega_org&
!    (dDdA_lambda,dDdA_chi,UMAT,lambda_mat,chi_mat)
!implicit none
!
!complex(kind(0d0)), intent(inout) :: dDdA_lambda&
!  (1:NMAT,1:NMAT,1:num_links, 1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)), intent(inout) :: dDdA_chi&
!  (1:NMAT,1:NMAT,1:num_faces, 1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)), intent(in) :: lambda_mat(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)), intent(in) :: chi_mat(1:NMAT,1:NMAT,1:num_necessary_faces)
!
!complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: Ufm(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: Cosinv(1:NMAT,1:NMAT)
!!complex(kind(0d0)) :: Sinmat(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: Omega_mat(1:NMAT,1:NMAT) ! traceless Omega
!complex(kind(0d0)) :: Omega_org(1:NMAT,1:NMAT) ! trace-ful Omega
!complex(kind(0d0)) :: dOmegadA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
!complex(kind(0d0)) :: dCosinvdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
!complex(kind(0d0)) :: Xmat(1:NMAT,1:NMAT),Ymat(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: UXmat(1:NMAT,1:NMAT),YUmat(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: dUfdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
!complex(kind(0d0)) :: tmp_diffmat(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
!complex(kind(0d0)) :: med_diffmat(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
!complex(kind(0d0)) :: diffmat(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
!complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: tmpmat3(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: dir_factor
!integer :: i,j,k,l,ll,f,a,b,l_place,ll_place,X_last
!complex(kind(0d0)) :: im_over_2
!
!im_over_2=(0d0,0.5d0)*dcmplx(dble(m_omega))
!
!
!do f=1,num_faces
!  call Make_face_variable(Uf(:,:),f,UMAT) 
!  if( m_omega == 1) then 
!    Ufm=Uf
!  else
!    call matrix_power(Ufm,Uf(:,:),m_omega)
!  endif
!  do i=1,NMAT
!    do j=1,NMAT
!      Cosinv(i,j) = Ufm(i,j) + dconjg(Ufm(j,i))
!      tmpmat1(i,j) =  Ufm(i,j) - dconjg(Ufm(j,i)) 
!    enddo
!  enddo
!  !! Cos^{-1}
!  !call HermitianMatrix_inverse(Cosinv)
!  call Matrix_inverse(Cosinv)
!  !! Omega
!  call ZHEMM('L','U',NMAT,NMAT, (0d0,-2d0)/dcmplx(dble(m_omega)), &
!  Cosinv, NMAT, &
!  tmpmat1, NMAT, &
!  (0d0,0d0), Omega_mat, NMAT)
!  if ( NMAT > 2 ) then 
!    Omega_org=Omega_mat
!    call make_matrix_traceless(Omega_mat)
!  endif
!
!  do l_place=1,links_in_f(f)%num_
!    l=links_in_f(f)%link_labels_(l_place)
!    dir_factor=dcmplx(dble(links_in_f(f)%link_dirs_(l_place)))&
!      *(0d0,2d0)/dcmplx(dble(m_omega))
!
!    call calc_Xmat(Xmat,f,l_place,UMAT)
!    call calc_Ymat(Ymat,f,l_place,UMAT)
!
!    if( links_in_f(f)%link_dirs_(l_place) == 1 ) then
!      X_last = l_place-1
!    else
!      X_last = l_place
!    endif
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ! For test
!    !tmpmat4=Cosinv
!    !call make_unit_matrix(Cosinv)
!    !call make_unit_matrix(Omega_mat)
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    !! dUf/dA,  d(Omega)/dA, d(Cos^{-1})/dA
!    do ll_place=1,links_in_f(f)%num_
!      ll=links_in_f(f)%link_labels_(ll_place)
!      !! dUf/dA
!      call calc_dUfdA(dUfdA,f,ll_place,UMAT)
!      !! d(Omega)/dA
!      call calc_dOmegadA_dCosinvdA(dOmegadA,dCosinvdA,dUfdA,Uf,Omega_org,Cosinv)
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      ! For test
!      !call calc_dOmegadA_dCosinvdA(dOmegadA,dCosinvdA,dUfdA,Uf,Omega_mat,tmpmat4)
!      !dCosinvdA=(0d0,0d0)
!      !dOmegadA=(0d0,0d0)
!      
!
!      do k=0,m_omega-1
!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        !! Preparation
!        !! UXmat
!        if(k==0) then
!          UXmat=Xmat
!        elseif(k==1) then
!          call matrix_product(UXmat,Uf,Xmat)
!        else
!          call matrix_power(tmpmat1,Uf,k)
!          call matrix_product(UXmat,tmpmat1,Xmat)
!        endif
!        !! YUmat
!        if(m_omega-k-1==0) then
!          YUmat=Ymat
!        elseif(m_omega-k-1==1) then
!          call matrix_product(YUmat,Ymat,Uf)
!        else
!          call matrix_power(tmpmat1,Uf,m_omega-k-1)
!          call matrix_product(YUmat,Ymat,tmpmat1)
!        endif
!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        !! d(Uk)/dA.Xmat + U^k d(X)/dA 
!        if( k==0 .and. ll_place <= X_last ) then
!          call calc_dXdA(diffmat,f,l_place,ll_place,UMAT)
!        elseif( k>0 ) then
!          call calc_dUfkdA(med_diffmat,k,dUfdA,Uf)
!          do b=1,NMAT
!            do a=1,NMAT
!              call matrix_product(diffmat(:,:,a,b),&
!                med_diffmat(:,:,a,b),Xmat)
!            enddo
!          enddo
!          if( ll_place <= X_last) then
!            call calc_dXdA(med_diffmat,f,l_place,ll_place,UMAT)
!            call matrix_power(tmpmat1,Uf,k)
!            do b=1,NMAT
!              do a=1,NMAT
!                call matrix_product(diffmat(:,:,a,b),&
!                  tmpmat1,med_diffmat(:,:,a,b),'N','N',(1d0,0d0),'ADD')
!              enddo
!            enddo
!          endif
!        endif
!
!        if( ll_place <= X_last .or. k>0 ) then
!          do b=1,NMAT
!            do a=1,NMAT
!              !!!!!!!!!!!!!!!!!!!!!!
!              call matrix_3_product(tmpmat1,diffmat(:,:,a,b),lambda_mat(:,:,l),YUmat)
!              call matrix_3_product(tmpmat2,YUmat,lambda_mat(:,:,l),diffmat(:,:,b,a),'C','N','C')
!
!              tmpmat3=tmpmat1+tmpmat2
!              call matrix_product(tmpmat3,&
!                tmpmat1-tmpmat2,&
!                Omega_mat,'N','N', im_over_2,'ADD')
!              call matrix_product(tmpmat1,Cosinv,tmpmat3)
!
!
!              dDdA_chi(:,:,f,a,b,ll)=dDdA_chi(:,:,f,a,b,ll)&
!                -dir_factor &
!                * dcmplx(alpha_f(f)*beta_f(f)) &
!                * tmpmat1
!                !* (tmpmat1+tmpmat2)
!     
!              !!!!!!!!!!!!!!!!!!!!!!
!              call matrix_product(tmpmat1, chi_mat(:,:,f),Cosinv)
!              call matrix_3_product(tmpmat3,YUmat,tmpmat1,diffmat(:,:,a,b))
!              call matrix_3_product(tmpmat3,diffmat(:,:,b,a),tmpmat1,YUmat,&
!                'C','N','C',(1d0,0d0),'ADD')
!
!              call matrix_3_product(tmpmat1, Omega_mat, chi_mat(:,:,f),Cosinv) 
!              call matrix_3_product(tmpmat3,YUmat,tmpmat1,diffmat(:,:,a,b)&
!                ,'N','N','N', im_over_2,'ADD')
!              call matrix_3_product(tmpmat3,diffmat(:,:,b,a),tmpmat1,YUmat&
!                ,'C','N','C', (-1d0,0d0)*im_over_2,'ADD')
!
!              dDdA_lambda(:,:,l,a,b,ll)=dDdA_lambda(:,:,l,a,b,ll)&
!                +dir_factor &
!                * dcmplx(alpha_f(f)*beta_f(f)) &
!                * tmpmat3
!                !* (tmpmat1+tmpmat2)
!            enddo
!          enddo
!        endif
!
!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        !! Y.d(U^{m-k-1})/dA + d(Ymat)/dA.U^{m-k-1}
!        if( m_omega-k-1==0 .and. ll_place > X_last) then
!          call calc_dYdA(diffmat,f,l_place,ll_place,UMAT)
!        elseif( m_omega-k-1 > 0 ) then
!          call calc_dUfkdA(med_diffmat,m_omega-k-1,dUfdA,Uf)
!          do b=1,NMAT
!            do a=1,NMAT
!              call matrix_product(diffmat(:,:,a,b),&
!                Ymat,med_diffmat(:,:,a,b))
!            enddo
!          enddo
!          if( ll_place > X_last) then
!            call calc_dYdA(med_diffmat,f,l_place,ll_place,UMAT)
!            call matrix_power(tmpmat1,Uf,m_omega-k-1)
!            do b=1,NMAT
!              do a=1,NMAT
!                call matrix_product(diffmat(:,:,a,b),&
!                  med_diffmat(:,:,a,b),tmpmat1,'N','N',(1d0,0d0),'ADD')
!              enddo
!            enddo
!          endif
!        endif
!
!        if( m_omega-k-1>0 .or. ll_place > X_last ) then
!          do b=1,NMAT
!            do a=1,NMAT
!              !!!!!!!!!!!!!!!!!!!!!!
!              call matrix_3_product(tmpmat1,UXmat,lambda_mat(:,:,l),diffmat(:,:,a,b))
!              call matrix_3_product(tmpmat2,diffmat(:,:,b,a),lambda_mat(:,:,l),UXmat,'C','N','C')
!              tmpmat3=tmpmat1+tmpmat2
!
!              call matrix_product(tmpmat3,tmpmat1-tmpmat2,Omega_mat,'N','N',&
!                im_over_2,'ADD')
!
!              call matrix_product(tmpmat1,Cosinv,tmpmat3)
!
!              dDdA_chi(:,:,f,a,b,ll)=dDdA_chi(:,:,f,a,b,ll)&
!                -dir_factor &
!                * dcmplx(alpha_f(f)*beta_f(f)) &
!                * tmpmat1
!                !* (tmpmat1+tmpmat2)
!    
!
!              !!!!!!!!!!!!!!!!!!!!!!!
!              call matrix_product(tmpmat1,chi_mat(:,:,f),Cosinv)
!              call matrix_3_product(tmpmat3,diffmat(:,:,a,b),tmpmat1,UXmat)
!              call matrix_3_product(tmpmat3,UXmat,tmpmat1,diffmat(:,:,b,a),&
!                'C','N','C',(1d0,0d0),'ADD')
!
!              call matrix_3_product(tmpmat1,Omega_mat,chi_mat(:,:,f),Cosinv) 
!              call matrix_3_product(tmpmat3,diffmat(:,:,a,b),tmpmat1,UXmat,&
!                'N','N','N', im_over_2,'ADD')
!              call matrix_3_product(tmpmat3,UXmat,tmpmat1,diffmat(:,:,b,a),&
!                'C','N','C', (-1d0,0d0)*im_over_2,'ADD')
!
!              dDdA_lambda(:,:,l,a,b,ll)=dDdA_lambda(:,:,l,a,b,ll)&
!                +dir_factor &
!                * dcmplx(alpha_f(f)*beta_f(f)) &
!                * tmpmat3
!                !* (tmpmat1+tmpmat2)
!            enddo
!          enddo
!        endif
!
!        !!!!!!!!!!!!
!        ! d(Omega)/dA and d(Cosinv)/dA
!        !!!!  DF_chi 
!        call matrix_3_product(tmpmat1,UXmat,lambda_mat(:,:,l),YUmat)
!        call matrix_3_product(tmpmat2,YUmat,lambda_mat(:,:,l),UXmat,'C','N','C')
!
!        tmpmat3=tmpmat1-tmpmat2 ! tmp1-tmp2
!        tmpmat1=tmpmat1+tmpmat2 ! tpm1+tmp2
!
!        ! tmpmat1 = UX.lambda.YU + YU^dag.lambda.XU^dag 
!        !           +im/2 ( UX.lambda.YU - YU^dag.lambda.XU^dag ).Omega
!        call matrix_product(tmpmat1,tmpmat3,Omega_mat, 'N','N',im_over_2,'ADD')
!
!        ! tmpmat2 = im/2 Cos^{-1}.(UX.lambda.YU - YU^dag.lambda.XU^dag)
!        call matrix_product(tmpmat2, Cosinv,tmpmat3, 'N','N',im_over_2)
!
!        do b=1,NMAT
!          do a=1,NMAT
!            ! Cos^{-1}.im/2( UX.lambda.YU-YU^dag.lambda.UX\dag ).dOmega/dA
!            call matrix_product(tmpmat3, tmpmat2, dOmegadA(:,:,a,b))
!            !!
!            call matrix_product(tmpmat3, dCosinvdA(:,:,a,b), tmpmat1,&
!              'N','N',(1d0,0d0),'ADD')
!
!            dDdA_chi(:,:,f,a,b,ll)=dDdA_chi(:,:,f,a,b,ll)&
!              -dir_factor &
!              * dcmplx(alpha_f(f)*beta_f(f)) &
!              * tmpmat3
!          enddo
!        enddo
!        
!        !!!!  DF_lambda
!        ! dOmega/dA 
!        call matrix_3_product(tmpmat1,chi_mat(:,:,f),Cosinv,UXmat)!,&
!          !'N','N','N',im_over_2)
!        call matrix_3_product(tmpmat2,chi_mat(:,:,f),Cosinv,YUmat, 'N','N','C')!,im_over_2)
!
!        do b=1,NMAT
!          do a=1,NMAT
!            call matrix_3_product(tmpmat3, &
!              YUmat,dOmegadA(:,:,a,b),tmpmat1,'N','N','N',im_over_2)
!            call matrix_3_product(tmpmat3, &
!              UXmat,dOmegadA(:,:,a,b),tmpmat2,'C','N','N',(-1d0,0d0)*im_over_2,'ADD')
!
!            dDdA_lambda(:,:,l,a,b,ll)=dDdA_lambda(:,:,l,a,b,ll)&
!              +dir_factor &
!              * dcmplx(alpha_f(f)*beta_f(f)) &
!              * tmpmat3 
!          enddo
!        enddo
!        ! dCosinv/dA
!        call matrix_product(tmpmat1,YUmat,chi_mat(:,:,f))
!        call matrix_3_product(tmpmat1,YUmat,Omega_mat,chi_mat(:,:,f),&
!          'N','N','N',im_over_2,'ADD')
!
!        call matrix_product(tmpmat2,UXmat,chi_mat(:,:,f),'C','N')
!        call matrix_3_product(tmpmat2,UXmat,Omega_mat,chi_mat(:,:,f),&
!          'C','N','N', (-1d0,0d0)*im_over_2,'ADD')
!
!        do b=1,NMAT
!          do a=1,NMAT
!            call matrix_3_product(tmpmat3,tmpmat1,dCosinvdA(:,:,a,b),UXmat)
!            call matrix_3_product(tmpmat3,tmpmat2,dCosinvdA(:,:,a,b),YUmat,&
!              'N','N','C',(1d0,0d0),'ADD')
!
!            dDdA_lambda(:,:,l,a,b,ll)=dDdA_lambda(:,:,l,a,b,ll)&
!              +dir_factor &
!              * dcmplx(alpha_f(f)*beta_f(f)) &
!              * tmpmat3 
!          enddo
!        enddo
!
!      enddo ! k
!    enddo ! ll_place
!  enddo
!enddo
!
!end subroutine calc_fermion_force_from_omega_org


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine calc_fermion_force_from_omega_m0_org&
!    (dDdA_lambda,dDdA_chi,UMAT,lambda_mat,chi_mat)
!implicit none
!
!complex(kind(0d0)), intent(inout) :: dDdA_lambda&
!  (1:NMAT,1:NMAT,1:num_links, 1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)), intent(inout) :: dDdA_chi&
!  (1:NMAT,1:NMAT,1:num_faces, 1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)), intent(in) :: lambda_mat(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)), intent(in) :: chi_mat(1:NMAT,1:NMAT,1:num_necessary_faces)
!
!complex(kind(0d0)) :: Xmat(1:NMAT,1:NMAT),Ymat(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: diffmat(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
!complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT),tmpmat2(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: dir_factor
!integer :: i,j,k,l,ll,f,a,b,l_place,ll_place,X_last
!
!
!do f=1,num_faces
!  do l_place=1,links_in_f(f)%num_
!    l=links_in_f(f)%link_labels_(l_place)
!    dir_factor=dcmplx(dble(links_in_f(f)%link_dirs_(l_place)))*(0d0,1d0)
!
!    call calc_Xmat(Xmat,f,l_place,UMAT)
!    call calc_Ymat(Ymat,f,l_place,UMAT)
!
!    if( links_in_f(f)%link_dirs_(l_place) == 1 ) then
!      X_last = l_place-1
!    else
!      X_last = l_place
!    endif
!
!    do ll_place=1,X_last
!      ll=links_in_f(f)%link_labels_(ll_place)
!      call calc_dXdA(diffmat,f,l_place,ll_place,UMAT)
!      
!      do b=1,NMAT
!        do a=1,NMAT
!          !!!!!!!!!!!!!!!!!!!!!!
!          call matrix_3_product(tmpmat1,diffmat(:,:,a,b),lambda_mat(:,:,l),Ymat)
!          call matrix_3_product(tmpmat2,Ymat,lambda_mat(:,:,l),diffmat(:,:,b,a),'C','N','C')
!          dDdA_chi(:,:,f,a,b,ll)=dDdA_chi(:,:,f,a,b,ll)&
!            -dir_factor * dcmplx(alpha_f(f)*beta_f(f)) &
!            * (tmpmat1+tmpmat2)
!
!          !!!!!!!!!!!!!!!!!!!!!!
!          call matrix_3_product(tmpmat1,Ymat,chi_mat(:,:,f),diffmat(:,:,a,b),'N','N','N')
!          call matrix_3_product(tmpmat2,diffmat(:,:,b,a),chi_mat(:,:,f),Ymat,'C','N','C')
!          dDdA_lambda(:,:,l,a,b,ll)=dDdA_lambda(:,:,l,a,b,ll)&
!            +dir_factor * dcmplx(alpha_f(f)*beta_f(f)) &
!            * (tmpmat1+tmpmat2)
!        enddo
!      enddo
!    enddo
!
!    do ll_place=X_last+1,links_in_f(f)%num_
!      ll=links_in_f(f)%link_labels_(ll_place)
!      call calc_dYdA(diffmat,f,l_place,ll_place,UMAT)
!
!      do b=1,NMAT
!        do a=1,NMAT
!          !!!!!!!!!!!!!!!!!!!!!!
!          call matrix_3_product(tmpmat1,Xmat,lambda_mat(:,:,l),diffmat(:,:,a,b))
!          call matrix_3_product(tmpmat2,diffmat(:,:,b,a),lambda_mat(:,:,l),Xmat,'C','N','C')
!          dDdA_chi(:,:,f,a,b,ll)=dDdA_chi(:,:,f,a,b,ll)&
!            -dir_factor * dcmplx(alpha_f(f)*beta_f(f)) &
!            * (tmpmat1+tmpmat2)
!
!          !!!!!!!!!!!!!!!!!!!!!!
!          call matrix_3_product(tmpmat1,diffmat(:,:,a,b),chi_mat(:,:,f),Xmat,'N','N','N')
!          call matrix_3_product(tmpmat2,Xmat,chi_mat(:,:,f),diffmat(:,:,b,a),'C','N','C')
!          dDdA_lambda(:,:,l,a,b,ll)=dDdA_lambda(:,:,l,a,b,ll)&
!            +dir_factor * dcmplx(alpha_f(f)*beta_f(f)) &
!            * (tmpmat1+tmpmat2)
!        enddo
!      enddo
!    enddo
!
!  enddo
!enddo
!
!end subroutine calc_fermion_force_from_omega_m0_org

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine calc_fermion_force_from_omega_adm_org&
!    (dDdA_lambda,dDdA_chi,UMAT,lambda_mat,chi_mat)
!implicit none
!
!complex(kind(0d0)), intent(inout) :: dDdA_lambda&
!  (1:NMAT,1:NMAT,1:num_links, 1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)), intent(inout) :: dDdA_chi&
!  (1:NMAT,1:NMAT,1:num_faces, 1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)), intent(in) :: lambda_mat(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)), intent(in) :: chi_mat(1:NMAT,1:NMAT,1:num_necessary_faces)
!
!complex(kind(0d0)) :: Xmat(1:NMAT,1:NMAT),Ymat(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: diffmat(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
!complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT),tmpmat2(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: tmpmat3(1:NMAT,1:NMAT),tmpmat4(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: dBdA(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: sinU(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: dUfdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
!complex(kind(0d0)) :: dsinUdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
!complex(kind(0d0)) :: dir_factor
!complex(kind(0d0)) :: Bval,trace,tmp
!integer :: i,j,k,l,ll,f,a,b,l_place,ll_place,X_last
!
!
!do f=1,num_faces
!  call Make_face_variable(tmpmat1,f,UMAT)
!  do j=1,NMAT
!    do i=1,NMAT
!      sinU(i,j)=tmpmat1(i,j)-dconjg( tmpmat1(j,i) )
!    enddo
!  enddo
!
!  Bval=(1d0,0d0)
!  do i=1,NMAT
!    Bval=Bval - (1d0,0d0)/(e_max*e_max) * dcmplx( 2d0 - 2d0*dble(tmpmat1(i,i)) )
!  enddo
!
!  do l_place=1,links_in_f(f)%num_
!    l=links_in_f(f)%link_labels_(l_place)
!    dir_factor=dcmplx(dble(links_in_f(f)%link_dirs_(l_place)))*(0d0,1d0)&
!       * alpha_f(f) * beta_f(f)
!
!    call calc_Xmat(Xmat,f,l_place,UMAT)
!    call calc_Ymat(Ymat,f,l_place,UMAT)
!
!
!    if( links_in_f(f)%link_dirs_(l_place) == 1 ) then
!      X_last = l_place-1
!    else
!      X_last = l_place
!    endif
!
!
!    do ll_place=1,links_in_f(f)%num_
!      ll=links_in_f(f)%link_labels_(ll_place)
!
!      !!!!!!!!!!!!!!!!!!
!      !! Prepare dB/dA, dUf/DA
!      call calc_dUfdA(dUfdA,f,ll_place,UMAT)
!      do j=1,NMAT
!        do i=1,NMAT
!          do b=1,NMAT
!            do a=1,NMAT
!              dSinUdA(i,j,a,b) = dUfdA(i,j,a,b) - dconjg(dUfdA(j,i,b,a))
!            enddo
!          enddo
!        enddo
!      enddo
!
!      dBdA=(0d0,0d0)
!      do b=1,NMAT
!        do a=1,NMAT
!          do i=1,NMAT
!            dBdA(a,b)=dBdA(a,b)+(dUfdA(i,i,a,b) + dconjg(dUfdA(i,i,b,a)))/(e_max*e_max)
!          enddo
!        enddo
!      enddo
!
!      !DF_chi(:,:,f)
!      ! dUf/dA
!      trace=(0d0,0d0)
!      call matrix_product(tmpmat2,lambda_mat(:,:,l),Ymat)
!      do j=1,NMAT
!        do i=1,NMAT
!          trace=trace+( Xmat(i,j) * tmpmat2(j,i) )
!        enddo
!      enddo
!      call matrix_product(tmpmat2,lambda_mat(:,:,l),Xmat,'N','C')
!      do j=1,NMAT
!        do i=1,NMAT
!          trace=trace-( dconjg(Ymat(j,i)) * tmpmat2(j,i) )
!        enddo
!      enddo
!      trace=trace/(Bval*Bval*e_max*e_max)
!
!      do b=1,NMAT
!        do a=1,NMAT
!          dDdA_chi(:,:,f,a,b,ll)=dDdA_chi(:,:,f,a,b,ll)&
!            + dir_factor * trace * dSinUdA(:,:,a,b)
!        enddo
!      enddo
!
!      ! dB/dA
!      tmp=(1d0,0d0)/(Bval*Bval)
!      call matrix_3_product(tmpmat1, Xmat,lambda_mat(:,:,l),Ymat,'N','N','N',tmp)
!      call matrix_3_product(tmpmat1, Ymat,lambda_mat(:,:,l),Xmat,'C','N','C',tmp,'ADD')
!      trace=trace*(-2d0,0d0)/Bval
!      tmpmat1=tmpmat1 + trace*SinU
!
!      do b=1,NMAT
!        do a=1,NMAT
!          dDdA_chi(:,:,f,a,b,ll)=dDdA_chi(:,:,f,a,b,ll)&
!            +dir_factor * tmpmat1*dBdA(a,b)
!        enddo
!      enddo
!
!      !DF_lambda(:,:,l)
!      ! dUf/dA
!      call matrix_product(tmpmat2, Ymat,Xmat,'N','N')!, tmp)
!      do j=1,NMAT
!        do i=1,NMAT
!          tmpmat1(i,j)=tmpmat2(i,j)-dconjg(tmpmat2(j,i))
!        enddo
!      enddo
!
!      do b=1,NMAT
!        do a=1,NMAT
!          trace=(0d0,0d0)
!          do i=1,NMAT
!            do j=1,NMAT
!              trace=trace+dSinUdA(i,j,a,b)*chi_mat(j,i,f)
!            enddo
!          enddo
!          trace=trace/(Bval*Bval*e_max*e_max)
!
!
!          dDdA_lambda(:,:,l,a,b,ll)=dDdA_lambda(:,:,l,a,b,ll)&
!            -dir_factor * trace*tmpmat1
!        enddo
!      enddo
!
!      ! At present
!      !  tmpmat1=(Y X-X^dag.Y^dag) 
!      ! dB/dA
!      trace=(0d0,0d0)
!      do i=1,NMAT
!        do j=1,NMAT
!          trace=trace+SinU(i,j)*chi_mat(j,i,f)
!        enddo
!      enddo
!      trace=trace*(-2d0,0d0)/(Bval*Bval*Bval*e_max*e_max)
!      tmpmat1=trace*tmpmat1
!
!      tmp=(1d0,0d0)/(Bval*Bval)
!      call matrix_3_product(tmpmat1,Ymat,chi_mat(:,:,f),Xmat,'N','N','N',tmp,'ADD')
!      call matrix_3_product(tmpmat1,Xmat,chi_mat(:,:,f),Ymat,'C','N','C',tmp,'ADD')
!
!      do b=1,NMAT
!        do a=1,NMAT
!          dDdA_lambda(:,:,l,a,b,ll)=dDdA_lambda(:,:,l,a,b,ll)&
!            -dir_factor * tmpmat1*dBdA(a,b)
!        enddo
!      enddo
!
!      !! dX/dA
!      if( ll_place <= X_last ) then
!        call calc_dXdA(diffmat,f,l_place,ll_place,UMAT)
!      
!        call matrix_product(tmpmat1,lambda_mat(:,:,l),Ymat,'N','N')
!        call matrix_product(tmpmat2,YMAT,lambda_mat(:,:,l),'C','N')
!
!        do b=1,NMAT
!          do a=1,NMAT
!            tmp=(-1d0,0d0)/Bval
!            call matrix_product(tmpmat3, &
!              diffmat(:,:,a,b),tmpmat1,'N','N',tmp)
!            call matrix_product(tmpmat3, &
!              tmpmat2,diffmat(:,:,b,a),'N','C',tmp,'ADD')
!
!            trace=(0d0,0d0)
!            do j=1,NMAT
!              do i=1,NMAT
!                trace=trace &
!                  + diffmat(i,j,a,b)*tmpmat1(j,i) &
!                  - tmpmat2(i,j)*dconjg(diffmat(i,j,b,a)) !HERE
!              enddo
!            enddo
!            trace=trace/(Bval*Bval*e_max*e_max)
!
!            tmpmat3=tmpmat3+trace*SinU
!
!            !call make_matrix_traceless(tmpmat3)
!
!            dDdA_chi(:,:,f,a,b,ll)=dDdA_chi(:,:,f,a,b,ll)&
!              +dir_factor * tmpmat3
!          enddo
!        enddo
!        
!
!        ! DF_lambda
!        trace=(0d0,0d0)
!        do j=1,NMAT
!          do i=1,NMAT
!            trace=trace+SinU(i,j)*chi_mat(j,i,f)
!          enddo
!        enddo
!        trace=trace/(Bval*Bval*e_max*e_max)
!
!        do b=1,NMAT
!          do a=1,NMAT
!            tmp=(-1d0,0d0)/Bval
!            call matrix_3_product(tmpmat3,&
!              Ymat,chi_mat(:,:,f),diffmat(:,:,a,b),'N','N','N',tmp)
!            call matrix_3_product(tmpmat3,&
!              diffmat(:,:,b,a),chi_mat(:,:,f),Ymat,'C','N','C',tmp,'ADD')
!
!            call matrix_product(tmpmat3,&
!              Ymat,diffmat(:,:,a,b),'N','N',trace,'ADD')
!            call matrix_product(tmpmat3,&
!              diffmat(:,:,b,a),Ymat,'C','C',(-1d0,0d0)*trace,'ADD')
!
!            dDdA_lambda(:,:,l,a,b,ll)=dDdA_lambda(:,:,l,a,b,ll)&
!              -dir_factor * tmpmat3
!          enddo
!        enddo
!
!      !! dY/dA
!      else
!        call calc_dYdA(diffmat,f,l_place,ll_place,UMAT)
!
!        call matrix_product(tmpmat1,Xmat,lambda_mat(:,:,l),'N','N')
!        call matrix_product(tmpmat2,lambda_mat(:,:,l),Xmat,'N','C')
!
!        do b=1,NMAT
!          do a=1,NMAT
!            tmp=(-1d0,0d0)/Bval
!            call matrix_product(tmpmat3, &
!              tmpmat1,diffmat(:,:,a,b),'N','N',tmp)
!            call matrix_product(tmpmat3, &
!              diffmat(:,:,b,a),tmpmat2,'C','N',tmp,'ADD')
!
!            trace=(0d0,0d0)
!            do j=1,NMAT
!              do i=1,NMAT
!                trace=trace&
!                  +tmpmat1(i,j)*diffmat(j,i,a,b) &
!                  -dconjg(diffmat(j,i,b,a))*tmpmat2(j,i) !HERE
!              enddo
!            enddo
!            trace=trace/(Bval*Bval*e_max*e_max)
!
!            tmpmat3=tmpmat3+trace*SinU
!
!            dDdA_chi(:,:,f,a,b,ll)=dDdA_chi(:,:,f,a,b,ll)&
!              +dir_factor * tmpmat3
!          enddo
!        enddo
!
!        ! DF_lambda
!        trace=(0d0,0d0)
!        do j=1,NMAT
!          do i=1,NMAT
!            trace=trace+SinU(i,j)*chi_mat(j,i,f)
!          enddo
!        enddo
!        trace=trace/(e_max*e_max*Bval*Bval)
!
!        do b=1,NMAT
!          do a=1,NMAT
!            tmp=(-1d0,0d0)/Bval
!            call matrix_3_product(tmpmat3, &
!              diffmat(:,:,a,b),chi_mat(:,:,f),Xmat,'N','N','N',tmp)
!            call matrix_3_product(tmpmat3, &
!              Xmat,chi_mat(:,:,f),diffmat(:,:,b,a),'C','N','C',tmp,'ADD')
!
!            call matrix_product(tmpmat3, &
!              diffmat(:,:,a,b),Xmat,'N','N',trace,'ADD')
!            call matrix_product(tmpmat3, &
!              Xmat,diffmat(:,:,b,a),'C','C',(-1d0,0d0)*trace,'ADD')
!
!            dDdA_lambda(:,:,l,a,b,ll)=dDdA_lambda(:,:,l,a,b,ll)&
!              -dir_factor * tmpmat3
!          enddo
!        enddo
!      endif
!
!    enddo
!  enddo
!enddo
!
!end subroutine calc_fermion_force_from_omega_adm_org

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine calc_fermion_force_from_omega_test&
!    (dDdA_lambda,dDdA_chi,UMAT,lambda_mat,chi_mat)
!implicit none
!
!complex(kind(0d0)), intent(inout) :: dDdA_lambda&
!  (1:NMAT,1:NMAT,1:num_links, 1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)), intent(inout) :: dDdA_chi&
!  (1:NMAT,1:NMAT,1:num_faces, 1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)), intent(in) :: lambda_mat(1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)), intent(in) :: chi_mat(1:NMAT,1:NMAT,1:num_faces)
!
!complex(kind(0d0)) :: UXmat(1:NMAT,1:NMAT),YUmat(1:NMAT,1:NMAT),Uf(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: Ufm(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT),tmpmat2(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: UnitMat(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: Cosinv(1:NMAT,1:NMAT),Omega_mat(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: dUfdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
!complex(kind(0d0)) :: dUfmdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
!complex(kind(0d0)) :: dUfminvdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
!complex(kind(0d0)) :: dOmegadA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
!complex(kind(0d0)) :: dCosinvdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
!complex(kind(0d0)) :: dSindA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
!complex(kind(0d0)) :: dMAT(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
!complex(kind(0d0)) :: dir_factor
!integer :: i,j,k,l,ll,f,a,b,l_place,ll_place,X_last
!integer :: power = 3
!
!UnitMat=(0d0,0d0)
!do i=1,NMAT
!  UnitMat(i,i)=(1d0,0d0)
!enddo
!
!
!do f=1,num_faces
!  do l_place=1,links_in_f(f)%num_
!    l=links_in_f(f)%link_labels_(l_place)
!    dir_factor&
!      =dcmplx(dble(links_in_f(f)%link_dirs_(l_place)))*(0d0,2d0)/dcmplx(dble(m_omega))
!    call Make_face_variable(Uf(:,:),f,UMAT) 
!    call matrix_power(Ufm,Uf(:,:),m_omega)
!    do i=1,NMAT
!      do j=1,NMAT
!        Cosinv(i,j) = Ufm(i,j) + dconjg(Ufm(j,i))
!        tmpmat1(i,j) = Ufm(i,j) - dconjg(Ufm(j,i))
!        tmpmat2(i,j) = Uf(i,j) - dconjg(Uf(j,i))
!      enddo
!    enddo
!    !! Cos^{-1}
!    !call HermitianMatrix_inverse(Cosinv)
!    call Matrix_inverse(Cosinv)
!    !! Omega
!    call matrix_product(Omega_mat,Cosinv,tmpmat1,&
!      'N','N', (0d0,-2d0)/dcmplx(dble(m_omega)))
!    !call ZHEMM('L','U',NMAT,NMAT,(0d0,-2d0)/dcmplx(dble(m_omega)), &
!    !  Cosinv(:,:), NMAT, &
!    !  tmpmat1, NMAT, &
!    !  (0d0,0d0),Omega_mat(:,:) , NMAT)
!
!    do ll_place=1,links_in_f(f)%num_
!      ll=links_in_f(f)%link_labels_(ll_place)
!
!      call calc_dUfdA(dUfdA,f,ll_place,Umat)
!      call calc_dUfkdA(dUfmdA,m_omega,dUfdA,Uf)
!      do j=1,NMAT
!        do i=1,NMAT
!          do b=1,NMAT
!            do a=1,NMAT
!              dUfminvdA(i,j,a,b)=dconjg(dUfmdA(j,i,b,a))
!              dSindA(i,j,a,b)=dUfdA(i,j,a,b)-dconjg(dUfdA(j,i,b,a))
!            enddo
!          enddo
!        enddo
!      enddo
!
!      call calc_dOmegadA_dCosinvda(dOmegadA,dCosinvdA,dUfdA,Uf,Omega_mat,Cosinv)
!    !if( NMAT > 2 ) then 
!    !  call make_matrix_traceless(Omega_mat)
!    !endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!dMAT=dSindA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!
!      do b=1,NMAT
!        do a=1,NMAT
!          !!!!!!!!!!!!!!!!!!!!!!
!          call matrix_3_product(tmpmat1,DMAT(:,:,a,b),lambda_mat(:,:,l),UnitMat,'N','N','N')
!          call matrix_3_product(tmpmat2,UnitMat,lambda_mat(:,:,l),DMAT(:,:,b,a),'C','N','C')
!          dDdA_chi(:,:,f,a,b,ll)=dDdA_chi(:,:,f,a,b,ll)&
!            -dir_factor * dcmplx(alpha_f(f)*beta_f(f)) &
!            * (tmpmat1+tmpmat2) 
!          !!!!!!!!!!!!!!!!!!!!!!
!          call matrix_3_product(tmpmat1,UnitMat,chi_mat(:,:,f),DMAT(:,:,a,b),'N','N','N')
!          call matrix_3_product(tmpmat2,DMAT(:,:,b,a),chi_mat(:,:,f),UnitMat,'C','N','C')
!          dDdA_lambda(:,:,l,a,b,ll)=dDdA_lambda(:,:,l,a,b,ll)&
!            +dir_factor * dcmplx(alpha_f(f)*beta_f(f)) &
!            * (tmpmat1+tmpmat2)
!        enddo
!      enddo
!    enddo
!
!  enddo
!enddo
!
!end subroutine calc_fermion_force_from_omega_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine calc_fermion_force_from_omega_org&
!    (dDdA_lambda,dDdA_chi,UMAT,lambda_mat,chi_mat)
!use Dirac_operator, only : calc_sinu_and_CosUinv, calc_Amat, Calc_Bmat
!implicit none 
!complex(kind(0d0)), intent(inout) :: dDdA_lambda&
!  (1:NMAT,1:NMAT,1:num_links, 1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)), intent(inout) :: dDdA_chi&
!  (1:NMAT,1:NMAT,1:num_faces, 1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)), intent(in) :: lambda_mat(1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)), intent(in) :: chi_mat(1:NMAT,1:NMAT,1:num_faces)
!
!complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT,1:num_faces)
!complex(kind(0d0)) :: Ufm(1:NMAT,1:NMAT,1:num_faces)
!complex(kind(0d0)) :: sinU(1:NMAT,1:NMAT,1:num_faces)
!complex(kind(0d0)) :: cosUinv(1:NMAT,1:NMAT,1:num_faces)
!complex(kind(0d0)) :: Amat(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: Bmat(1:NMAT,1:NMAT)
!
!complex(kind(0d0)) :: dCosUinvdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
!complex(kind(0d0)) :: dSinUdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
!complex(kind(0d0)) :: dAmatdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
!complex(kind(0d0)) :: dBmatdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
!
!complex(kind(0d0)) :: prodmat1(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: prodmat2(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: globalmat1(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: globalmat2(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: globalmat3(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: globalmat4(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: tmpmat3(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: tmpmat4(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: line(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
!integer :: l_label,ll_label,k,l,ll
!integer :: i,j,ii,jj,f
!complex(kind(0d0)) :: trace
!
!! label of face
!do f=1,num_faces
!  call Make_face_variable(Uf(:,:,f),f,UMAT) 
!  if (m_omega .ne. 0) then 
!    call matrix_power(Ufm(:,:,f),Uf(:,:,f),m_omega)
!    call calc_sinU_and_cosUinv(sinU(:,:,f),cosUinv(:,:,f),Ufm(:,:,f))
!  endif
!enddo
!
!! label of face
!do f=1,num_faces
!! label of link fermion \lambda_l
!do l_label=1,links_in_f(f)%num_
!  l=links_in_f(f)%link_labels_(l_label)
!
!  ! label of link to differentiate
!  do ll_label=1,links_in_f(f)%num_
!    ll=links_in_f(f)%link_labels_(ll_label)
!!subroutine calc_dCosUinvdA_dSinUdA(dCosUinvdA,dSinUdA,&
!!    cosUinv,sinU,Uf,UMAT,f,ll_label)
!
!    line=(0d0,0d0)
!    if( m_omega == 0 ) then
!      call calc_dABmatdA(dAmatdA,dBmatdA,&
!        Uf(:,:,f),UMAT,ll_label,f,l_label,1)
!      call calc_Amat(Amat,f,l_label,1,Uf(:,:,f),UMAT)
!      call calc_Bmat(Bmat,f,l_label,1,Uf(:,:,f),UMAT)
!
!      
!      do jj=1,NMAT
!      do ii=1,NMAT
!        call matrix_3_product(tmpmat1, &
!          dAmatdA(:,:,ii,jj),lambda_mat(:,:,l),BMAT)
!          !dAmatdA(:,:,ii,jj),lambda_mat(:,:,l),BMAT)
!        call matrix_3_product(tmpmat2, &
!          Amat,lambda_mat(:,:,l),dBmatdA(:,:,ii,jj))
!        call matrix_3_product(tmpmat3, &
!          dBmatdA(:,:,jj,ii),lambda_mat(:,:,l),Amat,'C','N','C')
!        call matrix_3_product(tmpmat4, &
!          Bmat,lambda_mat(:,:,l),dAmatdA(:,:,jj,ii),'C','N','C')
!        line(:,:,ii,jj)=tmpmat1+tmpmat2+tmpmat3+tmpmat4
!      enddo
!      enddo
!        !call matrix_product(tmpmat1,Amat,lambda_mat(:,:,l))
!        !call zgemm('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
!        !  tmpmat1, NMAT, &
!        !  dBmatdA(:,:,ii,jj), NMAT, &
!        !  (1d0,0d0), line(:,:,ii,jj), NMAT)
!
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      ! make traceless for (ii,jj)
!      tmpmat1=(0d0,0d0)
!      do ii=1,NMAT
!        tmpmat1=tmpmat1+line(:,:,ii,ii)
!      enddo
!      do ii=1,NMAT
!        line(:,:,ii,ii)=line(:,:,ii,ii)-tmpmat1/dcmplx(dble(NMAT))
!      enddo
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      do jj=1,NMAT
!        do ii=1,NMAT
!          dDdA_chi(:,:,f,ii,jj,ll)=dDdA_chi(:,:,f,ii,jj,ll) &
!            + dcmplx(dble(links_in_f(f)%link_dirs_(l_label) ))*(0d0,1d0) &
!             * dcmplx( alpha_f(f)*beta_f(f) ) &
!             * line(:,:,ii,jj) 
!        enddo
!      enddo
!    else
!      call calc_dCosUinvdA_dSinUdA(dCosUinvdA,dSinUdA,&
!        cosUinv(:,:,f),Uf(:,:,f),UMAT,f,ll_label)
!      do k=1,m_omega
!        call calc_dABmatdA(dAmatdA,dBmatdA,&
!          Uf(:,:,f),UMAT,ll_label,f,l_label,k)
!        call calc_Amat(Amat,f,l_label,k,Uf(:,:,f),UMAT)
!        call calc_Bmat(Bmat,f,l_label,k,Uf(:,:,f),UMAT)
!  
!        ! globalmat1 = A lambda B
!        call matrix_3_product(globalmat1,Amat,lambda_mat(:,:,l),Bmat)
!  
!        ! globalmat2 = B^dag lambda A^dag
!        call matrix_3_product(globalmat2,Bmat,lambda_mat(:,:,l),Amat,'C','N','C')
!  
!        do jj=1,NMAT
!        do ii=1,NMAT
!          ! prodat1 = ¥delta_A lambda B + A lambda ¥delta_B
!          call matrix_3_product(prodmat1, &
!            dAmatdA(:,:,ii,jj),lambda_mat(:,:,l),BMAT)
!          call matrix_product(tmpmat1,Amat,lambda_mat(:,:,l))
!          call zgemm('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
!            tmpmat1, NMAT, &
!            dBmatdA(:,:,ii,jj), NMAT, &
!            (1d0,0d0), prodmat1, NMAT)
!  
!          ! prodmat2 = ¥delta_B^dag lambda A^dag + B^dag lambda ¥delta_A^dag
!          call matrix_3_product(prodmat2,&
!            dBmatdA(:,:,jj,ii),lambda_mat(:,:,l),Amat,'C','N','C') 
!          call matrix_product(tmpmat1,Bmat,lambda_mat(:,:,l),'C','N')
!          call zgemm('N','C',NMAT,NMAT,NMAT,(1d0,0d0), &
!            tmpmat1, NMAT, &
!            dAmatdA(:,:,jj,ii), NMAT, & 
!            (1d0,0d0), prodmat2, NMAT)
!  
!          ! line1
!          call matrix_anticommutator(tmpmat1, prodmat1+prodmat2, CosUinv(:,:,f))
!          line(:,:,ii,jj)=line(:,:,ii,jj)+tmpmat1
!          ! line2
!          call matrix_anticommutator(tmpmat1,globalmat1+globalmat2, dCosUinvdA(:,:,ii,jj))
!          line(:,:,ii,jj) = line(:,:,ii,jj) + tmpmat1
!          ! line3
!          call matrix_product(tmpmat1,CosUinv(:,:,f),globalmat1-globalmat2)
!          call matrix_product(tmpmat2,tmpmat1,CosUinv(:,:,f)) 
!          call matrix_anticommutator(tmpmat3,dSinUdA(:,:,ii,jj),tmpmat2)
!          line(:,:,ii,jj) = line(:,:,ii,jj) - tmpmat3
!          ! line4
!           ! 4-1 collect in tmpmat2
!          call matrix_product(tmpmat3,dCosUinvdA(:,:,ii,jj),globalmat1-globalmat2)
!          call matrix_product(tmpmat2,tmpmat3,CosUinv(:,:,f))
!           ! 4-2 tmpmat1=CosUinv.(globalmat1-globalmat2) here
!          call zgemm('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
!            tmpmat1, NMAT, &
!            dCosUinvdA(:,:,ii,jj), NMAT, &
!            (1d0,0d0), tmpmat2, NMAT)
!           ! 4-3
!          call matrix_product(tmpmat3,CosUinv(:,:,f),prodmat1-prodmat2)
!          call zgemm('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
!            tmpmat3, NMAT, &
!            CosUinv(:,:,f), NMAT, &
!            (1d0,0d0), tmpmat2, NMAT)
!           ! take anti-commutator
!          call matrix_anticommutator(tmpmat1,tmpmat2,SinU(:,:,f))
!          line(:,:,ii,jj) = line(:,:,ii,jj) - tmpmat1
!  
!        enddo ! ii
!        enddo ! jj
!      enddo ! k
!      
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      ! make traceless for (ii,jj)
!      tmpmat1=(0d0,0d0)
!      do ii=1,NMAT
!        tmpmat1=tmpmat1+line(:,:,ii,ii)
!      enddo
!      do ii=1,NMAT
!        line(:,:,ii,ii)=line(:,:,ii,ii)-tmpmat1/dcmplx(dble(NMAT))
!      enddo
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      do jj=1,NMAT
!        do ii=1,NMAT
!          dDdA_chi(:,:,f,ii,jj,ll)=dDdA_chi(:,:,f,ii,jj,ll) &
!            + dcmplx(dble(links_in_f(f)%link_dirs_(l_label) ))*(0d0,1d0) &
!             * dcmplx( alpha_f(f)*beta_f(f)/dble(m_omega) ) &
!             * line(:,:,ii,jj)
!          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          ! make traceless for (i,j)
!          trace=(0d0,0d0)
!          do j=1,NMAT
!            trace=trace+dDdA_chi(j,j,f,ii,jj,ll)
!          enddo
!          do j=1,NMAT
!            dDdA_chi(j,j,f,ii,jj,ll)=dDdA_chi(j,j,f,ii,jj,ll)-trace/dcmplx(dble(NMAT))
!          enddo
!         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        enddo
!      enddo
!    endif
!  enddo ! ll
!enddo ! l
!enddo ! f
!
!do l=1,num_links
!do i=1,face_in_l(l)%num_       
!  f=face_in_l(l)%label_(i)
!  ! l_label: position of the link l in the face f
!  do l_label=1,links_in_f(f)%num_
!    if ( l == links_in_f(f)%link_labels_(l_label) ) exit
!  enddo
!
!  do ll_label=1,links_in_f(f)%num_
!    ll=links_in_f(f)%link_labels_(ll_label)
!
!    if (m_omega == 0) then 
!      call calc_dABmatdA(dAmatdA,dBmatdA,&
!        Uf(:,:,f),UMAT,ll_label,f,l_label,1)
!      call calc_Amat(Amat,f,l_label,1,Uf(:,:,f),UMAT)
!      call calc_Bmat(Bmat,f,l_label,1,Uf(:,:,f),UMAT)
!
!      do jj=1,NMAT
!        do ii=1,NMAT
!          call matrix_3_product(tmpmat1,dBmatdA(:,:,ii,jj),chi_mat(:,:,f),Amat)
!          call matrix_3_product(tmpmat2,Bmat,chi_mat(:,:,f),dAmatdA(:,:,ii,jj))
!          call matrix_3_product(tmpmat3,dAmatdA(:,:,jj,ii),chi_mat(:,:,f),Bmat,'C','N','C')
!          call matrix_3_product(tmpmat4,Amat,chi_mat(:,:,f),dBmatdA(:,:,jj,ii), 'C','N','C')
!          line(:,:,ii,jj)=tmpmat1+tmpmat2+tmpmat3+tmpmat4
!        enddo
!      enddo
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      ! make traceless for (ii,jj)
!      tmpmat1=(0d0,0d0)
!      do ii=1,NMAT
!        tmpmat1=tmpmat1+line(:,:,ii,ii)
!      enddo
!      do ii=1,NMAT
!        line(:,:,ii,ii)=line(:,:,ii,ii)-tmpmat1/dcmplx(dble(NMAT))
!      enddo
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  
!      do jj=1,NMAT
!        do ii=1,NMAT
!           dDdA_lambda(:,:,l,ii,jj,ll)=dDdA_lambda(:,:,l,ii,jj,ll) &
!              - dcmplx(dble(links_in_f(f)%link_dirs_(l_label) ))*(0d0,1d0) &
!               * dcmplx( alpha_f(f)*beta_f(f) ) &
!               * line(:,:,ii,jj)
!        enddo
!      enddo
!    else
!      call calc_dCosUinvdA_dSinUdA(dCosUinvdA,dSinUdA,&
!        cosUinv(:,:,f),Uf(:,:,f),UMAT,f,ll_label)
!  
!      globalmat3=(0d0,0d0)
!      globalmat4=(0d0,0d0)
!      ! globalmat3 = { CosUinv, chi_f }
!      call matrix_anticommutator(globalmat3,CosUinv(:,:,f),chi_mat(:,:,f))
!      ! globalmat4 = CosUinv.{ SinU, chi_f }.CosUinv
!      call matrix_anticommutator(tmpmat1,SinU(:,:,f),chi_mat(:,:,f))
!      call matrix_3_product(globalmat4,CosUinv(:,:,f),tmpmat1,CosUinv(:,:,f))
!  
!      line=(0d0,0d0)
!      do k=1,m_omega
!        call calc_dABmatdA(dAmatdA,dBmatdA,&
!          Uf(:,:,f),UMAT,ll_label,f,l_label,k)
!        call calc_Amat(Amat,f,l_label,k,Uf(:,:,f),UMAT)
!        call calc_Bmat(Bmat,f,l_label,k,Uf(:,:,f),UMAT)
!  
!        do jj=1,NMAT
!        do ii=1,NMAT
!          ! prodmat1
!          call matrix_anticommutator(prodmat1,dCosUinvdA(:,:,ii,jj),chi_mat(:,:,f))
!          ! prodmat2
!          ! 1st term
!          call matrix_anticommutator(tmpmat1,dSinUdA(:,:,ii,jj),chi_mat(:,:,f))
!          call matrix_3_product(prodmat2,CosUinv(:,:,f),tmpmat1,CosUinv(:,:,f))
!          ! 2nd term
!          call matrix_anticommutator(tmpmat1,SinU(:,:,f),chi_mat(:,:,f))
!          call matrix_product(tmpmat2,dCosUinvdA(:,:,ii,jj),tmpmat1)
!          call zgemm('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
!            tmpmat2, NMAT, &
!            CosUinv(:,:,f),NMAT, &
!            (1d0,0d0), prodmat2, NMAT)
!          ! 3rd term tmpmat1={ SinU, chi_mat } now
!          call matrix_product(tmpmat2,CosUinv(:,:,f),tmpmat1)
!          call zgemm('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
!            tmpmat2, NMAT, &
!            dCosUinvdA(:,:,ii,jj),NMAT, &
!            (1d0,0d0), prodmat2, NMAT)
!  
!          ! line1
!          call matrix_3_product(tmpmat1,&
!            dBmatdA(:,:,ii,jj),globalmat3-globalmat4,Amat)
!          line(:,:,ii,jj)=line(:,:,ii,jj)+tmpmat1
!          ! line2
!          call matrix_product(tmpmat1,Bmat,globalmat3-globalmat4)
!          call zgemm('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
!            tmpmat1, NMAT, &
!            dAmatdA(:,:,ii,jj),NMAT, &
!            (1d0,0d0), line(:,:,ii,jj), NMAT)
!          ! line3
!          call matrix_product(tmpmat1,dAmatdA(:,:,jj,ii),globalmat3+globalmat4,'C','N') 
!          call zgemm('N','C',NMAT,NMAT,NMAT,(1d0,0d0), &
!            tmpmat1, NMAT, &
!            Bmat,NMAT, &
!            (1d0,0d0), line(:,:,ii,jj), NMAT)
!          ! line4
!          call matrix_product(tmpmat1,Amat,globalmat3+globalmat4,'C','N')
!          call zgemm('N','C',NMAT,NMAT,NMAT,(1d0,0d0), &
!            tmpmat1, NMAT, &
!            dBmatdA(:,:,jj,ii),NMAT, & 
!            (1d0,0d0), line(:,:,ii,jj), NMAT)
!          ! line5
!          call matrix_product(tmpmat1,Bmat,prodmat1-prodmat2)
!          call zgemm('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
!            tmpmat1, NMAT, &
!            Amat,NMAT, &
!            (1d0,0d0), line(:,:,ii,jj), NMAT)
!          ! line6
!          call matrix_product(tmpmat1,Amat,prodmat1+prodmat2,'C','N')
!          call zgemm('N','C',NMAT,NMAT,NMAT,(1d0,0d0), &
!            tmpmat1, NMAT, &
!            Bmat,NMAT, &
!            (1d0,0d0), line(:,:,ii,jj), NMAT)
!        enddo ! ii
!        enddo ! jj
!      enddo ! k
!  
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      ! make traceless for (ii,jj)
!      tmpmat1=(0d0,0d0)
!      do ii=1,NMAT
!        tmpmat1=tmpmat1+line(:,:,ii,ii)
!      enddo
!      do ii=1,NMAT
!        line(:,:,ii,ii)=line(:,:,ii,ii)-tmpmat1/dcmplx(dble(NMAT))
!      enddo
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  
!      do jj=1,NMAT
!        do ii=1,NMAT
!           dDdA_lambda(:,:,l,ii,jj,ll)=dDdA_lambda(:,:,l,ii,jj,ll) &
!              - dcmplx(dble(links_in_f(f)%link_dirs_(l_label) ))*(0d0,1d0) &
!               * dcmplx( alpha_f(f)*beta_f(f)/dble(m_omega) ) &
!               * line(:,:,ii,jj)
!        enddo
!      enddo
!    endif 
!    do jj=1,NMAT
!      do ii=1,NMAT
!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        ! make traceless for (i,j)
!        trace=(0d0,0d0)
!        do j=1,NMAT
!          trace=trace+dDdA_lambda(j,j,l,ii,jj,ll)
!        enddo
!        do j=1,NMAT
!          dDdA_lambda(j,j,l,ii,jj,ll)=dDdA_lambda(j,j,l,ii,jj,ll)-trace/dcmplx(dble(NMAT))
!        enddo
!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      enddo
!    enddo
!  enddo ! ll
!enddo ! f
!enddo ! l
!
!end subroutine calc_fermion_force_from_omega_org






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MAT = UMAT or UMAT^\dagger in f
!subroutine UMAT_in_f( MAT, f, l_label, UMAT )
!implicit none
!
!complex(kind(0d0)), intent(out) :: MAT
!integer, intent(in) :: f,l_label
!complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
!
!integer :: i,j
!
!if ( links_in_f(f)%link_dirs_(l_label) == 1 ) then
!  MAT=UMAT(:,:,links_in_f(f)%link_labels_(l_label))
!else
!  do i=1,NMAT
!    do j=1,NMAT
!      MAT(i,j)=conjg(UMAT(j,i,links_in_f(f)%link_labels_(l_label)))
!    enddo
!  enddo
!endif
!
!end subroutine UMAT_in_f














!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to return product of d/dA D . vec
!subroutine prod_diff_Diracdagger_A(dDdagdA_chi,vec)
!implicit none
!
!complex(kind(0d0)), intent(in) :: vec(1:sizeD)
!complex(kind(0d0)), intent(inout) :: dDdagdA_chi(1:sizeD, 1:dimG,1:num_links)
!complex(kind(0d0)) :: conj_vec(1:sizeD)
!integer i,a,l
!
!do i=1,sizeD
!  conj_vec(i) = -conjg( vec(i) ) 
!enddo
!
!call prod_diff_Dirac_A(dDdagdA_chi,conj_vec)
!do l=1,num_links
!  do a=1,dimG
!    do i=1,sizeD
!      dDdagdA_chi(i,a,l) = dconjg( dDdagdA_chi(i,a,l) )
!    enddo
!  enddo
!enddo
!
!
!complex(kind(0d0)) :: eta_ele(1:dimG,1:num_sites)
!complex(kind(0d0)) :: eta_mat(1:NMAT,1:NMAT,1:num_sites)
!complex(kind(0d0)) :: lambda_ele(1:dimG,1:num_links)
!complex(kind(0d0)) :: lambda_mat(1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)) :: chi_ele(1:dimG,1:num_faces)
!complex(kind(0d0)) :: chi_mat(1:NMAT,1:NMAT,1:num_faces)
!
!complex(kind(0d0)) :: Uinv(1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)) :: Phimat(1:NMAT,1:NMAT,1:num_sites)
!complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT,1:num_faces)
!complex(kind(0d0)) :: Ufm(1:NMAT,1:NMAT,1:num_faces)
!complex(kind(0d0)) :: diffdiff_Omega(1:NMAT,1:NMAT,1:dimG,1:dimG)
!
!
!complex(kind(0d0)) :: trace,tmp
!complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT),tmpmat2(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: comm(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: FF_MAT(1:NMAT,1:NMAT) ! part of fermion force
!
!integer :: s,l,f
!integer :: i,j,k,nl,r
!integer :: a,b,c,d,e
!
!!! preparation
!dDdagdA_chi=(0d0,0d0)
!do s=1,num_sites
!  do a=1,dimG
!    eta_ele(a,s)=vec(site_index(a,s))
!  enddo
!  call make_traceless_matrix_from_modes(eta_mat(:,:,s),NMAT,eta_ele(:,s))
!enddo
!do l=1,num_links
!  do a=1,dimG
!    lambda_ele(a,l)=vec(link_index(a,s))
!  enddo
!  call make_traceless_matrix_from_modes(lambda_mat(:,:,l),NMAT,lambda_ele(:,l))
!enddo
!do f=1,num_faces
!  do a=1,dimG
!    chi_ele(a,f)=vec(face_index(a,f))
!  enddo
!  call make_traceless_matrix_from_modes(chi_mat(:,:,f),NMAT,chi_ele(:,f))
!enddo
!do s=1,num_sites
!  call make_traceless_matrix_from_modes(Phimat(:,:,s),NMAT,Phi)
!enddo
!
!! (1) Dirac from site
!!   no contribution
!
!!! (2) Dirac from link 1
!do s=1,num_sites
!  do k=1,linkorg_to_s(s)%num_
!    l=linkorg_to_s(s)%labels_(k)
!    do b=1,dimG
!      call MtimesT(tmpmat1,lambda_mat(:,:,l),b,NMAT)
!      call TtimesM(tmpmat2,lambda_mat(:,:,l),b,NMAT)
!      comm=tmpmat1-tmpmat2
!      call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
!        comm, NMAT, &
!        UMAT(:,:,l), NMAT, &
!        (0d0,0d0), tmpmat1, NMAT)
!      call ZGEMM('C','N',NMAT,NMAT,NMAT,dcmplx(alpha_l(l)), &
!        UMAT(:,:,l), NMAT, &
!        tmpmat1, NMAT, &
!        (0d0,0d0), FF_MAT, NMAT)
!      do a=1,dimG
!        call trace_MTa(trace,FF_MAT,a,NMAT)
!        dDdagdA_chi(site_index(a,s),b,l)= dDdagdA_chi(site_index(a,s),b,l) &
!          + trace
!      enddo
!    enddo
!  enddo
!enddo    
!
!
!do l=1,num_links
!  s=link_tip(l)
!  call ZGEMM('N','C',NMAT,NMAT,NMAT,(1d0,0d0), &
!    eta_mat(:,:,s), NMAT, &
!    UMAT(:,:,l), NMAT, &
!    (0d0,0d0), tmpmat1, NMAT)
!  call ZGEMM('N','N',NMAT,NMAT,NMAT,(0d0,-1d0)*dcmplx(alpha_l(l)), &
!    UMAT(:,:,l), NMAT, &
!    tmpmat1, NMAT, &
!    (0d0,0d0), FF_MAT, NMAT)
!  do r=1,NZF
!    a=NZF_index(1,r)
!    b=NZF_index(2,r)
!    c=NZF_index(3,r)
!    call trace_MTa(trace,FF_MAT,c,NMAT)
!    dDdagdA_chi(link_index(a,l),b,l)= dDdagdA_chi(link_index(a,l),b,l) &
!      + trace * NZF_value(r)
!  enddo
!enddo
!
!end subroutine prod_diff_Diracdagger_A

end module differential_Dirac

