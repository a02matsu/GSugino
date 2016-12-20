!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! module to compute forces
module forces
use global_parameters
use global_subroutines
implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! dS/Dphi 
!complex(kind(0d0)) function dSdPhi(a,s)
!subroutine Make_bosonic_force_Phi(dSdPhi_boson)
subroutine Make_force(dSdPhi,dSdA,UMAT,PhiMat,PF,info)
use SUN_generators, only : trace_mta
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(in) :: PF(1:sizeD)
complex(kind(0d0)), intent(out) :: dSdPhi(1:dimG,1:num_sites)
double precision, intent(out) :: dSdA(1:dimG,1:num_links)
integer, intent(inout) :: info

complex(kind(0d0)) :: dSdPhi_boson_mass_org(1:dimG,1:num_sites)
complex(kind(0d0)) :: dSdPhi_boson_site_org(1:dimG,1:num_sites)
complex(kind(0d0)) :: dSdPhi_boson_link_org(1:dimG,1:num_sites)
complex(kind(0d0)) :: dSdPhi_boson_mass(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: dSdPhi_boson_site(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: dSdPhi_boson_link(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: dSdPhi_fermion(1:dimG,1:num_sites)
!!!
!double precision :: dSdA_boson_test(1:dimG,1:num_links)
double precision :: dSdA_boson_link(1:dimG,1:num_links)
double precision :: dSdA_boson_face(1:dimG,1:num_links)
double precision :: dSdA_fermion(1:dimG,1:num_links)

integer :: s,a,l
do s=1,num_sites
  do a=1,dimG
    call trace_MTa(Phi(a,s),PhiMat(:,:,s),a,NMAT)
  enddo
enddo

dSdPhi=(0d0,0d0)
dSdPhi_boson_mass=(0d0,0d0)
dSdPhi_boson_site=(0d0,0d0)
dSdPhi_boson_link=(0d0,0d0)
dSdPhi_fermion=(0d0,0d0)
dSdA=0d0
dSdA_boson_link=0d0
dSdA_boson_face=0d0
dSdA_fermion=0d0
!! force for Phi from boson
call Make_bosonic_force_Phi_mass(dSdPhi_boson_mass,PhiMat)
call Make_bosonic_force_Phi_site(dSdPhi_boson_site,PhiMat)
call Make_bosonic_force_Phi_link(dSdPhi_boson_link,UMAT,PhiMat)

!! 戻す
do s=1,num_sites
  do a=1,dimG
    call trace_MTa(dSdPhi_boson_mass_org(a,s),dSdPhi_boson_mass(:,:,s),a,NMAT)
    call trace_MTa(dSdPhi_boson_site_org(a,s),dSdPhi_boson_site(:,:,s),a,NMAT)
    call trace_MTa(dSdPhi_boson_link_org(a,s),dSdPhi_boson_link(:,:,s),a,NMAT)
  enddo
enddo

!! force for A from boson
call Make_bosonic_force_A_link(dSdA_boson_link,UMAT,PhiMat)
call Make_bosonic_force_A_face(dSdA_boson_face,UMAT)
!call Make_bosonic_force_A_test(dSdA_boson_test,UMAT)

!! force from fermion
!write(*,*) "1"
call Make_fermionic_force(dSdPhi_fermion,dSdA_fermion,UMAT,PhiMat,PF,info)
!write(*,*) "2"
!write(*,*) dSdA_fermion

dSdPhi= dSdPhi_boson_mass_org  ! mass part
dSdPhi= dSdPhi + dSdPhi_boson_site_org  ! site part
dSdPhi= dSdPhi + dSdPhi_boson_link_org   ! link part
dSdPhi= dSdPhi + dSdPhi_fermion   ! link part

dSdA = dSdA_boson_link 
dSdA = dSdA + dSdA_boson_face 
dSdA = dSdA + dSdA_fermion
!dSdA = dSdA + dSdA_boson_test &


end subroutine Make_force

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! test
!subroutine Make_bosonic_force_A_test(dSdA_boson_test,UMAT)
!use matrix_functions, only : matrix_power
!implicit none
!
!complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
!double precision, intent(out) :: dSdA_boson_test(1:dimG,1:num_links)
!
!complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)) :: Ufm(1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)) :: diff_Ufm(1:NMAT,1:NMAT,1:dimG)
!complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT,1:num_faces)
!complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: tmp
!integer :: f,l,i,j,k,a
!integer :: m_tmp
!
!dSdA_boson_test=0d0
!
!m_tmp=m_omega
!
!do f=1,num_faces
!  call Make_face_variable(Uf(:,:,f),f,UMAT) 
!  call matrix_power(NMAT,Uf(:,:,f),m_tmp,Ufm(:,:,f))
!enddo
!
!tmpmat=(0d0,0d0)
!do f=1,num_faces
!  do i=1,NMAT
!    do j=1,NMAT
!      tmpmat(i,j,f)=(0d0,-1d0)*(Ufm(i,j,f)-dconjg( Ufm(j,i,f) ))
!    enddo
!  enddo
!enddo
!
!do l=1,num_links
!  do k=1,face_in_l(l)%num_
!    f=face_in_l(l)%label_(k)
!
!    if( m_tmp==1 ) then 
!      call calc_diff_Uf(diff_Ufm, Uf,f,l,UMAT)
!    else
!      call calc_diff_Ufm(diff_Ufm, Uf,f,l,UMAT)
!    endif
!    do a=1,dimG
!      do i=1,NMAT
!        do j=1,NMAT
!          tmpmat2(i,j)=(0d0,-1d0)*(diff_Ufm(i,j,a)-dconjg( diff_Ufm(j,i,a) ))
!        enddo
!      enddo
!      tmp=(0d0,0d0)
!      do i=1,NMAT
!        do j=1,NMAT
!          tmp=tmp+tmpmat(i,j,f)*tmpmat2(j,i)
!        enddo
!      enddo
!      dSdA_boson_test(a,l)=dSdA_boson_test(a,l)+dble(tmp)
!    enddo
!  enddo
!enddo
!
!end subroutine Make_bosonic_force_A_test


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Forces from mass terms, Sites, Links, Faces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! dS/dPhi from mass action
subroutine Make_bosonic_force_Phi_mass(dSdPhi_boson_mass,PhiMat)
implicit none

complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(out) :: dSdPhi_boson_mass(1:NMAT,1:NMAT,1:num_sites)
integer :: s,i,j

dSdPhi_boson_mass=(0d0,0d0)
do s=1,num_sites
  do j=1,NMAT
    do i=1,NMAT
  dSdPhi_boson_mass(i,j,s) = dcmplx(mass_square_phi) &!*  dconjg(Phi(a,s)) &
    * dconjg(PhiMat(j,i,s)) * dcmplx(overall_factor)
    !* cmplx(one_ov_2g2N)
    enddo
  enddo
enddo
end subroutine Make_bosonic_force_Phi_mass

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! dS/Dphi from site action
!complex(kind(0d0)) function dSdPhi_boson_site(a,s)
subroutine Make_bosonic_force_Phi_site(dSdPhi_boson_site,PhiMat)
use SUN_generators, only : make_traceless_matrix_from_modes,trace_MTa
use matrix_functions, only : matrix_commutator
implicit none

complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(out) :: dSdPhi_boson_site(1:NMAT,1:NMAT,1:num_sites)
integer :: a,s
complex(kind(0d0)) :: comm1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: comm2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmp(1:dimG), tmp1,tmp2
integer :: i,j

dSdPhi_boson_site=(0d0,0d0)

do s=1,num_sites
  call matrix_commutator(comm1,phimat(:,:,s),phimat(:,:,s),'N','C')
  call matrix_commutator(comm2,phimat(:,:,s),comm1,'C','N')
  !do a=1,dimG
  do j=1,NMAT
    do i=1,NMAT
    !call trace_MTa(trace,comm2,a,NMAT)
    dSdPhi_boson_site(i,j,s)=dcmplx(0.5d0*alpha_s(s))*comm2(i,j)& 
      * cmplx( overall_factor )
      !* cmplx( one_ov_2g2N )
    enddo
  enddo
enddo


!do s=1,num_sites
!do a=1,dimG
!  tmp=(0d0,0d0)
!  do i=1,NZF
!    b=NZF_index(1,i)
!    c=NZF_index(2,i)
!    d=NZF_index(3,i)
!    tmp(b)=tmp(b)+im_unit*NZF_value(i)*Phi(c,s)*dconjg(Phi(d,s))
!  enddo
!  
!  do i=1,NZF_a(a)%num_
!    b=NZF_a(a)%b_(i)
!    c=NZF_a(a)%c_(i)
!    dSdPhi_boson_site(a,s)=dSdPhi_boson_site(a,s) &
!        +im_unit*NZF_a(a)%value_(i)*tmp(b)*dconjg(Phi(c,s))
!  enddo
!  dSdPhi_boson_site(a,s) = dSdPhi_boson_site(a,s) * (-0.5d0)*alpha_s(s)
!
!enddo
!enddo
end subroutine Make_bosonic_force_Phi_site

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! dS/Dphi from link action
!complex(kind(0d0)) function dSdPhi_boson_link(a,s)
subroutine Make_bosonic_force_Phi_link(dSdPhi_boson_link,UMAT,PhiMat)
use SUN_generators, only : trace_MTa, make_traceless_matrix_from_modes
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(out) :: dSdPhi_boson_link(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: dPhi(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT), tmpmat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Ud_dPhibar_U(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: dPhibar(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: trace
integer :: a,s
integer :: l, l_i
integer :: i,j,k

dSdPhi_boson_link=(0d0,0d0)
do l=1,num_links
! diff( \Phi_s ) for the link l
  call Make_diff_PhiMat(dPhi,l,UMAT,PhiMat)
! U_l^\dagger diff( \bar\Phi_s ) U_l for all links 
  call ZGEMM('C','C',NMAT,NMAT,NMAT,(1d0,0d0), &
      UMAT(:,:,l), NMAT, &
      dPhi, NMAT, &
      (0d0,0d0), tmpmat, NMAT)
  call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
      tmpmat(:,:), NMAT, &
      UMAT(:,:,l), NMAT, &
      (0d0,0d0), Ud_dPhibar_U(:,:,l), NMAT)

  do j=1,NMAT
    do i=1,NMAT
      dPhibar(i,j,l)=conjg(dPhi(j,i))
    enddo
  enddo
! Tr( Ta U_l^\dagger diff( \bar\Phi_s ) U_l ) for all links 
    !call Trace_MTa(Tr_Ta_Ud_dPhibar_U(a,l),tmpmat2,a,NMAT)
! Tr( Ta diff( \bar\Phi_s ) )
    !call Trace_MTa(trace,dPhi,a,NMAT)
    !Tr_Ta_dPhibar(a,l) = dconjg(trace)
enddo

dSdPhi_boson_link=(0d0,0d0)
do s=1,num_sites
  do i=1,linkorg_to_s(s)%num_
    l_i=linkorg_to_s(s)%Labels_(i)
    do k=1,NMAT
      do j=1,NMAT
        dSdPhi_boson_link(j,k,s)=dSdPhi_boson_link(j,k,s) &
          + alpha_l(l_i) *  Ud_dPhibar_U(j,k,l_i)
      enddo
    enddo
  enddo

  do i=1,linktip_from_s(s)%num_
    l_i=linktip_from_s(s)%Labels_(i)
    do k=1,NMAT
      do j=1,NMAT
        dSdPhi_boson_link(j,k,s)=dSdPhi_boson_link(j,k,s) - alpha_l(l_i) &
          * dPhibar(j,k,l_i)
      enddo
    enddo
  enddo
enddo

dSdPhi_boson_link = dSdPhi_boson_link  &
  * cmplx( overall_factor )

end subroutine Make_bosonic_force_Phi_link

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! dS/DA from link action
subroutine Make_bosonic_force_A_link(dSdA_boson_link,UMAT,PhiMat)
use SUN_generators, only : trace_MTa, make_traceless_matrix_from_modes
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
double precision, intent(out) :: dSdA_boson_link(1:dimG,1:num_links)
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: dPhi(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Phi_tip(1:NMAT,1:NMAT)
complex(kind(0d0)) :: UPhiUinv(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Comm(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmp
integer :: a,l,i,j

!integer :: b,c,d,i
!double precision :: f_cad
!complex(kind(0d0)) :: trace

dSdA_boson_link=0d0
do l=1,num_links
  call Make_diff_PhiMAT(dPhi, l,UMAT,PhiMat)

  Phi_tip=PhiMat(:,:,link_tip(l))
  !call Make_traceless_matrix_from_modes(Phi_tip,NMAT,Phi(:,link_tip(l)))

  ! U_l.Phi_tip
  call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
    UMAT(:,:,l), NMAT, &
    Phi_tip, NMAT, &
    (0d0,0d0), tmpmat, NMAT)
  ! U_l.Phi_tip.U_l^\dagger
  call ZGEMM('N','C',NMAT,NMAT,NMAT,(1d0,0d0), &
    tmpmat, NMAT, &
    UMAT(:,:,l), NMAT, &
    (0d0,0d0), UPhiUinv, NMAT)

  ! i[ UPhiUinv, dSPhi^\dagger ]
  call ZGEMM('N','C',NMAT,NMAT,NMAT,(0d0,1d0), &
    UPhiUinv, NMAT, &
    dPhi, NMAT, &
    (0d0,0d0), Comm, NMAT)
  call ZGEMM('C','N',NMAT,NMAT,NMAT,(0d0,-1d0), &
    dPhi, NMAT, &
    UPhiUinv, NMAT, &
    (1d0,0d0), Comm, NMAT)

  do a=1,dimG
    call Trace_MTa( tmp, Comm, a, NMAT )
    dSdA_boson_link(a,l) = alpha_l(l) * dble( tmp + dconjg(tmp) )
  enddo
enddo

dSdA_boson_link=dSdA_boson_link * overall_factor !*one_ov_2g2N

end subroutine Make_bosonic_force_A_link


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! dS/DA from face action
subroutine Make_bosonic_force_A_face(dSdA_boson_face,UMAT)
use matrix_functions, only : matrix_power
implicit none
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
double precision, intent(out) :: dSdA_boson_face(1:dimG,1:num_links)
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT,1:num_faces), Ufm(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: Omega(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: diff_Omega(1:NMAT,1:NMAT,1:dimG)
integer :: FaceSize
integer, allocatable :: sites(:),link_labels(:),link_dirs(:)
integer, allocatable :: faces_l(:)
complex(kind(0d0)) :: trace
integer :: a,l,f
integer :: l_f
integer :: i,j,k

dSdA_boson_face=0d0
do f=1,num_faces
  call Make_face_variable(Uf(:,:,f),f,UMAT) 
  call matrix_power(Ufm(:,:,f),Uf(:,:,f),m_omega)
  call Make_moment_map(Omega(:,:,f),Ufm(:,:,f))
enddo

do l=1,num_links
  !call get_faces_in_link_sc(sc,l,faces_l)
  do k=1,face_in_l(l)%num_
    f=face_in_l(l)%label_(k)

    call calc_diff_omega(diff_Omega(:,:,:),Uf(:,:,f),Ufm(:,:,f),f,l,UMAT)
    do a=1,dimG
      trace=(0d0,0d0)
      do i=1,NMAT
        do j=1,NMAT
          trace=trace+Omega(i,j,f)*diff_Omega(j,i,a)
        enddo
      enddo
      dSdA_boson_face(a,l)=dSdA_boson_face(a,l) & 
        + 0.5d0*alpha_f(f)*beta_f(f)*beta_f(f)*dble(trace)
    enddo
  enddo
enddo

dSdA_boson_face=dSdA_boson_face * overall_factor !*one_ov_2g2N

end subroutine Make_bosonic_force_A_face

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! dS/dPhi and dS/dA from fermion part
subroutine Make_fermionic_force(dSdPhi_fermion,dSdA_fermion,UMAT,Phimat,PF,info)
use rational_algorithm
use Dirac_operator
use differential_Dirac
implicit none 

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(in) :: PF(1:sizeD)
complex(kind(0d0)), intent(out) :: dSdPhi_fermion(1:dimG,1:num_sites)
double precision, intent(out) :: dSdA_fermion(1:dimG,1:num_links)
integer, intent(inout) :: info

! chi_r = (D\dagger D + \beta_r)^{-1} F
complex(kind(0d0)) :: chi(1:sizeD,1:N_Remez4)
! Dchi_r = D.chi_r
complex(kind(0d0)) :: Dchi(1:sizeD,1:N_Remez4)
! [\frac{dD^\dagger D}{dPhi_s^a}]_{ij} \chi_j
!complex(kind(0d0)) :: dDdagD_dPhi_chi(1:sizeD, 1:dimG,1:num_sites,1:N_Remez4)
! [\frac{dD^\dagger D}{dA_l^a}]_{ij} \chi_j
!complex(kind(0d0)) :: dDdagD_dA_chi(1:sizeD, 1:dimG,1:num_links,1:N_Remez4)
! 
complex(kind(0d0)) :: dDdPhi_chi(1:sizeD,1:dimG,1:num_sites,1:N_Remez4)
complex(kind(0d0)) :: dDdbPhi_chi(1:sizeD,1:dimG,1:num_sites,1:N_Remez4)
complex(kind(0d0)) :: dDdA_chi(1:sizeD,1:dimG,1:num_links,1:N_Remez4)

complex(kind(0d0)) :: tmp
integer :: CGite
integer :: r,i,s,l,a,b

do s=1,num_sites
  do a=1,dimG
    call trace_MTa(Phi(a,s),PhiMat(:,:,s),a,NMAT)
  enddo
enddo

!! for test
!complex(kind(0d0)) :: DdagD(1:sizeD,1:sizeD)
!complex(kind(0d0)) :: tmpmat(1:sizeD,1:sizeD)
!complex(kind(0d0)) :: chi_direct(1:sizeD,1:N_Remez4)
!double precision :: distance
!integer :: j


dSdPhi_fermion=(0d0,0d0)
dSdA_fermion=0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! comutation by CG
! compute chi_r(i,r) =  (D^\dagger D + \beta_r)^{-1} F(i)
!write(*,*) "   mmBiCG start"
!call check_Dirac(UMAT,Phi)
call mmBiCG( chi, PF, Remez_beta4, sizeD, N_Remez4, epsilon, &
             CG_max, info, CGite, UMAT, PhiMat, Prod_DdagD )
!write(*,*) "   mmBiCG end", cgite
do r=1,N_Remez4
  call prod_Dirac(Dchi(:,r),chi(:,r),sizeD,UMAT,PhiMat)
  
  call prod_dDdPhi(dDdPhi_chi(:,:,:,r),chi(:,r),UMAT,Phi)
  call prod_dDdbPhi(dDdbPhi_chi(:,:,:,r),chi(:,r),UMAT,Phi)
  
  call prod_dDdA(dDdA_chi(:,:,:,r),chi(:,r),UMAT,Phi)
enddo

do s=1,num_sites
  do a=1,dimG
    do r=1,N_Remez4
      do i=1,sizeD
        dSdPhi_fermion(a,s)=dSdPhi_fermion(a,s) & 
          - dcmplx(Remez_alpha4(r)) &
            * ( dconjg( Dchi(i,r) ) * dDdPhi_chi(i,a,s,r) &
               + dconjg( dDdbPhi_chi(i,a,s,r) ) * Dchi(i,r) )
      enddo
    enddo
  enddo
enddo

do l=1,num_links
  do a=1,dimG
    do r=1,N_Remez4
      do i=1,sizeD
        dSdA_fermion(a,l)=dSdA_fermion(a,l) & 
          - Remez_alpha4(r)   &
            *dble( dconjg( Dchi(i,r) ) * dDdA_chi(i,a,l,r) &
                   + dconjg( dDdA_chi(i,a,l,r) ) * Dchi(i,r) ) 
      !write(*,*) "real?", &
      !       dconjg( Dchi(i,r) ) * dDdA_chi(i,a,l,r) &
      !       + dconjg( dDdA_chi(i,a,l,r) ) * Dchi(i,r) 
      enddo
    enddo
  enddo
enddo

! changed 2015/11/08
!dSdPhi_fermion = dSdPhi_fermion / cmplx(overall_factor) ! / cmplx(one_ov_2g2N)
!dSdA_fermion = dSdA_fermion / cmplx(overall_factor) ! / cmplx(one_ov_2g2N)
!write(*,*) "   end: fermionic force"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! direct computation of chi
!chi_direct=(0d0,0d0)
!call make_DdagD(DdagD)
!do r=1,N_Remez4
!  tmpmat=DdagD
!  do i=1,sizeD
!    tmpmat(i,i)=tmpmat(i,i)+Remez_beta4(r)
!  enddo
!  call matrix_inverse(sizeD,tmpmat)
!  do i=1,sizeD
!    do j=1,sizeD
!      chi_direct(i,r)=chi_direct(i,r)+tmpmat(i,j)*PF(j)
!    enddo
!  enddo
!enddo
!
!distance=0d0
!do r=1,N_Remez4
!  do i=1,sizeD
!    tmp=chi(i,r)-chi_direct(i,r)
!    distance=distance+dble( tmp * dconjg(tmp) )
!enddo
!enddo
!
!write(*,*) distance
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  


  

end subroutine Make_fermionic_force

end module forces
