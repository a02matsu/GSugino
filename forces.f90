!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! module to compute forces
module forces
use global_parameters
use global_subroutines
implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! dS/Dphi 
subroutine Make_force(dSdPhi,dSdA,UMAT,PhiMat,PF,info)
use SUN_generators, only : trace_mta, make_traceless_matrix_from_modes
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(in) :: PF(1:sizeD)
complex(kind(0d0)), intent(out) :: dSdPhi(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(out) :: dSdA(1:NMAT,1:NMAT,1:num_links)
integer, intent(inout) :: info

complex(kind(0d0)) :: dSdPhi_boson_mass(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: dSdPhi_boson_site(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: dSdPhi_boson_link(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: dSdPhi_fermion(1:NMAT,1:NMAT,1:num_sites)
!!!
complex(kind(0d0)) :: dSdA_boson_link(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: dSdA_boson_face(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: dSdA_fermion(1:NMAT,1:NMAT,1:num_links)

integer :: s,a,ii,jj,l
complex(kind(0d0)) :: tmp

dSdPhi=(0d0,0d0)
dSdPhi_boson_mass=(0d0,0d0)
dSdPhi_boson_site=(0d0,0d0)
dSdPhi_boson_link=(0d0,0d0)
dSdPhi_fermion=(0d0,0d0)
dSdA=(0d0,0d0)
dSdA_boson_face=(0d0,0d0)
dSdA_boson_link=(0d0,0d0)
dSdA_fermion=(0d0,0d0)
!! force for Phi from boson
call Make_bosonic_force_Phi_mass(dSdPhi_boson_mass,PhiMat)
call Make_bosonic_force_Phi_site(dSdPhi_boson_site,PhiMat)
call Make_bosonic_force_Phi_link(dSdPhi_boson_link,UMAT,PhiMat)
!! force from fermion
call Make_fermionic_force(dSdPhi_fermion,dSdA_fermion,UMAT,PhiMat,PF,info)

!! force for A from boson
call Make_bosonic_force_A_link(dSdA_boson_link,UMAT,PhiMat)
call Make_bosonic_force_A_face(dSdA_boson_face,UMAT)
!dSdA_boson_face=(0d0,0d0)

dSdPhi= dSdPhi_boson_mass &   ! mass part
+ dSdPhi_boson_site  & ! site part
+ dSdPhi_boson_link  &  ! link part
+ dSdPhi_fermion   ! link part

dSdA = dSdA_boson_link  &
+ dSdA_boson_face &
+ dSdA_fermion

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
    dSdPhi_boson_site(i,j,s)=dcmplx(0.5d0*alpha_s(s))*comm2(i,j)& 
      * cmplx( overall_factor )
      !* cmplx( one_ov_2g2N )
    enddo
  enddo
enddo

end subroutine Make_bosonic_force_Phi_site

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! dS/Dphi from link action
!complex(kind(0d0)) function dSdPhi_boson_link(a,s)
subroutine Make_bosonic_force_Phi_link(dSdPhi_boson_link,UMAT,PhiMat)
use matrix_functions, only : matrix_3_product
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
  call matrix_3_product(Ud_dPhibar_U(:,:,l),UMAT(:,:,l),dPhi,UMAT(:,:,l),'C','C','N')

  do j=1,NMAT
    do i=1,NMAT
      dPhibar(i,j,l)=conjg(dPhi(j,i))
    enddo
  enddo
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
use matrix_functions, only : matrix_3_product, matrix_commutator
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
!double precision, intent(out) :: dSdA_boson_link_org(1:dimG,1:num_links)
complex(kind(0d0)), intent(out) :: dSdA_boson_link(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: dPhi(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Phi_tip(1:NMAT,1:NMAT)
complex(kind(0d0)) :: UPhiUinv(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Comm(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmp
integer :: a,l,i,j

dSdA_boson_link=(0d0,0d0)
do l=1,num_links
  call Make_diff_PhiMAT(dPhi, l,UMAT,PhiMat)

  Phi_tip=(0d0,1d0)*PhiMat(:,:,link_tip(l))

  ! U_l.Phi_tip
!  ! i*U_l.Phi_tip.U_l^\dagger
  call matrix_3_product(UPhiUinv,UMAT(:,:,l),Phi_tip,UMAT(:,:,l),'N','N','C')

  ! i[ UPhiUinv, dSPhi^\dagger ]
  call matrix_commutator(Comm,UPhiUinv,dPhi,'N','C')

  do j=1,NMAT
    do i=1,NMAT
      dSdA_boson_link(i,j,l) = alpha_l(l) * ( Comm(i,j) + conjg( Comm(j,i)) )
    enddo
  enddo
enddo

dSdA_boson_link=dSdA_boson_link * overall_factor !*one_ov_2g2N

end subroutine Make_bosonic_force_A_link


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! dS/DA from face action
subroutine Make_bosonic_force_A_face(dSdA_boson_face,UMAT)
use SUN_generators, only : trace_mta, make_traceless_matrix_from_modes
use matrix_functions, only : matrix_power,matrix_product
!use Dirac_operator, only : calc_Sinu_and_CosUinv
!use diferential_Dirac, only :  calc_dCosUinvdA_dSinUdA
implicit none
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
double precision dSdA_boson_face_org(1:dimG,1:num_links)
complex(kind(0d0)), intent(out) :: dSdA_boson_face(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT,1:num_faces), Ufm(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: Omega(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: diff_Omega_org(1:NMAT,1:NMAT,1:dimG)
!!!!!!!!!
complex(kind(0d0)) :: diff_Omega(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
complex(kind(0d0)) :: sinU(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: CosUinv(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: dCosUinvdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
complex(kind(0d0)) :: dSinUdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
!!!!!!!!!
!integer :: FaceSize
!integer, allocatable :: sites(:),link_labels(:),link_dirs(:)
!integer, allocatable :: faces_l(:)
complex(kind(0d0)) :: trace,tmp
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: tmp2(1:dimG)
!complex(kind(0d0)) :: comp(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:num_faces,1:num_links)
integer :: a,l,f!,ll
integer :: l_label
integer :: i,j,k,ii,jj

!Omega=(0d0,0d0)
sinU=(0d0,0d0)
cosUinv=(0d0,0d0)
do f=1,num_faces
  call Make_face_variable(Uf(:,:,f),f,UMAT) 
  call matrix_power(Ufm(:,:,f),Uf(:,:,f),m_omega)
  call Make_moment_map(Omega(:,:,f),Ufm(:,:,f))
  call calc_sinU_and_cosUinv(sinU(:,:,f),cosUinv(:,:,f),Ufm(:,:,f))
enddo

dSdA_boson_face=(0d0,0d0)
do l=1,num_links
  do k=1,face_in_l(l)%num_
    f=face_in_l(l)%label_(k)
    do l_label=1,links_in_f(f)%num_
      if ( l == links_in_f(f)%link_labels_(l_label) ) exit
    enddo

    dCosUinvdA=(0d0,0d0)
    dSinUdA=(0d0,0d0)
    call calc_dCosUinvdA_dSinUdA(&
      dCosUinvdA(:,:,:,:),dSinUdA(:,:,:,:),&
      cosUinv(:,:,f),Uf(:,:,f),UMAT,f,l_label)

    diff_Omega=(0d0,0d0)
    tmpmat=(0d0,0d0)
    do jj=1,NMAT
      do ii=1,NMAT
         call matrix_product(&
           diff_Omega(:,:,ii,jj),dSinUdA(:,:,ii,jj),cosUinv(:,:,f))
         call zgemm('N','N',NMAT,NMAT,NMAT,(0d0,-1d0), &
           SinU(:,:,f), NMAT, &
           dCosUinvdA(:,:,ii,jj), NMAT, &
           (0d0,-1d0), diff_Omega(:,:,ii,jj), NMAT)
         call zgemm('N','N',NMAT,NMAT,NMAT,(0d0,-1d0), &
           dCosUinvdA(:,:,ii,jj), NMAT, &
           SinU(:,:,f), NMAT, &
           (1d0,0d0), diff_Omega(:,:,ii,jj), NMAT)
         call zgemm('N','N',NMAT,NMAT,NMAT,(0d0,-1d0), &
           CosUinv(:,:,f), NMAT, &
           dSinUdA(:,:,ii,jj), NMAT, &
           (1d0,0d0), diff_Omega(:,:,ii,jj), NMAT)
         if( ii==jj ) then
           tmpmat=tmpmat+diff_Omega(:,:,ii,ii)
         endif
       enddo
     enddo
     ! make diff_Omega traceless
     do ii=1,NMAT
       diff_Omega(:,:,ii,ii)=diff_Omega(:,:,ii,ii)-tmpmat/cmplx(dble( NMAT ))
     enddo
     do jj=1,NMAT
       do ii=1,NMAT
        !!!!!!!!!!!!!!!!!!!!!
        ! make traceless
        if( NMAT > 2 ) then
          trace=(0d0,0d0)
          do i=1,NMAT
            trace=trace+diff_Omega(i,i,ii,jj)
          enddo
          do i=1,NMAT
            diff_Omega(i,i,ii,jj)=diff_Omega(i,i,ii,jj) &
              - trace / cmplx(dble( NMAT ))
          enddo
        endif
        !!!!!!!!!!!!!!!!!!!!!
      enddo
    enddo
    
    do jj=1,NMAT
      do ii=1,NMAT
        do j=1,NMAT
          do i=1,NMAT
            dSdA_boson_face(ii,jj,l)=dSdA_boson_face(ii,jj,l) &
              + (0.5d0,0d0) &
                *cmplx( overall_factor * alpha_f(f)*beta_f(f)*beta_f(f) / dble(m_omega)) &
                *diff_Omega(i,j,ii,jj) *  Omega(j,i,f) 
                    !+  Omega(i,j,f) *  diff_Omega(j,i,ii,jj) )
            !comp(i,j,ii,jj,f,l)=diff_Omega(i,j,ii,jj)
          enddo
        enddo
      enddo
    enddo
  enddo
enddo

!!!!!!!!!!! OLD !!!!!!!!!!!
!dSdA_boson_face_org=0d0
!do l=1,num_links
!  !call get_faces_in_link_sc(sc,l,faces_l)
!  do k=1,face_in_l(l)%num_
!    f=face_in_l(l)%label_(k)
!    do l_label=1,links_in_f(f)%num_
!      if ( l == links_in_f(f)%link_labels_(l_label) ) exit
!    enddo
!
!    call calc_diff_omega(diff_Omega_org(:,:,:),Uf(:,:,f),Ufm(:,:,f),f,l,UMAT)
!    !!!!!!!!!!!!!!!!!!
!    do j=1,NMAT
!      do i=1,NMAT
!        do a=1,dimG
!          tmp2(a)=diff_omega_org(i,j,a)
!        enddo
!        call make_traceless_matrix_from_modes(tmpmat,NMAT,tmp2)
!        do jj=1,NMAT
!          do ii=1,NMAT
!            write(*,*) f,l,i,j,ii,jj,tmpmat(ii,jj)-comp(i,j,ii,jj,f,l)
!          enddo
!        enddo
!      enddo
!    enddo
!    !!!!!!!!!!!!!!!!!!
!    do a=1,dimG
!    trace=(0d0,0d0)
!    do i=1,NMAT
!      do j=1,NMAT
!        trace=trace+Omega(i,j,f)*diff_Omega_org(j,i,a)
!      enddo
!    enddo
!    dSdA_boson_face_org(a,l)=dSdA_boson_face_org(a,l) & 
!      + 0.5d0*alpha_f(f)*beta_f(f)*beta_f(f)*dble(trace)
!    enddo
!  enddo
!enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!
!stop
!
!dSdA_boson_face_org=dSdA_boson_face_org * overall_factor !*one_ov_2g2N
!do l=1,num_links
!  !do a=1,dimG
!    !call trace_MTa(tmp,dSdA_boson_face(:,:,l),a,NMAT)
!    !write(*,*) tmp
!  !enddo
!  call make_traceless_matrix_from_modes(dSdA_boson_face(:,:,l),NMAT,dcmplx( dSdA_boson_face_org(:,l)) )
!enddo


end subroutine Make_bosonic_force_A_face

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! dS/dPhi and dS/dA from fermion part
subroutine Make_fermionic_force(dSdPhi_fermion,dSdA_fermion,UMAT,Phimat,PF,info)
use rational_algorithm
use Dirac_operator
use differential_Dirac
implicit none 

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(in) :: PF(1:sizeD)
complex(kind(0d0)), intent(out) :: dSdPhi_fermion(1:NMAT,1:NMAT,1:num_sites)
!complex(kind(0d0)) :: dSdPhi_fermion_org(1:NMAT,1:NMAT,1:num_sites)
!double precision, intent(out) :: dSdA_fermion_org(1:dimG,1:num_links)
complex(kind(0d0)), intent(out) :: dSdA_fermion(1:NMAT,1:NMAT,1:num_links)
integer, intent(inout) :: info

! vec_r = (D\dagger D + \beta_r)^{-1} F
complex(kind(0d0)) :: vec(1:sizeD,1:N_Remez4)
! Dvec_r = D.vec_r
complex(kind(0d0)) :: Dvec(1:sizeD,1:N_Remez4)
complex(kind(0d0)) :: dvec_eta(1:NMAT,1:NMAT,1:num_sites,1:N_Remez4)
complex(kind(0d0)) :: dvec_lambda(1:NMAT,1:NMAT,1:num_links,1:N_Remez4)
complex(kind(0d0)) :: dvec_chi(1:NMAT,1:NMAT,1:num_faces,1:N_Remez4)
! [\frac{dD^\dagger D}{dPhi_s^a}]_{ij} \vec_j
!complex(kind(0d0)) :: dDdagD_dPhi_vec(1:sizeD, 1:dimG,1:num_sites,1:N_Remez4)
! [\frac{dD^\dagger D}{dA_l^a}]_{ij} \vec_j
!complex(kind(0d0)) :: dDdagD_dA_vec(1:sizeD, 1:dimG,1:num_links,1:N_Remez4)
! 
complex(kind(0d0)) :: dDdPhi_vec(1:sizeD,1:dimG,1:num_sites,1:N_Remez4)
complex(kind(0d0)) :: dDdbPhi_vec(1:sizeD,1:dimG,1:num_sites,1:N_Remez4)
complex(kind(0d0)) :: dDdA_vec(1:sizeD,1:dimG,1:num_links,1:N_Remez4)
!!
complex(kind(0d0)) :: dDdPhi_eta(1:NMAT,1:NMAT,1:num_sites,1:NMAT,1:NMAT,1:num_sites,1:N_Remez4)
complex(kind(0d0)) :: dDdPhi_lambda(1:NMAT,1:NMAT,1:num_links,1:NMAT,1:NMAT,1:num_sites,1:N_Remez4)
complex(kind(0d0)) :: dDdPhi_chi(1:NMAT,1:NMAT,1:num_faces,1:NMAT,1:NMAT,1:num_sites,1:N_Remez4)
complex(kind(0d0)) :: dDdbPhi_lambda(1:NMAT,1:NMAT,1:num_links,1:NMAT,1:NMAT,1:num_sites,1:N_Remez4)
complex(kind(0d0)) :: dDdA_eta(1:NMAT,1:NMAT,1:num_sites,1:NMAT,1:NMAT,1:num_links,1:N_Remez4)
complex(kind(0d0)) :: dDdA_lambda(1:NMAT,1:NMAT,1:num_links,1:NMAT,1:NMAT,1:num_links,1:N_Remez4)
complex(kind(0d0)) :: dDdA_chi(1:NMAT,1:NMAT,1:num_faces,1:NMAT,1:NMAT,1:num_links,1:N_Remez4)

complex(kind(0d0)) :: eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: chi(1:NMAT,1:NMAT,1:num_faces)

complex(kind(0d0)) :: tmp
double precision :: rtmp
integer :: CGite
integer :: r,i,j,s,l,a,b,ss,f,ll,ii,jj


!! for test
!complex(kind(0d0)) :: DdagD(1:sizeD,1:sizeD)
!complex(kind(0d0)) :: tmpmat(1:sizeD,1:sizeD)
!complex(kind(0d0)) :: vec_direct(1:sizeD,1:N_Remez4)
!double precision :: distance
!integer :: j


dSdPhi_fermion=(0d0,0d0)
!dSdPhi_fermion_org=(0d0,0d0)
!dSdA_fermion_org=0d0
dSdA_fermion=(0d0,0d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! comutation by CG
! compute vec_r(i,r) =  (D^\dagger D + \beta_r)^{-1} F(i)
!write(*,*) "   mmBiCG start"
!call check_Dirac(UMAT,Phi)
call mmBiCG( vec, PF, Remez_beta4, sizeD, N_Remez4, epsilon, &
             CG_max, info, CGite, UMAT, PhiMat, Prod_DdagD )
!write(*,*) "   mmBiCG end", cgite
do r=1,N_Remez4
  call vec_to_mat(eta,lambda,chi,vec(:,r))
  call prod_Dirac(Dvec(:,r),vec(:,r),sizeD,UMAT,PhiMat)
  
  !call prod_dDdPhi(dDdPhi_vec(:,:,:,r),vec(:,r),UMAT)
  call prod_dDdPhi(&
    dDdPhi_eta(:,:,:,:,:,:,r),&
    dDdPhi_chi(:,:,:,:,:,:,r), &
    eta,chi,UMAT)
  call prod_dDdbPhi(&
    dDdbPhi_lambda(:,:,:,:,:,:,r),&
    lambda,UMAT)
  !call prod_dDdA(dDdA_vec(:,:,:,r),eta,lambda,chi,UMAT,PhiMat)
  call prod_dDdA(&
    dDdA_eta(:,:,:,:,:,:,r),&
    dDdA_lambda(:,:,:,:,:,:,r), &
    dDdA_chi(:,:,:,:,:,:,r),&
    eta,lambda,chi,UMAT,PhiMat)
enddo

do r=1,N_Remez4
  call vec_to_mat( &
    Dvec_eta(:,:,:,r),&
    Dvec_lambda(:,:,:,r),&
    Dvec_chi(:,:,:,r), &
    Dvec(:,r) )
enddo

do s=1,num_sites
  do ii=1,NMAT
    do jj=1,NMAT
      do r=1,N_Remez4
        tmp=(0d0,0d0)
        do ss=1,num_sites
          do i=1,NMAT
            do j=1,NMAT
              tmp=tmp&
                +conjg( Dvec_eta(i,j,ss,r) )*dDdPhi_eta(i,j,ss,ii,jj,s,r) !&
                !+conjg( dDdbPhi_eta(i,j,ss,ii,jj,s,r) )*Dvec_eta(i,j,ss,r)
            enddo
          enddo
        enddo
        do l=1,num_links
          do i=1,NMAT
            do j=1,NMAT
              tmp=tmp&
                !+conjg( Dvec_lambda(i,j,l,r) )*dDdPhi_lambda(i,j,l,ii,jj,s,r) &
                +conjg( dDdbPhi_lambda(i,j,l,jj,ii,s,r) )*Dvec_lambda(i,j,l,r)
            enddo
          enddo
        enddo
        do f=1,num_faces
          do i=1,NMAT
            do j=1,NMAT
              tmp=tmp&
                +conjg( Dvec_chi(i,j,f,r) )*dDdPhi_chi(i,j,f,ii,jj,s,r)! &
                !+conjg( dDdbPhi_chi(i,j,f,ii,jj,s,r) )*Dvec_chi(i,j,f,r)
            enddo
          enddo
        enddo
        dSdPhi_fermion(ii,jj,s)=dSdPhi_fermion(ii,jj,s) &
          - dcmplx(Remez_alpha4(r))*tmp
      enddo
    enddo
  enddo
enddo



do ll=1,num_links
  do ii=1,NMAT
    do jj=1,NMAT
      do r=1,N_Remez4
        tmp=(0d0,0d0)
        do ss=1,num_sites
          do i=1,NMAT
            do j=1,NMAT
              tmp=tmp &
                + conjg( Dvec_eta(i,j,ss,r) ) * dDdA_eta(i,j,ss,ii,jj,ll,r)&
                + conjg( dDdA_eta(i,j,ss,jj,ii,ll,r) ) * Dvec_eta(i,j,ss,r) 
            enddo
          enddo
        enddo
        do l=1,num_links
          do i=1,NMAT
            do j=1,NMAT
              tmp=tmp &
                +conjg( Dvec_lambda(i,j,l,r) ) * dDdA_lambda(i,j,l,ii,jj,ll,r)&
                + conjg( dDdA_lambda(i,j,l,jj,ii,ll,r) ) * Dvec_lambda(i,j,l,r) 
            enddo
          enddo
        enddo
        do f=1,num_faces
          do i=1,NMAT
            do j=1,NMAT
              tmp=tmp &
                + conjg( Dvec_chi(i,j,f,r) ) * dDdA_chi(i,j,f,ii,jj,ll,r) &
                + conjg( dDdA_chi(i,j,f,jj,ii,ll,r) ) * Dvec_chi(i,j,f,r) 
            enddo
          enddo
        enddo
        dSdA_fermion(ii,jj,ll)=dSdA_fermion(ii,jj,ll) & 
          - Remez_alpha4(r) * tmp 
      enddo
    enddo
  enddo
enddo

!do l=1,num_links
  !do a=1,dimG
    !call trace_MTa(tmp,dSdA_fermion(:,:,l),a,NMAT)
    !dSdA_fermion_org(a,l)=dble(tmp)
  !enddo
!enddo
end subroutine Make_fermionic_force

end module forces
