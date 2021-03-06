!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Make_bosonic_force(dSdPhi,dSdA,UMAT,PhiMat)
use SUN_generators, only : trace_mta, make_traceless_matrix_from_modes
#ifdef PARALLEL
use parallel
#endif
implicit none

complex(kind(0d0)), intent(out) :: dSdPhi(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(out) :: dSdA(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) :: dSdPhi_boson_mass(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: dSdPhi_boson_site(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: dSdPhi_boson_link(1:NMAT,1:NMAT,1:num_sites)
!!!
complex(kind(0d0)) :: dSdA_boson_link(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: dSdA_boson_face(1:NMAT,1:NMAT,1:num_links)

integer :: i,j,s,l,f
complex(kind(0d0)) :: tmp
double precision :: tmp_bosonic_force_Phi
double precision :: tmp_bosonic_force_A
double precision :: rtmp1, rtmp2, rtmp3

dSdPhi=(0d0,0d0)
dSdPhi_boson_mass=(0d0,0d0)
dSdPhi_boson_site=(0d0,0d0)
dSdPhi_boson_link=(0d0,0d0)
dSdA=(0d0,0d0)
dSdA_boson_face=(0d0,0d0)
dSdA_boson_link=(0d0,0d0)

!! force for Phi from boson
if(pb_mass==0) call Make_bosonic_force_Phi_mass(dSdPhi_boson_mass,PhiMat)
if(pb_site==0) call Make_bosonic_force_Phi_site(dSdPhi_boson_site,PhiMat)
if(pb_link==0) call Make_bosonic_force_Phi_link(dSdPhi_boson_link,UMAT,PhiMat)

!! force for A from boson
if(pb_link==0) call Make_bosonic_force_A_link(dSdA_boson_link,UMAT,PhiMat)
if(pb_face==0) call Make_bosonic_force_A_face(dSdA_boson_face,UMAT)

dSdPhi= dSdPhi_boson_mass &   ! mass part
+ dSdPhi_boson_site  & ! site part
+ dSdPhi_boson_link  

dSdA = dSdA_boson_link  &
+ dSdA_boson_face 

if ( force_measurement == 1 ) then 
  tmp_bosonic_force_Phi=0d0
  rtmp2=0d0
  do s=1,num_sites
    rtmp1=0d0
    tmp=(0d0,0d0)
    do j=1,NMAT
      do i=1,NMAT
        tmp=dSdPhi_boson_mass(i,j,s) &
          + dSdPhi_boson_site(i,j,s) &
          + dSdPhi_boson_link(i,j,s)
        rtmp1 = rtmp1 + dble(tmp * dconjg(tmp))
      enddo
    enddo
    !tmp_bosonic_force_Phi= tmp_bosonic_force_Phi + dsqrt(rtmp1)
    rtmp2 = rtmp2 + dsqrt(rtmp1) 
  enddo
#ifdef PARALLEL
  call MPI_REDUCE(rtmp2,tmp_bosonic_force_Phi,1,MPI_DOUBLE_PRECISION, &
    MPI_SUM,0,MPI_COMM_WORLD,IERR)
#endif
  tmp_bosonic_force_Phi= tmp_bosonic_force_Phi / dble(global_num_sites)
  

  tmp_bosonic_force_A=0d0
  rtmp2=0d0
  do l=1,num_links
    rtmp1=0d0
    tmp=(0d0,0d0)
    do j=1,NMAT
      do i=1,NMAT
        tmp=dSdA_boson_link(i,j,l) + dSdA_boson_face(i,j,l)
        rtmp1 = rtmp1 + dble(tmp * dconjg(tmp))
      enddo
    enddo
    rtmp2=rtmp2+dsqrt(rtmp1) 
    !tmp_bosonic_force_A= tmp_bosonic_force_A + dsqrt(rtmp1)
  enddo
#ifdef PARALLEL
  call MPI_REDUCE(rtmp2,tmp_bosonic_force_A,1,MPI_DOUBLE_PRECISION, &
    MPI_SUM,0,MPI_COMM_WORLD,IERR)
#endif
  tmp_bosonic_force_A= tmp_bosonic_force_A / dble(global_num_links)

bosonic_force_Phi = bosonic_force_Phi + tmp_bosonic_force_Phi 
bosonic_force_A = bosonic_force_A + tmp_bosonic_force_A 
b_phi_count=b_phi_count+1
b_A_count=b_A_count+1
endif
end subroutine Make_bosonic_force

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Make_bosonic_force_nomass(dSdPhi,dSdA,UMAT,PhiMat)
use SUN_generators, only : trace_mta, make_traceless_matrix_from_modes
#ifdef PARALLEL
use parallel
#endif
implicit none

complex(kind(0d0)), intent(out) :: dSdPhi(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(out) :: dSdA(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) :: dSdPhi_boson_site(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: dSdPhi_boson_link(1:NMAT,1:NMAT,1:num_sites)
!!!
complex(kind(0d0)) :: dSdA_boson_link(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: dSdA_boson_face(1:NMAT,1:NMAT,1:num_links)

integer :: i,j,s,l,f
complex(kind(0d0)) :: tmp
double precision :: tmp_bosonic_force_Phi
double precision :: tmp_bosonic_force_A
double precision :: rtmp1, rtmp2, rtmp3

dSdPhi=(0d0,0d0)
dSdPhi_boson_site=(0d0,0d0)
dSdPhi_boson_link=(0d0,0d0)
dSdA=(0d0,0d0)
dSdA_boson_face=(0d0,0d0)
dSdA_boson_link=(0d0,0d0)
!! force for Phi from boson
if(pb_site==0) call Make_bosonic_force_Phi_site(dSdPhi_boson_site,PhiMat)
if(pb_link==0) call Make_bosonic_force_Phi_link(dSdPhi_boson_link,UMAT,PhiMat)

!! force for A from boson
if(pb_link==0) call Make_bosonic_force_A_link(dSdA_boson_link,UMAT,PhiMat)
if(pb_face==0) call Make_bosonic_force_A_face(dSdA_boson_face,UMAT)

dSdPhi= + dSdPhi_boson_site  & ! site part
+ dSdPhi_boson_link  

dSdA = dSdA_boson_link  &
+ dSdA_boson_face 

end subroutine Make_bosonic_force_nomass

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Forces from mass terms, Sites, Links, Faces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! dS/dPhi from mass action
subroutine Make_bosonic_force_Phi_mass(dSdPhi_boson_mass,PhiMat)
implicit none

complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(out) :: dSdPhi_boson_mass(1:NMAT,1:NMAT,1:num_sites)
integer :: s,i,j

dSdPhi_boson_mass=(0d0,0d0)
do s=1,num_sites
  do j=1,NMAT
    do i=1,NMAT
  dSdPhi_boson_mass(i,j,s) = dcmplx(mass_square_phi*0.5d0) &!*  dconjg(Phi(a,s)) &
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

complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(out) :: dSdPhi_boson_site(1:NMAT,1:NMAT,1:num_sites)
integer :: a,s
complex(kind(0d0)) :: comm1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: comm2(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: tmp(1:dimG), tmp1,tmp2
integer :: i,j

dSdPhi_boson_site=(0d0,0d0)

do s=1,num_sites
  call matrix_commutator(comm1,phimat(:,:,s),phimat(:,:,s),'N','C')
  call matrix_commutator(comm2,phimat(:,:,s),comm1,'C','N')
  !do a=1,dimG
  do j=1,NMAT
    do i=1,NMAT
    dSdPhi_boson_site(i,j,s)=dcmplx(0.5d0*alpha_s(s))*comm2(i,j)& 
      * dcmplx( overall_factor )
      !* cmplx( one_ov_2g2N )
    enddo
  enddo
enddo

end subroutine Make_bosonic_force_Phi_site

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! dS/Dphi from link action
!complex(kind(0d0)) function dSdPhi_boson_link(a,s)
subroutine Make_bosonic_force_Phi_link(dSdPhi_boson_link,UMAT,PhiMat)
use matrix_functions, only : matrix_3_product, make_matrix_traceless
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(out) :: dSdPhi_boson_link(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: dPhi(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT), tmpmat2(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: Ud_dPhibar_U(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)) :: dPhibar(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: Ud_dPhibar_U(1:NMAT,1:NMAT)
complex(kind(0d0)) :: dPhibar(1:NMAT,1:NMAT)
complex(kind(0d0)) :: trace
integer :: a,s
integer :: l_i
integer :: i,j,k,l

dSdPhi_boson_link=(0d0,0d0)
do s=1,num_sites
  do i=1,linkorg_to_s(s)%num_
    l_i=linkorg_to_s(s)%Labels_(i)
    !write(*,*) "org",s,l_i
    call Make_diff_PhiMat(dPhi,l_i,UMAT,PhiMat)
    call matrix_3_product(Ud_dPhibar_U,UMAT(:,:,l_i),dPhi,UMAT(:,:,l_i),'C','C','N')
    dSdPhi_boson_link(:,:,s)=dSdPhi_boson_link(:,:,s) &
      + alpha_l(l_i) *  Ud_dPhibar_U &
      * U1Rfactor_link(l_i)**2d0 ! * U1R_ratio(l_i)**2d0
  enddo

  !write(*,*) "tip",linktip_from_s(s)%num_
  !write(*,*) "org",linkorg_to_s(s)%num_
  do i=1,linktip_from_s(s)%num_
    l_i=linktip_from_s(s)%Labels_(i)
    call Make_diff_PhiMat(dPhi,l_i,UMAT,PhiMat)
    do k=1,NMAT
      do j=1,NMAT
        dPhibar(j,k)=dconjg(dPhi(k,j))
      enddo
    enddo
    dSdPhi_boson_link(:,:,s)=dSdPhi_boson_link(:,:,s) - alpha_l(l_i) &
      * dPhibar
  enddo
enddo

dSdPhi_boson_link = dSdPhi_boson_link * dcmplx( overall_factor )
end subroutine Make_bosonic_force_Phi_link

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! dS/DA from link action
subroutine Make_bosonic_force_A_link(dSdA_boson_link,UMAT,PhiMat)
use SUN_generators, only : trace_MTa, make_traceless_matrix_from_modes
use matrix_functions, only : matrix_3_product, matrix_commutator,make_matrix_traceless
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(out) :: dSdA_boson_link(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: dPhi(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Phi_tip(1:NMAT,1:NMAT)
complex(kind(0d0)) :: UPhiUinv(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Comm(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmp,trace
integer :: a,l,i,j,k

dSdA_boson_link=(0d0,0d0)
do l=1,num_links
  call Make_diff_PhiMAT(dPhi, l,UMAT,PhiMat)

  Phi_tip=(0d0,1d0)*PhiMat(:,:,link_tip(l))

  ! U_l.Phi_tip
!  ! i*U_l.Phi_tip.U_l^\dagger
  call matrix_3_product(UPhiUinv,UMAT(:,:,l),Phi_tip,UMAT(:,:,l),'N','N','C')

  ! i[ UPhiUinv, dPhi^\dagger ]
  call matrix_commutator(Comm,UPhiUinv,dPhi,'N','C')
  Comm=Comm*U1Rfactor_link(l)**2d0 ! * U1R_ratio(l)**2d0

  do j=1,NMAT
    do i=1,NMAT
      dSdA_boson_link(i,j,l) = dcmplx(alpha_l(l)) * ( Comm(i,j) + dconjg( Comm(j,i)) )
    enddo
  enddo
  call make_matrix_traceless(dSdA_boson_link(:,:,l))
enddo

dSdA_boson_link=dSdA_boson_link * dcmplx(overall_factor)

end subroutine Make_bosonic_force_A_link


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! dS/DA from face action
subroutine Make_bosonic_force_A_face(dSdA_boson_face,UMAT)
use SUN_generators, only : trace_mta, make_traceless_matrix_from_modes
use matrix_functions, only : matrix_power,matrix_product,matrix_3_product,make_matrix_traceless
implicit none
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(out) :: dSdA_boson_face(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT), Ufm(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Omega(1:NMAT,1:NMAT)
!!!!!!!!!
complex(kind(0d0)) :: diff_Omega(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
complex(kind(0d0)) :: dUfdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
complex(kind(0d0)) :: sinU(1:NMAT,1:NMAT)
complex(kind(0d0)) :: CosUinv(1:NMAT,1:NMAT)
complex(kind(0d0)) :: dCosUinvdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
complex(kind(0d0)) :: dSinUdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
complex(kind(0d0)) :: dBdA(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Mae(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Ushiro(1:NMAT,1:NMAT)
complex(kind(0d0)) :: l_factor
!!!!!!!!!
complex(kind(0d0)) :: trace,tmp
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Bval
integer :: a,l,f!,ll
integer :: l_label
integer :: i,j,k,ii,jj

dSdA_boson_face=(0d0,0d0)
do l=1,num_links
  do k=1,face_in_l(l)%num_
    f=face_in_l(l)%label_(k)
    sinU=(0d0,0d0)
    cosUinv=(0d0,0d0)
    call Make_face_variable(Uf,f,UMAT)
    if( m_omega == 0 ) then
      call Make_moment_map0(Omega(:,:),Uf(:,:))
    elseif( m_omega == -1 ) then
      call Make_moment_map_adm(Omega(:,:),Uf(:,:))
    else
      call matrix_power(Ufm(:,:),Uf(:,:),m_omega)
      call Make_moment_map(Omega(:,:),Ufm(:,:))
      call calc_sinU_and_cosUinv(sinU(:,:),cosUinv(:,:),Ufm(:,:))
    endif
    do l_label=1,links_in_f(f)%num_
      if ( l == links_in_f(f)%link_labels_(l_label) ) exit
    enddo

  if( m_omega == 0 ) then 
    dSinUdA=(0d0,0d0)
    call calc_dSinUdA(dSinUdA(:,:,:,:),Uf(:,:),UMAT,f,l_label)
    diff_Omega=(0d0,-1d0)*dSinUdA
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! make diff_Omega traceless
    tmpmat=(0d0,0d0)
    do ii=1,NMAT
      tmpmat=tmpmat+diff_Omega(:,:,ii,ii)
    enddo
    do ii=1,NMAT
      diff_Omega(:,:,ii,ii)=diff_Omega(:,:,ii,ii)-tmpmat/dcmplx(dble( NMAT ))
    enddo
    !!!!!!!!!!!!!!!!!!!!!
    ! make traceless
    if( NMAT > 2 ) then
      do jj=1,NMAT
        do ii=1,NMAT
          call make_matrix_traceless(diff_Omega(:,:,ii,jj))
        enddo
      enddo
    endif
    !!!!!!!!!!!!!!!!!!!!!
      
    do jj=1,NMAT
      do ii=1,NMAT
        do j=1,NMAT
          do i=1,NMAT
            dSdA_boson_face(ii,jj,l)=dSdA_boson_face(ii,jj,l) &
              + (0.5d0,0d0) &
                *dcmplx( overall_factor * alpha_f(f)*beta_f(f)*beta_f(f) )&
                *diff_Omega(i,j,ii,jj) *  Omega(j,i) 
          enddo
        enddo
      enddo
    enddo
 
  elseif( m_omega == -1 ) then
    call calc_dUfdA_dBdA(Mae,Ushiro,l_factor,dBdA,Umat,f,l_label)
    Bval=(1d0,0d0)
    do i=1,NMAT
      Bval=Bval - (1d0,0d0)/(e_max*e_max) * dcmplx( 2d0 - 2d0*dble(Uf(i,i)) )
    enddo

    !! dBdA contribution
    trace=(0d0,0d0)
    do i=1,NMAT
      do j=1,NMAT
        trace=trace+Omega(i,j)*Omega(j,i)
      enddo
    enddo
    trace=trace&
      *(-0.5d0,0d0)*dcmplx( overall_factor * alpha_f(f)*beta_f(f)*beta_f(f) )/Bval
    !!!!!!!!!!!!!!!!!
    dSdA_boson_face(:,:,l) = dSdA_boson_face(:,:,l) + trace * dBdA
    !!!!!!!!!!!!!!!!!

    !! dOmega/dA contribution
    call matrix_3_product(tmpmat,Ushiro,Omega,Mae)
    call matrix_3_product(tmpmat,Mae,Omega,Ushiro,'C','N','C',(1d0,0d0),'ADD')
    tmp=(0.5d0,0d0)*dcmplx( overall_factor * alpha_f(f)*beta_f(f)*beta_f(f) )&
      *l_factor*(0d0,-1d0)/Bval
    call make_matrix_traceless(tmpmat)
    !!!!!!!!!!!!!!!!!
    dSdA_boson_face(:,:,l) = dSdA_boson_face(:,:,l) + tmp * tmpmat
    !!!!!!!!!!!!!!!!!

  else
    dCosUinvdA=(0d0,0d0)
    dSinUdA=(0d0,0d0)
    call calc_dCosUinvdA_dSinUdA(&
      dCosUinvdA(:,:,:,:),dSinUdA(:,:,:,:),&
      cosUinv(:,:),Uf(:,:),UMAT,f,l_label)

    diff_Omega=(0d0,0d0)
    tmpmat=(0d0,0d0)
    do jj=1,NMAT
      do ii=1,NMAT
         call matrix_product(&
           diff_Omega(:,:,ii,jj),dSinUdA(:,:,ii,jj),cosUinv(:,:))
         call zgemm('N','N',NMAT,NMAT,NMAT,(0d0,-1d0), &
           SinU(:,:), NMAT, &
           dCosUinvdA(:,:,ii,jj), NMAT, &
           (0d0,-1d0), diff_Omega(:,:,ii,jj), NMAT)
         call zgemm('N','N',NMAT,NMAT,NMAT,(0d0,-1d0), &
           dCosUinvdA(:,:,ii,jj), NMAT, &
           SinU(:,:), NMAT, &
           (1d0,0d0), diff_Omega(:,:,ii,jj), NMAT)
         call zgemm('N','N',NMAT,NMAT,NMAT,(0d0,-1d0), &
           CosUinv(:,:), NMAT, &
           dSinUdA(:,:,ii,jj), NMAT, &
           (1d0,0d0), diff_Omega(:,:,ii,jj), NMAT)
         if( ii==jj ) then
           tmpmat=tmpmat+diff_Omega(:,:,ii,ii)
         endif
       enddo
     enddo
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! make diff_Omega traceless
     do ii=1,NMAT
       diff_Omega(:,:,ii,ii)=diff_Omega(:,:,ii,ii)-tmpmat/dcmplx(dble( NMAT ))
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
              - trace / dcmplx(dble( NMAT ))
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
                  *dcmplx( overall_factor * alpha_f(f)*beta_f(f)*beta_f(f) / dble(m_omega)) &
                  *diff_Omega(i,j,ii,jj) *  Omega(j,i) 
            enddo
          enddo
        enddo
      enddo
    endif
  enddo
  !call make_matrix_traceless(dSdA_boson_face(:,:,l))
enddo

end subroutine Make_bosonic_force_A_face

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! dS/dPhi and dS/dA from fermion part
subroutine Make_fermionic_force(dSdPhi_fermion,dSdA_fermion,UMAT,Phimat,PF_eta,PF_lambda,PF_chi,info)
use rational_algorithm
use Dirac_operator
use differential_Dirac
#ifdef PARALLEL
use parallel
#endif
implicit none 

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
!complex(kind(0d0)) :: PF(1:sizeD)
complex(kind(0d0)), intent(in) :: PF_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: PF_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PF_chi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), intent(out) :: dSdPhi_fermion(1:NMAT,1:NMAT,1:num_sites)
!complex(kind(0d0)) :: dSdPhi_fermion_org(1:NMAT,1:NMAT,1:num_sites)
!double precision, intent(out) :: dSdA_fermion_org(1:dimG,1:num_links)
complex(kind(0d0)), intent(out) :: dSdA_fermion(1:NMAT,1:NMAT,1:num_links)
integer, intent(inout) :: info

complex(kind(0d0)) :: Deta(1:NMAT,1:NMAT,1:num_necessary_sites,1:N_Remez4)
complex(kind(0d0)) :: Dlambda(1:NMAT,1:NMAT,1:num_necessary_links,1:N_Remez4)
complex(kind(0d0)) :: Dchi(1:NMAT,1:NMAT,1:num_necessary_faces,1:N_Remez4)

complex(kind(0d0)) :: eta(1:NMAT,1:NMAT,1:num_necessary_sites,N_Remez4)
complex(kind(0d0)) :: lambda(1:NMAT,1:NMAT,1:num_necessary_links,N_Remez4)
complex(kind(0d0)) :: chi(1:NMAT,1:NMAT,1:num_necessary_faces,N_Remez4)
complex(kind(0d0)) :: tmp_eta(1:NMAT,1:NMAT,1:num_sites,N_Remez4)
complex(kind(0d0)) :: tmp_lambda(1:NMAT,1:NMAT,1:num_links,N_Remez4)
complex(kind(0d0)) :: tmp_chi(1:NMAT,1:NMAT,1:num_faces,N_Remez4)
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT,1:num_necessary_faces) 
complex(kind(0d0)) :: Ufm(1:NMAT,1:NMAT,1:num_necessary_faces) 

complex(kind(0d0)) :: tmp
double precision :: rtmp1,rtmp2
double precision :: tmp_fermionic_force_Phi
double precision :: tmp_fermionic_force_A
integer :: CGite
integer :: r,i,j,s,l,a,b,ss,f,ll,ii,jj
integer :: k,l_label



dSdPhi_fermion=(0d0,0d0)
dSdA_fermion=(0d0,0d0)
Deta=(0d0,0d0)
Dlambda=(0d0,0d0)
Dchi=(0d0,0d0)
chi=(0d0,0d0)
lambda=(0d0,0d0)
chi=(0d0,0d0)
tmp_chi=(0d0,0d0)
tmp_lambda=(0d0,0d0)
tmp_chi=(0d0,0d0)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! comutation by CG
! compute vec_r(i,r) =  (D^\dagger D + \beta_r)^{-1} F(i)
!write(*,*) "   mmBiCG start"
!call check_Dirac(UMAT,Phi)
!call mat_to_vec(PF,PF_eta,PF_lambda,PF_chi)

!! preparation
call mmBiCG( &
  tmp_eta, tmp_lambda, tmp_chi, &
  PF_eta, PF_lambda, PF_chi, &
  Remez_beta4, epsilon, CG_max, info, CGite, UMAT, PhiMat, Prod_DdagD)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! if CG is failed, we should reject
if( info==1 ) return
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


do r=1,N_Remez4
  do s=1,num_sites
    eta(:,:,s,r)=tmp_eta(:,:,s,r)
  enddo
  do l=1,num_links
    lambda(:,:,l,r)=tmp_lambda(:,:,l,r)
  enddo
  do f=1,num_faces
    chi(:,:,f,r)=tmp_chi(:,:,f,r)
  enddo
#ifdef PARALLEL
  call syncronize_sites(eta(:,:,:,r))
  call syncronize_links(lambda(:,:,:,r))
  call syncronize_faces(chi(:,:,:,r))
#endif
  call prod_Dirac( &
    Deta(:,:,1:num_sites,r),Dlambda(:,:,1:num_links,r),Dchi(:,:,1:num_faces,r),&
    eta(:,:,:,r),lambda(:,:,:,r),chi(:,:,:,r), &
    UMAT,PhiMat)
#ifdef PARALLEL
  call syncronize_sites(Deta(:,:,:,r))
  call syncronize_links(Dlambda(:,:,:,r))
  call syncronize_faces(Dchi(:,:,:,r))
#endif
enddo

! \Phi(ss,ii,jj)での微分
do ss=1,num_sites
  call calc_force_from_dDdPhi(dSdPhi_fermion(:,:,ss),eta,lambda,chi,Deta,Dlambda,Dchi,Umat,ss)

  call make_matrix_traceless(dSdPhi_fermion(:,:,ss))
enddo

! A(ll,ii,jj)で微分
if( p5==0 ) then 
  do f=1,num_necessary_faces
    call make_face_variable(Uf(:,:,f),f,Umat)
    if( m_omega==1 ) then
      Ufm(:,:,f)=Uf(:,:,f)
    elseif( m_omega>=2 ) then
      call matrix_power(Ufm(:,:,f),Uf(:,:,f),m_omega)
    endif
  enddo
endif
do ll=1,num_links
  call calc_force_from_dDdA(dSdA_fermion(:,:,ll),&
      eta,lambda,chi,Deta,Dlambda,Dchi,UMat,PhiMat,Uf,Ufm,ll)
  call make_matrix_traceless(dSdA_fermion(:,:,ll))
enddo


if ( force_measurement == 1 ) then 
  tmp_fermionic_force_Phi=0d0
  rtmp2=0d0
  do s=1,num_sites
    rtmp1=0d0
    tmp=(0d0,0d0)
    do j=1,NMAT
      do i=1,NMAT
        rtmp1 = rtmp1 + dble( dSdPhi_fermion(i,j,s) * dconjg( dSdPhi_fermion(i,j,s) ) )
      enddo
    enddo
    rtmp2=rtmp2+dsqrt(rtmp1)
    !tmp_fermionic_force_Phi= tmp_fermionic_force_Phi + dsqrt(rtmp1)
  enddo
#ifdef PARALLEL
  call MPI_REDUCE(rtmp2,tmp_fermionic_force_Phi,1,MPI_DOUBLE_PRECISION, &
    MPI_SUM,0,MPI_COMM_WORLD,IERR)
#endif
  tmp_fermionic_force_Phi= tmp_fermionic_force_Phi / dble(global_num_sites)
  
  tmp_fermionic_force_A=0d0
  rtmp2=0d0
  do l=1,num_links
    rtmp1=0d0
    do j=1,NMAT
      do i=1,NMAT
        rtmp1 = rtmp1 + dble(dSdA_fermion(i,j,l) * dconjg(dSdA_fermion(i,j,l)) )
      enddo
    enddo
    rtmp2=rtmp2+dsqrt(rtmp1)
    !tmp_fermionic_force_A= tmp_fermionic_force_A + dsqrt(rtmp1)
  enddo
#ifdef PARALLEL
  call MPI_REDUCE(rtmp2,tmp_fermionic_force_A,1,MPI_DOUBLE_PRECISION, &
    MPI_SUM,0,MPI_COMM_WORLD,IERR)
#endif
  tmp_fermionic_force_A= tmp_fermionic_force_A / dble(global_num_links)

  fermionic_force_Phi = fermionic_force_Phi + tmp_fermionic_force_Phi 
  fermionic_force_A = fermionic_force_A + tmp_fermionic_force_A 
  f_phi_count=f_phi_count+1
  f_A_count=f_A_count+1
endif
end subroutine Make_fermionic_force

