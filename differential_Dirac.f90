!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! module for fermion matrix
module differential_Dirac
use global_parameters
use global_subroutines
use SUN_generators
use matrix_functions
implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to return product of d/d\Phi D
subroutine prod_dDdPhi(dDdPhi_chi,vec,UMAT,Phi)
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links) 
complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)), intent(in) :: vec(1:sizeD)
! dDdPhi_chi(a,s,i) = \frac{d}{d\Phi_{s,a}} ( D \Psi )_i
complex(kind(0d0)), intent(out) :: dDdPhi_chi(1:sizeD, 1:dimG,1:num_sites)

!complex(kind(0d0)) :: eta_ele(1:dimG,1:num_sites)
complex(kind(0d0)) :: chi_ele(1:dimG,1:num_faces)

integer :: s,f
integer :: i,j,r
integer :: a,b,c

!! site: (s,a) <--> r=a+dimG*(s-1)
!! link: (l,a) <--> r=a+dimG*(num_sites + l -1)
!! face: (f,a) <--> r=a+dimG*(num_sites + num_links + f -1 )
!! preparation
!dDdPhi_chi=(0d0,0d0)
!do s=1,num_sites
!  do i=1,dimG
!    eta_ele(i,s)=vec(i+dimG*(s-1))
!  enddo
!enddo
do f=1,num_faces
  do a=1,dimG
    chi_ele(a,f)=vec(face_index(a,f))
  enddo
enddo

dDdPhi_chi=(0d0,0d0)

if( p1 == 0 ) then
    !write(*,*) p1
!! (1) Dirac from site
do s=1,num_sites
  do r=1,NZF
    a=NZF_index(1,r)
    b=NZF_index(2,r)
    c=NZF_index(3,r)
    j=a+dimG*(s-1) ! j <=> (a,s)
    dDdPhi_chi(j,b,s)=dDdPhi_chi(j,b,s)& 
      +(0d0,-0.5d0)*dcmplx(alpha_s(s))*NZF_value(r) &
      * vec(c+dimG*(s-1)) * overall_factor ! \eta_{s,c}
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
  do r=1,NZF
    a=NZF_index(1,r)
    b=NZF_index(2,r)
    c=NZF_index(3,r)
    dDdPhi_chi(face_index(a,f),b,s)= dDdPhi_chi(face_index(a,f),b,s) &
      +(0d0,-2d0)*dcmplx(alpha_f(f))*NZF_value(r)*chi_ele(c,f)*overall_factor
  enddo
enddo
endif

!dDdPhi_chi=dDdPhi_chi*overall_factor
  

end subroutine prod_dDdPhi


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to return product of d/d\Phi D^\dagger
subroutine prod_dDdbPhi(dDdbPhi_chi,vec,UMAT,Phi)
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links) 
complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)), intent(in) :: vec(1:sizeD)
complex(kind(0d0)), intent(out) :: dDdbPhi_chi(1:sizeD,1:dimG,1:num_sites)

complex(kind(0d0)) :: eta_ele(1:dimG,1:num_sites)
complex(kind(0d0)) :: eta_mat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: lambda_ele(1:dimG,1:num_links)
complex(kind(0d0)) :: lambda_mat(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: chi_ele(1:dimG,1:num_faces)
complex(kind(0d0)) :: chi_mat(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: Phi_mat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: bPhi_mat(1:NMAT,1:NMAT,1:num_sites)


complex(kind(0d0)) :: Udag(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: comm(1:NMAT,1:NMAT)


complex(kind(0d0)) :: trace
complex(kind(0d0)) :: T(1:NMAT,1:NMAT,1:dimG)

integer :: s,f,r
integer :: i,j,k,l
integer :: a,b,c,d

!! preparation
dDdbPhi_chi=(0d0,0d0)

do s=1,num_sites
  do a=1,dimG
    eta_ele(a,s)=vec(site_index(a,s))
  enddo
  call make_traceless_matrix_from_modes(eta_mat(:,:,s),NMAT,eta_ele(:,s))
enddo
do l=1,num_links
  do a=1,dimG    
    lambda_ele(a,l)=vec(link_index(a,l))
  enddo          
  call make_traceless_matrix_from_modes(lambda_mat(:,:,l),NMAT,lambda_ele(:,l))
enddo            
do f=1,num_faces 
  do a=1,dimG    
    chi_ele(a,f)=vec(face_index(a,f))
  enddo          
  call make_traceless_matrix_from_modes(chi_mat(:,:,f),NMAT,chi_ele(:,f))
enddo            
do s=1,num_sites 
  call make_traceless_matrix_from_modes(Phi_mat(:,:,s),NMAT,Phi(:,s))
  call make_traceless_matrix_from_modes(bPhi_mat(:,:,s),NMAT,dconjg(Phi(:,s)))
enddo            


! (1) Dirac from site
! no contribution

!! (2) Dirac from link 1
! no contribution

if ( p3 == 0 ) then
!! (3) Dirac from link 2
call make_SUN_generators(T,NMAT)
do l=1,num_links
  do b=1,dimG
    call MTMd(tmpmat1,UMAT(:,:,l),b,NMAT)
    call Matrix_Commutator(comm,tmpmat1, lambda_mat(:,:,l)) !,'N','N')
    do a=1,dimG
      call trace_MTa(trace,comm,a,NMAT)
      dDdbPhi_chi(link_index(a,l),b,link_tip(l)) = &
        dDdbPhi_chi(link_index(a,l),b,link_tip(l)) &
        + (1d0,0d0) * dcmplx( alpha_l(l) ) * trace*overall_factor
    enddo
  enddo
  !!!
  s=link_org(l)
  do r=1,NZF
    a=NZF_index(1,r)
    b=NZF_index(2,r)
    c=NZF_index(3,r)
    dDdbPhi_chi(link_index(a,l),b,s) = dDdbPhi_chi(link_index(a,l),b,s) &
      + dcmplx( alpha_l(l) ) * (0d0,1d0) * NZF_value(r) * lambda_ele(c,l)&
      * overall_factor
  enddo
enddo
endif

!! (4) Dirac from face 1
!   no contribution

!dDdbPhi_chi=dDdbPhi_chi*overall_factor
end subroutine prod_dDdbPhi


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to return product of d/dA D . vec
subroutine prod_dDdA(dDdA_vec,vec,UMAT,PhiMat)
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links) 
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(in) :: vec(1:sizeD)
complex(kind(0d0)), intent(out) :: dDdA_vec(1:sizeD, 1:dimG,1:num_links)

complex(kind(0d0)) :: eta_ele(1:dimG,1:num_sites)
complex(kind(0d0)) :: eta_mat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: lambda_ele(1:dimG,1:num_links)
complex(kind(0d0)) :: lambda_mat(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: chi_ele(1:dimG,1:num_faces)
complex(kind(0d0)) :: chi_mat(1:NMAT,1:NMAT,1:num_faces)

! d/dA_{ll,ii,jj} (D\Psi)_{s,i,j}=dDdA_eta(i,j,s,ii,jj,ll)
complex(kind(0d0)) :: dDdA_eta(1:NMAT,1:NMAT,1:num_sites,1:NMAT,1:NMAT,1:num_links)
! d/dA_{ll,ii,jj} (D\Psi)_{l,i,j}=dDdA_lambda(i,j,l,ii,jj,ll)
complex(kind(0d0)) :: dDdA_lambda(1:NMAT,1:NMAT,1:num_links,1:NMAT,1:NMAT,1:num_links)
! d/dA_{ll,ii,jj} (D\Psi)_{f,i,j}=dDdA_chi(i,j,f,ii,jj,ll)
complex(kind(0d0)) :: dDdA_chi(1:NMAT,1:NMAT,1:num_faces,1:NMAT,1:NMAT,1:num_links)

! d/dA_{ll,ii,jj} (Ta)_{ji} (D\Psi)_{s,i,j}=tmpdDdA_eta(ii,jj,ll,a,s)
complex(kind(0d0)) :: tmpdDdA_eta(1:NMAT,1:NMAT,1:num_links,1:dimG,1:num_sites)
! d/dA_{ll,ii,jj} (Ta)_{ji} (D\Psi)_{l,i,j}=tmpdDdA_lambda(ii,jj,ll,a,l)
complex(kind(0d0)) :: tmpdDdA_lambda(1:NMAT,1:NMAT,1:num_links,1:dimG,1:num_links)
! d/dA_{ll,ii,jj} (Ta)_{ji} (D\Psi)_{f,i,j}=tmpdDdA_chi(ii,jj,ll,a,f)
complex(kind(0d0)) :: tmpdDdA_chi(1:NMAT,1:NMAT,1:num_links,1:dimG,1:num_faces)

complex(kind(0d0)) :: Uinv(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: Ufm(1:NMAT,1:NMAT,1:num_faces)
!complex(kind(0d0)) :: diffdiff_Omega(1:NMAT,1:NMAT,1:dimG,1:dimG)
complex(kind(0d0)) :: diffdiff_Omega_lambda(1:NMAT,1:NMAT,1:dimG) 

complex(kind(0d0)) :: trace,tmp
complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT),tmpmat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat3(1:NMAT,1:NMAT),tmpmat4(1:NMAT,1:NMAT)
complex(kind(0d0)) :: comm(1:NMAT,1:NMAT)
complex(kind(0d0)) :: FF_MAT(1:NMAT,1:NMAT) ! part of fermion force

!for test
!complex(kind(0d0)) :: T(1:NMAT,1:NMAT,1:dimG)
!complex(kind(0d0)) :: tmp_diffdiff_Omega(1:NMAT,1:NMAT,1:dimG,1:dimG)
!complex(kind(0d0)) :: tmp_diffdiff_Omega2(1:NMAT,1:NMAT,1:dimG,1:dimG)
!complex(kind(0d0)) :: tmp_diff_Omega(1:NMAT,1:NMAT,1:dimG)
!complex(kind(0d0)) :: tmp_diff_Omega2(1:NMAT,1:NMAT,1:dimG)


integer :: s,l,f,ll
integer :: i,j,k,nl,r,ii,jj
integer :: a,b,c,d,e

type diffdiff_by_linkvals_in_face
  complex(kind(0d0)), allocatable :: val(:,:,:,:,:,:)
end type diffdiff_by_linkvals_in_face

type(diffdiff_by_linkvals_in_face) :: diffdiff_Omega2(1:num_faces)

do f=1,num_faces
  allocate( diffdiff_Omega2(f)%val(1:NMAT,1:NMAT,1:dimG,1:dimG,&
    1:links_in_f(f)%num_, 1:links_in_f(f)%num_) )
enddo

!! preparation
dDdA_vec=(0d0,0d0)
dDdA_eta=(0d0,0d0)
dDdA_lambda=(0d0,0d0)
dDdA_chi=(0d0,0d0)
do s=1,num_sites
  do a=1,dimG
    eta_ele(a,s)=vec(site_index(a,s))
  enddo
  call make_traceless_matrix_from_modes(eta_mat(:,:,s),NMAT,eta_ele(:,s))
enddo
do l=1,num_links
  do a=1,dimG
    lambda_ele(a,l)=vec(link_index(a,l))
  enddo
  call make_traceless_matrix_from_modes(lambda_mat(:,:,l),NMAT,lambda_ele(:,l))
enddo
do f=1,num_faces
  do a=1,dimG
    chi_ele(a,f)=vec(face_index(a,f))
  enddo
  call make_traceless_matrix_from_modes(chi_mat(:,:,f),NMAT,chi_ele(:,f))
enddo

! for test
!call make_SUN_generators(T,NMAT)

!! (1) Dirac from site
!   no contribution

if ( p2==0 ) then
!! (2) Dirac from link 1
do s=1,num_sites
  do k=1,linkorg_to_s(s)%num_
    l=linkorg_to_s(s)%labels_(k)
    ! tmpmat1 = lambda_l.U_l
    call matrix_product(tmpmat1,lambda_mat(:,:,l),UMAT(:,:,l))
    ! tmpmat2 = Ul^{-1}.lambda_l
    call matrix_product(tmpmat2,UMAT(:,:,l),lambda_mat(:,:,l),'C','N')
    do jj=1,NMAT
      do ii=1,NMAT
        do j=1,NMAT
          do i=1,NMAT
            dDdA_eta(i,j,s,ii,jj,l)= dDdA_eta(i,j,s,ii,jj,l) &
              - cmplx(alpha_l(l))*conjg(UMAT(jj,i,l))*tmpmat1(ii,j) &
              + cmplx(alpha_l(l))*tmpmat2(i,jj)*UMAT(ii,j,l)
          enddo
        enddo
      enddo
    enddo
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
          - cmplx(alpha_l(l)) * tmpmat2(ii,j)
        dDdA_lambda(i,k,l,k,jj,l) = dDdA_lambda(i,k,l,k,jj,l) & 
          + cmplx(alpha_l(l)) * tmpmat2(i,jj)
      enddo
    enddo
  enddo
enddo
endif 

if ( p3==0 ) then
!! (3) Dirac from link 2
do l=1,num_links
  s=link_tip(l)
! d/dA_{ll,ii,jj} (D\Psi)_{s,i,j}=dDdA_eta(i,j,s,ii,jj,ll)
! d/dA_{ll,ii,jj} (D\Psi)_{l,i,j}=dDdA_lambda(i,j,l,ii,jj,ll)
! d/dA_{ll,ii,jj} (D\Psi)_{f,i,j}=dDdA_chi(i,j,f,ii,jj,ll)

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
              +(0d0,1d0)*cmplx(alpha_l(l))*tmpmat2(ii,j)
          endif
          if ( j==ii ) then
            dDdA_lambda(i,j,l,ii,jj,l) = dDdA_lambda(i,j,l,ii,jj,l) & 
              +(0d0,1d0)*cmplx(alpha_l(l))*tmpmat3(i,jj)
          endif
          dDdA_lambda(i,j,l,ii,jj,l) = dDdA_lambda(i,j,l,ii,jj,l) & 
            -(0d0,1d0)*cmplx(alpha_l(l))*lambda_mat(i,jj,l)*tmpmat1(ii,j) &
            -(0d0,1d0)*cmplx(alpha_l(l))*lambda_mat(ii,j,l)*tmpmat1(i,jj)
        enddo
      enddo
    enddo
  enddo
enddo
endif

      
! d/dA_{ll,ii,jj} (Ta)_{ji} (D\Psi)_{s,i,j}=tmpdDdA_eta(ii,jj,ll,a,s)
! d/dA_{ll,ii,jj} (Ta)_{ji} (D\Psi)_{l,i,j}=tmpdDdA_lambda(ii,jj,ll,a,l)
! d/dA_{ll,ii,jj} (Ta)_{ji} (D\Psi)_{f,i,j}=tmpdDdA_chi(ii,jj,ll,a,f)


      
!! (4) Dirac from face 1
!   no contribution

if( p5 == 0 ) then 
!! (5) Dirac from face 2
! preparation
do f=1,num_faces
  call Make_face_variable(Uf(:,:,f),f,UMAT) 
  call matrix_power(Ufm(:,:,f),Uf(:,:,f),m_omega)
  !call calc_sinU_and_cosUinv(sinU(:,:,f),cosUinv(:,:,f),Ufm(:,:,f))
enddo


do f=1,num_faces
  do j=1,links_in_f(f)%num_
    ll=links_in_f(f)%link_labels_(j)
    do i=1,links_in_f(f)%num_
      l=links_in_f(f)%link_labels_(i)
      call calc_diffdiff_Omega(&
        diffdiff_Omega2(f)%val(:,:,:,:,i,j),&
        Uf(:,:,f),Ufm(:,:,f),f,l,ll,UMAT)
    enddo
  enddo
enddo

do f=1,num_faces
  do i=1,links_in_f(f)%num_
    l=links_in_f(f)%link_labels_(i)

    diffdiff_Omega_lambda=(0d0,0d0)
    do j=1,links_in_f(f)%num_
      ll=links_in_f(f)%link_labels_(j)

      do c=1,dimG
        diffdiff_Omega_lambda(:,:,:) = diffdiff_Omega_lambda(:,:,:) &
          +diffdiff_Omega2(f)%val(:,:,:,c,i,j) * lambda_ele(c,ll)
       enddo
    enddo

    do b=1,dimG
      do a=1,dimG
        call trace_MTa(trace,diffdiff_Omega_lambda(:,:,b),a,NMAT)
        dDdA_vec(face_index(a,f),b,l)= dDdA_vec(face_index(a,f),b,l) &
          + (0d0,1d0)*dcmplx(alpha_f(f)*beta_f(f))*trace*overall_factor
      enddo
    enddo

  enddo
enddo


do l=1,num_links
  do i=1,face_in_l(l)%num_
    f=face_in_l(l)%label_(i)
    do k=1,links_in_f(f)%num_
      if ( links_in_f(f)%link_labels_(k) == l ) then 
        exit
      endif
    enddo
    do j=1,links_in_f(f)%num_
      ll=links_in_f(f)%link_labels_(j)

      do a=1,dimG
        do b=1,dimG
          tmp=(0d0,0d0)
          do ii=1,NMAT
            do jj=1,NMAT
              tmp=tmp+chi_mat(ii,jj,f)&
                *diffdiff_Omega2(f)%val(jj,ii,b,a,j,k)
            enddo
          enddo
          dDdA_vec(link_index(a,l),b,ll)= dDdA_vec(link_index(a,l),b,ll) &
              -(0d0,1d0)*dcmplx(alpha_f(f)*beta_f(f))*tmp*overall_factor
        enddo
      enddo

    enddo
  enddo
enddo
endif

do ll=1,num_links
  do ii=1,NMAT
    do jj=1,NMAT
      !!!!!!!!!!!
      do s=1,num_sites
        do a=1,dimG
          call trace_MTa(tmpdDdA_eta(ii,jj,ll,a,s),dDdA_eta(:,:,s,ii,jj,ll),a,NMAT)
        enddo
      enddo
      !!!!!!!!!!!
      do l=1,num_links
        do a=1,dimG
          call trace_MTa(tmpdDdA_lambda(ii,jj,ll,a,l),dDdA_lambda(:,:,l,ii,jj,ll),a,NMAT)
        enddo
      enddo
      !!!!!!!!!!!
    enddo
  enddo
enddo
do s=1,num_sites
  do a=1,dimG
    do ll=1,num_links
      do b=1,dimG
        call trace_MTa(tmp,tmpdDdA_eta(:,:,ll,a,s),b,NMAT)
        dDdA_vec(site_index(a,s),b,ll) &
          = dDdA_vec(site_index(a,s),b,ll) + overall_factor * tmp
      enddo
    enddo
  enddo
enddo
do l=1,num_links
  do a=1,dimG
    do ll=1,num_links
      do b=1,dimG
        call trace_MTa(tmp,tmpdDdA_lambda(:,:,ll,a,l),b,NMAT)
        dDdA_vec(link_index(a,l),b,ll) &
          = dDdA_vec(link_index(a,l),b,ll) + overall_factor * tmp
      enddo
    enddo
  enddo
enddo

!dDdA_chi=dDdA_chi*overall_factor

!! (T1) test action
!do f=1,num_faces
!  call Make_face_variable(Uf(:,:,f),f,UMAT) 
!  call matrix_power(NMAT,Uf(:,:,f),m_omega,Ufm(:,:,f))
!enddo
!
!do f=1,num_faces
!  do i=1,links_in_f(f)%num_
!    l=links_in_f(f)%link_labels_(i)
!
!    diffdiff_Omega_lambda=(0d0,0d0)
!    do j=1,links_in_f(f)%num_
!      ll=links_in_f(f)%link_labels_(j)
!
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      !! change here
!      !call calc_diffdiff_Uf(tmp_diffdiff_Omega,Uf(:,:,f),f,l,ll,UMAT) ! OK
!      call calc_diffdiff_Ufm(tmp_diffdiff_Omega,Uf(:,:,f),f,l,ll,UMAT) ! OK
!      !call calc_diffdiff_SandC(tmp_diffdiff_Omega,tmp_diffdiff_Omega2,Uf(:,:,f),f,l,ll,UMAT) 
!      !diffdiff_Omega=tmp_diffdiff_Omega2
!      do a=1,dimG
!        do b=1,dimG
!          do ii=1,NMAT
!            do jj=1,NMAT
!              diffdiff_Omega(ii,jj,a,b) &
!                =(tmp_diffdiff_Omega(ii,jj,a,b)-conjg(tmp_diffdiff_Omega(jj,ii,a,b)))
!            enddo
!          enddo
!        enddo
!      enddo
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      do c=1,dimG
!        diffdiff_Omega_lambda(:,:,:) = diffdiff_Omega_lambda(:,:,:) &
!          +diffdiff_Omega(:,:,:,c) * lambda_ele(c,ll)
!      enddo
!    enddo
!
!    do b=1,dimG
!      do a=1,dimG
!        call trace_MTa(trace,diffdiff_Omega_lambda(:,:,b),a,NMAT)
!        dDdA_chi(face_index(a,f),b,l) = dDdA_chi(face_index(a,f),b,l) &
!          + trace
!      enddo
!    enddo
!
!  enddo
!enddo
!
!
!do l=1,num_links
!  do i=1,face_in_l(l)%num_
!    f=face_in_l(l)%label_(i)
!    do j=1,links_in_f(f)%num_
!      ll=links_in_f(f)%link_labels_(j)
!
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      !! change here
!      !call calc_diffdiff_Uf(tmp_diffdiff_Omega,Uf(:,:,f),f,ll,l,UMAT) ! OK 
!      call calc_diffdiff_Ufm(tmp_diffdiff_Omega,Uf(:,:,f),f,ll,l,UMAT) 
!      !call calc_diffdiff_SandC(tmp_diffdiff_Omega,tmp_diffdiff_Omega2,Uf(:,:,f),f,ll,l,UMAT) 
!      !diffdiff_Omega=tmp_diffdiff_Omega2
!
!      do a=1,dimG
!        do b=1,dimG
!          do ii=1,NMAT
!            do jj=1,NMAT
!              diffdiff_Omega(ii,jj,a,b) &
!                =(tmp_diffdiff_Omega(ii,jj,a,b)-conjg(tmp_diffdiff_Omega(jj,ii,a,b)))
!            enddo
!          enddo
!        enddo
!      enddo
!
!      !call calc_diff_Ufm(tmp_diff_Omega,Uf(:,:,f),f,ll,UMAT) 
!      !call calc_diff_Ufm(tmp_diff_Omega2,Uf(:,:,f),f,l,UMAT) 
!      !write(*,*) "====",ll,l,"====="
!      !do a=1,dimG
!      !  do b=1,dimG
!      !    call ZGEMM('N','C',NMAT,NMAT,NMAT,(1d0,0d0), &
!      !      Ufm(:,:,f), NMAT, &
!      !      tmp_diffdiff_Omega(:,:,a,b), NMAT, &
!      !      (0d0,0d0), tmpmat1, NMAT)
!      !    call ZGEMM('N','C',NMAT,NMAT,NMAT,(1d0,0d0), &
!      !      tmp_diffdiff_Omega(:,:,a,b), NMAT, &
!      !      Ufm(:,:,f), NMAT, &
!      !      (1d0,0d0), tmpmat1, NMAT)
!      !    call ZGEMM('N','C',NMAT,NMAT,NMAT,(1d0,0d0), &
!      !      tmp_diff_Omega(:,:,a), NMAT, &
!      !      tmp_diff_Omega2(:,:,b), NMAT, &
!      !      (1d0,0d0), tmpmat1, NMAT)
!      !    call ZGEMM('N','C',NMAT,NMAT,NMAT,(1d0,0d0), &
!      !      tmp_diff_Omega2(:,:,b), NMAT, &
!      !      tmp_diff_Omega(:,:,a), NMAT, &
!      !      (1d0,0d0), tmpmat1, NMAT)
!      !    tmp=(0d0,0d0)
!      !    do ii=1,NMAT
!      !      do jj=1,NMAT
!      !        tmp=tmp+( tmpmat1(ii,jj)*dconjg(tmpmat1(ii,jj)) ) 
!      !      enddo
!      !    enddo
!      !    write(*,*) a,b,dble(tmp)
!      !  enddo
!      !enddo
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      
!      do a=1,dimG
!        do b=1,dimG
!          tmp=(0d0,0d0)
!          do ii=1,NMAT
!            do jj=1,NMAT
!              tmp=tmp+chi_mat(ii,jj,f)*diffdiff_Omega(jj,ii,b,a)
!            enddo
!          enddo
!          dDdA_chi(link_index(a,l),b,ll) = dDdA_chi(link_index(a,l),b,ll)-tmp
!        enddo
!      enddo
!
!    enddo
!  enddo
!enddo

end subroutine prod_dDdA



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to return product of d/dA D . vec
subroutine calc_dCosUinvdA_dSinUdA(dCosUinvdA,dSinUdA,&
    cosUinv,sinU,Uf,UMAT,f,ll_label)
use Dirac_operator, only : calc_Amat, calc_Bmat
implicit none

! for given f and ll
! d cosUinv(i,j,f) / dA_{ii,jj,ll)
complex(kind(0d0)), intent(out) :: dCosUinvdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
complex(kind(0d0)), intent(out) :: dSinUdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: sinU(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: cosUinv(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: UMat(1:NMAT,1:NMAT,1:num_links)
integer, intent(in) :: f,ll_label

complex(kind(0d0)) :: Amat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Bmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: UinvA(1:NMAT,1:NMAT), BUinv(1:NMAT,1:NMAT)

integer :: k,i,j,ii,jj

dCosUinvdA=(0d0,0d0)
do k=1,m_omega
  call calc_Amat(Amat,f,ll_label,k,Uf,Umat)
  call calc_Bmat(Bmat,f,ll_label,k,Uf,Umat)

  call matrix_product(UinvA,cosUinv,Amat)
  call matrix_product(BUinv,Bmat,cosUinv)

  do jj=1,NMAT
    do ii=1,NMAT
      do j=1,NMAT
        do i=1,NMAT
          dCosUinvdA(i,j,ii,jj)=dCosUinvdA(i,j,ii,jj) &
            + UinvA(i,jj)*BUinv(ii,j) &
            - conjg( BUinv(jj,i) ) * conjg( UinvA(j,ii) )
          dSinUdA(i,j,ii,jj)=dSinUdA(i,j,ii,jj) &
            + Amat(i,jj)*Bmat(ii,j) &
            + conjg(Bmat(jj,i))*conjg(Amat(j,ii))
        enddo
      enddo
    enddo
  enddo
enddo
dCosUinvdA=dCosUinvdA*(-links_in_f(f)%link_dirs_(ll_label))*(0d0,1d0)
dSinUdA=dSinUdA * links_in_f(f)%link_dirs_(ll_label) * (0d0,1d0)

end subroutine calc_dCosUinvdA_dSinUdA



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to return product of d/dA D . vec
subroutine calc_dABmatdA(dAmatdA,dBmatdA,Uf,UMAT,ll_label,f,l_label,k)
implicit none

! for given f,l,k and ll
! d Amat(i,j,f,l,k) / dA_{ii,jj,ll)
complex(kind(0d0)), intent(out) :: dAmatdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
complex(kind(0d0)), intent(out) :: dBmatdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: UMat(1:NMAT,1:NMAT,1:num_links)
integer, intent(in) :: f,ll_label,l_label,k

complex(kind(0d0)) :: UU_1_to_ll(1:NMAT,1:NMAT) ! 1..ll
complex(kind(0d0)) :: UU_ll_to_n(1:NMAT,1:NMAT) ! ll..n
complex(kind(0d0)) :: UU_ll_to_l(1:NMAT,1:NMAT) ! ll..l
complex(kind(0d0)) :: UU_1_to_l(1:NMAT,1:NMAT) ! 1..l
complex(kind(0d0)) :: UU_l_to_n(1:NMAT,1:NMAT) ! l..n
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
if ( links_in_f(f)%link_labels_(ll_label) == 1 ) then
  label2 = ll_label - 1
else
  label2 = ll_label
endif 
call calc_prodUl_from_n1_to_n2_in_Uf(UU_1_to_ll,f,1,label2,UMAT)
!!!!!!!!!!!!!!
if ( links_in_f(f)%link_labels_(ll_label) == 1 ) then
  label1 = ll_label 
else
  label1 = ll_label+1
endif 
call calc_prodUl_from_n1_to_n2_in_Uf(UU_ll_to_n,f,label1,n,UMAT)
!!!!!!!!!!!!!!
if ( links_in_f(f)%link_labels_(ll_label) == 1 ) then
  label1 = ll_label 
else
  label1 = ll_label+1
endif 
if ( links_in_f(f)%link_labels_(l_label) == 1 ) then
  label2 = l_label-1
else
  label2 = l_label
endif 
call calc_prodUl_from_n1_to_n2_in_Uf(UU_ll_to_l,f,label1,label2,UMAT)
!!!!!!!!!!!!!!
if ( links_in_f(f)%link_labels_(l_label) == 1 ) then
  label2 = l_label-1
else
  label2 = l_label
endif 
call calc_prodUl_from_n1_to_n2_in_Uf(UU_1_to_l,f,1,label2,UMAT)
!!!!!!!!!!!!!!
if ( links_in_f(f)%link_labels_(l_label) == 1 ) then
  label1 = l_label
else
  label1 = l_label+1
endif 
call calc_prodUl_from_n1_to_n2_in_Uf(UU_l_to_n,f,label1,n,UMAT)
!!!!!!!!!!!!!!
if ( links_in_f(f)%link_labels_(l_label) == 1 ) then
  label1 = l_label 
else
  label1 = l_label+1
endif 
call calc_prodUl_from_n1_to_n2_in_Uf(UU_l_to_n,f,label1,n,UMAT)
!!!!!!!!!!!!!!
if ( links_in_f(f)%link_labels_(l_label) == 1 ) then
  label1 = l_label 
else
  label1 = l_label+1
endif 
if ( links_in_f(f)%link_labels_(ll_label) == 1 ) then
  label2 = ll_label-1
else
  label2 = ll_label
endif 
call calc_prodUl_from_n1_to_n2_in_Uf(UU_l_to_ll,f,label1,label2,UMAT)

        
!!!!!!!!!!!!!! dAmatdA !!!!!!!!!!!!!!!!!
if ( k > 1 ) then
  do kk=1,k-1
    call matrix_power(tmpmat2,Uf,kk-1)
    call matrix_product(tmpmat1,tmpmat2,UU_1_to_ll)

    call matrix_power(tmpmat2,Uf,k-kk-1)
    call matrix_product(tmpmat3,UU_ll_to_n,tmpmat2)
    call matrix_product(tmpmat2,tmpmat3,UU_1_to_l)

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
  call matrix_power(tmpmat2,Uf,k-1)
  call matrix_product(tmpmat1,tmpmat2,UU_1_to_ll)

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
    call matrix_power(tmpmat2,Uf,k-1)
    call matrix_product(tmpmat1, tmpmat2,UU_1_to_l)
    do jj=1,NMAT
      do ii=1,NMAT
        j=ii
        do i=1,NMAT
          dAmatdA(i,j,ii,jj)=dAmatdA(i,j,ii,jj) &
            + tmpmat1(i,jj) * epsilon_r
        enddo
      enddo
    enddo
  endif
endif

!!!!!!!!!!! dBmatdA !!!!!!!!!!!!!!!
if (m_omega-k > 0) then
  do kk=1,m_omega-k
    call matrix_power(tmpmat1, Uf, kk-1)
    call matrix_product(tmpmat2,UU_l_to_n,tmpmat1)
    call matrix_product(tmpmat1,tmpmat2,UU_1_to_ll)

    call matrix_power(tmpmat3,Uf,m_omega-k-kk)
    call matrix_product(tmpmat2,UU_ll_to_n,tmpmat3)

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
  call matrix_power(tmpmat1,Uf,m_omega-k)
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
    call matrix_power(tmpmat2,Uf,m_omega-k)
    call matrix_product(tmpmat1,UU_ll_to_n,tmpmat2)
    do jj=1,NMAT
      do ii=1,NMAT
        i=jj
        do j=1,NMAT
          dBmatdA(i,j,ii,jj)=dBmatdA(i,j,ii,jj) &
            + tmpmat1(ii,j) * epsilon_r
        enddo
      enddo
    enddo
  endif
endif

end subroutine calc_dABmatdA


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

