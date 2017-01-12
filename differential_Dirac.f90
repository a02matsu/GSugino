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
subroutine prod_dDdPhi(dDdPhi_eta,dDdPhi_chi,eta_mat,chi_mat,UMAT)
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links) 
complex(kind(0d0)), intent(in) :: eta_mat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(in) :: chi_mat(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)), intent(out) :: dDdPhi_eta(1:NMAT,1:NMAT,1:num_sites,1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(out) :: dDdPhi_chi(1:NMAT,1:NMAT,1:num_faces,1:NMAT,1:NMAT,1:num_sites)


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
    !write(*,*) p1
!! (1) Dirac from site
do s=1,num_sites
  do jj=1,NMAT
    i=jj
    do ii=1,NMAT
      do j=1,NMAT 
        dDdPhi_eta(i,j,s,ii,jj,s)= dDdPhi_eta(i,j,s,ii,jj,s) &
          + (-0.5d0,0d0)*cmplx( alpha_s(s)*overall_factor ) &
            * eta_mat(ii,j,s)
      enddo
    enddo
    !!
    do ii=1,NMAT
      j=ii
      do i=1,NMAT
        dDdPhi_eta(i,j,s,ii,jj,s)= dDdPhi_eta(i,j,s,ii,jj,s) &
          - (-0.5d0,0d0)*cmplx( alpha_s(s)*overall_factor ) &
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
          + (-2d0,0d0)*cmplx( alpha_f(f)*overall_factor ) &
            * chi_mat(ii,j,f)
      enddo
    enddo
    !!
    do ii=1,NMAT
      j=ii
      do i=1,NMAT
        dDdPhi_chi(i,j,f,ii,jj,s)= dDdPhi_chi(i,j,f,ii,jj,s) &
          - (-2d0,0d0)*cmplx( alpha_f(f)*overall_factor ) &
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

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links) 
complex(kind(0d0)), intent(in) :: lambda_mat(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(out) :: dDdbPhi_lambda(1:NMAT,1:NMAT,1:num_links,1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: dDdbPhi_vec(1:sizeD,1:dimG,1:num_sites)


complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmp

integer :: s,l
integer :: i,j,ii,jj
integer :: a,b

!! preparation
dDdbPhi_vec=(0d0,0d0)
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
             + cmplx( alpha_l(l) * overall_factor ) &
               * ( UMAT(i,jj,l) * tmpmat1(ii,j) &
                   - tmpmat2(i,jj) * conjg( UMAT(j,ii,l) ))
          if ( i == jj ) then
            dDdbPhi_lambda(i,j,l,ii,jj,link_org(l)) &
             = dDdbPhi_lambda(i,j,l,ii,jj,link_org(l)) &
               + cmplx( alpha_l(l) * overall_factor ) &
                 * lambda_mat(ii,j,l)
          endif
          if ( ii == j ) then
            dDdbPhi_lambda(i,j,l,ii,jj,link_org(l)) &
             = dDdbPhi_lambda(i,j,l,ii,jj,link_org(l)) &
               - cmplx( alpha_l(l) * overall_factor ) &
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
complex(kind(0d0)), intent(inout) :: dDdA_eta(1:NMAT,1:NMAT,1:num_sites,1:NMAT,1:NMAT,1:num_links)
! d/dA_{ll,ii,jj} (D\Psi)_{l,i,j}=dDdA_lambda(i,j,l,ii,jj,ll)
complex(kind(0d0)), intent(inout) :: dDdA_lambda(1:NMAT,1:NMAT,1:num_links,1:NMAT,1:NMAT,1:num_links)
! d/dA_{ll,ii,jj} (D\Psi)_{f,i,j}=dDdA_chi(i,j,f,ii,jj,ll)
complex(kind(0d0)),intent(inout) :: dDdA_chi(1:NMAT,1:NMAT,1:num_faces,1:NMAT,1:NMAT,1:num_links)

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links) 
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(in) :: eta_mat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(in) :: lambda_mat(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: chi_mat(1:NMAT,1:NMAT,1:num_faces)

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

!type diffdiff_by_linkvals_in_face
  !complex(kind(0d0)), allocatable :: val(:,:,:,:,:,:)
!end type diffdiff_by_linkvals_in_face

!type(diffdiff_by_linkvals_in_face) :: diffdiff_Omega2(1:num_faces)

!do f=1,num_faces
  !allocate( diffdiff_Omega2(f)%val(1:NMAT,1:NMAT,1:dimG,1:dimG,&
    !1:links_in_f(f)%num_, 1:links_in_f(f)%num_) )
!enddo

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

!! (4) Dirac from face 1
!   no contribution

if( p5 == 0 ) then 
!! (5) Dirac from face 2
! preparation
  call calc_fermion_force_from_omega&
      (dDdA_lambda,dDdA_chi,UMAT,lambda_mat,chi_mat)
endif 
dDdA_eta=dDdA_eta*cmplx(overall_factor)
dDdA_lambda=dDdA_lambda*cmplx(overall_factor)
dDdA_chi=dDdA_chi*cmplx(overall_factor)
!do l=1,num_links
!  do j=1,NMAT
!    do i=1,NMAT
!      do s=1,num_sites
!        tmp=(0d0,0d0)
!        do ii=1,NMAT
!          tmp=tmp+dDdA_eta(ii,ii,s,i,j,l)
!        enddo
!        do ii=1,NMAT
!          dDdA_eta(ii,ii,s,i,j,l)=dDdA_eta(ii,ii,s,i,j,l)-tmp/cmplx(dble(NMAT))
!        enddo
!      enddo
!      do ll=1,num_links
!        tmp=(0d0,0d0)
!        do ii=1,NMAT
!          tmp=tmp+dDdA_lambda(ii,ii,ll,i,j,l)
!        enddo
!        do ii=1,NMAT
!          dDdA_lambda(ii,ii,ll,i,j,l)=dDdA_lambda(ii,ii,ll,i,j,l)-tmp/cmplx(dble(NMAT))
!        enddo
!      enddo
!      do f=1,num_faces
!        tmp=(0d0,0d0)
!        do ii=1,NMAT
!          tmp=tmp+dDdA_chi(ii,ii,f,i,j,l)
!        enddo
!        do ii=1,NMAT
!          dDdA_chi(ii,ii,f,i,j,l)=dDdA_chi(ii,ii,f,i,j,l)-tmp/cmplx(dble(NMAT))
!        enddo
!      enddo
!    enddo
!  enddo
!enddo


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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_fermion_force_from_omega&
    (dDdA_lambda,dDdA_chi,UMAT,lambda_mat,chi_mat)
use Dirac_operator, only : calc_sinu_and_CosUinv, calc_Amat, Calc_Bmat
implicit none 
complex(kind(0d0)), intent(inout) :: dDdA_lambda&
  (1:NMAT,1:NMAT,1:num_links, 1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(inout) :: dDdA_chi&
  (1:NMAT,1:NMAT,1:num_faces, 1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: lambda_mat(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: chi_mat(1:NMAT,1:NMAT,1:num_faces)

complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: Ufm(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: sinU(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: cosUinv(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: Amat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Bmat(1:NMAT,1:NMAT)

complex(kind(0d0)) :: dCosUinvdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
complex(kind(0d0)) :: dSinUdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
complex(kind(0d0)) :: dAmatdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
complex(kind(0d0)) :: dBmatdA(1:NMAT,1:NMAT,1:NMAT,1:NMAT)

complex(kind(0d0)) :: prodmat1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: prodmat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: globalmat1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: globalmat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: globalmat3(1:NMAT,1:NMAT)
complex(kind(0d0)) :: globalmat4(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat3(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat4(1:NMAT,1:NMAT)
complex(kind(0d0)) :: line(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
integer :: l_label,ll_label,k,l,ll
integer :: i,j,ii,jj,f
complex(kind(0d0)) :: trace

! label of face
do f=1,num_faces
  call Make_face_variable(Uf(:,:,f),f,UMAT) 
  if (m_omega .ne. 0) then 
    call matrix_power(Ufm(:,:,f),Uf(:,:,f),m_omega)
    call calc_sinU_and_cosUinv(sinU(:,:,f),cosUinv(:,:,f),Ufm(:,:,f))
  endif
enddo

! label of face
do f=1,num_faces
! label of link fermion \lambda_l
do l_label=1,links_in_f(f)%num_
  l=links_in_f(f)%link_labels_(l_label)

  ! label of link to differentiate
  do ll_label=1,links_in_f(f)%num_
    ll=links_in_f(f)%link_labels_(ll_label)
!subroutine calc_dCosUinvdA_dSinUdA(dCosUinvdA,dSinUdA,&
!    cosUinv,sinU,Uf,UMAT,f,ll_label)

    line=(0d0,0d0)
    if( m_omega == 0 ) then
      call calc_dABmatdA(dAmatdA,dBmatdA,&
        Uf(:,:,f),UMAT,ll_label,f,l_label,1)
      call calc_Amat(Amat,f,l_label,1,Uf(:,:,f),UMAT)
      call calc_Bmat(Bmat,f,l_label,1,Uf(:,:,f),UMAT)

      
      do jj=1,NMAT
      do ii=1,NMAT
        call matrix_3_product(prodmat1, &
          dAmatdA(:,:,ii,jj),lambda_mat(:,:,l),BMAT)
        call matrix_product(tmpmat1,Amat,lambda_mat(:,:,l))
        call zgemm('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
          tmpmat1, NMAT, &
          dBmatdA(:,:,ii,jj), NMAT, &
          (1d0,0d0), line(:,:,ii,jj), NMAT)
      enddo
      enddo

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! make traceless for (ii,jj)
      tmpmat1=(0d0,0d0)
      do ii=1,NMAT
        tmpmat1=tmpmat1+line(:,:,ii,ii)
      enddo
      do ii=1,NMAT
        line(:,:,ii,ii)=line(:,:,ii,ii)-tmpmat1/cmplx(dble(NMAT))
      enddo
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do jj=1,NMAT
        do ii=1,NMAT
          dDdA_chi(:,:,f,ii,jj,ll)=dDdA_chi(:,:,f,ii,jj,ll) &
            + cmplx(dble(links_in_f(f)%link_dirs_(l_label) ))*(0d0,1d0) &
             * cmplx( alpha_f(f)*beta_f(f) ) &
             * line(:,:,ii,jj)
        enddo
      enddo
    else
      call calc_dCosUinvdA_dSinUdA(dCosUinvdA,dSinUdA,&
        cosUinv(:,:,f),Uf(:,:,f),UMAT,f,ll_label)
      do k=1,m_omega
        call calc_dABmatdA(dAmatdA,dBmatdA,&
          Uf(:,:,f),UMAT,ll_label,f,l_label,k)
        call calc_Amat(Amat,f,l_label,k,Uf(:,:,f),UMAT)
        call calc_Bmat(Bmat,f,l_label,k,Uf(:,:,f),UMAT)
  
        ! globalmat1 = A lambda B
        call matrix_3_product(globalmat1,Amat,lambda_mat(:,:,l),Bmat)
  
        ! globalmat2 = B^dag lambda A^dag
        call matrix_3_product(globalmat2,Bmat,lambda_mat(:,:,l),Amat,'C','N','C')
  
        do jj=1,NMAT
        do ii=1,NMAT
          ! prodat1 = ¥delta_A lambda B + A lambda ¥delta_B
          call matrix_3_product(prodmat1, &
            dAmatdA(:,:,ii,jj),lambda_mat(:,:,l),BMAT)
          call matrix_product(tmpmat1,Amat,lambda_mat(:,:,l))
          call zgemm('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat1, NMAT, &
            dBmatdA(:,:,ii,jj), NMAT, &
            (1d0,0d0), prodmat1, NMAT)
  
          ! prodmat2 = ¥delta_B^dag lambda A^dag + B^dag lambda ¥delta_A^dag
          call matrix_3_product(prodmat2,&
            dBmatdA(:,:,jj,ii),lambda_mat(:,:,l),Amat,'C','N','C') 
          call matrix_product(tmpmat1,Bmat,lambda_mat(:,:,l),'C','N')
          call zgemm('N','C',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat1, NMAT, &
            dAmatdA(:,:,jj,ii), NMAT, & 
            (1d0,0d0), prodmat2, NMAT)
  
          ! line1
          call matrix_anticommutator(tmpmat1, prodmat1+prodmat2, CosUinv(:,:,f))
          line(:,:,ii,jj)=line(:,:,ii,jj)+tmpmat1
          ! line2
          call matrix_anticommutator(tmpmat1,globalmat1+globalmat2, dCosUinvdA(:,:,ii,jj))
          line(:,:,ii,jj) = line(:,:,ii,jj) + tmpmat1
          ! line3
          call matrix_product(tmpmat1,CosUinv(:,:,f),globalmat1-globalmat2)
          call matrix_product(tmpmat2,tmpmat1,CosUinv(:,:,f)) 
          call matrix_anticommutator(tmpmat3,dSinUdA(:,:,ii,jj),tmpmat2)
          line(:,:,ii,jj) = line(:,:,ii,jj) - tmpmat3
          ! line4
           ! 4-1 collect in tmpmat2
          call matrix_product(tmpmat3,dCosUinvdA(:,:,ii,jj),globalmat1-globalmat2)
          call matrix_product(tmpmat2,tmpmat3,CosUinv(:,:,f))
           ! 4-2 tmpmat1=CosUinv.(globalmat1-globalmat2) here
          call zgemm('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat1, NMAT, &
            dCosUinvdA(:,:,ii,jj), NMAT, &
            (1d0,0d0), tmpmat2, NMAT)
           ! 4-3
          call matrix_product(tmpmat3,CosUinv(:,:,f),prodmat1-prodmat2)
          call zgemm('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat3, NMAT, &
            CosUinv(:,:,f), NMAT, &
            (1d0,0d0), tmpmat2, NMAT)
           ! take anti-commutator
          call matrix_anticommutator(tmpmat1,tmpmat2,SinU(:,:,f))
          line(:,:,ii,jj) = line(:,:,ii,jj) - tmpmat1
  
        enddo ! ii
        enddo ! jj
      enddo ! k
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! make traceless for (ii,jj)
      tmpmat1=(0d0,0d0)
      do ii=1,NMAT
        tmpmat1=tmpmat1+line(:,:,ii,ii)
      enddo
      do ii=1,NMAT
        line(:,:,ii,ii)=line(:,:,ii,ii)-tmpmat1/cmplx(dble(NMAT))
      enddo
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do jj=1,NMAT
        do ii=1,NMAT
          dDdA_chi(:,:,f,ii,jj,ll)=dDdA_chi(:,:,f,ii,jj,ll) &
            + cmplx(dble(links_in_f(f)%link_dirs_(l_label) ))*(0d0,1d0) &
             * cmplx( alpha_f(f)*beta_f(f)/dble(m_omega) ) &
             * line(:,:,ii,jj)
        enddo
      enddo 
    endif 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! make traceless for (i,j)
    trace=(0d0,0d0)
    do j=1,NMAT
      trace=trace+dDdA_chi(j,j,f,ii,jj,ll)
    enddo
    do j=1,NMAT
      dDdA_chi(j,j,f,ii,jj,ll)=dDdA_chi(j,j,f,ii,jj,ll)-trace/cmplx(dble(NMAT))
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  enddo ! ll
enddo ! l
enddo ! f

do l=1,num_links
do i=1,face_in_l(l)%num_       
  f=face_in_l(l)%label_(i)
  ! j: position of the link l in the face l
  do l_label=1,links_in_f(f)%num_
    if ( l == links_in_f(f)%link_labels_(l_label) ) exit
  enddo

  do ll_label=1,links_in_f(f)%num_
    ll=links_in_f(f)%link_labels_(ll_label)

    if (m_omega == 0) then 
      call calc_dABmatdA(dAmatdA,dBmatdA,&
        Uf(:,:,f),UMAT,ll_label,f,l_label,1)
      call calc_Amat(Amat,f,l_label,1,Uf(:,:,f),UMAT)
      call calc_Bmat(Bmat,f,l_label,1,Uf(:,:,f),UMAT)

      do jj=1,NMAT
        do ii=1,NMAT
          call matrix_3_product(tmpmat1,dBmatdA(:,:,ii,jj),chi_mat(:,:,f),Amat)
          call matrix_3_product(tmpmat2,Bmat,chi_mat(:,:,f),dAmatdA(:,:,ii,jj))
          call matrix_3_product(tmpmat3,dAmatdA(:,:,jj,ii),chi_mat(:,:,f),Bmat,'C','M','C')
          call matrix_3_product(tmpmat4,Amat,chi_mat(:,:,f),dBmatdA(:,:,jj,ii), 'C','N','C')
          line(:,:,ii,jj)=tmpmat1+tmpmat2+tmpmat3+tmpmat4
        enddo
      enddo
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! make traceless for (ii,jj)
      tmpmat1=(0d0,0d0)
      do ii=1,NMAT
        tmpmat1=tmpmat1+line(:,:,ii,ii)
      enddo
      do ii=1,NMAT
        line(:,:,ii,ii)=line(:,:,ii,ii)-tmpmat1/cmplx(dble(NMAT))
      enddo
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
      do jj=1,NMAT
        do ii=1,NMAT
           dDdA_lambda(:,:,l,ii,jj,ll)=dDdA_lambda(:,:,l,ii,jj,ll) &
              - cmplx(dble(links_in_f(f)%link_dirs_(l_label) ))*(0d0,1d0) &
               * cmplx( alpha_f(f)*beta_f(f) ) &
               * line(:,:,ii,jj)
        enddo
      enddo
    else
      call calc_dCosUinvdA_dSinUdA(dCosUinvdA,dSinUdA,&
        cosUinv(:,:,f),Uf(:,:,f),UMAT,f,ll_label)
  
      globalmat3=(0d0,0d0)
      globalmat4=(0d0,0d0)
      ! globalmat3 = { CosUinv, chi_f }
      call matrix_anticommutator(globalmat3,CosUinv(:,:,f),chi_mat(:,:,f))
      ! globalmat4 = CosUinv.{ SinU, chi_f }.CosUinv
      call matrix_anticommutator(tmpmat1,SinU(:,:,f),chi_mat(:,:,f))
      call matrix_3_product(globalmat4,CosUinv(:,:,f),tmpmat1,CosUinv(:,:,f))
  
      line=(0d0,0d0)
      do k=1,m_omega
        call calc_dABmatdA(dAmatdA,dBmatdA,&
          Uf(:,:,f),UMAT,ll_label,f,l_label,k)
        call calc_Amat(Amat,f,l_label,k,Uf(:,:,f),UMAT)
        call calc_Bmat(Bmat,f,l_label,k,Uf(:,:,f),UMAT)
  
        do jj=1,NMAT
        do ii=1,NMAT
          ! prodmat1
          call matrix_anticommutator(prodmat1,dCosUinvdA(:,:,ii,jj),chi_mat(:,:,f))
          ! prodmat2
          ! 1st term
          call matrix_anticommutator(tmpmat1,dSinUdA(:,:,ii,jj),chi_mat(:,:,f))
          call matrix_3_product(prodmat2,CosUinv(:,:,f),tmpmat1,CosUinv(:,:,f))
          ! 2nd term
          call matrix_anticommutator(tmpmat1,SinU(:,:,f),chi_mat(:,:,f))
          call matrix_product(tmpmat2,dCosUinvdA(:,:,ii,jj),tmpmat1)
          call zgemm('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat2, NMAT, &
            CosUinv(:,:,f),NMAT, &
            (1d0,0d0), prodmat2, NMAT)
          ! 3rd term tmpmat1={ SinU, chi_mat } now
          call matrix_product(tmpmat2,CosUinv(:,:,f),tmpmat1)
          call zgemm('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat2, NMAT, &
            dCosUinvdA(:,:,ii,jj),NMAT, &
            (1d0,0d0), prodmat2, NMAT)
  
          ! line1
          call matrix_3_product(tmpmat1,&
            dBmatdA(:,:,ii,jj),globalmat3-globalmat4,Amat)
          line(:,:,ii,jj)=line(:,:,ii,jj)+tmpmat1
          ! line2
          call matrix_product(tmpmat1,Bmat,globalmat3-globalmat4)
          call zgemm('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat1, NMAT, &
            dAmatdA(:,:,ii,jj),NMAT, &
            (1d0,0d0), line(:,:,ii,jj), NMAT)
          ! line3
          call matrix_product(tmpmat1,dAmatdA(:,:,jj,ii),globalmat3+globalmat4,'C','N') 
          call zgemm('N','C',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat1, NMAT, &
            Bmat,NMAT, &
            (1d0,0d0), line(:,:,ii,jj), NMAT)
          ! line4
          call matrix_product(tmpmat1,Amat,globalmat3+globalmat4,'C','N')
          call zgemm('N','C',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat1, NMAT, &
            dBmatdA(:,:,jj,ii),NMAT, & 
            (1d0,0d0), line(:,:,ii,jj), NMAT)
          ! line5
          call matrix_product(tmpmat1,Bmat,prodmat1-prodmat2)
          call zgemm('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat1, NMAT, &
            Amat,NMAT, &
            (1d0,0d0), line(:,:,ii,jj), NMAT)
          ! line6
          call matrix_product(tmpmat1,Amat,prodmat1+prodmat2,'C','N')
          call zgemm('N','C',NMAT,NMAT,NMAT,(1d0,0d0), &
            tmpmat1, NMAT, &
            Bmat,NMAT, &
            (1d0,0d0), line(:,:,ii,jj), NMAT)
        enddo ! ii
        enddo ! jj
      enddo ! k
  
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! make traceless for (ii,jj)
      tmpmat1=(0d0,0d0)
      do ii=1,NMAT
        tmpmat1=tmpmat1+line(:,:,ii,ii)
      enddo
      do ii=1,NMAT
        line(:,:,ii,ii)=line(:,:,ii,ii)-tmpmat1/cmplx(dble(NMAT))
      enddo
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
      do jj=1,NMAT
        do ii=1,NMAT
           dDdA_lambda(:,:,l,ii,jj,ll)=dDdA_lambda(:,:,l,ii,jj,ll) &
              - cmplx(dble(links_in_f(f)%link_dirs_(l_label) ))*(0d0,1d0) &
               * cmplx( alpha_f(f)*beta_f(f)/dble(m_omega) ) &
               * line(:,:,ii,jj)
        enddo
      enddo
    endif 
    do jj=1,NMAT
      do ii=1,NMAT
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! make traceless for (i,j)
        trace=(0d0,0d0)
        do j=1,NMAT
          trace=trace+dDdA_lambda(j,j,l,ii,jj,ll)
        enddo
        do j=1,NMAT
          dDdA_lambda(j,j,l,ii,jj,ll)=dDdA_lambda(j,j,l,ii,jj,ll)-trace/cmplx(dble(NMAT))
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      enddo
    enddo
  enddo ! ll
enddo ! f
enddo ! l

end subroutine calc_fermion_force_from_omega






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

