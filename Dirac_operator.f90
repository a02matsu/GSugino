!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! module for fermion matrix
module Dirac_operator
use global_parameters
use global_subroutines
use SUN_generators
use matrix_functions
#ifdef PARALLEL
use parallel
#endif
implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compute D.vec
!!  S_f = 1/2 \Psi^T D \Psi
subroutine Prod_Dirac(&
    DF_eta, DF_lambda, DF_chi, &
    eta_mat,lambda_mat,chi_mat, &
    UMAT,PhiMat)
use matrix_functions, only : matrix_3_product
#ifdef PARALLEL
use parallel
#endif
implicit none

!integer, intent(in) :: Vsize !! vecsize must be sizeD
!complex(kind(0d0)), intent(in) :: vec(1:Vsize)
!complex(kind(0d0)), intent(inout) :: D_vec(1:Vsize)
complex(kind(0d0)),intent(in) :: eta_mat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)),intent(in) :: lambda_mat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)),intent(in) :: chi_mat(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)), intent(out) :: DF_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(out) :: DF_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(out) :: DF_chi(1:NMAT,1:NMAT,1:num_faces)

complex(kind(0d0)) :: bPhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)


!complex(kind(0d0)) :: Uinv(1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT,1:num_faces),Ufm(1:NMAT,1:NMAT,1:num_faces)
!complex(kind(0d0)) :: diff_Omega(1:NMAT,1:NMAT,1:dimG)

complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: comm(1:NMAT,1:NMAT)
complex(kind(0d0)) :: trace,tmp
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT,1:num_necessary_faces)


!! for fermion mass term
complex(kind(0d0)) :: tmp_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: tmp_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: tmp_chi(1:NMAT,1:NMAT,1:num_necessary_faces)
!complex(kind(0d0)) :: tmp_vec(1:sizeD),tmp_Dvec(1:sizeD)
complex(kind(0d0)), allocatable :: tmp_vec(:),tmp_Dvec(:)
integer :: sizeD

integer :: s,l,f
integer :: i,j,k
integer :: a,b,c
integer :: r
integer :: label

!! for test
!integer :: ii,jj

!! preparation
!D_vec=(0d0,0d0)
!r_site=(0d0,0d0)
!r_link=(0d0,0d0)
!r_face=(0d0,0d0)

!call vec_to_mat(eta_mat,lambda_mat,chi_mat,vec)

if( p_mass == 0 ) then
#ifdef PARALLEL
  sizeD = dimG*(global_num_sites+global_num_links+global_num_faces)
#else
  sizeD = dimG*(num_sites+num_links+num_faces)
#endif
  allocate( tmp_vec(1:sizeD), tmp_Dvec(1:sizeD) )
endif

do s=1,num_necessary_sites
  do i=1,NMAT
    do j=1,NMAT
  !call make_traceless_matrix_from_modes(PhiMat(:,:,s),NMAT,Phi(:,s))
  bPhiMat(i,j,s)=dconjg(PhiMat(j,i,s))
  !call make_traceless_matrix_from_modes(bPhiMat(:,:,s),NMAT,dconjg(Phi(:,s)))
    enddo
  enddo
enddo

DF_eta=(0d0,0d0)
DF_lambda=(0d0,0d0)
DF_chi=(0d0,0d0)

if ( p1 == 0 ) then 
    !write(*,*) p1
!! (1) site action 
do s=1,num_sites
  call matrix_commutator(tmpmat1,PhiMat(:,:,s),eta_mat(:,:,s))
  DF_eta(:,:,s)=DF_eta(:,:,s) &
      +dcmplx(alpha_s(s))*(-0.5d0,0d0) *overall_factor*tmpmat1
enddo
endif

if (p2==0) then
!! (2) link action 1
do s=1,num_sites
  !tmpmat2 = -i \sum_{l\in\org(s)}( \alpha_l U_l^{-1}.\lambda_l.U_l ) 
  !          +i \sum_{l\in\tip(s)}( \alpha_l \lambda_l )
  tmpmat2=(0d0,0d0)
  do k=1,linkorg_to_s(s)%num_
    l=linkorg_to_s(s)%labels_(k)
    call matrix_3_product(tmpmat2,Umat(:,:,l),lambda_mat(:,:,l),Umat(:,:,l),&
      'C','N','N',&
      (0d0,-1d0)*dcmplx(alpha_l(l))*dconjg(U1Rfactor_link(l)*U1R_ratio(l)),'ADD')
  enddo
  do k=1,linktip_from_s(s)%num_
    l=linktip_from_s(s)%labels_(k)
    tmpmat2=tmpmat2 + (0d0,1d0)*dcmplx(alpha_l(l)) * lambda_mat(:,:,l)
  enddo
  DF_eta(:,:,s)=DF_eta(:,:,s) &
    + dcmplx(overall_factor) * tmpmat2
enddo

do l=1,num_links
  !tmpmat2=i alpha_l U_l \lambda_l U_l^{-1}
  tmpmat2=(0d0,-1d0)*eta_mat(:,:,link_org(l))
  call matrix_3_product(tmpmat2,Umat(:,:,l),eta_mat(:,:,link_tip(l)),Umat(:,:,l),&
    'N','N','C',(0d0,1d0)*dconjg(U1Rfactor_link(l)*U1R_ratio(l)),'ADD')
  DF_lambda(:,:,l) = DF_lambda(:,:,l) &
      + dcmplx(alpha_l(l) * overall_factor) * tmpmat2 
enddo
endif

    
if(p3==0) then
!! (3) link action 2 
do l=1,num_links
  tmpmat1=bPhiMat(:,:,link_org(l))
  call matrix_3_product(tmpmat1,Umat(:,:,l),PhiMat(:,:,link_tip(l)),Umat(:,:,l),&
    'N','C','C',dconjg(U1Rfactor_link(l)**2d0 * U1R_ratio(l)**2d0),'ADD')
    !'N','C','C',(1d0,0d0),'ADD')
  call matrix_commutator(tmpmat2,tmpmat1,lambda_mat(:,:,l))
  DF_lambda(:,:,l)=DF_lambda(:,:,l) &
      + dcmplx( alpha_l(l) * overall_factor ) * tmpmat2
enddo
endif


if(p4==0) then
!! (4) face action 1
do f=1,num_faces
  s=sites_in_f(f)%label_(1)
  call matrix_commutator(comm,PhiMat(:,:,s),chi_mat(:,:,f))
  DF_chi(:,:,f)=DF_chi(:,:,f) &
      + (-2d0,0d0) * dcmplx(alpha_f(f) * overall_factor) * comm
enddo
endif

!! (5) face action 2
if( p5==0 ) then
  if( m_omega == 0 ) then
    call Dirac_Omega_m0(DF_chi,DF_lambda,Umat,lambda_mat,chi_mat)
  elseif( m_omega == -1 ) then
    call Dirac_Omega_adm(DF_chi,DF_lambda,Umat,lambda_mat,chi_mat)
  else
    !if( p5_test == 0 ) then 
      call Dirac_Omega(DF_chi,DF_lambda,Umat,lambda_mat,chi_mat)
    !else
    !  call Dirac_Omega_test(DF_chi,DF_lambda,Umat,lambda_mat,chi_mat)
    !endif
  endif
endif

!! (mass) fermion mass term
if( p_mass == 0 ) then
!  call Dirac_mass(DF_chi,DF_lambda,DF_eta,Umat,eta_mat,lambda_mat,chi_mat)

#ifdef PARALLEL
  call localmat_to_globalvec(tmp_vec,eta_mat,lambda_mat,chi_mat)
  if(MYRANK==0) then 
    tmp_Dvec=(0d0,0d0)
    do i=1,sizeD/2
      tmp_Dvec(2*i-1)=tmp_Dvec(2*i-1) + overall_factor * dcmplx(mass_f)*tmp_vec(2*i) 
      tmp_Dvec(2*i)=tmp_Dvec(2*i) - overall_factor * dcmplx(mass_f)*tmp_vec(2*i-1) 
    enddo
  endif
  call MPI_BCAST(tmp_Dvec,sizeD,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
  call globalvec_to_localmat(tmp_eta,tmp_lambda,tmp_chi,tmp_Dvec)
#else
  call mat_to_vec(tmp_vec,eta_mat,lambda_mat,chi_mat)
  tmp_Dvec=(0d0,0d0)
  do i=1,sizeD/2
    tmp_Dvec(2*i-1)=tmp_Dvec(2*i-1) + overall_factor * dcmplx(mass_f)*tmp_vec(2*i) 
    tmp_Dvec(2*i)=tmp_Dvec(2*i) - overall_factor * dcmplx(mass_f)*tmp_vec(2*i-1) 
  enddo
  call vec_to_mat(tmp_eta,tmp_lambda,tmp_chi,tmp_Dvec)
#endif

  DF_eta=DF_eta+tmp_eta
  DF_lambda=DF_lambda+tmp_lambda
  DF_chi=DF_chi+tmp_chi
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! make traceless for (i,j)
do f=1,num_faces
  trace=(0d0,0d0)
  do j=1,NMAT
    trace=trace+DF_chi(j,j,f)
  enddo
  do j=1,NMAT
    DF_chi(j,j,f)=DF_chi(j,j,f)-trace/dcmplx(dble(NMAT))
  enddo
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do l=1,num_links
  trace=(0d0,0d0)
  do j=1,NMAT
    trace=trace+DF_lambda(j,j,l)
  enddo
  do j=1,NMAT
    DF_lambda(j,j,l)=DF_lambda(j,j,l)-trace/dcmplx(dble(NMAT))
  enddo
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do s=1,num_sites
  trace=(0d0,0d0)
  do j=1,NMAT
    trace=trace+DF_eta(j,j,s)
  enddo
  do j=1,NMAT
    DF_eta(j,j,s)=DF_eta(j,j,s)-trace/dcmplx(dble(NMAT))
  enddo
enddo
end subroutine Prod_Dirac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
subroutine Dirac_Omega(DF_chi,DF_lambda,Umat,lambda_mat,chi_mat)
use matrix_functions, only : make_matrix_traceless
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: lambda_mat(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: chi_mat(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)), intent(out) :: DF_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(out) :: DF_chi(1:NMAT,1:NMAT,1:num_faces)

complex(kind(0d0)) :: Ufm(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Uf0tom(1:NMAT,1:NMAT,0:m_omega-1)
complex(kind(0d0)) :: Xmat(1:NMAT,1:NMAT),Ymat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Cosinv(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Sinmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Omega(1:NMAT,1:NMAT)
complex(kind(0d0)) :: UXmat(1:NMAT,1:NMAT),YUmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat3(1:NMAT,1:NMAT)
integer :: i,j,k,l,f,a, l_place
complex(kind(0d0)) :: dir_factor
complex(kind(0d0)) :: im_over_2

integer :: ll,last_place
complex(kind(0d0)) :: U1Rfactor_fl


do f=1,num_necessary_faces
!! preparation( Cos^{-1} and Omega )
  call Make_face_variable(Uf(:,:),f,UMAT) 
  if( m_omega == 1) then 
    Ufm=Uf
  else
    call matrix_power(Ufm,Uf(:,:),m_omega)
  endif
  !! Sinmat and Cos^{-1}
  do i=1,NMAT
    do j=1,NMAT
      Cosinv(i,j) = Ufm(i,j) + dconjg(Ufm(j,i))
      Sinmat(i,j) = Ufm(i,j) - dconjg(Ufm(j,i))
    enddo
  enddo
  call Matrix_inverse(Cosinv)

  !! Omega = Cosinv . Sinmat
  call matrix_product(Omega,Cosinv,Sinmat)
  !call make_matrix_traceless(Omega)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! for development
  !call make_unit_matrix(Cosinv)
  !call make_unit_matrix(Sinmat)
  !call make_unit_matrix(Omega)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do l_place=1,links_in_f(f)%num_
    l=links_in_f(f)%link_labels_(l_place)
    !! U1Rfacor
    !call calc_U1Rfactor_fl(U1Rfactor_fl,f,l)
    call calc_U1Rfactor_fl_by_route(U1Rfactor_fl,f,l_place)

    dir_factor=&
      (0d0,-2d0)/dcmplx(m_omega)*overall_factor &
      *dcmplx(links_in_f(f)%link_dirs_(l_place))&
      *dcmplx(alpha_f(f)*beta_f(f)) &
      *U1Rfactor_fl


    if( f<=num_faces .or. l<=num_links ) then 

      call calc_XYmat(Xmat,Ymat,f,l_place,UMAT)
      !!
      do k=0,m_omega-1
        call matrix_power(Uf0tom(:,:,k),Uf(:,:),k)
      enddo

      do k=0,m_omega-1
        !! UXmat
        call matrix_product(UXmat,Uf0tom(:,:,k),Xmat)
        !! YUmat
        call matrix_product(YUmat,Ymat,Uf0tom(:,:,m_omega-k-1))
  
        !!!!  DF_chi
        if( f<=num_faces ) then 
          tmpmat1=(0d0,0d0)
          tmpmat2=(0d0,0d0)
          ! term 1
          ! tmpmat1 = UX.lambda.YU 
          call matrix_3_product(tmpmat1,UXmat,lambda_mat(:,:,l),YUmat)
          ! term 2
          ! tmpmat2 = YU^dag.lambda.UX^dag
          call matrix_3_product(tmpmat2,YUmat,lambda_mat(:,:,l),UXmat,'C','N','C')
  
          tmpmat3=tmpmat1+tmpmat2
          ! term 3 and term 4
          call matrix_product(tmpmat3,&
            tmpmat1-tmpmat2,Omega, &
            'N','N',(-1d0,0d0),'ADD')
  
          call matrix_product(tmpmat1,Cosinv,tmpmat3)

          DF_chi(:,:,f)=DF_chi(:,:,f) +  dir_factor * tmpmat1 
        endif
  
        !!!!  DF_lambda
        if( l<=num_links ) then 
          ! tmpmat1 = chi.Cosinv
          call matrix_product(tmpmat1,chi_mat(:,:,f),Cosinv)

          ! term 5 
          ! tmpmat3 = -YU.(chi.Cinv).UX - UX^dag.(chi.Cinv).(YU)^dag
          call matrix_3_product(tmpmat3,YUmat,-tmpmat1,UXmat)
          ! term 6
          call matrix_3_product(tmpmat3,UXmat,-tmpmat1,YUmat,'C','N','C',(1d0,0d0),'ADD')
  
          ! tmpmat2 = Sin.Cosinv.chi.Cosinv
          call matrix_product(tmpmat2,Omega,tmpmat1)
          
          call matrix_3_product(tmpmat3,YUmat,tmpmat2,UXmat,'N','N','N',&
            (1d0,0d0),'ADD')
          call matrix_3_product(tmpmat3,UXmat,tmpmat2,YUmat,'C','N','C',&
            (-1d0,0d0),'ADD')
  
          DF_lambda(:,:,l)=DF_lambda(:,:,l) + dir_factor * tmpmat3 
        endif
      enddo
    endif
  enddo
enddo


end subroutine Dirac_Omega



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
subroutine Dirac_Omega_m0(DF_chi,DF_lambda,Umat,lambda_mat,chi_mat)
use matrix_functions, only : make_unit_matrix
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: lambda_mat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: chi_mat(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), intent(out) :: DF_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(out) :: DF_chi(1:NMAT,1:NMAT,1:num_faces)


complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT,1:num_faces)
!complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Xmat(1:NMAT,1:NMAT),Ymat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
integer :: i,j,k,l,f,l_place
complex(kind(0d0)) :: dir_factor
complex(kind(0d0)) :: U1Rfactor_fl

!! m_omega=0 case
do f=1,num_necessary_faces
  do l_place=1,links_in_f(f)%num_
    l=links_in_f(f)%link_labels_(l_place)

    if( f <= num_faces .or. l <= num_links ) then 
      call calc_XYmat(Xmat,Ymat,f,l_place,UMAT)
      !call make_unit_matrix(Ymat)
       
      !call calc_U1Rfactor_fl(U1Rfactor_fl,f,l)
      call calc_U1Rfactor_fl_by_route(U1Rfactor_fl,f,l_place)

      !dir_factor=dcmplx(dble(links_in_f(f)%link_dirs_(l_place)) &
        !* alpha_f(f)*beta_f(f)*overall_factor )*(0d0,1d0)
      dir_factor=&
        (0d0,-1d0)*overall_factor &
        *dcmplx(links_in_f(f)%link_dirs_(l_place))&
        *dcmplx(alpha_f(f)*beta_f(f)) &
        *U1Rfactor_fl

      if( f <= num_faces ) then 
        call matrix_3_product(tmpmat1,Xmat,lambda_mat(:,:,l),Ymat)
        call matrix_3_product(tmpmat1,Ymat,lambda_mat(:,:,l),Xmat,'C','N','C',(1d0,0d0),'ADD')
        DF_chi(:,:,f)=DF_chi(:,:,f) + dir_factor * tmpmat1
      endif
  
      if( l <= num_links ) then 
        call matrix_3_product(tmpmat1,Ymat,chi_mat(:,:,f),Xmat)
        call matrix_3_product(tmpmat1,Xmat,chi_mat(:,:,f),Ymat,'C','N','C',(1d0,0d0),'ADD')
        DF_lambda(:,:,l)=DF_lambda(:,:,l) - dir_factor * tmpmat1
      endif
    endif
  enddo
enddo
end subroutine Dirac_Omega_m0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
subroutine Dirac_Omega_adm(DF_chi,DF_lambda,Umat,lambda_mat,chi_mat)
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: lambda_mat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: chi_mat(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), intent(inout) :: DF_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(inout) :: DF_chi(1:NMAT,1:NMAT,1:num_faces)

complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)) :: SinU(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Xmat(1:NMAT,1:NMAT),Ymat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
integer :: i,j,l,f,l_place
complex(kind(0d0)) :: dir_factor
complex(kind(0d0)) :: Bval
complex(kind(0d0)) :: trace

do f=1,num_necessary_faces
  !!!!!!!!!!!!!!!
  !! Uf and SinU
  call make_face_variable(Uf(:,:),f,Umat)
  do i=1,NMAT
    do j=1,NMAT
      sinU(i,j) = Uf(i,j) - dconjg( Uf(j,i) )
    enddo
  enddo
  !!!!!!!!!!!!!!!
  !! Bval
  Bval=(1d0,0d0)
  do i=1,NMAT
    Bval=Bval - ((2d0,0d0)-Uf(i,i)-dconjg(Uf(i,i)))/(e_max*e_max) 
  enddo

  do l_place=1,links_in_f(f)%num_
    l=links_in_f(f)%link_labels_(l_place)

    if( f <= num_faces .or. l <= num_links ) then 
      dir_factor&
        =(0d0,-1d0)*dcmplx(&
          dble(links_in_f(f)%link_dirs_(l_place)) &
          * alpha_f(f) * beta_f(f) * overall_factor)

      !!!!!!!!!!!!!!!!!!!!!!!!
      !! X_(f,l) and Y_(f,l)
      call calc_XYmat(Xmat,Ymat,f,l_place,UMAT)
    endif

    if( f<= num_faces ) then
      !!!!!!!!!!!!!!!!!!!!!!!!
      !! tmpmat1 = X.lambda.Y  
      !! tmpmat2 = Y^\dag .lambda. X^\dag
      call matrix_3_product(tmpmat1,Xmat,lambda_mat(:,:,l),Ymat)
      call matrix_3_product(tmpmat2,Ymat,lambda_mat(:,:,l),Xmat, 'C','N','C')
      trace=(0d0,0d0)
      do i=1,NMAT
        trace=trace+tmpmat1(i,i)-tmpmat2(i,i)
      enddo

      Df_chi(:,:,f)=Df_chi(:,:,f) &
        + dir_factor/Bval * (tmpmat1+tmpmat2) &
        - dir_factor/(Bval*Bval*e_max*e_max) * trace * SinU
    endif

    if( l<= num_links ) then
      !!!!!!!!!!!!!!!!!!!!!!!!
      !! tmpmat1 = Y.chi.X + X^\dag .chi. Y^\dag
      call matrix_3_product(tmpmat1,Ymat,chi_mat(:,:,f),Xmat)
      call matrix_3_product(tmpmat1,Xmat,chi_mat(:,:,f),Ymat,&
        'C','N','C',(1d0,0d0),'ADD')
      !!!!!!!!!!!!!!!!!!!!!!!!
      !! tmpmat2 = Y.X - X^dag.Y^dag
      call matrix_product(tmpmat2,Ymat,Xmat)
      call matrix_product(tmpmat2,Xmat,Ymat,'C','C',(-1d0,0d0),'ADD')

      trace=(0d0,0d0)
      do i=1,NMAT
        do j=1,NMAT
          trace=trace+chi_mat(i,j,f)*SinU(j,i)
        enddo
      enddo

      Df_lambda(:,:,l)=Df_lambda(:,:,l) - dir_factor/Bval * tmpmat1
      Df_lambda(:,:,l)=Df_lambda(:,:,l) &
        + dir_factor/(Bval*Bval*e_max*e_max) * trace * tmpmat2
    endif
      
  enddo
enddo

end subroutine Dirac_Omega_adm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
subroutine Dirac_mass(DF_chi,DF_lambda,DF_eta,Umat,eta_mat,lambda_mat,chi_mat)
implicit none
complex(kind(0d0)), intent(out) :: DF_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(out) :: DF_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(out) :: DF_chi(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: eta_mat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: lambda_mat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: chi_mat(1:NMAT,1:NMAT,1:num_necessary_faces)

complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Xmat(1:NMAT,1:NMAT),Ymat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
integer :: i,j,k,l,s,f,l_place
complex(kind(0d0)) :: dir_factor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! lambda-eta part
! modify "\lambda_l D_l\eta" term
do s=1,num_sites
  !tmpmat2 = -i (\sum_{l\in\org(s)}( \alpha_l U_l^{-1}.\lambda_l.U_l ) 
  !              + \sum_{l\in\tip(s)}( \alpha_l \lambda_l ) )
  tmpmat2=(0d0,0d0)
  do k=1,linkorg_to_s(s)%num_
    l=linkorg_to_s(s)%labels_(k)
    call matrix_3_product(tmpmat2,Umat(:,:,l),lambda_mat(:,:,l),Umat(:,:,l),&
      'C','N','N',(0d0,-1d0)*dcmplx(alpha_l(l)),'ADD')
  enddo
  do k=1,linktip_from_s(s)%num_
    l=linktip_from_s(s)%labels_(k)
    tmpmat2=tmpmat2 + (0d0,-1d0)*dcmplx(alpha_l(l)) * lambda_mat(:,:,l)
  enddo
  DF_eta(:,:,s)=DF_eta(:,:,s) &
    + dcmplx(overall_factor*mass_f) * tmpmat2
enddo

do l=1,num_links
  !tmpmat2=i alpha_l U_l \lambda_l U_l^{-1}
  tmpmat2=(0d0,1d0)*eta_mat(:,:,link_org(l))
  call matrix_3_product(tmpmat2,Umat(:,:,l),eta_mat(:,:,link_tip(l)),Umat(:,:,l),&
    'N','N','C',(0d0,1d0),'ADD')
  DF_lambda(:,:,l) = DF_lambda(:,:,l) &
      + dcmplx(alpha_l(l) * overall_factor * mass_f) * tmpmat2
enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! lambda-chi part
! modify "\chi Q(Uf-Uf^-1) "
do f=1,num_necessary_faces
  dir_factor=dcmplx(mass_f*overall_factor*alpha_f(f)*beta_f(f))*(0d0,-1d0)
  do l_place=1,links_in_f(f)%num_
    l=links_in_f(f)%link_labels_(l_place)

    if( f <= num_faces .or. l <= num_links ) then 
      call calc_XYmat(Xmat,Ymat,f,l_place,UMAT)
    endif

    if( f <= num_faces ) then 
      call matrix_3_product(tmpmat1,Xmat,lambda_mat(:,:,l),Ymat)
      call matrix_3_product(tmpmat1,Ymat,lambda_mat(:,:,l),Xmat,'C','N','C',(1d0,0d0),'ADD')
      DF_chi(:,:,f)=DF_chi(:,:,f) &
        + dir_factor * tmpmat1
    endif

    if( l <= num_links ) then 
      call matrix_3_product(tmpmat1,Ymat,chi_mat(:,:,f),Xmat)
      call matrix_3_product(tmpmat1,Xmat,chi_mat(:,:,f),Ymat,'C','N','C',(1d0,0d0),'ADD')
      DF_lambda(:,:,l)=DF_lambda(:,:,l) &
        - dir_factor * tmpmat1
    endif
  enddo
enddo

end subroutine Dirac_mass

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
subroutine Dirac_Omega_test(DF_chi,DF_lambda,Umat,lambda_mat,chi_mat)
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: lambda_mat(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: chi_mat(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)), intent(out) :: DF_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(out) :: DF_chi(1:NMAT,1:NMAT,1:num_faces)


complex(kind(0d0)) :: Ufm(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Ufminv(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)) :: UnitMat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Mat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Xmat(1:NMAT,1:NMAT),Ymat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Cosinv(1:NMAT,1:NMAT),Omega_mat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: UXmat(1:NMAT,1:NMAT),YUmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat3(1:NMAT,1:NMAT)
integer :: i,j,k,l,f,a, l_place
complex(kind(0d0)) :: dir_factor

UnitMat=(0d0,0d0)
do i=1,NMAT
  UnitMat(i,i)=(1d0,0d0)
enddo

do f=1,num_faces
!! preparation( Cos^{-1} and Omega )
  call Make_face_variable(Uf(:,:),f,UMAT) 
  call matrix_power(Ufm,Uf(:,:),m_omega)
  Ufminv=Ufm
  call Matrix_inverse(Ufminv)
  do j=1,NMAT
    do i=1,NMAT
      Cosinv(i,j) = Ufm(i,j) + dconjg(Ufm(j,i))
      tmpmat1(i,j) = Ufm(i,j) - dconjg(Ufm(j,i))
      tmpmat2(i,j) = Uf(i,j) - dconjg(Uf(j,i))
    enddo
  enddo
  !Cosinv = Ufm + Ufminv
  !tmpmat1 = Ufm - Ufminv
  !! Cos^{-1}
  !call HermitianMatrix_inverse(Cosinv)
  call Matrix_inverse(Cosinv)
  !! Omega
  call matrix_product(Omega_mat,Cosinv,tmpmat1,&
    'N','N', (0d0,-2d0)/dcmplx(dble(m_omega)))
  !call ZHEMM('L','U',NMAT,NMAT,(0d0,-2d0)/dcmplx(dble(m_omega)), &
  !  Cosinv(:,:), NMAT, &
  !  tmpmat1, NMAT, &
  !  (0d0,0d0),Omega_mat , NMAT)
  if( NMAT > 2 ) then 
    call make_matrix_traceless(Omega_mat)
  endif
  !Ufminv=Ufminv+Ufm
  !call make_matrix_traceless(Ufminv)

!!!!!!!!!!!!!!!!!!!!!!!!!
  MAT=tmpmat2
!!!!!!!!!!!!!!!!!!!!!!!!!

  do l_place=1,links_in_f(f)%num_
    l=links_in_f(f)%link_labels_(l_place)
    dir_factor=&
      dcmplx(dble(links_in_f(f)%link_dirs_(l_place)))*(0d0,2d0)/dcmplx(dble(m_omega)) 


      !!!!  DF_chi
      call matrix_3_product(tmpmat1,MAT,lambda_mat(:,:,l),UnitMat,'N','N','N')
      call matrix_3_product(tmpmat2,UnitMat,lambda_mat(:,:,l),MAT,'C','N','C')
      DF_chi(:,:,f)=DF_chi(:,:,f) &
        - overall_factor &
          * dcmplx(alpha_f(f)*beta_f(f)) &
          * dir_factor * (tmpmat1+tmpmat2)
      
      !!!!  DF_lambda
      call matrix_3_product(tmpmat1,UnitMat,chi_mat(:,:,f),MAT,'N','N','N')
      call matrix_3_product(tmpmat2,MAT,chi_mat(:,:,f),UnitMat,'C','N','C')
      DF_lambda(:,:,l)=DF_lambda(:,:,l) &
        + overall_factor &
          * dcmplx(alpha_f(f)*beta_f(f)) &
          * dir_factor * (tmpmat1+tmpmat2)
  enddo

enddo 
end subroutine Dirac_Omega_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
subroutine Dirac_Omega_org(DF_chi,DF_lambda,Umat,lambda_mat,chi_mat)

implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: lambda_mat(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: chi_mat(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)), intent(out) :: DF_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(out) :: DF_chi(1:NMAT,1:NMAT,1:num_faces)


complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat3(1:NMAT,1:NMAT)
complex(kind(0d0)) :: comm(1:NMAT,1:NMAT)
complex(kind(0d0)) :: acomm(1:NMAT,1:NMAT)
complex(kind(0d0)) :: line1(1:NMAT,1:NMAT)
complex(kind(0d0)) :: line2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: trace,tmp

complex(kind(0d0)) :: AMAT(1:NMAT,1:NMAT)
complex(kind(0d0)) :: BMAT(1:NMAT,1:NMAT)
complex(kind(0d0)) :: sinU(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: cosUinv(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT,1:num_faces)

complex(kind(0d0)) :: Ufm(1:NMAT,1:NMAT,1:num_faces)
integer :: l,f,i,j,k

  do f=1,num_faces
      call Make_face_variable(Uf(:,:,f),f,UMAT) 
    if( m_omega .ne. 0) then 
      call matrix_power(Ufm(:,:,f),Uf(:,:,f),m_omega)
      call calc_sinU_and_cosUinv(sinU(:,:,f),cosUinv(:,:,f),Ufm(:,:,f))
    endif
  enddo



  do f=1,num_faces
    do i=1,links_in_f(f)%num_
      l=links_in_f(f)%link_labels_(i)
  
      if( m_omega == 0 ) then
        call calc_Amat(Amat,f,i,1,Uf(:,:,f),UMAT)
        call calc_Bmat(Bmat,f,i,1,Uf(:,:,f),UMAT)
  
        call matrix_3_product(tmpmat1,Amat,lambda_mat(:,:,l),Bmat)
        call matrix_3_product(tmpmat2,Bmat,lambda_mat(:,:,l),Amat,'C','N','C')
  
        DF_chi(:,:,f)=DF_chi(:,:,f)&
          + dcmplx(dble(links_in_f(f)%link_dirs_(i)))*(0d0,1d0)&
            * dcmplx(overall_factor) * dcmplx(alpha_f(f)*beta_f(f)) & 
            * (tmpmat1 + tmpmat2)
      else
        line1=(0d0,0d0)
        line2=(0d0,0d0)
        do k=1,m_omega
          call calc_Amat(Amat,f,i,k,Uf(:,:,f),UMAT)
          call calc_Bmat(Bmat,f,i,k,Uf(:,:,f),UMAT)
    
          ! tmpmat2=A.lambda.B
          call matrix_3_product(tmpmat2,Amat,lambda_mat(:,:,l),Bmat)
          ! tmpmat3=B^dag.lambda.A^dag 
          call matrix_3_product(tmpmat3,Bmat,lambda_mat(:,:,l),Amat,'C','N','C')
    
          ! line1 = A.lambda.B+Bdag.lambda.Adag
          line1=line1+tmpmat2+tmpmat3
          ! line2 = cosUinv.(A.lambda.B-Bdag.lambda.Adag).cosUinv
          call matrix_3_product(tmpmat1,cosUinv(:,:,f),tmpmat2-tmpmat3,cosUinv(:,:,f))
          line2=line2+tmpmat1
        enddo
    
        ! line1 = {A.lambda.B+B^dag.lambda.A^dag , (Um+Uminv)^{-1}
        tmpmat1=line1
        call matrix_AntiCommutator(line1,tmpmat1,cosUinv(:,:,f))
        ! line2 = {cosUinv.(A.lambda.B-Bdag.lambda.Adag).cosUinv, sinU}
        tmpmat1=line2
        call matrix_AntiCommutator(line2,tmpmat1,sinU(:,:,f))
        
        DF_chi(:,:,f)=DF_chi(:,:,f)&
          + dcmplx(dble(links_in_f(f)%link_dirs_(i)))*(0d0,1d0)&
            * dcmplx(overall_factor) / dcmplx(dble(m_omega))&
            * dcmplx(alpha_f(f)*beta_f(f)) * (line1-line2)
      endif
    enddo
  enddo
  
  do l=1,num_links
    do i=1,face_in_l(l)%num_       
      f=face_in_l(l)%label_(i)
      ! j: position of the link l in the face l
      do j=1,links_in_f(f)%num_
        if ( l == links_in_f(f)%link_labels_(j) ) exit
      enddo
  
      if( m_omega == 0 ) then 
        call calc_Amat(Amat,f,j,1,Uf(:,:,f),UMAT)
        call calc_Bmat(Bmat,f,j,1,Uf(:,:,f),UMAT)
  
        call matrix_3_product(tmpmat1,Bmat,chi_mat(:,:,f),Amat)
        call matrix_3_product(tmpmat2,Amat,chi_mat(:,:,f),Bmat,'C','N','C')
  
        DF_lambda(:,:,l)=DF_lambda(:,:,l) &
          + dcmplx(dble(links_in_f(f)%link_dirs_(j)))*(0d0,-1d0)&
            * dcmplx(overall_factor) * dcmplx(alpha_f(f)*beta_f(f)) &
            * (tmpmat1+tmpmat2)
      else
        ! tmpmat1= { (Uf+Ufinv)^{-1} , \chi_f }
        call matrix_anticommutator(tmpmat1,cosUinv(:,:,f),chi_mat(:,:,f))
    
        ! tmpmat2= (Uf+Ufinv)^{-1}.{ Uf-Ufinv , \chi_f }.(Uf+Ufinv)^{-1} 
        call matrix_anticommutator(tmpmat3,sinU(:,:,f),chi_mat(:,:,f))
        call matrix_3_product(tmpmat2,cosUinv(:,:,f),tmpmat3,cosUinv(:,:,f))
    
        line1=tmpmat1-tmpmat2
        line2=tmpmat1+tmpmat2
    
        acomm=(0d0,0d0)
        do k=1,m_omega
          call calc_Amat(Amat,f,j,k,Uf(:,:,f),UMAT)
          call calc_Bmat(Bmat,f,j,k,Uf(:,:,f),UMAT)
          ! tmpmat1= B.(tmpmat1-tmpmat2).A
          call matrix_3_product(tmpmat1,BMAT,line1,AMAT)
          ! tmpmat2 = Adag.(tmpmat1+tmpmat2).Bdag
          call matrix_3_product(tmpmat2, AMAT,line2,BMAT,'C','N','C')
  
          acomm=acomm+tmpmat1+tmpmat2
        enddo
    
        DF_lambda(:,:,l)=DF_lambda(:,:,l) &
          + dcmplx(dble(links_in_f(f)%link_dirs_(j)))*(0d0,1d0)&
            * dcmplx(-overall_factor) / dcmplx(dble(m_omega))&
            * dcmplx(alpha_f(f)*beta_f(f)) * acomm
      endif
    enddo
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     ! make traceless for (i,j)
!     trace=(0d0,0d0)
!     do j=1,NMAT
!       trace=trace+DF_lambda(j,j,l)
!     enddo
!     do j=1,NMAT
!       DF_lambda(j,j,l)=DF_lambda(j,j,l)-trace/cmplx(dble(NMAT))
!     enddo
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  enddo

end subroutine Dirac_Omega_org

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! check anti-symmetricity of D
subroutine check_Dirac(UMAT,PhiMat)
implicit none 

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)

!complex(kind(0d0)) :: PF(1:sizeD)
!complex(kind(0d0)) :: Dirac(1:sizeD,1:sizeD), Det
complex(kind(0d0)), allocatable :: PF(:)
complex(kind(0d0)), allocatable :: Dirac(:,:)
!complex(kind(0d0)), allocatable :: mat(:,:)
!complex(kind(0d0)), allocatable :: eigen(:)
complex(kind(0d0)) :: argDet
double precision :: logDet
integer :: sizeD
double precision :: rtmp,dtmp
integer :: seed
integer :: i,j
integer :: s,l,f,a,b
integer :: ss,ll,ff

#ifdef PARALLEL
  sizeD = dimG*(global_num_sites+global_num_links+global_num_faces)
#else
  sizeD = dimG*(num_sites+num_links+num_faces)
#endif
allocate( PF(1:sizeD), Dirac(1:sizeD,1:sizeD) )
!allocate( eigen(1:sizeD), mat(1:sizeD,1:sizeD) )

call make_Dirac(Dirac,UMAT,PhiMat)

!mat=Dirac
!call matrix_eigenvalues(eigen,mat)
!do i=1,sizeD
!  write(*,*) "eigenvalue",i,eigen(i)
!enddo


if( Dirac_write == 1 ) then 
#ifdef PARALLEL
if(MYRANK==0) then
#endif 
write(*,*) "## PhiMat"
do l=1,num_sites
  do i=1,NMAT
    do j=1,NMAT
      write(*,*) l,i,j,PhiMat(i,j,l)
    enddo
  enddo
enddo
write(*,*) "## UMAT"
do l=1,num_links
  do i=1,NMAT
    do j=1,NMAT
      write(*,*) l,i,j,UMAT(i,j,l)
    enddo
  enddo
enddo
write(*,*) "###################"




do s=1,global_num_sites
  do a=1,dimG
    i=site_index(a,s,NMAT)

    do ss=1,global_num_sites
      do b=1,dimG
        j=site_index(b,ss,NMAT)
        if( abs(Dirac(i,j)) > epsilon ) then 
          write(*,*) "site",s,"site",ss,abs(Dirac(i,j))
        endif
      enddo
    enddo
    do ll=1,global_num_links
      do b=1,dimG
        j=link_index(b,ll,NMAT,global_num_sites)
        if( abs(Dirac(i,j)) > epsilon ) then 
          write(*,*) "site",s,"link",ll,abs(Dirac(i,j))
        endif
      enddo
    enddo
    do ff=1,global_num_faces
      do b=1,dimG
        j=face_index(b,ss,NMAT,global_num_sites,global_num_links)
        if( abs(Dirac(i,j)) > epsilon ) then 
          write(*,*) "site",s,"face",ff,abs(Dirac(i,j))
        endif
      enddo
    enddo
  enddo
enddo

do l=1,global_num_links
  do a=1,dimG
    i=link_index(a,l,NMAT,global_num_sites)

    do ss=1,global_num_sites
      do b=1,dimG
        j=site_index(b,ss,NMAT)
        if( abs(Dirac(i,j)) > epsilon ) then 
          write(*,*) "link",l,"site",ss,abs(Dirac(i,j))
        endif
      enddo
    enddo
    do ll=1,global_num_links
      do b=1,dimG
        j=link_index(b,ll,NMAT,global_num_sites)
        if( abs(Dirac(i,j)) > epsilon ) then 
          write(*,*) "link",l,"link",ll,abs(Dirac(i,j))
        endif
      enddo
    enddo
    do ff=1,global_num_faces
      do b=1,dimG
        j=face_index(b,ss,NMAT,global_num_sites,global_num_links)
        if( abs(Dirac(i,j)) > epsilon ) then 
          write(*,*) "link",l,"face",ff,abs(Dirac(i,j))
        endif
      enddo
    enddo
  enddo
enddo

do f=1,global_num_faces
  do a=1,dimG
    i=face_index(a,f,NMAT,global_num_sites,global_num_links)

    do ss=1,global_num_sites
      do b=1,dimG
        j=site_index(b,ss,NMAT)
        if( abs(Dirac(i,j)) > epsilon ) then 
          write(*,*) "face",f,"site",ss,abs(Dirac(i,j))
        endif
      enddo
    enddo
    do ll=1,global_num_links
      do b=1,dimG
        j=link_index(b,ll,NMAT,global_num_sites)
        if( abs(Dirac(i,j)) > epsilon ) then 
          write(*,*) "face",f,"link",ll,abs(Dirac(i,j))
        endif
      enddo
    enddo
    do ff=1,global_num_faces
      do b=1,dimG
        j=face_index(b,ss,NMAT,global_num_sites,global_num_links)
        if( abs(Dirac(i,j)) > epsilon ) then 
          write(*,*) "face",f,"face",ff,abs(Dirac(i,j))
        endif
      enddo
    enddo
  enddo
enddo

!do i=1,sizeD
!  do j=1,sizeD
!    if( abs(Dirac(i,j)) > epsilon ) then 
!      write(*,*) i,j,abs(Dirac(i,j))
!    endif
!  enddo
!enddo
#ifdef PARALLEL
endif
#endif

endif ! if ( Dirac_write == 1 ) then

#ifdef PARALLEL
if( MYRANK == 0 ) then 
#endif
  rtmp=0d0
  dtmp=0d0
  do i=1,sizeD-1
    rtmp=rtmp+dble(Dirac(i,i)*dconjg(Dirac(i,i)))
    do j=i+1,sizeD
      rtmp=rtmp+dble( (Dirac(i,j)+Dirac(j,i)) * dconjg(Dirac(i,j)+Dirac(j,i)) ) 
    enddo
  enddo
  write(*,*) "# D's diagonal elements?", dtmp+dble(Dirac(sizeD,sizeD)*dconjg(Dirac(sizeD,sizeD)))
  write(*,*) "# Is D anti-symmetric?", rtmp
  
#ifdef PARALLEL
endif
#endif
  call make_DdagD(Dirac,UMAT,PhiMat)
#ifdef PARALLEL
if( MYRANK == 0 ) then 
#endif
  rtmp=0d0
  do i=1,sizeD
    do j=i,sizeD
      rtmp=rtmp&
          +dble ( &
            ( Dirac(i,j) - dconjg( Dirac(j,i) ) ) &
            *dconjg( Dirac(i,j) - dconjg( Dirac(j,i) ) ) )
    enddo
  enddo
  write(*,*) "# Is D^\dag D hermitian?", rtmp
  
  call Matrix_Determinant(logDet,argDet,Dirac)
  write(*,*) "# log|Det(D)|=", logDet
#ifdef PARALLEL
endif
#endif

end subroutine check_Dirac


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! make Dirac matrix
subroutine make_Dirac(Dirac,UMAT,PhiMat)
use SUN_generators, only : trace_MTa
implicit none

complex(kind(0d0)), intent(inout) :: Dirac(:,:)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
!complex(kind(0d0)) :: Phi(1:dimG,1:num_sites)

complex(kind(0d0)) :: DF_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: DF_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: DF_chi(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: eta_mat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: lambda_mat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: chi_mat(1:NMAT,1:NMAT,1:num_necessary_faces)

complex(kind(0d0)), allocatable :: unitvec(:)
integer i,j
integer s,a
integer :: sizeD

#ifdef PARALLEL
sizeD=(NMAT*NMAT-1)*(global_num_sites+global_num_links+global_num_faces)
#else
sizeD=(NMAT*NMAT-1)*(num_sites+num_links+num_faces)
#endif


if( size(Dirac,1) /= sizeD ) then
  write(*,*) "size of Dirac is not ",sizeD,"but",size(Dirac,1)
  call stop_for_test
endif

allocate( unitvec(1:sizeD) )

do i=1,sizeD
  unitvec=(0d0,0d0)
  unitvec(i)=(1d0,0d0)
  !write(*,*) MYRANK,"===========",i,"========="
  !call prod_Dirac_vec(Dirac(:,i),unitvec,sizeD,UMAT,PhiMat)

#ifdef PARALLEL
  call globalvec_to_localmat(eta_mat,lambda_mat,chi_mat,unitvec)
#else
  call vec_to_mat(eta_mat,lambda_mat,chi_mat,unitvec)
#ENDIF

  call Prod_Dirac(&
    DF_eta, DF_lambda, DF_chi, &
    eta_mat,lambda_mat,chi_mat, &
    UMAT,Phimat)

#ifdef PARALLEL
  ! only Rank0 possesses Dirac
  call localmat_to_globalvec(Dirac(:,i),DF_eta,DF_lambda,DF_chi)
#else
  call mat_to_vec(Dirac(:,i),DF_eta,DF_lambda,DF_chi)
#endif

enddo
!call MPI_FINALIZE(IERR)
!stop
!do i=1,sizeD
!  do j=1,sizeD
!    if (abs(Dirac(i,j)) > 1d-5) then
!      write(*,*) MYRANK,i,j,abs(Dirac(i,j))
!    endif
!  enddo
!enddo


end subroutine make_Dirac


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! make D^\dagger D matrix
subroutine make_DdagD(DdagD,UMAT,PhiMat)
implicit none

!complex(kind(0d0)), intent(inout) :: DdagD(1:sizeD,1:sizeD)
complex(kind(0d0)), intent(inout) :: DdagD(:,:)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) :: DF_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: DF_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: DF_chi(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: eta_mat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: lambda_mat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: chi_mat(1:NMAT,1:NMAT,1:num_necessary_faces)

!complex(kind(0d0)) :: unitvec(1:sizeD),dirac_vec(1:sizeD)
complex(kind(0d0)), allocatable :: unitvec(:)
integer i,j,sizeD



#ifdef PARALLEL
sizeD=(NMAT*NMAT-1)*(global_num_sites+global_num_links+global_num_faces)
#else
sizeD=(NMAT*NMAT-1)*(num_sites+num_links+num_faces)
#endif

if( size(DdagD,1) /= sizeD ) then
  write(*,*) "size of DdagD is not ",sizeD,"but",size(DdagD,1)
  call stop_for_test
endif

allocate( unitvec(1:sizeD) )

do i=1,sizeD
  unitvec=(0d0,0d0)
  unitvec(i)=(1d0,0d0)
  !call prod_DdagD_vec(DdagD(:,i),unitvec,sizeD,UMAT,Phimat)

  !write(*,*) i
#ifdef PARALLEL
  call globalvec_to_localmat(eta_mat,lambda_mat,chi_mat,unitvec)
#else
  call vec_to_mat(eta_mat,lambda_mat,chi_mat,unitvec)
#ENDIF

  call Prod_DdagD(&
    DF_eta, DF_lambda, DF_chi, &
    eta_mat,lambda_mat,chi_mat, &
    UMAT,Phimat)

#ifdef PARALLEL
  ! only Rank0 possesses Dirac
  call localmat_to_globalvec(DdagD(:,i),DF_eta,DF_lambda,DF_chi)
#else

  call mat_to_vec(DdagD(:,i),DF_eta,DF_lambda,DF_chi)
#endif

enddo
end subroutine make_DdagD


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compute D^dagger.D.vec
!!   1) D_{ij}=-D_{ji}
!!   2) D^\dag)_{ij} = (D^*)_{ji} = - D^*_{ij}
!! ==>
!!   [D^\dagger D v]_i 
!!     = conjg(D_{ji}) D_{jk} v_k 
!!     = -conjg(D_{ij}) D_{jk} v_k                 
!! w=Dv => D^\dagger w = -( D.conjg(w) )^\dagger 
!!                     = -conjg[ D.conjg( D v ) ]
!subroutine Prod_DdagD(DdagD_vec, vec, Vsize,UMAT,Phimat)
subroutine Prod_DdagD(&
    DF_eta, DF_lambda, DF_chi, &
    eta_mat,lambda_mat,chi_mat, &
    UMAT,Phimat)
implicit none

!integer, intent(in) :: Vsize !! vecsize must be sizeD
complex(kind(0d0)), intent(out) :: DF_eta(:,:,:)!(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(out) :: DF_lambda(:,:,:)!(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(out) :: DF_chi(:,:,:)!(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)),intent(in) :: eta_mat(:,:,:)!(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)),intent(in) :: lambda_mat(:,:,:)!(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)),intent(in) :: chi_mat(:,:,:)!(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)),intent(in) :: UMAT(:,:,:)!(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)),intent(in) :: PhiMat(:,:,:)!(1:NMAT,1:NMAT,1:num_necessary_sites)
!complex(kind(0d0)),intent(in) :: vec(:)
!complex(kind(0d0)),intent(out) :: DdagD_vec(:)

complex(kind(0d0)) :: tmp_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: tmp_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: tmp_chi(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: tmp2_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: tmp2_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: tmp2_chi(1:NMAT,1:NMAT,1:num_necessary_faces)

integer :: i,j,s,l,f
double precision :: rtmp

call Prod_Dirac(&
    tmp_eta, tmp_lambda, tmp_chi, &
    eta_mat,lambda_mat,chi_mat, &
    UMAT,PhiMat)

do s=1,num_sites
  do j=1,NMAT
    do i=1,NMAT
      tmp2_eta(i,j,s) = -dconjg( tmp_eta(j,i,s) )
    enddo
  enddo
enddo
do l=1,num_links
  do j=1,NMAT
    do i=1,NMAT
      tmp2_lambda(i,j,l) = -dconjg( tmp_lambda(j,i,l) )
    enddo
  enddo
enddo
do f=1,num_faces
  do j=1,NMAT
    do i=1,NMAT
      tmp2_chi(i,j,f) = -dconjg( tmp_chi(j,i,f) )
    enddo
  enddo
enddo
#ifdef PARALLEL
call syncronize_sites(tmp2_eta)
!write(*,*) MYRANK,"testS"
call syncronize_faces(tmp2_chi)
!write(*,*) MYRANK,"testF"
call syncronize_links(tmp2_lambda)
!write(*,*) MYRANK,"testL"
#endif

call Prod_Dirac(&
    tmp_eta, tmp_lambda, tmp_chi, &
    tmp2_eta,tmp2_lambda,tmp2_chi, &
    UMAT,PhiMat)

do s=1,num_sites
  do j=1,NMAT
    do i=1,NMAT
      DF_eta(i,j,s) = dconjg( tmp_eta(j,i,s) )
    enddo
  enddo
enddo
do l=1,num_links
  do j=1,NMAT
    do i=1,NMAT
      DF_lambda(i,j,l) = dconjg( tmp_lambda(j,i,l) )
    enddo
  enddo
enddo
do f=1,num_faces
  do j=1,NMAT
    do i=1,NMAT
      DF_chi(i,j,f) = dconjg( tmp_chi(j,i,f) )
    enddo
  enddo
enddo

end subroutine Prod_DdagD


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
subroutine Prod_DiracDag(&
    DF_eta, DF_lambda, DF_chi, &
    eta_mat,lambda_mat,chi_mat, &
    UMAT,Phimat)
implicit none

!integer, intent(in) :: Vsize !! vecsize must be sizeD
complex(kind(0d0)), intent(out) :: DF_eta(:,:,:)!(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(out) :: DF_lambda(:,:,:)!(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(out) :: DF_chi(:,:,:)!(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)),intent(in) :: eta_mat(:,:,:)!(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)),intent(in) :: lambda_mat(:,:,:)!(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)),intent(in) :: chi_mat(:,:,:)!(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)),intent(in) :: UMAT(:,:,:)!(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)),intent(in) :: PhiMat(:,:,:)!(1:NMAT,1:NMAT,1:num_necessary_sites)
!complex(kind(0d0)),intent(in) :: vec(:)
!complex(kind(0d0)),intent(out) :: DdagD_vec(:)

complex(kind(0d0)) :: tmp_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: tmp_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: tmp_chi(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: tmp2_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: tmp2_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: tmp2_chi(1:NMAT,1:NMAT,1:num_necessary_faces)

integer :: i,j,s,l,f
double precision :: rtmp


do s=1,num_necessary_sites
  do j=1,NMAT
    do i=1,NMAT
      tmp2_eta(i,j,s) = dconjg( eta_mat(j,i,s) )
    enddo
  enddo
enddo
do l=1,num_necessary_links
  do j=1,NMAT
    do i=1,NMAT
      tmp2_lambda(i,j,l) = dconjg( lambda_mat(j,i,l) )
    enddo
  enddo
enddo
do f=1,num_necessary_faces
  do j=1,NMAT
    do i=1,NMAT
      tmp2_chi(i,j,f) = dconjg( chi_mat(j,i,f) )
    enddo
  enddo
enddo

call Prod_Dirac(&
    tmp_eta, tmp_lambda, tmp_chi, &
    tmp2_eta,tmp2_lambda,tmp2_chi, &
    UMAT,PhiMat)

do s=1,num_sites
  do j=1,NMAT
    do i=1,NMAT
      DF_eta(i,j,s) = -dconjg( tmp_eta(j,i,s) )
    enddo
  enddo
enddo
do l=1,num_links
  do j=1,NMAT
    do i=1,NMAT
      DF_lambda(i,j,l) = -dconjg( tmp_lambda(j,i,l) )
    enddo
  enddo
enddo
do f=1,num_faces
  do j=1,NMAT
    do i=1,NMAT
      DF_chi(i,j,f) = -dconjg( tmp_chi(j,i,f) )
    enddo
  enddo
enddo
end subroutine Prod_DiracDag






end module Dirac_operator
