!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! module for fermion matrix
module Dirac_operator
use global_parameters
use global_subroutines
use SUN_generators
use matrix_functions
implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! make Dirac matrix
subroutine make_Dirac(Dirac,UMAT,PhiMat)
use SUN_generators, only : trace_MTa
implicit none

complex(kind(0d0)), intent(inout) :: Dirac(1:sizeD,1:sizeD)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
!complex(kind(0d0)) :: Phi(1:dimG,1:num_sites)

complex(kind(0d0)) :: unitvec(1:sizeD),dirac_vec(1:sizeD)
integer i,j
integer s,a

!do s=1,num_sites
  !do a=1,dimG
    !call trace_MTa(Phi(a,s),PhiMat(:,:,s),a,NMAT)
  !enddo
!enddo
do i=1,sizeD
  unitvec=(0d0,0d0)
  unitvec(i)=(1d0,0d0)
  call prod_Dirac(Dirac(:,i),unitvec,sizeD,UMAT,PhiMat)
enddo

end subroutine make_Dirac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! make D^\dagger D matrix
subroutine make_DdagD(DdagD,UMAT,PhiMat)
implicit none

complex(kind(0d0)), intent(inout) :: DdagD(1:sizeD,1:sizeD)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)


complex(kind(0d0)) :: unitvec(1:sizeD),dirac_vec(1:sizeD)
integer i,j

do i=1,sizeD
  unitvec=(0d0,0d0)
  unitvec(i)=(1d0,0d0)
  call prod_DdagD(DdagD(:,i),unitvec,sizeD,UMAT,Phimat)
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
subroutine Prod_DdagD(DdagD_vec, vec, Vsize,UMAT,Phimat)
implicit none

integer, intent(in) :: Vsize !! vecsize must be sizeD
complex(kind(0d0)), intent(in) :: vec(1:Vsize)
complex(kind(0d0)), intent(inout) :: DdagD_vec(1:Vsize)
complex(kind(0d0)), intent(in) :: UMAT(:,:,:),PhiMat(:,:,:)

complex(kind(0d0)) :: tmpvec(1:Vsize),tmpvec2(1:Vsize)

integer :: i



call Prod_Dirac(tmpvec,vec,Vsize,UMAT,PhiMat)
do i=1,Vsize
  tmpvec2(i)=-dconjg(tmpvec(i))
enddo

call Prod_Dirac(tmpvec,tmpvec2,Vsize,UMAT,PhiMat)

do i=1,Vsize
  DdagD_vec(i)=dconjg(tmpvec(i))!*0.25d0 !! HERE!!
enddo


end subroutine Prod_DdagD



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compute D.vec
!!  S_f = 1/2 \Psi^T D \Psi
subroutine Prod_Dirac(D_vec, vec, Vsize,UMAT,PhiMat)
implicit none

integer, intent(in) :: Vsize !! vecsize must be sizeD
complex(kind(0d0)), intent(in) :: vec(1:Vsize)
complex(kind(0d0)), intent(inout) :: D_vec(1:Vsize)
complex(kind(0d0)), intent(in) :: UMAT(:,:,:)
complex(kind(0d0)), intent(in) :: PhiMat(:,:,:)

!complex(kind(0d0)), parameter :: mass_f=1d0

complex(kind(0d0)) :: eta_mat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: eta_ele(1:dimG,1:num_sites)
complex(kind(0d0)) :: lambda_mat(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: lambda_ele(1:dimG,1:num_links)
complex(kind(0d0)) :: chi_mat(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: chi_ele(1:dimG,1:num_faces)

complex(kind(0d0)) :: DF_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: DF_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: DF_chi(1:NMAT,1:NMAT,1:num_faces)



!complex(kind(0d0)) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: bPhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: r_site(1:dimG,1:num_sites)
complex(kind(0d0)) :: r_link(1:dimG,1:num_links)
complex(kind(0d0)) :: r_face(1:dimG,1:num_faces)


complex(kind(0d0)) :: Uinv(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT,1:num_faces),Ufm(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: diff_Omega(1:NMAT,1:NMAT,1:dimG)

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


integer :: s,l,f
integer :: i,j,k
integer :: a,b,c
integer :: r

!! for test
complex(kind(0d0)) :: tmp_diff_Omega(1:NMAT,1:NMAT,1:dimG)
complex(kind(0d0)) :: tmp_diff_Omega2(1:NMAT,1:NMAT,1:dimG)
integer :: ii,jj

!! preparation
D_vec=(0d0,0d0)
r_site=(0d0,0d0)
r_link=(0d0,0d0)
r_face=(0d0,0d0)

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
  do i=1,NMAT
    do j=1,NMAT
  !call make_traceless_matrix_from_modes(PhiMat(:,:,s),NMAT,Phi(:,s))
  bPhiMat(i,j,s)=conjg(PhiMat(j,i,s))
  !call make_traceless_matrix_from_modes(bPhiMat(:,:,s),NMAT,dconjg(Phi(:,s)))
    enddo
  enddo
enddo

DF_eta=(0d0,0d0)
DF_lambda=(0d0,0d0)
DF_chi=(0d0,0d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! (0) fermion mass
do l=1,num_links
  do a=1,dimG
    r_link(a,l)=r_link(a,l) + (0d0,0.5d0)*mass_f & 
      * ( eta_ele( a, link_tip(l) ) - eta_ele( a, link_org(l) ) )
  enddo
enddo 

do s=1,num_sites
  do a=1,dimG
    do i=1,linkorg_to_s(s)%num_
      l=linkorg_to_s(s)%labels_(i)
      r_site(a,s)=r_site(a,s) - (0d0,0.5d0)*mass_f &
        * lambda_ele(a,l)
    enddo
    do i=1,linktip_from_s(s)%num_
      l=linktip_from_s(s)%labels_(i)
      r_site(a,s)=r_site(a,s) + (0d0,0.5d0)*mass_f &
        * lambda_ele(a,l)
    enddo
  enddo
enddo

do f=1,num_faces
  do a=1,dimG
    do i=1,links_in_f(f)%num_
      l=links_in_f(f)%link_labels_(i)
      r_face(a,f) = r_face(a,f) + (0d0,0.5d0)*mass_f &
        * lambda_ele(a,l)
    enddo
  enddo
enddo

do l=1,num_links
  do a=1,dimG
    do i=1,face_in_l(l)%num_       
      f=face_in_l(l)%label_(i)
      r_link(a,l) = r_link(a,l) - (0d0,0.5d0)*mass_f &
        * chi_ele(a,f) 
    enddo
  enddo
enddo

!! preparation
!do l=1,num_links
!  do i=1,NMAT
!    do j=1,NMAT
!      Uinv(i,j,l)=dconjg(UMAT(j,i,l))
!    enddo
!  enddo
!enddo


if ( p1 == 0 ) then 
    !write(*,*) p1
!! (1) site action 
do s=1,num_sites
  call matrix_commutator(tmpmat1,PhiMat(:,:,s),eta_mat(:,:,s))
  DF_eta(:,:,s)=DF_eta(:,:,s) &
      +dcmplx(alpha_s(s))*(-0.5d0,0d0) *overall_factor*tmpmat1
!  do r=1,NZF
!    a=NZF_index(1,r)
!    b=NZF_index(2,r)
!    c=NZF_index(3,r)
!    D_vec(site_index(a,s))=D_vec(site_index(a,s)) & 
!      +dcmplx(alpha_s(s))*(0d0,-0.5d0)*NZF_value(r) &
!        *Phi(b,s)*vec(site_index(c,s))*overall_factor
!  enddo
enddo

endif



if (p2==0) then
    !write(*,*) p2
!! (2) link action 1
do s=1,num_sites
  tmpmat2=(0d0,0d0)
  do k=1,linkorg_to_s(s)%num_
    ! tmpmat2 = -i \sum_l( \alpha_l U_l^{-1}.\lambda_l.U_l )
    l=linkorg_to_s(s)%labels_(k)
    call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
      lambda_mat(:,:,l), NMAT, &
      UMAT(:,:,l), NMAT, &
      (0d0,0d0), tmpmat1, NMAT)
    call ZGEMM('C','N',NMAT,NMAT,NMAT,(0d0,-1d0)*dcmplx(alpha_l(l)), &
      UMAT(:,:,l), NMAT, &
      tmpmat1, NMAT, &
      (1d0,0d0), tmpmat2, NMAT)  
  enddo
  DF_eta(:,:,s)=DF_eta(:,:,s) &
    + overall_factor * tmpmat2
  !do a=1,dimG
    !call trace_MTa(trace,tmpmat2,a,NMAT)
    !D_vec(site_index(a,s))=D_vec(site_index(a,s)) &
      !+ trace *overall_factor
  !enddo
  !!!
  do k=1,linktip_from_s(s)%num_
    l=linktip_from_s(s)%labels_(k)
    DF_eta(:,:,s)=DF_eta(:,:,s) &
      + (0d0,1d0) * dcmplx(alpha_l(l)) * overall_factor * lambda_mat(:,:,l)
  enddo
enddo

do l=1,num_links
  !tmpmat2=i alpha_l U_l \lambda_l U_l^{-1}
  s=link_tip(l)
  call ZGEMM('N','C',NMAT,NMAT,NMAT,(1d0,0d0), &
    eta_mat(:,:,s), NMAT, &
    UMAT(:,:,l), NMAT, &
    (0d0,0d0), tmpmat1, NMAT)
  call ZGEMM('N','N',NMAT,NMAT,NMAT,(0d0,1d0)*dcmplx(alpha_l(l)), &
    UMAT(:,:,l), NMAT, &
    tmpmat1, NMAT, &
    (0d0,0d0), tmpmat2, NMAT)
  DF_lambda(:,:,l) = DF_lambda(:,:,l) &
    + overall_factor * tmpmat2
  s=link_org(l)
  DF_lambda(:,:,l) = DF_lambda(:,:,l) &
      - (0d0,1d0) * dcmplx(alpha_l(l)) * overall_factor * eta_mat(:,:,s)
enddo
endif

    
if(p3==0) then
!! (3) link action 2 
do l=1,num_links
  ! compute Ul.bPhi_tip(l).Ul^\dagger + bPhi_org(l)
  tmpmat2=bPhiMat(:,:,link_org(l))
  s=link_tip(l) 
  call ZGEMM('N','C',NMAT,NMAT,NMAT,(1d0,0d0), &
    bPhiMat(:,:,s), NMAT, &
    UMAT(:,:,l), NMAT, &
    (0d0,0d0), tmpmat1, NMAT)
  call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
    UMAT(:,:,l), NMAT, &
    tmpmat1, NMAT, &
    (1d0,0d0), tmpmat2, NMAT)
  !!!
  ! compute the commutator 
  call Matrix_Commutator(comm,tmpmat2,lambda_mat(:,:,l))
  !!
  DF_lambda(:,:,l)=DF_lambda(:,:,l) &
      + dcmplx(alpha_l(l)) * overall_factor * comm
enddo
endif

if(p4==0) then
!! (4) face action 1
do f=1,num_faces
  s=sites_in_f(f)%label_(1)
  call matrix_commutator(comm,PhiMat(:,:,s),chi_mat(:,:,f))
  DF_chi(:,:,f)=DF_chi(:,:,f) &
      + (-2d0,0d0) * dcmplx(alpha_f(f)) * overall_factor * comm
enddo
endif

if( p5==0 ) then
!! (5) face action 2
do f=1,num_faces
  call Make_face_variable(Uf(:,:,f),f,UMAT) 
  call matrix_power(Ufm(:,:,f),Uf(:,:,f),m_omega)
  call calc_sinU_and_cosUinv(sinU(:,:,f),cosUinv(:,:,f),Ufm(:,:,f))
enddo

do f=1,num_faces
  do i=1,links_in_f(f)%num_
    l=links_in_f(f)%link_labels_(i)
    line1=(0d0,0d0)
    line2=(0d0,0d0)
    do k=1,m_omega
      call calc_Amat(Amat,f,i,k,Uf(:,:,f),UMAT)
      call calc_Bmat(Bmat,f,i,k,Uf(:,:,f),UMAT)

      ! tmpmat2=A.lambda.B
      call matrix_product(tmpmat1,Amat,lambda_mat(:,:,l))
      call matrix_product(tmpmat2,tmpmat1,Bmat)
      ! tmpmat3=B^dag.lambda.A^dag 
      call matrix_product(tmpmat1,Bmat,lambda_mat(:,:,l),'C','N')
      call matrix_product(tmpmat3,tmpmat1,Amat,'N','C') 

      line1=line1+tmpmat2+tmpmat3
      line2=line2+tmpmat2-tmpmat3
    enddo

    ! line2 = cosUinv.(A.lambda.B-Bdag.lambda.Adag).cosUinv
    call matrix_product(tmpmat1,cosUinv(:,:,f),line2)
    call matrix_product(line2,tmpmat1,cosUinv(:,:,f))

    ! line1 = {A.lambda.B+B^dag.lambda.A^dag , (Um+Uminv)^{-1}
    tmpmat1=line1
    call matrix_AntiCommutator(line1,tmpmat1,cosUinv(:,:,f))
    ! line2 = {cosUinv.(A.lambda.B-Bdag.lambda.Adag).cosUinv, sinU}
    tmpmat1=line2
    call matrix_AntiCommutator(line2,tmpmat1,sinU(:,:,f))
    
    DF_chi(:,:,f)=DF_chi(:,:,f)&
      + cmplx(dble(links_in_f(f)%link_dirs_(i)))*(0d0,1d0)&
        * cmplx(overall_factor) / cmplx(dble(m_omega))&
        * cmplx(alpha_f(f)*beta_f(f)) * (line1-line2)
  enddo
enddo
!  do a=1,dimG
!    call trace_MTa(trace,DF_chi(:,:,f),a,NMAT)
!    D_vec(face_index(a,f))=D_vec(face_index(a,f))+trace
!  enddo
!enddo

!do f=1,num_faces
!  tmpmat1=(0d0,0d0)
!  do i=1,links_in_f(f)%num_
!    l=links_in_f(f)%link_labels_(i)
!    call calc_diff_omega(diff_Omega(:,:,:),Uf(:,:,f),Ufm(:,:,f),f,l,UMAT)
!    
!    do b=1,dimG
!      tmpmat1=tmpmat1+lambda_ele(b,l)*diff_Omega(:,:,b)
!    enddo
!  enddo  
!  do a=1,dimG
!    call trace_MTa(trace,tmpmat1,a,NMAT)
!    D_vec(face_index(a,f))=D_vec(face_index(a,f)) &
!      +(0d0,1d0)*dcmplx(alpha_f(f)*beta_f(f))*trace*overall_factor
!  enddo
!enddo

do l=1,num_links
  do i=1,face_in_l(l)%num_       
    f=face_in_l(l)%label_(i)
    ! j: position of the link l in the face l
    do j=1,links_in_f(f)%num_
      if ( l == links_in_f(f)%link_labels_(j) ) exit
    enddo

    ! tmpmat1= { (Uf+Ufinv)^{-1} , \chi_f }
    call matrix_anticommutator(tmpmat1,cosUinv(:,:,f),chi_mat(:,:,f))

    ! tmpmat2= (Uf+Ufinv)^{-1}.{ Uf-Ufinv , \chi_f }.(Uf+Ufinv)^{-1} 
    call matrix_anticommutator(tmpmat2,sinU(:,:,f),chi_mat(:,:,f))
    call matrix_product(tmpmat3,cosUinv(:,:,f),tmpmat2)
    call matrix_product(tmpmat2,tmpmat3,cosUinv(:,:,f))
    line1=tmpmat1-tmpmat2
    line2=tmpmat1+tmpmat2

    acomm=(0d0,0d0)
    do k=1,m_omega
      call calc_Amat(Amat,f,j,k,Uf(:,:,f),UMAT)
      call calc_Bmat(Bmat,f,j,k,Uf(:,:,f),UMAT)
      ! tmpmat1= B.(tmpmat1-tmpmat2).A
      call matrix_product(tmpmat3,BMAT,line1)
      call matrix_product(tmpmat1,tmpmat3,AMAT)
      ! tmpmat2 = Adag.(tmpmat1+tmpmat2).Bdag
      call matrix_product(tmpmat3,AMAT,line2,'C','N')
      call matrix_product(tmpmat2,tmpmat3,BMAT,'N','C')

      acomm=acomm+tmpmat1+tmpmat2
    enddo

    DF_lambda(:,:,l)=DF_lambda(:,:,l) &
      + cmplx(dble(links_in_f(f)%link_dirs_(j)))*(0d0,1d0)&
        * cmplx(-overall_factor) / cmplx(dble(m_omega))&
        * cmplx(alpha_f(f)*beta_f(f)) * acomm
  enddo
enddo


!    call calc_diff_omega(diff_Omega(:,:,:),Uf(:,:,f),Ufm(:,:,f),f,l,UMAT)
!    do a=1,dimG
!      tmp=(0d0,0d0)
!      do j=1,NMAT
!        do k=1,NMAT
!          tmp=tmp+chi_mat(j,k,f)*diff_Omega(k,j,a)
!        enddo
!      enddo
!
!      D_vec(link_index(a,l))=D_vec(link_index(a,l)) &
!        +(0d0,-1d0)*dcmplx(alpha_f(f)*beta_f(f))*tmp*overall_factor
!    enddo
  !enddo
!enddo
endif

do s=1,num_sites
  do a=1,dimG
    call trace_MTa(trace,DF_eta(:,:,s),a,NMAT)
    D_vec(site_index(a,s))=D_vec(site_index(a,s))+trace
  enddo
enddo
do l=1,num_links
  do a=1,dimG
    call trace_MTa(trace,DF_lambda(:,:,l),a,NMAT)
    D_vec(link_index(a,l))=D_vec(link_index(a,l))+trace
  enddo
enddo
do f=1,num_faces
  do a=1,dimG
    call trace_MTa(trace,DF_chi(:,:,f),a,NMAT)
    D_vec(face_index(a,f))=D_vec(face_index(a,f))+trace
  enddo
enddo
!! (T1) test action
!do f=1,num_faces
!  call Make_face_variable(Uf(:,:,f),f,UMAT) 
!  call matrix_power(NMAT,Uf(:,:,f),m_omega,Ufm(:,:,f))
!enddo
!
!do f=1,num_faces
!  tmpmat1=(0d0,0d0)
!  do i=1,links_in_f(f)%num_
!    l=links_in_f(f)%link_labels_(i)
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    !! change here 
!    !call calc_diff_Uf(tmp_diff_Omega,Uf(:,:,f),f,l,UMAT) ! OK
!    call calc_diff_Ufm(tmp_diff_Omega,Uf(:,:,f),f,l,UMAT) 
!    !call calc_diff_SandC(tmp_diff_Omega,tmp_diff_Omega2,Uf(:,:,f),f,l,UMAT) 
!    !diff_Omega=tmp_diff_Omega
!    do a=1,dimG
!      do ii=1,NMAT
!        do jj=1,NMAT
!          diff_Omega(ii,jj,a)&
!            =(tmp_diff_Omega(ii,jj,a)-conjg(tmp_diff_Omega(jj,ii,a)))
!        enddo
!      enddo
!    enddo
!
!    !do a=1,dimG
!    !  call ZGEMM('N','C',NMAT,NMAT,NMAT,(1d0,0d0), &
!    !    Ufm(:,:,f), NMAT, &
!    !    tmp_diff_Omega(:,:,a), NMAT, &
!    !    (0d0,0d0), tmpmat1, NMAT)
!    !  call ZGEMM('N','C',NMAT,NMAT,NMAT,(1d0,0d0), &
!    !    tmp_diff_Omega(:,:,a), NMAT, &
!    !    Ufm(:,:,f), NMAT, &
!    !    (1d0,0d0), tmpmat1, NMAT)
!    !  tmp=(0d0,0d0)
!    !  do ii=1,NMAT
!    !    do jj=1,NMAT
!    !      tmp=tmp+( tmpmat1(ii,jj)*dconjg(tmpmat1(ii,jj)) ) 
!    !    enddo
!    !  enddo
!    !  write(*,*) tmp
!    !enddo
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    do b=1,dimG
!      tmpmat1=tmpmat1+lambda_ele(b,l)*diff_Omega(:,:,b)
!    enddo
!  enddo  
!  do a=1,dimG
!    call trace_MTa(trace,tmpmat1,a,NMAT)
!    D_vec(face_index(a,f))=D_vec(face_index(a,f)) + trace
!  enddo
!enddo
!
!do l=1,num_links
!  do i=1,face_in_l(l)%num_       
!    f=face_in_l(l)%label_(i)
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    !! change here 
!    !call calc_diff_Uf(tmp_diff_Omega,Uf(:,:,f),f,l,UMAT) ! OK 
!    call calc_diff_Ufm(tmp_diff_Omega,Uf(:,:,f),f,l,UMAT)  
!    !call calc_diff_SandC(tmp_diff_Omega,tmp_diff_Omega2,Uf(:,:,f),f,l,UMAT) 
!    !diff_Omega=tmp_diff_Omega2
!    do a=1,dimG
!      do ii=1,NMAT
!        do jj=1,NMAT
!          diff_Omega(ii,jj,a)&
!            =(tmp_diff_Omega(ii,jj,a))-conjg(tmp_diff_Omega(jj,ii,a))
!        enddo
!      enddo
!    enddo
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    do a=1,dimG
!      tmp=(0d0,0d0)
!      do j=1,NMAT
!        do k=1,NMAT
!          tmp=tmp+chi_mat(j,k,f)*diff_Omega(k,j,a)
!        enddo
!      enddo
!
!      D_vec(link_index(a,l))=D_vec(link_index(a,l)) - tmp
!    enddo
!  enddo
!enddo
!!! regularize
!!   If there is no hole, sizeD is even. 
do i=1,sizeD,2
  D_vec(i)=D_vec(i) + dcmplx(mass_f)*vec(i+1)
  D_vec(i+1)=D_vec(i+1) - dcmplx(mass_f)*vec(i)
enddo
!D_vec=D_vec*overall_factor
!write(*,*) D_vec
end subroutine Prod_Dirac


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! check anti-symmetricity of D
subroutine check_Dirac(UMAT,PhiMat)
implicit none 

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)

complex(kind(0d0)) :: PF(1:sizeD)
complex(kind(0d0)) :: Dirac(1:sizeD,1:sizeD), Det
double precision :: rtmp,dtmp
integer :: seed
integer :: i,j
integer :: s,l,a,b


!seed=12345
!call genrand_init( put=seed )
! produce pseudo-fermion
!call make_pseudo_fermion(PF,UMAT,Phi)
!! calculate Hamiltonian 

call make_Dirac(Dirac,UMAT,PhiMat)

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

call make_DdagD(Dirac,UMAT,PhiMat)
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

end subroutine check_Dirac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate Amat used in fermionic action and force
subroutine calc_Amat(Amat,f,l_label,k,Uf,UMAT)
use matrix_functions, only : matrix_power,matrix_product
implicit none

complex(kind(0d0)), intent(out) :: Amat(1:NMAT,1:NMAT)
integer, intent(in) :: f,l_label,k
complex(kind(0d0)), intent(in) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)

complex(kind(0d0)):: tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)):: tmpmat2(1:NMAT,1:NMAT)
character :: C1
integer :: l,i

l=links_in_f(f)%link_labels_(l_label)

! tmpmat=U_f^{k-1}
call matrix_power(Amat,Uf,k-1)


! tmpmat=U_f^{k-1} U_{l1}^{e1}...U_{l-1}^{el-1}
if ( l_label .ne. 1 ) then
  do i=1,l_label-1
    if ( links_in_f(f)%link_dirs_(i)==1) then
      C1='N'
    else
      C1='C'
    endif
    call matrix_product(tmpmat,Amat,UMAT(:,:,links_in_f(f)%link_labels_(i)),'N',c1)
    Amat=tmpmat
  enddo
endif

if ( links_in_f(f)%link_dirs_(l_label)==-1 ) then
  call matrix_product(tmpmat,Amat,UMAT(:,:,l),'N','C')
  Amat=tmpmat
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

complex(kind(0d0)):: tmpmat(1:NMAT,1:NMAT)
character :: C1
integer :: l,i

l=links_in_f(f)%link_labels_(l_label)

! tmpmat=U_f^{m-k-1}
call matrix_power(Bmat,Uf,m_omega-k)


! tmpmat=U_f^{k-1} U_{l1}^{e1}...U_{l-1}^{el-1}
if ( l_label .ne. links_in_f(f)%num_ ) then
  do i=links_in_f(f)%num_ , l_label+1, -1
    if ( links_in_f(f)%link_dirs_(i)==1) then
      C1='N'
    else
      C1='C'
    endif
    call matrix_product(tmpmat,UMAT(:,:,links_in_f(f)%link_labels_(i)),BMAT,C1,'N')
    Bmat=tmpmat
  enddo
endif

if ( links_in_f(f)%link_dirs_(l_label)==1 ) then
  call matrix_product(tmpmat,UMAT(:,:,l),BMAT,'N','N')
  Bmat=tmpmat
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

end module Dirac_operator
