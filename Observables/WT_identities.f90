!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compute \Xi Tr(Phi eta)
subroutine calc_XiPhiEta(XiPhiEta, &
    Geta_eta, Glambda_eta, Gchi_eta, &
    Geta_lambda, Glambda_lambda, Gchi_lambda, &
    Geta_chi, Glambda_chi, Gchi_chi, &
    Umat,PhiMat)
use parallel
use matrix_functions, only : trace_MM
implicit none

complex(kind(0d0)), intent(out) :: XiPhiEta
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) ::  Geta_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_sites) 
complex(kind(0d0)) , intent(in) :: Glambda_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_sites) 
complex(kind(0d0)) , intent(in) :: Gchi_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_sites) 
!!!
complex(kind(0d0)) , intent(in) :: Geta_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_links) 
complex(kind(0d0)) , intent(in) :: Glambda_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_links) 
complex(kind(0d0)) , intent(in) :: Gchi_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_links) 
!!!
complex(kind(0d0)) , intent(in) :: Geta_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_faces)  
complex(kind(0d0)) , intent(in) :: Glambda_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_faces) 
complex(kind(0d0)) , intent(in) :: Gchi_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_faces) 

complex(kind(0d0)) Xi_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) Xi_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) Xi_chi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)) DinvPhi_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) DinvPhi_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) DinvPhi_chi(1:NMAT,1:NMAT,1:num_faces)

complex(kind(0d0)) eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) chi(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) tmp,trace
integer :: ls, ll, lf

call make_XiVec(Xi_eta,Xi_lambda,Xi_chi,Umat,PhiMat)

call DinvPF_direct(&
  DinvPhi_eta, DinvPhi_lambda, DinvPhi_chi, &
  Xi_eta(:,:,1:num_sites),Xi_lambda(:,:,1:num_links),Xi_chi(:,:,1:num_faces),&
  Geta_eta, Glambda_eta, Gchi_eta, &
  Geta_lambda, Glambda_lambda, Gchi_lambda, &
  Geta_chi, Glambda_chi, Gchi_chi, &
  Umat,PhiMat)

tmp=(0d0,0d0)
do ls=1,num_sites
  call trace_MM(trace, PhiMat(:,:,ls),DinvPhi_eta(:,:,ls)) 
  tmp=tmp-trace
enddo

XiPHiEta=(0d0,0d0)
call MPI_REDUCE(tmp,XiPhiEta,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)

end subroutine calc_XiPhiEta


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! compute \Xi Tr(Phi eta)
!subroutine calc_XiPhiEta_org(XiPhiEta,Umat,PhiMat,triger)
!#ifdef PARALLEL
!use parallel
!#endif
!implicit none
!
!complex(kind(0d0)), intent(out) :: XiPhiEta
!complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
!integer, intent(in) :: triger
!
!complex(kind(0d0)) Xi_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
!complex(kind(0d0)) Xi_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)) Xi_chi(1:NMAT,1:NMAT,1:num_necessary_faces)
!complex(kind(0d0)) DinvXi_eta(1:NMAT,1:NMAT,1:num_sites)
!complex(kind(0d0)) DinvXi_lambda(1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)) DinvXi_chi(1:NMAT,1:NMAT,1:num_faces)
!
!complex(kind(0d0)) tmp
!integer :: i,j,s
!integer :: info
!
!call make_XiVec(Xi_eta,Xi_lambda,Xi_chi,Umat,PhiMat)
!info=0
!
!!call test_calc_DinvF(Xi_eta,Xi_lambda,Xi_chi,Umat,PhiMat)
!!stop
!if( triger == 0) then 
!  call calc_dinvf(dinvxi_eta,dinvxi_lambda,dinvxi_chi,xi_eta,xi_lambda,xi_chi,umat,phimat,info)
!else
!  call calc_DinvF_direct(DinvXi_eta,DinvXi_lambda,DinvXi_chi,Xi_eta,Xi_lambda,Xi_chi,Umat,PhiMat)
!endif
!
!XiPhiEta=(0d0,0d0)
!if( info==1 ) then
!#ifdef PARALLEL
!  if(MYRANK==0) then
!#endif
!  write(*,*) "# Dirac is singular"
!#ifdef PARALLEL
!  endif
!#endif
!  return
!endif
!
!tmp=(0d0,0d0)
!do s=1,num_sites
!  do j=1,NMAT
!    do i=1,NMAT
!      tmp=tmp + DinvXi_eta(i,j,s)*PhiMat(j,i,s)
!    enddo
!  enddo
!enddo
!#ifdef PARALLEL
!call MPI_REDUCE(tmp,XiPhiEta,1,MPI_DOUBLE_COMPLEX, &
!  MPI_SUM,0,MPI_COMM_WORLD,IERR)
!#else
!XiPhiEta=tmp
!#endif
!
!end subroutine calc_XiPhiEta_org 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to make S_eta, S_lambda, S_chi of
!!  \Zia = Tr(eta S_eta) + Tr(lambda S_lambda) + Tr(chi S_chi)
subroutine make_XiVec(Xi_eta,Xi_lambda,Xi_chi,Umat,PhiMat)
use matrix_functions, only : matrix_commutator,matrix_3_product
#ifdef PARALLEL
use parallel
#endif
implicit none

complex(kind(0d0)), intent(out) :: Xi_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(out) :: Xi_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(out) :: Xi_chi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT)
integer :: s,l,f
integer :: i,j

do s=1,num_sites
  call matrix_commutator(tmpmat,PhiMat(:,:,s),Phimat(:,:,s),'N','C')
  Xi_eta(:,:,s)=dcmplx(alpha_s(s)*0.25d0*overall_factor)*tmpmat
enddo

do l=1,num_links
  do j=1,NMAT
    do i=1,NMAT
      tmpmat(i,j)=-dconjg( PhiMat(j,i,link_org(l)) )
    enddo
  enddo
  call matrix_3_product(tmpmat,&
    Umat(:,:,l),PhiMat(:,:,link_tip(l)),Umat(:,:,l),&
    'N','C','C',(1d0,0d0),'ADD') 
  Xi_lambda(:,:,l)=(0d0,-1d0)*dcmplx(alpha_l(l)*overall_factor)*tmpmat
enddo

do f=1,num_faces
  call make_face_variable(Uf,f,Umat)
  if(m_omega==0) call make_moment_map0(tmpmat,Uf)
  if(m_omega==-1) call make_moment_map_adm(tmpmat,Uf)
  Xi_chi(:,:,f)=(0d0,-0.5d0)*dcmplx(alpha_f(f)*beta_f(f)*overall_factor)*tmpmat
enddo

#ifdef PARALLEL
  call syncronize_sites(Xi_eta)
  call syncronize_links(Xi_lambda)
  call syncronize_faces(Xi_chi)
#endif
end subroutine make_XiVec



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! (D D†~)^{-1} F
subroutine calc_DdagDinvF(X_eta,X_lambda,X_chi,eta,lambda,chi,Umat,Phimat,info)
use Dirac_operator, only :prod_DdagD
#ifdef PARALLEL
use parallel
use global_subroutines, only : syncronize_sites, syncronize_links, syncronize_faces
#endif
implicit none 

complex(kind(0d0)), intent(out) :: X_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(out) :: X_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(out) :: X_chi(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)), intent(in) :: eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: chi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
integer, intent(out) :: info

complex(kind(0d0)) :: r_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: r_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: r_chi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)) :: p_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: p_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: p_chi(1:NMAT,1:NMAT,1:num_necessary_faces)
!complex(kind(0d0)) :: s_eta(1:NMAT,1:NMAT,1:num_sites)
!complex(kind(0d0)) :: s_lambda(1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)) :: s_chi(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: tmp_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: tmp_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: tmp_chi(1:NMAT,1:NMAT,1:num_faces)

complex(kind(0d0)) :: rkrk,alpha_k,tmp_c1,beta_k
integer :: MaxIte=10000

integer :: ite,s,l,f


!! 初期配置
!Xvec=(0d0,0d0)
X_eta=(0d0,0d0)
X_lambda=(0d0,0d0)
X_chi=(0d0,0d0)
!r_vec = Bvec
r_eta=eta
r_lambda=lambda
r_chi=chi

p_eta=r_eta
p_lambda=r_lambda
p_chi=r_chi

info=0
!! iteration start
do ite=1,MaxIte
    ! (1) construct \alpha_k
  ! rkrk = (r_k, r_k)
  call InnerProd(rkrk, &
    r_eta, r_lambda, r_chi, &
    r_eta, r_lambda, r_chi)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!== DEPEND ON ProdMat YOUR MADE ==
  call Prod_DdagD(&
    tmp_eta,tmp_lambda,tmp_chi, &
    p_eta,p_lambda,p_chi, &
    Umat,PhiMat)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call InnerProd(tmp_c1, &
    p_eta,p_lambda,p_chi, &
    tmp_eta, tmp_lambda, tmp_chi)
  alpha_k = rkrk / dcmplx(dble(tmp_c1))
 
 ! (2) update Xvec and r_vec
  X_eta = X_eta + dcmplx(alpha_k)*p_eta(:,:,1:num_sites)
  X_lambda = X_lambda + dcmplx(alpha_k)*p_lambda(:,:,1:num_links)
  X_chi = X_chi + dcmplx(alpha_k)*p_chi(:,:,1:num_faces)

  r_eta(:,:,1:num_sites) = r_eta(:,:,1:num_sites) - dcmplx(alpha_k)*tmp_eta
  r_lambda(:,:,1:num_links) = r_lambda(:,:,1:num_links) - dcmplx(alpha_k)*tmp_lambda
  r_chi(:,:,1:num_faces) = r_chi(:,:,1:num_faces) - dcmplx(alpha_k)*tmp_chi
#ifdef PARALLEL
  call syncronize_sites(r_eta)
  call syncronize_faces(r_chi)
  call syncronize_links(r_lambda)
#endif
    
 ! (3) conservation check
  call InnerProd( tmp_c1, &
    r_eta(:,:,1:num_sites), r_lambda(:,:,1:num_links), r_chi(:,:,1:num_faces), &
    r_eta(:,:,1:num_sites), r_lambda(:,:,1:num_links), r_chi(:,:,1:num_faces))
  if ( dsqrt( dble(tmp_c1) ) < epsilon ) then
    return
  endif

! (4) update p_k --> p_{k+1}
!construct beta_k
  !call InnerProd(tmp_c1, r_vec, r_vec)
  beta_k = tmp_c1 / rkrk
  p_eta(:,:,1:num_sites) = r_eta(:,:,1:num_sites) + dcmplx(beta_k) * p_eta(:,:,1:num_sites)
  p_lambda(:,:,1:num_links) = r_lambda(:,:,1:num_links) + dcmplx(beta_k) * p_lambda(:,:,1:num_links)
  p_chi(:,:,1:num_faces) = r_chi(:,:,1:num_faces) + dcmplx(beta_k) * p_chi(:,:,1:num_faces)
#ifdef PARALLEL
  call syncronize_sites(p_eta)
  call syncronize_faces(p_chi)
  call syncronize_links(p_lambda)
#endif
enddo

info = 1
return

end subroutine calc_DdagDinvF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! check calc_DinvF
subroutine test_calc_DinvF(eta,lambda,chi,Umat,Phimat)
use Dirac_operator
use matrix_functions
implicit none

complex(kind(0d0)), intent(in) :: eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: chi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

integer :: info
complex(kind(0d0)) :: X_eta1(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: X_lambda1(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: X_chi1(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)) :: X_eta2(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: X_lambda2(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: X_chi2(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), allocatable :: Dinv(:,:),vec(:),DinvV(:)
integer :: sizeD,i,j

sizeD=(NMAT*NMAT-1)*(global_num_sites+global_num_links+global_num_faces)
allocate( Dinv(1:sizeD,1:sizeD) )
allocate( vec(1:sizeD) )
allocate( DinvV(1:sizeD) )
call make_Dirac(Dinv,Umat,PhiMat)
call Matrix_Inverse(Dinv)
call MPI_BCAST(Dinv, sizeD*sizeD, MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)

call localmat_to_globalvec(vec, eta,lambda,chi)
call MPI_BCAST(vec, sizeD, MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)

DinvV=(0d0,0d0)
do i=1,sizeD
  do j=1,sizeD
    DinvV(j) = DinvV(j) + Dinv(j,i)*vec(i)
  enddo
enddo
call globalvec_to_localmat(X_eta1,X_lambda1,X_chi1,DinvV)

call calc_DinvF(X_eta2,X_lambda2,X_chi2,eta,lambda,chi,Umat,Phimat,info)

do i=1, num_sites
  write(*,*) X_eta1(:,:,i)-X_eta2(:,:,i)
enddo
do i=1, num_links
  write(*,*) X_lambda1(:,:,i)-X_lambda2(:,:,i)
enddo
do i=1, num_faces
  write(*,*) X_chi1(:,:,i)-X_chi2(:,:,i)
enddo


end subroutine test_calc_DinvF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Dirac inverse by direct computation
subroutine calc_DinvF_direct(X_eta,X_lambda,X_chi,eta,lambda,chi,Umat,Phimat)
use Dirac_operator
use matrix_functions
implicit none

complex(kind(0d0)), intent(out) :: X_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(out) :: X_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(out) :: X_chi(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)), intent(in) :: eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: chi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)), allocatable :: Dinv(:,:),vec(:),DinvV(:)
complex(kind(0d0)) :: tmpX_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: tmpX_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: tmpX_chi(1:NMAT,1:NMAT,1:num_necessary_faces)
integer :: sizeD,i,j

sizeD=(NMAT*NMAT-1)*(global_num_sites+global_num_links+global_num_faces)
allocate( Dinv(1:sizeD,1:sizeD) )
allocate( vec(1:sizeD) )
allocate( DinvV(1:sizeD) )

call make_Dirac(Dinv,Umat,PhiMat)
call localmat_to_globalvec(vec, eta,lambda,chi)

if( MYRANK == 0 ) then 
  call Matrix_Inverse(Dinv)
  DinvV=(0d0,0d0)
  do j=1,sizeD
    do i=1,sizeD
      DinvV(j) = DinvV(j) + Dinv(j,i)*vec(i)
    enddo
  enddo
endif
call MPI_BCAST(DinvV, sizeD, MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
call globalvec_to_localmat(tmpX_eta,tmpX_lambda,tmpX_chi,DinvV)

do i=1,num_sites
  X_eta(:,:,i)=tmpX_eta(:,:,i)
enddo
do i=1,num_links
  X_lambda(:,:,i)=tmpX_lambda(:,:,i)
enddo
do i=1,num_faces
  X_chi(:,:,i)=tmpX_chi(:,:,i)
enddo


deallocate( Dinv,vec,DinvV)

end subroutine calc_DinvF_direct

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Dirac inverse by CG
subroutine calc_DinvF(X_eta,X_lambda,X_chi,eta,lambda,chi,Umat,Phimat,info)
use Dirac_operator, only : prod_Dirac,prod_DiracDag
#ifdef PARALLEL
use parallel
use global_subroutines, only : syncronize_sites, syncronize_links, syncronize_faces
#endif
implicit none 

complex(kind(0d0)), intent(out) :: X_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(out) :: X_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(out) :: X_chi(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)), intent(in) :: eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: chi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
integer, intent(out) :: info

complex(kind(0d0)) :: r_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: r_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: r_chi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)) :: p_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: p_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: p_chi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)) :: s_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: s_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: s_chi(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: Dp_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: Dp_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: Dp_chi(1:NMAT,1:NMAT,1:num_faces)

complex(kind(0d0)) :: rkrk,sksk,sksk2,DpDp,lambda_k,mu_k
integer :: MaxIte=50000

integer :: ite,s,l,f

!! 初期配置
!Xvec=(0d0,0d0)
X_eta=(0d0,0d0)
X_lambda=(0d0,0d0)
X_chi=(0d0,0d0)
!r_vec = Bvec
r_eta=eta
r_lambda=lambda
r_chi=chi
! p_vec = s_vec = D^\dag r_vec
call prod_DiracDag(p_eta(:,:,1:num_sites),p_lambda(:,:,1:num_links),p_chi(:,:,1:num_faces),r_eta,r_lambda,r_chi,Umat,Phimat)
#ifdef PARALLEL
  call syncronize_sites(p_eta)
  call syncronize_faces(p_chi)
  call syncronize_links(p_lambda)
#endif
s_eta=p_eta(:,:,1:num_sites)
s_lambda=p_lambda(:,:,1:num_links)
s_chi=p_chi(:,:,1:num_faces)

!! preparation
call InnerProd(sksk, &
 s_eta, s_lambda, s_chi, &
 s_eta, s_lambda, s_chi)

info=0
!! iteration start
do ite=1,MaxIte
  ! (1) construct lambda_k 
  ! rkrk = (r_k, r_k)
  !call InnerProd(rkrk, r_vec, r_vec)
  call Prod_Dirac(Dp_eta,Dp_lambda,Dp_chi, &
    p_eta,p_lambda,p_chi,Umat,PhiMat)
  call InnerProd(DpDp, &
    Dp_eta, Dp_lambda, Dp_chi, &
    Dp_eta, Dp_lambda, Dp_chi)
  lambda_k = dexp( dlog(dble(sksk)) - dlog(dble(DpDp)) ) 
   
  X_eta=X_eta+dcmplx(lambda_k)*p_eta(:,:,1:num_sites)
  X_lambda=X_lambda+dcmplx(lambda_k)*p_lambda(:,:,1:num_links)
  X_chi=X_chi+dcmplx(lambda_k)*p_chi(:,:,1:num_faces)

 
 ! (2) update r_vec
  !r_vec = r_vec - dcmplx(alpha_k)*tmp_vec
  r_eta(:,:,1:num_sites)=r_eta(:,:,1:num_sites)-dcmplx(lambda_k)*Dp_eta
  r_lambda(:,:,1:num_links)=r_lambda(:,:,1:num_links)-dcmplx(lambda_k)*Dp_lambda
  r_chi(:,:,1:num_faces)=r_chi(:,:,1:num_faces)-dcmplx(lambda_k)*Dp_chi
#ifdef PARALLEL
  call syncronize_sites(r_eta)
  call syncronize_faces(r_chi)
  call syncronize_links(r_lambda)
#endif
    
 ! (3) conservation check
  !call InnerProd( tmp_c1, r_vec, r_vec)
  call InnerProd(rkrk, &
    r_eta(:,:,1:num_sites),r_lambda(:,:,1:num_links),r_chi(:,:,1:num_faces),&
    r_eta(:,:,1:num_sites),r_lambda(:,:,1:num_links),r_chi(:,:,1:num_faces))

  if ( dsqrt( dble(rkrk) ) < 1d-6 ) then
    return
  endif

! (4) update p_k --> p_{k+1}
!construct beta_k
  !call InnerProd(tmp_c1, r_vec, r_vec)
  !beta_k = dble(tmp_c1) / dble(rkrk)
  !p_vec = r_vec + dcmplx(beta_k) * p_vec

  call prod_DiracDag(s_eta,s_lambda,s_chi,r_eta,r_lambda,r_chi,Umat,Phimat)

  call InnerProd(sksk2, &
   s_eta, s_lambda, s_chi, &
   s_eta, s_lambda, s_chi)
  mu_k = dcmplx( dexp( dlog(dble(sksk2)) - dlog(dble(sksk)) ) ) 

  !! skskを更新
  sksk=sksk2
  !! p_vecを更新
  p_eta(:,:,1:num_sites) = s_eta(:,:,1:num_sites) + dcmplx(mu_k) * p_eta(:,:,1:num_sites)
  p_lambda(:,:,1:num_links) = s_lambda(:,:,1:num_links) + dcmplx(mu_k) * p_lambda(:,:,1:num_links)
  p_chi(:,:,1:num_faces) = s_chi(:,:,1:num_faces) + dcmplx(mu_k) * p_chi(:,:,1:num_faces)
#ifdef PARALLEL
  call syncronize_sites(p_eta)
  call syncronize_faces(p_chi)
  call syncronize_links(p_lambda)
#endif
enddo

info = 1
return
end subroutine calc_DinvF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate WT identity 1 and 2 and take average over the same latitude
!! *** ASSUMING MISUMI-BUNKATSU ***
subroutine calc_WT12_average(WT1,WT2,Umat,PhiMat,fop,sizeM,sizeN)
implicit none

complex(kind(0d0)), intent(out) :: WT1(:), WT2(:)
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
integer, intent(in) :: fop ! global face where \phi^{N_c^2-1} sits
integer, intent(in) :: sizeM, sizeN
complex(kind(0d0)) :: tmp1, tmp2

integer :: m,n,gf

if( size(WT1,1) .ne. sizeN/2-1 ) then 
  write(*,*) "WT1 must be WT1(1:sizeN/2-1)"
  call stop_for_test
endif
if( size(WT2,1) .ne. sizeN/2-1 ) then 
  write(*,*) "WT2 must be WT1(1:sizeN/2-1)"
  call stop_for_test
endif

WT1=(0d0,0d0)
WT2=(0d0,0d0)
! (n-1)*M+2 ... n*M+1
do n=1, sizeN/2-1
  do m=1, sizeM
    gf=(n-1)*sizeM + m + 1
    call calc_WT12(tmp1,tmp2,Umat,PhiMat,fop,gf)
    WT1(n)=WT1(n)+tmp1
    WT2(n)=WT2(n)+tmp2
  enddo
  WT1(n)=WT1(n)/dcmplx(dble(sizeM))
  WT2(n)=WT2(n)/dcmplx(dble(sizeM))
enddo


end subroutine calc_WT12_average


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
subroutine calc_WT12(WT1,WT2,Umat,PhiMat,fop,fcurr)
use global_subroutines, only : &
  make_face_variable, &
  make_moment_map0, &
  make_diff_PhiMat
use matrix_functions, only : &
  matrix_commutator, &
  matrix_power, &
  trace_MM
use parallel
implicit none

complex(kind(0d0)), intent(out) :: WT1 ! f: position of current
complex(kind(0d0)), intent(out) :: WT2 ! f: position of current
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
integer, intent(in) :: fop ! global face where \phi^{N_c^2-1} sits
integer, intent(in) :: fcurr ! global face where the current sits

complex(kind(0d0)) :: F12(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)) :: DPhi(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: CommPhi(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)) :: op(1:NMAT,1:NMAT)

complex(kind(0d0)) :: eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: chi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)) :: WTeta1(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: WTlambda1(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: WTchi1(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)) :: WTeta2(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: WTlambda2(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: WTchi2(1:NMAT,1:NMAT,1:num_necessary_faces)

complex(kind(0d0)) :: op_in(1:NMAT,1:NMAT) ! local operator at fop
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT) 
complex(kind(0d0)) :: tmp,trace

integer :: i,j,s,l,f
integer :: fop_local,fop_rank,s0_local,f_local,f_rank,info
complex(kind(0d0)) :: tmp1, tmp2

WT1=(0d0,0d0)
WT2=(0d0,0d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Preparation
do f=1,num_faces
  call make_face_variable(Uf,f,Umat)
  call make_moment_map0(F12(:,:,f),Uf)
  F12(:,:,f)=F12(:,:,f)*(0d0,-0.5d0)*beta_f(f)
enddo

do l=1,num_links
  call Make_diff_PhiMat(DPhi(:,:,l),l,Umat,PhiMat)
enddo

do s=1,num_sites
  call matrix_commutator(CommPhi(:,:,s),PhiMat(:,:,s),PhiMat(:,:,s),'N','C')
enddo

call syncronize_faces(F12)
call syncronize_links(DPhi)
call syncronize_sites(CommPhi)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fop_local=local_face_of_global(fop)%label_
fop_rank=local_face_of_global(fop)%rank_
s0_local=sites_in_f(fop_local)%label_(1)

op_in=(0d0,0d0)
if( MYRANK==fop_rank ) then
  do i=1,NMAT
    do j=1,NMAT
      tmpmat(i,j)=dconjg( PhiMat(j,i,s0_local) )
    enddo
  enddo
  call matrix_power(op_in,tmpmat,NMAT*NMAT-1)
endif

call make_div34rot12(eta,lambda,chi,F12,Dphi,CommPhi,Umat,PhiMat,fcurr)
call calc_DinvF(WTeta1,WTlambda1,WTchi1,eta,lambda,chi,Umat,PhiMat,info)

call make_div12rot34(eta,lambda,chi,F12,Dphi,CommPhi,Umat,PhiMat,fcurr)
call calc_DinvF(WTeta2,WTlambda2,WTchi2,eta,lambda,chi,Umat,PhiMat,info)

!!!!!!!!!!!!!!!!!!!!!!!!!!!
tmp=(0d0,0d0)
if( MYRANK == fop_rank ) then
  call trace_MM(trace, op_in, WTchi1(:,:,fop_local))
  tmp=tmp + trace

  call trace_MM(trace, op_in, WTeta2(:,:,s0_local))
  tmp=tmp + (0.5d0,0d0)*trace
endif
call MPI_REDUCE(tmp,WT1,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)

!!!!!!!!!!!!!!!!!!!!!!!!!!!
tmp=(0d0,0d0)
if( MYRANK == fop_rank ) then
  call trace_MM(trace, op_in, WTchi2(:,:,fop_local))
  tmp=tmp + trace
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fop_local=local_face_of_global(fop)%label_
fop_rank=local_face_of_global(fop)%rank_
s0_local=sites_in_f(fop_local)%label_(1)

op_in=(0d0,0d0)
if( MYRANK==fop_rank ) then
  do i=1,NMAT
    do j=1,NMAT
      tmpmat(i,j)=dconjg( PhiMat(j,i,s0_local) )
    enddo
  enddo
  call matrix_power(op_in,tmpmat,NMAT*NMAT-1)
endif

call make_div34rot12(eta,lambda,chi,F12,Dphi,CommPhi,Umat,PhiMat,fcurr)
call calc_DinvF(WTeta1,WTlambda1,WTchi1,eta,lambda,chi,Umat,PhiMat,info)

call make_div12rot34(eta,lambda,chi,F12,Dphi,CommPhi,Umat,PhiMat,fcurr)
call calc_DinvF(WTeta2,WTlambda2,WTchi2,eta,lambda,chi,Umat,PhiMat,info)


  call trace_MM(trace, op_in, WTeta1(:,:,s0_local))
  tmp=tmp - (0.5d0,0d0)*trace
endif
call MPI_REDUCE(tmp,WT2,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)


end subroutine calc_WT12

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!subroutine test_WT1(WT1,Umat,PhiMat,fop,fcurr,option)
!use global_subroutines, only : &
!  make_face_variable, &
!  make_moment_map0, &
!  make_diff_PhiMat
!use matrix_functions, only : &
!  matrix_commutator, &
!  matrix_power, &
!  trace_MM
!use parallel
!implicit none
!
!complex(kind(0d0)), intent(out) :: WT1 ! f: position of current
!complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
!integer, intent(in) :: fop ! global face where \phi^{N_c^2-1} sits
!integer, intent(in) :: fcurr ! global face where the current sits
!integer, intent(in) :: option
!
!complex(kind(0d0)) :: F12(1:NMAT,1:NMAT,1:num_necessary_faces)
!complex(kind(0d0)) :: DPhi(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)) :: CommPhi(1:NMAT,1:NMAT,1:num_necessary_sites)
!complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: op(1:NMAT,1:NMAT)
!
!complex(kind(0d0)) :: eta(1:NMAT,1:NMAT,1:num_necessary_sites)
!complex(kind(0d0)) :: lambda(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)) :: chi(1:NMAT,1:NMAT,1:num_necessary_faces)
!complex(kind(0d0)) :: WTeta1(1:NMAT,1:NMAT,1:num_necessary_sites)
!complex(kind(0d0)) :: WTlambda1(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)) :: WTchi1(1:NMAT,1:NMAT,1:num_necessary_faces)
!
!complex(kind(0d0)) :: op_in(1:NMAT,1:NMAT) ! local operator at fop
!complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT) 
!complex(kind(0d0)) :: tmp,trace
!
!integer :: i,j,s,l,f
!integer :: fop_local,fop_rank,s0_local,f_local,f_rank,info
!complex(kind(0d0)) :: tmp1, tmp2
!
!WT1=(0d0,0d0)
!!! Preparation
!do f=1,num_faces
!  call make_face_variable(Uf,f,Umat)
!  call make_moment_map0(F12(:,:,f),Uf)
!  F12(:,:,f)=F12(:,:,f)*(0d0,-0.5d0)*beta_f(f)
!enddo
!
!do l=1,num_links
!  call Make_diff_PhiMat(DPhi(:,:,l),l,Umat,PhiMat)
!enddo
!
!do s=1,num_sites
!  call matrix_commutator(CommPhi(:,:,s),PhiMat(:,:,s),PhiMat(:,:,s),'N','C')
!enddo
!
!call syncronize_faces(F12)
!call syncronize_links(DPhi)
!call syncronize_sites(CommPhi)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!fop_local=local_face_of_global(fop)%label_
!fop_rank=local_face_of_global(fop)%rank_
!s0_local=sites_in_f(fop_local)%label_(1)
!
!op_in=(0d0,0d0)
!if( MYRANK==fop_rank ) then
!  do i=1,NMAT
!    do j=1,NMAT
!      tmpmat(i,j)=dconjg( PhiMat(j,i,s0_local) )
!    enddo
!  enddo
!  call matrix_power(op_in,tmpmat,NMAT*NMAT-1)
!endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!eta=(0d0,0d0)
!lambda=(0d0,0d0)
!chi=(0d0,0d0)
!!!
!if( option== 1) call make_divK3(eta,lambda,chi,F12,Dphi,CommPhi,Umat,PhiMat,fcurr)
!if( option== 2) call make_divK4(eta,lambda,chi,F12,Dphi,CommPhi,Umat,PhiMat,fcurr)
!if( option== 3) call make_rotK1(eta,lambda,chi,F12,Dphi,CommPhi,Umat,PhiMat,fcurr)
!if( option== 4) call make_rotK2(eta,lambda,chi,F12,Dphi,CommPhi,Umat,PhiMat,fcurr)
!if( option== 5) call make_phichi(eta,lambda,chi,F12,Dphi,CommPhi,Umat,PhiMat,fcurr)
!!!
!call syncronize_sites(eta)
!call syncronize_links(lambda)
!call syncronize_faces(chi)
!!!
!call calc_DinvF(WTeta1,WTlambda1,WTchi1,eta,lambda,chi,Umat,PhiMat,info)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!tmp=(0d0,0d0)
!if( MYRANK == fop_rank ) then
!  call trace_MM(trace, op_in, WTchi1(:,:,fop_local))
!  tmp=tmp + trace
!endif
!call MPI_REDUCE(tmp,WT1,1,MPI_DOUBLE_COMPLEX, &
!  MPI_SUM,0,MPI_COMM_WORLD,IERR)
!
!end subroutine test_WT1



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! parts of WT identity in
!!  (div(K3+K4) + rot(K1+K2))
subroutine make_div34rot12(eta,lambda,chi,F12,Dphi,CommPhi,Umat,PhiMat,fcurr)
use parallel
implicit none

complex(kind(0d0)), intent(out) :: eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(out) :: lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(out) :: chi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), intent(in) :: F12(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), intent(in) :: DPhi(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: CommPhi(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: Phimat(1:NMAT,1:NMAT,1:num_necessary_sites)
integer, intent(in) :: fcurr ! global face

integer :: dir, Nfaces
integer :: s,l,f,i,j,k
integer :: rank
integer :: gs,gl,gf
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)

eta=(0d0,0d0)
lambda=(0d0,0d0)
chi=(0d0,0d0)

!! divergence part
do i=1, global_sites_in_f(fcurr)%num_
  gs=global_sites_in_f(fcurr)%label_(i)
  Nfaces=global_face_in_s(gs)%num_
  do j=1, global_linktip_from_s(gs)%num_
    gl=global_linktip_from_s(gs)%labels_(j)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    rank=local_link_of_global(gl)%rank_
    l=local_link_of_global(gl)%label_
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    if( MYRANK == rank ) then
      do k=1, face_in_l(l)%num_
        f=face_in_l(l)%label_(k)
        !write(*,*) MYRANK,"test1"
        call carry_MAT(tmpmat,F12(:,:,f),sites_in_f(f)%label_(1),link_org(l),f,Umat)
        lambda(:,:,l)=lambda(:,:,l) &
          + dcmplx(0.5d0 * alpha_l(l) / dble(Nfaces)) * tmpmat
      enddo
    endif
    do k=1, global_face_in_l(gl)%num_
      gf=global_face_in_l(gl)%label_(k)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      rank=local_face_of_global(gf)%rank_
      f=local_face_of_global(gf)%label_
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if( MYRANK == rank ) then
        do l=1,num_necessary_links
          if( global_link_of_local(l) == gl ) exit
        enddo
        !write(*,*) MYRANK,"test2"
        call carry_MAT(tmpmat,DPhi(:,:,l),link_org(l),sites_in_f(f)%label_(1),f,Umat)
        chi(:,:,f)=chi(:,:,f) &
          - dcmplx( 0.5d0 * alpha_l(l) / dble(Nfaces) ) * tmpmat
      endif
    enddo
  enddo
  !!!
  do j=1, global_linkorg_to_s(gs)%num_
    gl=global_linkorg_to_s(gs)%labels_(j)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    rank=local_link_of_global(gl)%rank_
    l=local_link_of_global(gl)%label_
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    if( MYRANK == rank ) then
      do k=1, face_in_l(l)%num_
        f=face_in_l(l)%label_(k)
        !write(*,*) MYRANK,"test3"
        call carry_MAT(tmpmat,F12(:,:,f),sites_in_f(f)%label_(1),link_org(l),f,Umat)
        lambda(:,:,l)=lambda(:,:,l) &
          - dcmplx(0.5d0 * alpha_l(l) / dble(Nfaces)) * tmpmat
      enddo
    endif
    do k=1, global_face_in_l(gl)%num_
      gf=global_face_in_l(gl)%label_(k)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      rank=local_face_of_global(gf)%rank_
      f=local_face_of_global(gf)%label_
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if( MYRANK == rank ) then
        do l=1,num_necessary_links
          if( global_link_of_local(l) == gl ) exit
        enddo
        !write(*,*) MYRANK,"test4"
        call carry_MAT(tmpmat,DPhi(:,:,l),link_org(l),sites_in_f(f)%label_(1),f,Umat)
        chi(:,:,f)=chi(:,:,f) &
          + dcmplx( 0.5d0 * alpha_l(l) / dble(Nfaces) ) * tmpmat
      endif
    enddo
  enddo
enddo

!! rotation part
do i=1, global_links_in_f(fcurr)%num_
  gl=global_links_in_f(fcurr)%link_labels_(i)
  dir=global_links_in_f(fcurr)%link_dirs_(i)
  if( dir==1 ) then 
    gs=global_link_org(gl)
  else
    gs=global_link_tip(gl)
  endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  rank=local_site_of_global(gs)%rank_
  s=local_site_of_global(gs)%label_
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if( MYRANK == rank ) then 
    do l=1,num_necessary_links
      if( global_link_of_local(l) == gl ) exit
    enddo
    eta(:,:,s)=eta(:,:,s) &
      + dcmplx( 0.5d0 * beta_f(f) * dble(dir) ) * DPhi(:,:,l)
  endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  rank=local_link_of_global(gl)%rank_
  l=local_link_of_global(gl)%label_
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if( MYRANK == rank ) then
    s=link_org(l)
    lambda(:,:,l) = lambda(:,:,l) & 
      + dcmplx( -0.5d0 * beta_f(f) * dble(dir) ) * (0d0,1d0) * CommPhi(:,:,s)
  endif
enddo

!! mass part
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
rank=local_face_of_global(fcurr)%rank_
f=local_face_of_global(fcurr)%label_
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if( MYRANK == rank ) then
  s=sites_in_f(f)%label_(1)
  chi(:,:,f)=chi(:,:,f) + dcmplx(mass_square_phi) * PhiMat(:,:,s)
endif

!!!!!!!!!!!!!!!!
call syncronize_sites(eta)
call syncronize_links(lambda)
call syncronize_faces(chi)

end subroutine make_div34rot12

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! parts of WT identity in
!!  div(-K1+K2) + rot(-K3+K4)
subroutine make_div12rot34(eta,lambda,chi,F12,Dphi,CommPhi,Umat,PhiMat,fcurr)
use parallel
implicit none

complex(kind(0d0)), intent(out) :: eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(out) :: lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(out) :: chi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), intent(in) :: F12(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), intent(in) :: DPhi(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: CommPhi(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: Phimat(1:NMAT,1:NMAT,1:num_necessary_sites)
integer :: fcurr 

integer :: dir, Nfaces, rank
integer :: s,l,f,i,j
integer :: gs,gl,gf
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)

eta=(0d0,0d0)
lambda=(0d0,0d0)
chi=(0d0,0d0)

!! divergence part
do i=1, global_sites_in_f(fcurr)%num_
  gs=global_sites_in_f(fcurr)%label_(i)
  Nfaces=global_face_in_s(gs)%num_
  do j=1, global_linktip_from_s(gs)%num_
    gl=global_linktip_from_s(gs)%labels_(j)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    rank=local_link_of_global(gl)%rank_
    l=local_link_of_global(gl)%label_
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if( MYRANK == rank ) then
      eta(:,:,link_org(l)) = eta(:,:,link_org(l)) &
        - dcmplx( 0.5d0 * alpha_l(l) / dble(Nfaces) ) * DPhi(:,:,l)
      lambda(:,:,l) = lambda(:,:,l) &
        - dcmplx( 0.5d0 * alpha_l(l) / dble(Nfaces) ) * (0d0,1d0) &
          * CommPhi(:,:,link_org(l))
    endif
  enddo
  !!!
  do j=1, global_linkorg_to_s(gs)%num_
    gl=global_linkorg_to_s(gs)%labels_(j)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    rank=local_link_of_global(gl)%rank_
    l=local_link_of_global(gl)%label_
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if( MYRANK == rank ) then
      eta(:,:,link_org(l)) = eta(:,:,link_org(l)) &
        + dcmplx( 0.5d0 * alpha_l(l) / dble(Nfaces) ) * DPhi(:,:,l)
      lambda(:,:,l) = lambda(:,:,l) &
        + dcmplx( 0.5d0 * alpha_l(l) / dble(Nfaces) ) * (0d0,1d0) &
          * CommPhi(:,:,link_org(l))
    endif
  enddo
enddo

!! rotation part
do i=1, global_links_in_f(fcurr)%num_
  gl=global_links_in_f(fcurr)%link_labels_(i)
  dir=global_links_in_f(fcurr)%link_dirs_(i)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  rank=local_link_of_global(gl)%rank_
  l=local_link_of_global(gl)%label_
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if( MYRANK == rank ) then
    do j=1,face_in_l(l)%num_
      f=face_in_l(l)%label_(j)
      call carry_MAT(tmpmat,F12(:,:,f),sites_in_f(f)%label_(1),link_org(l),f,Umat)
      lambda(:,:,l) = lambda(:,:,l) &
        - dcmplx( 0.5d0 * beta_f(f) * dble(dir) ) * tmpmat
    enddo
  endif
  do j=1, global_face_in_l(gl)%num_
    gf=global_face_in_l(gl)%label_(j)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    rank=local_face_of_global(gf)%rank_
    f=local_face_of_global(gf)%label_
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if( MYRANK == rank ) then
      do l=1, num_necessary_links
        if( global_link_of_local(l) == gl ) exit
      enddo
      call carry_MAT(tmpmat,DPhi(:,:,link_org(l)),link_org(l),sites_in_f(f)%label_(1),f,Umat)
      chi(:,:,f) = chi(:,:,f) & 
        - dcmplx( 0.5d0 * beta_f(f) * dble(dir) ) * tmpmat
    endif
  enddo
enddo

!! mass part
gs=global_sites_in_f(fcurr)%label_(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
rank=local_site_of_global(gs)%rank_
s=local_site_of_global(gs)%label_
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if( MYRANK == rank ) then
  eta(:,:,s)=eta(:,:,s) &
    + dcmplx(0.5d0 * mass_square_phi) * PhiMat(:,:,s)
endif

!!!!!!!!!!!!!!!!
call syncronize_sites(eta)
call syncronize_links(lambda)
call syncronize_faces(chi)

end subroutine make_div12rot34

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compute U_{s_carry...s_mat} MAT(s_mat) U_{s_mat...s_carry} in a face f
subroutine carry_MAT(MAT2,MAT1,s_mat,s_carry,face,Umat)
use global_subroutines, only : calc_prodUl_from_n1_to_n2_in_Uf
use matrix_functions, only : matrix_3_product
use parallel
implicit none

complex(kind(0d0)), intent(out) :: MAT2(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: MAT1(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
integer, intent(in) :: s_mat, s_carry, face

integer :: s_place, t_place
complex(kind(0d0)) :: Ucarry(1:NMAT,1:NMAT)


call find_s_in_f(s_place, s_carry, face)
call find_s_in_f(t_place, s_mat, face)

call calc_prodUl_from_n1_to_n2_in_Uf(Ucarry,face,s_place,t_place-1,Umat)
call matrix_3_product(MAT2,Ucarry,MAT1,Ucarry,'N','N','C')

end subroutine carry_MAT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! find place of site in face
subroutine find_s_in_f(s_place,site,face)
use parallel
implicit none

integer, intent(out) :: s_place
integer, intent(in) :: site, face

integer :: i

!write(*,*) face, site, ":", sites_in_f(face)%label_
do i=1, sites_in_f(face)%num_
  if( sites_in_f(face)%label_(i) == site ) then 
    s_place=i
    return
  endif
enddo

write(*,*) "[find_s_in_f] site is not in the face"
call stop_for_test

end subroutine find_s_in_f


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate Tr( div(V_l) )
subroutine calc_trace_div(TrDiv, Vec)
use parallel
implicit none

complex(kind(0d0)), intent(out) :: TrDiv(1:num_faces)
complex(kind(0d0)), intent(in) :: Vec(1:NMAT,1:NMAT,1:num_necessary_links)

complex(kind(0d0)) :: trace
integer :: s,l,f,i,j,a

TrDiv=(0d0,0d0)
do f=1,num_faces
  do i=1,sites_in_f(f)%num_
    s=sites_in_f(f)%label_(i)
    trace=(0d0,0d0)
    do j=1,linktip_from_s(s)%num_
      l=linktip_from_s(s)%labels_(j)
      do a=1,NMAT
        trace=trace +  alpha_l(l)*Vec(a,a,l)
      enddo
    enddo
    do j=1,linkorg_to_s(s)%num_
      l=linkorg_to_s(s)%labels_(j)
      do a=1,NMAT
        trace=trace -  alpha_l(l)*Vec(a,a,l)
      enddo
    enddo
    TrDiv(f)=TrDiv(f) + trace / dcmplx(dble(num_faces_in_s(s)))
  enddo
enddo

end subroutine calc_trace_div

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate Tr( rot(V_l) )
subroutine calc_trace_rot(TrRot, Vec)
use parallel
implicit none

complex(kind(0d0)), intent(out) :: TrRot(1:num_faces)
complex(kind(0d0)), intent(in) :: Vec(1:NMAT,1:NMAT,1:num_necessary_links)

complex(kind(0d0)) :: trace
integer :: s,l,f,i,j,a,b

TrRot=(0d0,0d0)
do f=1,num_faces
  do i=1,links_in_f(f)%num_
    l=links_in_f(f)%link_labels_(i)
    trace=(0d0,0d0)
    do a=1,NMAT
      trace=trace+Vec(a,a,l)
    enddo
    TrRot=TrRot + links_in_f(f)%link_dirs_(i) * trace
  enddo
  TrRot = TrRot * beta_f(f)
enddo

end subroutine calc_trace_rot


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine simplest_op(op,PhiMat,sop)
use matrix_functions, only : &
  matrix_power, &
  trace_MM
use parallel
implicit none

complex(kind(0d0)), intent(out) :: op ! f: position of current
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
integer, intent(in) :: sop ! global face where \phi^{N_c^2-1} sits

complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmp
integer :: sop_local,sop_rank,i,j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
sop_local=local_site_of_global(sop)%label_
sop_rank=local_site_of_global(sop)%rank_

op=(0d0,0d0)
tmp=(0d0,0d0)
if( MYRANK==sop_rank ) then
  do i=1,NMAT
    do j=1,NMAT
      tmpmat(i,j)=dconjg( PhiMat(j,i,sop_local) )
    enddo
  enddo
  call matrix_power(tmpmat2,tmpmat,NMAT*NMAT-1)
  do i=1,NMAT
    tmp=tmp+tmpmat2(i,i)
  enddo
endif

call MPI_REDUCE(tmp,op,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)


end subroutine simplest_op


