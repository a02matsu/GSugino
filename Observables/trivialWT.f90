!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! trivial WT identity
!!  <Sb> + \mu^2/2 \Xi \sum_s( \Phi_s \eta_s ) - (N^2-1)/2 (Ns+Nl)
subroutine calc_trivialWT(WT,Geta_eta,Geta_lambda,Geta_chi,Umat,PhiMat)
use matrix_functions, only : trace_MM
use parallel
implicit none

complex(kind(0d0)), intent(out) :: WT
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

double precision :: Sbr
complex(kind(0d0)) :: Sb


complex(kind(0d0)), intent(in) ::  Geta_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_sites) 
complex(kind(0d0)) , intent(in) :: Geta_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_links) 
complex(kind(0d0)) , intent(in) :: Geta_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_faces) 

complex(kind(0d0)) Xi_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) Xi_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) Xi_chi(1:NMAT,1:NMAT,1:num_necessary_faces)

complex(kind(0d0)) mass_cont
complex(kind(0d0)) Fnumber

complex(kind(0d0)) tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) DinvXi(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) trace
integer gs,ll,lf,ls,ls2,rank
integer k,l

!! (1) bosonic action
call calc_bosonic_action(Sbr,Umat,PhiMat)
Sb=dcmplx(Sbr)

!! (2) mass contribution
call mass_contribution(mass_cont,Geta_eta,Geta_lambda,Geta_chi,Umat,PhiMat)

!! (3) contribution from fermion number 
Fnumber=dcmplx(0.5d0*dble( (NMAT*NMAT-1)*(global_num_sites+global_num_links) ) )

if( MYRANK==0 ) then
  WT=Sb+mass_cont-Fnumber
endif

end subroutine calc_trivialWT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! mass contribution to the trivial WT identity 
!!  -1/2g^2 mu^2/2  Tr( \phi_s \eta_s) \Xi 
subroutine mass_contribution(mass_cont,Geta_eta,Geta_lambda,Geta_chi,Umat,PhiMat)
implicit none

complex(kind(0d0)), intent(out) :: mass_cont
complex(kind(0d0)), intent(in) ::  Geta_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_sites) 
complex(kind(0d0)) , intent(in) :: Geta_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_links) 
complex(kind(0d0)) , intent(in) :: Geta_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_faces) 
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)


complex(kind(0d0)) Xi_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) Xi_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) Xi_chi(1:NMAT,1:NMAT,1:num_necessary_faces)

complex(kind(0d0)) tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) DinvXi(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) trace
integer gs,ll,lf,ls,ls2,rank
integer k,l


mass_cont=(0d0,0d0)
call make_XiVec(Xi_eta,Xi_lambda,Xi_chi,Umat,PhiMat)
DinvXi=(0d0,0d0)
do gs=1,global_num_sites
  ls2=local_site_of_global(gs)%label_
  rank=local_site_of_global(gs)%rank_
  tmpmat=(0d0,0d0)
  do ls=1,num_sites
    do k=1,NMAT
      do l=1,NMAT
        tmpmat=tmpmat+Geta_eta(:,:,l,k,gs,ls)*Xi_eta(k,l,ls)!& *site_U1Rfactor(ls)
      enddo
    enddo
  enddo
  !!!
  do ll=1,num_links
    do k=1,NMAT
      do l=1,NMAT
        tmpmat=tmpmat+Geta_lambda(:,:,l,k,gs,ll)*Xi_lambda(k,l,ll)! & *site_U1Rfactor(link_org(ls))
      enddo
    enddo
  enddo
  !!!
  do lf=1,num_faces
    do k=1,NMAT
      do l=1,NMAT
        tmpmat=tmpmat+Geta_chi(:,:,l,k,gs,lf)*Xi_chi(k,l,lf) ! & *site_U1Rfactor(sites_in_f(lf)%label_(1))
      enddo
    enddo
  enddo
  call MPI_REDUCE(tmpmat,DinvXi(:,:,ls2),NMAT*NMAT,MPI_DOUBLE_COMPLEX, &
    MPI_SUM,rank,MPI_COMM_WORLD,IERR)
enddo
trace=(0d0,0d0)
do ls=1,num_sites
  do k=1,NMAT
    do l=1,NMAT
      trace=trace - DinvXi(k,l,ls)*PhiMat(l,k,ls)
    enddo
  enddo
enddo
call MPI_REDUCE(trace,mass_cont,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)
mass_cont = mass_cont * dcmplx( 0.5d0*mass_square_phi*overall_factor )

end subroutine mass_contribution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Q( C_{\bar\phi} ) \Xi term 
!!  r/NS \sum_s ( B_s^{r-1} 2/N Tr(\bar\phi_s \eta_s) \Xi
!! where
!!  B_s=1/N Tr(\bar\phi^2)
!!  r = \chi dimG / 4
subroutine calc_QC_Xi(QC_Xi,Geta_eta,Geta_lambda,Geta_chi,Umat,PhiMat)
implicit none

complex(kind(0d0)), intent(out) :: QC_Xi
complex(kind(0d0)), intent(in) ::  Geta_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_sites) 
complex(kind(0d0)) , intent(in) :: Geta_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_links) 
complex(kind(0d0)) , intent(in) :: Geta_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_faces) 
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)


complex(kind(0d0)) Xi_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) Xi_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) Xi_chi(1:NMAT,1:NMAT,1:num_necessary_faces)

complex(kind(0d0)) tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) DinvXi(1:NMAT,1:NMAT)
complex(kind(0d0)) trace
complex(kind(0d0)) tmp_QC_XI
integer gs,ll,lf,ls,ls2,rank
integer i,j,k,l

complex(kind(0d0)) tr_phibar2
double precision :: ratio,eular
double precision :: radius, phase

eular=global_num_sites-global_num_links+global_num_faces 
ratio=dble((NMAT*NMAT-1)*eular)/4d0 


QC_Xi=(0d0,0d0)
tmp_QC_Xi=(0d0,0d0)
call make_XiVec(Xi_eta,Xi_lambda,Xi_chi,Umat,PhiMat)
do gs=1,global_num_sites
  ls2=local_site_of_global(gs)%label_
  rank=local_site_of_global(gs)%rank_

  !! tr_phibar2 = (1/N Tr( \bar\phi^2 ))^{r-1} at ls2
  if( MYRANK == rank ) then 
    trace=(0d0,0d0)
    do i=1,NMAT
      do j=1,NMAT
        trace = trace+PhiMat(j,i,ls2)*PhiMat(i,j,ls2)
      enddo
    enddo
    trace = dconjg(trace)/dcmplx(dble(NMAT))
  
    radius=cdabs(trace)
    phase=atan2(dble(trace),dble(trace*(0d0,-1d0)))
    tr_phibar2 &
      = dcmplx(radius**(ratio-1d0)) * cdexp( (0d0,1d0)*dcmplx(phase*(ratio-1d0)) )
  endif

  !! 2/N Tr(\bar\phi \eta) \Xi at ls2
  tmpmat=(0d0,0d0)
  do ls=1,num_sites
    do k=1,NMAT
      do l=1,NMAT
        tmpmat=tmpmat+Geta_eta(:,:,l,k,gs,ls)*Xi_eta(k,l,ls)
      enddo
    enddo
  enddo
  !!!
  do ll=1,num_links
    do k=1,NMAT
      do l=1,NMAT
        tmpmat=tmpmat+Geta_lambda(:,:,l,k,gs,ll)*Xi_lambda(k,l,ll)
      enddo
    enddo
  enddo
  !!!
  do lf=1,num_faces
    do k=1,NMAT
      do l=1,NMAT
        tmpmat=tmpmat+Geta_chi(:,:,l,k,gs,lf)*Xi_chi(k,l,lf)
      enddo
    enddo
  enddo
  DinvXi=(0d0,0d0)
  call MPI_REDUCE(tmpmat,DinvXi(:,:),NMAT*NMAT,MPI_DOUBLE_COMPLEX, &
    MPI_SUM,rank,MPI_COMM_WORLD,IERR)
  
  !! trace = 2/N Tr( \bar\phi \eta ) \Xi at ls2
  if( MYRANK == rank ) then
    trace=(0d0,0d0)
    do i=1,NMAT
      do j=1,NMAT
        trace = trace + DinvXi(i,j)*dconjg(PhiMat(i,j,ls2))
      enddo
    enddo
    trace = trace * dcmplx( 2d0/dble(Nmat) )
  endif


  !! tmp_QC_Xi = tr_phibar2 * trace at ls2
  if( MYRANK == rank ) then
    tmp_QC_Xi = tmp_QC_Xi + tr_phibar2 * trace
  endif
enddo

call MPI_REDUCE(tmp_QC_Xi,QC_Xi,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)
QC_Xi = QC_Xi * dcmplx( ratio / dble(global_num_sites) )

end subroutine calc_QC_Xi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to make S_eta, S_lambda, S_chi of
!!  \Zia = Tr(eta S_eta) + Tr(lambda S_lambda) + Tr(chi S_chi)
!subroutine make_XiVec(Xi_eta,Xi_lambda,Xi_chi,Umat,PhiMat)
!use matrix_functions, only : matrix_commutator,matrix_3_product
!use parallel
!implicit none
!
!complex(kind(0d0)), intent(out) :: Xi_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
!complex(kind(0d0)), intent(out) :: Xi_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)), intent(out) :: Xi_chi(1:NMAT,1:NMAT,1:num_necessary_faces)
!complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
!
!complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT)
!integer :: s,l,f
!integer :: i,j
!
!do s=1,num_sites
!  call matrix_commutator(tmpmat,PhiMat(:,:,s),Phimat(:,:,s),'N','C')
!  Xi_eta(:,:,s)=dcmplx(alpha_s(s)*0.25d0*overall_factor)*tmpmat
!enddo
!
!do l=1,num_links
!  do j=1,NMAT
!    do i=1,NMAT
!      tmpmat(i,j)=-dconjg( PhiMat(j,i,link_org(l)) )
!    enddo
!  enddo
!  call matrix_3_product(tmpmat,&
!    Umat(:,:,l),PhiMat(:,:,link_tip(l)),Umat(:,:,l),&
!    'N','C','C',(1d0,0d0),'ADD') 
!  Xi_lambda(:,:,l)=(0d0,-1d0)*dcmplx(alpha_l(l)*overall_factor)*tmpmat
!enddo
!
!do f=1,num_faces
!  call make_face_variable(Uf,f,Umat)
!  if(m_omega==0) call make_moment_map0(tmpmat,Uf)
!  if(m_omega==-1) call make_moment_map_adm(tmpmat,Uf)
!  Xi_chi(:,:,f)=(0d0,-0.5d0)*dcmplx(alpha_f(f)*beta_f(f)*overall_factor)*tmpmat
!enddo
!
!  call syncronize_sites(Xi_eta)
!  call syncronize_links(Xi_lambda)
!  call syncronize_faces(Xi_chi)
!end subroutine make_XiVec



