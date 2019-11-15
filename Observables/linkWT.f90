!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! trivial WT identity
!!  <Sb^L> + <Sf^L> + \Xi^L QS_b
subroutine calc_linkWT(WT,Glambda_eta,Geta_lambda,Glambda_lambda,PhiMat,Umat)
use matrix_functions, only : trace_MM
use parallel
implicit none

complex(kind(0d0)), intent(out) :: WT
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: UMat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) ::  Glambda_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_sites) 
complex(kind(0d0)), intent(in) ::  Geta_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_links) 
complex(kind(0d0)), intent(in) ::  Glambda_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_links) 

double precision :: Sblink
complex(kind(0d0)) :: mass_cont
complex(kind(0d0)) :: Sflink
complex(kind(0d0)) Xi_lambda(1:NMAT,1:NMAT,1:num_necessary_links)

!! (0) preparation (making \Xi)
call make_XiVec_link(Xi_lambda,Umat,Phimat)

!! (1) bosonic action
call calc_bosonic_action_link(Sblink,Umat,PhiMat)

!! (2) mass contribution
call mass_contribution_link(mass_cont,Glambda_eta,Xi_lambda,Umat,PhiMat)

!! (3) contribution from fermion number 
call fermion_action_link(Sflink,Glambda_eta,Glambda_lambda,Umat,PhiMat)

if( MYRANK==0 ) then
  WT=dcmplx(Sblink)+mass_cont+Sflink  
endif

end subroutine calc_linkWT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! mass contribution to the trivial WT identity 
!!  -1/2g^2 mu^2/2  Tr( \phi_s \eta_s) \Xi_L
subroutine mass_contribution_link(mass_cont,Glambda_eta,Xi_lambda,Umat,PhiMat)
implicit none

complex(kind(0d0)), intent(out) :: mass_cont
complex(kind(0d0)), intent(in) :: Glambda_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_sites) 
complex(kind(0d0)), intent(in) :: Xi_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: UMat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) DinvXi(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) trace
integer gl,ll, rank
integer ls
integer k,l

mass_cont=(0d0,0d0)

DinvXi=(0d0,0d0)
do gl=1,global_num_links
  ll=local_link_of_global(gl)%label_
  rank=local_link_of_global(gl)%rank_
  tmpmat=(0d0,0d0)
  !! D^{-1} \Phi
  do ls=1,num_sites
    do k=1,NMAT
      do l=1,NMAT
        tmpmat=tmpmat+Glambda_eta(:,:,l,k,gl,ls)*Phimat(k,l,ls)
      enddo
    enddo
  enddo
  call MPI_REDUCE(tmpmat,DinvXi(:,:,ll),NMAT*NMAT,MPI_DOUBLE_COMPLEX, &
    MPI_SUM,rank,MPI_COMM_WORLD,IERR)
enddo
trace=(0d0,0d0)
do ll=1,num_links
  do k=1,NMAT
    do l=1,NMAT
      trace=trace + DinvXi(k,l,ll)*Xi_lambda(l,k,ll)
    enddo
  enddo
enddo
call MPI_REDUCE(trace,mass_cont,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)
mass_cont = mass_cont * dcmplx( 0.5d0*mass_square_phi*overall_factor )

end subroutine mass_contribution_link

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! mass contribution to the trivial WT identity 
!!  -1/2g^2 mu^2/2  Tr( \phi_s \eta_s) \Xi 
subroutine fermion_action_link(Sflink,Geta_lambda,Glambda_lambda,Umat,PhiMat)
use matrix_functions, only : matrix_3_product
implicit none

complex(kind(0d0)), intent(out) :: Sflink
complex(kind(0d0)), intent(in) ::  Geta_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_links) 
complex(kind(0d0)), intent(in) ::  Glambda_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_links) 
complex(kind(0d0)), intent(in) :: UMat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) trace,tmp
complex(kind(0d0)) tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) tmpmat2(1:NMAT,1:NMAT)
integer ll,gl,ls,gs
integer i,j,k,l

Sflink=(0d0,0d0)
trace=(0d0,0d0)
do ll=1,num_links
  gl=global_link_of_local(ll)

!!! lambda-eta
  tmp=(0d0,0d0)
  ls=link_tip(ll)
  gs=global_site_of_local(ls)
  do l=1,NMAT
    do k=1,NMAT
      do j=1,NMAT
        do i=1,NMAT
          tmp=tmp&
            - Geta_lambda(k,l,i,j,gs,ll)*Umat(j,k,ll)*dconjg(Umat(i,l,ll))
            !+ Glambda_eta(i,j,k,l,gl,ls)*Umat(j,k,ll)*dconjg(Umat(i,l,ll))
        enddo
      enddo
    enddo
  enddo
  ls=link_org(ll)
  gs=global_site_of_local(ls)
  do j=1,NMAT
    do i=1,NMAT
      tmp=tmp&
        + Geta_lambda(i,j,j,i,gs,ll)
    enddo
  enddo
  trace = trace + tmp * (0d0,1d0) * dcmplx( alpha_l(ll) )
!!! lambda-lambda
  call matrix_3_product(tmpmat,&
    Umat(:,:,ll),Phimat(:,:,link_tip(ll)),Umat(:,:,ll), 'N','C','C')
  tmp=(0d0,0d0)
  do i=1,NMAT
    do j=1,NMAT
      do k=1,NMAT
        tmp = tmp + Glambda_lambda(i,k,k,j,gl,ll)*tmpmat(j,i)
      enddo
    enddo
  enddo
  trace = trace + tmp * (-2d0,0d0) * dcmplx( alpha_l(ll) )
enddo
call MPI_REDUCE(trace,Sflink,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)

Sflink=Sflink * dcmplx( overall_factor )


end subroutine fermion_action_link

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Q( C_{\bar\phi} ) \Xi term 
!!  r/NS \sum_s ( B_s^{r-1} 2/N Tr(\bar\phi_s \eta_s) \Xi
!! where
!!  B_s=1/N Tr(\bar\phi^2)
!!  r = \chi dimG / 4
!subroutine calc_QC_Xi(QC_Xi,Geta_eta,Geta_lambda,Geta_chi,Umat,PhiMat)
!implicit none
!
!complex(kind(0d0)), intent(out) :: QC_Xi
!complex(kind(0d0)), intent(in) ::  Geta_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_sites) 
!complex(kind(0d0)) , intent(in) :: Geta_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_links) 
!complex(kind(0d0)) , intent(in) :: Geta_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_faces) 
!complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
!
!
!complex(kind(0d0)) Xi_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
!complex(kind(0d0)) Xi_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)) Xi_chi(1:NMAT,1:NMAT,1:num_necessary_faces)
!
!complex(kind(0d0)) tmpmat(1:NMAT,1:NMAT)
!complex(kind(0d0)) DinvXi(1:NMAT,1:NMAT)
!complex(kind(0d0)) trace
!complex(kind(0d0)) tmp_QC_XI
!integer gs,ll,lf,ls,ls2,rank
!integer i,j,k,l
!
!complex(kind(0d0)) tr_phibar2
!double precision :: ratio,eular
!double precision :: radius, phase
!
!eular=global_num_sites-global_num_links+global_num_faces 
!ratio=dble((NMAT*NMAT-1)*eular)/4d0 
!
!
!QC_Xi=(0d0,0d0)
!tmp_QC_Xi=(0d0,0d0)
!call make_XiVec(Xi_eta,Xi_lambda,Xi_chi,Umat,PhiMat)
!do gs=1,global_num_sites
!  ls2=local_site_of_global(gs)%label_
!  rank=local_site_of_global(gs)%rank_
!
!  !! tr_phibar2 = (1/N Tr( \bar\phi^2 ))^{r-1} at ls2
!  if( MYRANK == rank ) then 
!    trace=(0d0,0d0)
!    do i=1,NMAT
!      do j=1,NMAT
!        trace = trace+PhiMat(j,i,ls2)*PhiMat(i,j,ls2)
!      enddo
!    enddo
!    trace = dconjg(trace)/dcmplx(dble(NMAT))
!  
!    radius=cdabs(trace)
!    phase=atan2(dble(trace),dble(trace*(0d0,-1d0)))
!    tr_phibar2 &
!      = dcmplx(radius**(ratio-1d0)) * cdexp( (0d0,1d0)*dcmplx(phase*(ratio-1d0)) )
!  endif
!
!  !! 2/N Tr(\bar\phi \eta) \Xi at ls2
!  tmpmat=(0d0,0d0)
!  do ls=1,num_sites
!    do k=1,NMAT
!      do l=1,NMAT
!        tmpmat=tmpmat+Geta_eta(:,:,l,k,gs,ls)*Xi_eta(k,l,ls)
!      enddo
!    enddo
!  enddo
!  !!!
!  do ll=1,num_links
!    do k=1,NMAT
!      do l=1,NMAT
!        tmpmat=tmpmat+Geta_lambda(:,:,l,k,gs,ll)*Xi_lambda(k,l,ll)
!      enddo
!    enddo
!  enddo
!  !!!
!  do lf=1,num_faces
!    do k=1,NMAT
!      do l=1,NMAT
!        tmpmat=tmpmat+Geta_chi(:,:,l,k,gs,lf)*Xi_chi(k,l,lf)
!      enddo
!    enddo
!  enddo
!  DinvXi=(0d0,0d0)
!  call MPI_REDUCE(tmpmat,DinvXi(:,:),NMAT*NMAT,MPI_DOUBLE_COMPLEX, &
!    MPI_SUM,rank,MPI_COMM_WORLD,IERR)
!  
!  !! trace = 2/N Tr( \bar\phi \eta ) \Xi at ls2
!  if( MYRANK == rank ) then
!    trace=(0d0,0d0)
!    do i=1,NMAT
!      do j=1,NMAT
!        trace = trace + DinvXi(i,j)*dconjg(PhiMat(i,j,ls2))
!      enddo
!    enddo
!    trace = trace * dcmplx( 2d0/dble(Nmat) )
!  endif
!
!
!  !! tmp_QC_Xi = tr_phibar2 * trace at ls2
!  if( MYRANK == rank ) then
!    tmp_QC_Xi = tmp_QC_Xi + tr_phibar2 * trace
!  endif
!enddo
!
!call MPI_REDUCE(tmp_QC_Xi,QC_Xi,1,MPI_DOUBLE_COMPLEX, &
!  MPI_SUM,0,MPI_COMM_WORLD,IERR)
!QC_Xi = QC_Xi * dcmplx( ratio / dble(global_num_sites) )
!
!end subroutine calc_QC_Xi

