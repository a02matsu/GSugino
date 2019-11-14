!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! trivial WT identity
!!  <Sb> + \mu^2/2 \Xi \sum_s( \Phi_s \eta_s ) - (N^2-1)/2 (Ns+Nl)
subroutine calc_siteWT(WT,Geta_eta,PhiMat)
use matrix_functions, only : trace_MM
use parallel
implicit none

complex(kind(0d0)), intent(out) :: WT
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) ::  Geta_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_sites) 

double precision :: Sbsite
complex(kind(0d0)) :: mass_cont
complex(kind(0d0)) :: Sfsite
complex(kind(0d0)) Xi_eta(1:NMAT,1:NMAT,1:num_necessary_sites)

!! (0) preparation (making \Xi)
call make_XiVec_site(Xi_eta,Phimat)

!! (1) bosonic action
call calc_bosonic_action_site(Sbsite,PhiMat)

!! (2) mass contribution
call mass_contribution_site(mass_cont,Geta_eta,Xi_eta,PhiMat)

!! (3) contribution from fermion number 
call fermion_action_site(Sfsite,Geta_eta,PhiMat)

if( MYRANK==0 ) then
  WT=dcmplx(Sbsite)+mass_cont+Sfsite
endif

end subroutine calc_siteWT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! mass contribution to the trivial WT identity 
!!  -1/2g^2 mu^2/2  Tr( \phi_s \eta_s) \Xi 
subroutine mass_contribution_site(mass_cont,Geta_eta,Xi_eta,PhiMat)
implicit none

complex(kind(0d0)), intent(out) :: mass_cont
complex(kind(0d0)), intent(in) ::  Geta_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_sites) 
complex(kind(0d0)), intent(in) :: Xi_eta(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) DinvXi(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) trace
integer gt,lt, rank
integer ls
integer k,l

mass_cont=(0d0,0d0)

DinvXi=(0d0,0d0)
do gt=1,global_num_sites
  lt=local_site_of_global(gt)%label_
  rank=local_site_of_global(gt)%rank_
  tmpmat=(0d0,0d0)
  do ls=1,num_sites
    do k=1,NMAT
      do l=1,NMAT
        tmpmat=tmpmat+Geta_eta(:,:,l,k,gt,ls)*Xi_eta(k,l,ls)
      enddo
    enddo
  enddo
  call MPI_REDUCE(tmpmat,DinvXi(:,:,lt),NMAT*NMAT,MPI_DOUBLE_COMPLEX, &
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

end subroutine mass_contribution_site

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! mass contribution to the trivial WT identity 
!!  -1/2g^2 mu^2/2  Tr( \phi_s \eta_s) \Xi 
subroutine fermion_action_site(Sfsite,Geta_eta,PhiMat)
implicit none

complex(kind(0d0)), intent(out) :: Sfsite
complex(kind(0d0)), intent(in) ::  Geta_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_sites) 
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) trace,tmp
integer ls,gs
integer i,j,k

Sfsite=(0d0,0d0)
trace=(0d0,0d0)
do ls=1,num_sites
  gs=global_site_of_local(ls)
  tmp=(0d0,0d0)
  do k=1,NMAT
    do j=1,NMAT
      do i=1,NMAT
        tmp=tmp&
          +Geta_eta(k,i,j,k,gs,ls)*Phimat(i,j,ls) &
          -Geta_eta(i,k,k,j,gs,ls)*PhiMat(j,i,ls)
      enddo
    enddo
  enddo
  trace = trace + tmp * dcmplx( -0.25d0 * alpha_s(ls) )
enddo
call MPI_REDUCE(trace,Sfsite,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)

Sfsite=Sfsite * dcmplx( -overall_factor )


end subroutine fermion_action_site

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

