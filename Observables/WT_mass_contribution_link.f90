!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! mass contribution to the trivial WT identity 
!!  -1/2g^2 mu^2/2  Tr( \phi_s \eta_s) \Xi_L
subroutine mass_contribution_link(mass_cont,Glambda_eta,Umat,PhiMat)
implicit none

complex(kind(0d0)), intent(out) :: mass_cont
complex(kind(0d0)), intent(in) :: Glambda_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_sites) 
complex(kind(0d0)), intent(in) :: UMat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)):: Xi_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) DinvPhi(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) trace
integer gl,ll, rank
integer ls
integer k,l

!! (0) preparation (making \Xi)
call make_XiVec_link(Xi_lambda,Umat,Phimat)

mass_cont=(0d0,0d0)

DinvPhi=(0d0,0d0)
do gl=1,global_num_links
  ll=local_link_of_global(gl)%label_
  rank=local_link_of_global(gl)%rank_
  tmpmat=(0d0,0d0)
  !! D^{-1} \Phi
  do ls=1,num_sites
    do k=1,NMAT
      do l=1,NMAT
        tmpmat=tmpmat+Glambda_eta(:,:,l,k,gl,ls)*Phimat(k,l,ls)*U1Rfactor_site(ls)
      enddo
    enddo
  enddo
  call MPI_REDUCE(tmpmat,DinvPhi(:,:,ll),NMAT*NMAT,MPI_DOUBLE_COMPLEX, &
    MPI_SUM,rank,MPI_COMM_WORLD,IERR)
enddo
trace=(0d0,0d0)
do ll=1,num_links
  do k=1,NMAT
    do l=1,NMAT
      trace=trace + DinvPhi(k,l,ll)*Xi_lambda(l,k,ll)
    enddo
  enddo
enddo
call MPI_REDUCE(trace,mass_cont,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)
mass_cont = mass_cont * dcmplx( 0.5d0*mass_square_phi*overall_factor )

end subroutine mass_contribution_link

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! trivial WT identity
!!!  <Sb^L> + <Sf^L> + \Xi^L QS_b
!subroutine calc_linkWT(WT,Glambda_eta,Geta_lambda,Glambda_lambda,PhiMat,Umat)
!use matrix_functions, only : trace_MM
!use parallel
!implicit none
!
!complex(kind(0d0)), intent(out) :: WT
!complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
!complex(kind(0d0)), intent(in) :: UMat(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)), intent(in) ::  Glambda_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_sites) 
!complex(kind(0d0)), intent(in) ::  Geta_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_links) 
!complex(kind(0d0)), intent(in) ::  Glambda_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_links) 
!
!double precision :: Sblink
!complex(kind(0d0)) :: mass_cont
!complex(kind(0d0)) :: Sflink1, Sflink2
!!complex(kind(0d0)) Xi_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
!
!
!!! (1) bosonic action
!call calc_bosonic_action_link(Sblink,Umat,PhiMat)
!
!!! (2) contribution from fermion number 
!call calc_Sf_link1(Sflink1,Glambda_eta,Umat,PhiMat)
!call calc_Sf_link2(Sflink2,PhiMat,Umat,Glambda_eta)
!
!!! (3) mass contribution
!call mass_contribution_link(mass_cont,Glambda_eta,Umat,PhiMat)
!
!
!if( MYRANK==0 ) then
!  WT=dcmplx(Sblink)+Sflink1+Sflink2+mass_cont
!endif
!
!end subroutine calc_linkWT



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! fermion action link 1
!subroutine calc_Sf_link1(Sflink,Geta_lambda,Umat,PhiMat)
!use matrix_functions, only : matrix_3_product
!implicit none
!
!complex(kind(0d0)), intent(out) :: Sflink
!complex(kind(0d0)), intent(in) ::  Geta_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_links) 
!complex(kind(0d0)), intent(in) :: UMat(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
!
!complex(kind(0d0)) trace,tmp
!complex(kind(0d0)) tmpmat(1:NMAT,1:NMAT)
!complex(kind(0d0)) tmpmat2(1:NMAT,1:NMAT)
!integer ll,gl,ls,gs
!integer i,j,k,l
!
!Sflink=(0d0,0d0)
!trace=(0d0,0d0)
!do ll=1,num_links
!  gl=global_link_of_local(ll)
!
!!!! lambda-eta
!  ls=link_tip(ll)
!  gs=global_site_of_local(ls)
!  tmp=(0d0,0d0)
!  do l=1,NMAT
!    do k=1,NMAT
!      do j=1,NMAT
!        do i=1,NMAT
!          tmp=tmp&
!            - Geta_lambda(i,j,k,l,gs,ll)*dconjg(Umat(k,j,ll))*Umat(l,i,ll)
!        enddo
!      enddo
!    enddo
!  enddo
!  ls=link_org(ll)
!  gs=global_site_of_local(ls)
!  do j=1,NMAT
!    do i=1,NMAT
!      tmp = tmp + Geta_lambda(i,j,j,i,gs,ll)
!    enddo
!  enddo
!  trace = trace + tmp * (0d0,1d0) * dcmplx( alpha_l(ll) )
!enddo
!
!call MPI_REDUCE(trace,Sflink,1,MPI_DOUBLE_COMPLEX, &
!  MPI_SUM,0,MPI_COMM_WORLD,IERR)
!
!Sflink=Sflink * dcmplx( overall_factor ) 
!
!
!end subroutine calc_Sf_link1

