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
call calc_Sf_site(Sfsite,Geta_eta,PhiMat)

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
complex(kind(0d0)), intent(in) :: Xi_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
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
!subroutine fermion_action_site(Sfsite,Geta_eta,PhiMat)
!implicit none
!
!complex(kind(0d0)), intent(out) :: Sfsite
!complex(kind(0d0)), intent(in) ::  Geta_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_sites) 
!complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
!
!complex(kind(0d0)) trace,tmp
!integer ls,gs
!integer i,j,k
!
!Sfsite=(0d0,0d0)
!trace=(0d0,0d0)
!do ls=1,num_sites
!  gs=global_site_of_local(ls)
!  tmp=(0d0,0d0)
!  do k=1,NMAT
!    do j=1,NMAT
!      do i=1,NMAT
!        tmp=tmp&
!          +Geta_eta(k,i,j,k,gs,ls)*Phimat(i,j,ls) &
!          -Geta_eta(i,k,k,j,gs,ls)*PhiMat(j,i,ls)
!      enddo
!    enddo
!  enddo
!  trace = trace + tmp * dcmplx( -0.25d0 * alpha_s(ls) )
!enddo
!call MPI_REDUCE(trace,Sfsite,1,MPI_DOUBLE_COMPLEX, &
!  MPI_SUM,0,MPI_COMM_WORLD,IERR)
!
!Sfsite=Sfsite * dcmplx( overall_factor )
!
!
!end subroutine fermion_action_site

