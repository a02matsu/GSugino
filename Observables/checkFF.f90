
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! check Dinv.PF
subroutine check_DinvPF(&
    Geta_eta, Glambda_eta, Gchi_eta, &
    Geta_lambda, Glambda_lambda, Gchi_lambda, &
    Geta_chi, Glambda_chi, Gchi_chi, &
    Umat,PhiMat)
use parallel
use SUN_generators, only : make_SUN_generators
use Dirac_operator, only : Prod_Dirac
use matrix_functions, only : make_random_matrix, make_matrix_traceless
implicit none

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

complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) :: eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: chi(1:NMAT,1:NMAT,1:num_necessary_faces)

complex(kind(0d0)) :: DinvPF_eta1(1:NMAT,1:NMAT,1:num_necessary_sites) 
complex(kind(0d0)) :: DinvPF_lambda1(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: DinvPF_chi1(1:NMAT,1:NMAT,1:num_necessary_faces)

complex(kind(0d0)) :: DinvPF_eta2(1:NMAT,1:NMAT,1:num_necessary_sites) 
complex(kind(0d0)) :: DinvPF_lambda2(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: DinvPF_chi2(1:NMAT,1:NMAT,1:num_necessary_faces)

complex(kind(0d0)) :: DDeta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: DDlambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: DDchi(1:NMAT,1:NMAT,1:num_faces)

complex(kind(0d0)) :: T(1:NMAT,1:NMAT,1:dimG) ! SU(N) generators

integer :: gs,gl,gf,ls,ll,lf
integer :: i,j
integer :: rank
integer :: order
integer :: switch

! SU(N) generators
call make_SUN_generators(T,NMAT)

eta=(0d0,0d0)
lambda=(0d0,0d0)
chi=(0d0,0d0)

do gs=1,global_num_sites
  rank=local_site_of_global(gs)%rank_
  ls=local_site_of_global(gs)%label_
  if( MYRANK==rank ) then
    call make_random_matrix(eta(:,:,ls))
    call make_matrix_traceless(eta(:,:,ls))
  endif
enddo
do gl=1,global_num_links
  rank=local_link_of_global(gl)%rank_
  ll=local_link_of_global(gl)%label_
  if( MYRANK==rank ) then
    call make_random_matrix(lambda(:,:,ll))
    call make_matrix_traceless(lambda(:,:,ll))
  endif
enddo
do gf=1,global_num_faces
  rank=local_face_of_global(gf)%rank_
  lf=local_face_of_global(gf)%label_
  if( MYRANK==rank ) then
    call make_random_matrix(chi(:,:,lf))
    call make_matrix_traceless(chi(:,:,lf))
  endif
enddo
call syncronize_sites(eta)
call syncronize_links(lambda)
call syncronize_faces(chi)


switch=1
if( switch == 0 ) then
  call DinvPF_CG(&
      DinvPF_eta1(:,:,1:num_sites), &
      DinvPF_lambda1(:,:,1:num_links), &
      DinvPF_chi1(:,:,1:num_faces), &
      eta,lambda,chi,Umat,PhiMat)
  call syncronize_sites(DinvPF_eta1)
  call syncronize_links(DinvPF_lambda1)
  call syncronize_faces(DinvPF_chi1)
  call Prod_Dirac(&
        DDeta, DDlambda, DDchi, &
        DinvPF_eta1,DinvPF_lambda1,DinvPF_chi1, UMAT,PhiMat)
else
  call DinvPF_direct( &
      DinvPF_eta2(:,:,1:num_sites), &
      DinvPF_lambda2(:,:,1:num_links), &
      DinvPF_chi2(:,:,1:num_faces), &
      eta,lambda,chi,&
      Geta_eta, Glambda_eta, Gchi_eta, &
      Geta_lambda, Glambda_lambda, Gchi_lambda, &
      Geta_chi, Glambda_chi, Gchi_chi, &
      Umat,PhiMat)
  call syncronize_sites(DinvPF_eta2)
  call syncronize_links(DinvPF_lambda2)
  call syncronize_faces(DinvPF_chi2)
  call Prod_Dirac(&
        DDeta, DDlambda, DDchi, &
        DinvPF_eta2,DinvPF_lambda2,DinvPF_chi2, UMAT,PhiMat)
endif

do order=0,NPROCS-1
  if( MYRANK==order ) then
    do ls=1,num_sites
      do i=1,NMAT
        do j=1,NMAT
          if( cdabs( DDeta(i,j,ls)-eta(i,j,ls) ) > 1d-3 ) then
            write(*,*) global_site_of_local(ls), i,j
            write(*,*) DDeta(i,j,ls)
            write(*,*) eta(i,j,ls)
          endif
        enddo
      enddo
    enddo
    call MPI_BARRIER(MPI_COMM_WORLD,IERR)
  endif
enddo
do order=0,NPROCS-1
  if( MYRANK==order ) then
    do ll=1,num_links
      do i=1,NMAT
        do j=1,NMAT
          if( cdabs( DDlambda(i,j,ll)-lambda(i,j,ll) ) > 1d-3 ) then
            write(*,*) global_link_of_local(ll), i,j
            write(*,*) DDlambda(i,j,ll)
            write(*,*) lambda(i,j,ll)
          endif
        enddo
      enddo
    enddo
    call MPI_BARRIER(MPI_COMM_WORLD,IERR)
  endif
enddo
do order=0,NPROCS-1
  if( MYRANK==order ) then
    do lf=1,num_faces
      do i=1,NMAT
        do j=1,NMAT
          if( cdabs( DDchi(i,j,lf)-chi(i,j,lf) ) > 1d-3 ) then
            write(*,*) global_link_of_local(lf), i,j
            write(*,*) DDchi(i,j,lf)
            write(*,*) chi(i,j,lf)
          endif
        enddo
      enddo
    enddo
    call MPI_BARRIER(MPI_COMM_WORLD,IERR)
  endif
enddo

stop


call DinvPF_direct( &
    DinvPF_eta2(:,:,1:num_sites), &
    DinvPF_lambda2(:,:,1:num_links), &
    DinvPF_chi2(:,:,1:num_faces), &
    eta,lambda,chi,&
    Geta_eta, Glambda_eta, Gchi_eta, &
    Geta_lambda, Glambda_lambda, Gchi_lambda, &
    Geta_chi, Glambda_chi, Gchi_chi, &
    Umat,PhiMat)

do order=0,NPROCS-1
  if( MYRANK==order ) then
    do ls=1,num_sites
      write(*,*) global_site_of_local(ls), DinvPF_eta1(:,:,ls)-DinvPF_eta2(:,:,ls)
    enddo
    call MPI_BARRIER(MPI_COMM_WORLD,IERR)
  endif
enddo


end subroutine check_DinvPF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Subroutine to compute Dinv.PF by CG
subroutine DinvPF_CG(DinvPF_eta, DinvPF_lambda, DinvPF_chi, &
    eta,lambda,chi,Umat,PhiMat)
use global_parameters
use parallel
implicit none

complex(kind(0d0)), intent(in) :: eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: chi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)), intent(out) :: DinvPF_eta(1:NMAT,1:NMAT,1:num_sites) 
complex(kind(0d0)), intent(out) :: DinvPF_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(out) :: DinvPF_chi(1:NMAT,1:NMAT,1:num_faces)

integer info

call calc_dinvf(DinvPF_eta,DinvPF_lambda,DinvPF_chi,&
 eta,lambda,chi,umat,phimat,info)


end subroutine DinvPF_CG

