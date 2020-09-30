!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to test Q-transformation of S
subroutine check_QS(Umat,PhiMat)
use matrix_functions, only : matrix_commutator, matrix_3_product
use Dirac_operator, only : Prod_Dirac
!use simulation, only : make_bosonic_force_nomass
#ifdef PARALLEL
use parallel
#endif
implicit none

complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) :: Qeta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: Qlambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: Qchi(1:NMAT,1:NMAT,1:num_necessary_faces)

complex(kind(0d0)) :: DQeta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: DQlambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: DQchi(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: Omega(1:NMAT,1:NMAT,1:num_faces)

complex(kind(0d0)) :: Bforce_s(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: Bforce_l(1:NMAT,1:NMAT,1:num_links)

complex(kind(0d0)) :: QS_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: QS_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: QS_chi(1:NMAT,1:NMAT,1:num_faces)

complex(kind(0d0)) :: ctmp
double precision :: tmp,QS,tmp2
integer :: info,s,l,f,i,j,triger

!call Prod_Dirac(DQeta,DQlambda,DQchi,Qeta,Qlambda,Qchi,UMAT,Phimat)

!! boson part
!call make_bosonic_force_nomass(Bforce_s,Bforce_l,Umat,PhiMat)
!do s=1,num_sites
  !Bforce_s(:,:,s)=Bforce_s(:,:,s)*dconjg(U1Rfactor_site(s))
!enddo
!do l=1,num_links
  !Bforce_l(:,:,l)=Bforce_l(:,:,l)*U1Rfactor_site(link_org(l))
!enddo

!! Qeta, Qlambda, Qchi
call make_Qfermion(Qeta,Qlambda,Qchi,Omega,Umat,PhiMat)

if(MYRANK==0) write(*,*) "# QS = 0 ?"
!!!!!!!!!!!!!!!!!!!!!
!! QS-site
call calc_QS_site(QS_eta,PhiMat,Qeta)
tmp=0d0
tmp2=0d0
do s=1,num_sites
  do i=1,NMAT
    do j=1,NMAT
      tmp=tmp+dble(QS_eta(i,j,s)*dconjg(QS_eta(i,j,s)))
    enddo
  enddo
enddo
call MPI_REDUCE(tmp,tmp2,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
if(MYRANK==0) write(*,*) "#   S-site:",tmp2

!!!!!!!!!!!!!!!!!!!!!
!! QS-link
call calc_QS_link(QS_eta,QS_lambda,PhiMat,Umat,Qeta,Qlambda)
tmp=0d0
tmp2=0d0
do s=1,num_sites
  do i=1,NMAT
    do j=1,NMAT
      tmp=tmp+dble(QS_eta(i,j,s)*dconjg(QS_eta(i,j,s)))
    enddo
  enddo
enddo
call MPI_REDUCE(tmp,tmp2,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
if(MYRANK==0) write(*,*) "#   L-site:",tmp2
!!
tmp=0d0
tmp2=0d0
do l=1,num_links
  do i=1,NMAT
    do j=1,NMAT
      tmp=tmp+dble(QS_lambda(i,j,l)*dconjg(QS_lambda(i,j,l)))
    enddo
  enddo
enddo
call MPI_REDUCE(tmp,tmp2,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
if(MYRANK==0) write(*,*) "#   L-link:",tmp2

!!!!!!!!!!!!!!!!!!!!!
!! QS-face
call calc_QS_face(QS_lambda,QS_chi,PhiMat,Umat,Qlambda,Qchi)
tmp=0d0
tmp2=0d0
do l=1,num_links
  do i=1,NMAT
    do j=1,NMAT
      tmp=tmp+dble(QS_lambda(i,j,l)*dconjg(QS_lambda(i,j,l)))
    enddo
  enddo
enddo
call MPI_REDUCE(tmp,tmp2,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
if(MYRANK==0) write(*,*) "#   F-link:",tmp2
!!
tmp=0d0
tmp2=0d0
do f=1,num_faces
  do i=1,NMAT
    do j=1,NMAT
      tmp=tmp+dble(QS_chi(i,j,f)*dconjg(QS_chi(i,j,f)))
    enddo
  enddo
enddo
call MPI_REDUCE(tmp,tmp2,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
if(MYRANK==0) write(*,*) "#   F-face:",tmp2


! Q^2 \Omega を care する
!do f=1,num_faces
!  call matrix_commutator(tmpmat,PhiMat(:,:,sites_in_f(f)%label_(1)),Omega(:,:,f))
!  DQchi(:,:,f)=DQchi(:,:,f)+(0d0,1d0)*dcmplx(beta_f(f))*tmpmat
!enddo

!if(MYRANK==0) write(*,*) "# QS = 0 ?"
!write(*,*) DQchi
!QS=0d0
!tmp=0d0
!!call make_bosonic_force_phi_site(Bforce_s,PhiMat)
!do s=1,num_sites
!  tmp=0d0
!  do i=1,NMAT
!    do j=1,NMAT
!      ctmp=&
!        - DQeta(i,j,s) &
!        + dconjg(Bforce_s(j,i,s))!*U1Rfactor_site(s)
!      tmp=tmp+dble( ctmp*dconjg(ctmp) )
!    enddo
!    !write(*,*) global_site_of_local(s), cdabs(ctmp), &
!    !  dble(cdlog(site_U1Rfactor(s))/LatticeSpacing*(0d0,-1d0))
!  enddo
!enddo
!call MPI_REDUCE(tmp,QS,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
!if(MYRANK==0) write(*,*) "#   site:",QS
!
!QS=0d0
!tmp=0d0
!do l=1,num_links
!  do i=1,NMAT
!    do j=1,NMAT
!      !tmpmat(i,j) = &
!      ctmp = &
!        - DQlambda(i,j,l) & !*site_U1Rfactor(link_org(l))&
!        + Bforce_l(i,j,l) !*site_U1Rfactor(link_org(l))
!      tmp=tmp+dble( ctmp*dconjg(ctmp) )
!    enddo
!  enddo
!  !write(*,*) "(L)",global_link_of_local(l), cdabs(ctmp)
!enddo
!call MPI_REDUCE(tmp,QS,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
!if(MYRANK==0) write(*,*) "#   link:",QS
!
!QS=0d0
!tmp=0d0
!do f=1,num_faces
!  !tmp=0d0
!  do i=1,NMAT
!    do j=1,NMAT
!      ctmp=-DQchi(j,i,f) !*site_U1Rfactor(sites_in_f(f)%label_(1))
!      tmp=tmp+dble( ctmp*dconjg(ctmp) )
!    enddo
!  enddo
!  !write(*,*) "(F)",global_face_of_local(f), cdabs(ctmp)
!enddo
!call MPI_REDUCE(tmp,QS,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
!if(MYRANK==0) write(*,*) "#   face:",QS

end subroutine check_QS

 

