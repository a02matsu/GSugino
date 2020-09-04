!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! A = 1/NS sum_s ( 1/N Tr(\bar(\phi)^2) )^{dimG*\chi/4}
subroutine calc_phibar_compensator(Acomp,PhiMat)
use parallel
use global_parameters
use matrix_functions, only : matrix_product
implicit none

complex(kind(0d0)), intent(out) :: Acomp
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
double precision :: ratio,eular
double precision :: radius, phase

complex(kind(0d0)) :: tmp_Acomp, tmp
integer :: ls
integer :: i,j

eular=global_num_sites-global_num_links+global_num_faces 
ratio=dble((NMAT*NMAT-1)*eular)/4d0 
Acomp=(0d0,0d0)
tmp_Acomp=(0d0,0d0)
do ls=1,num_sites
  tmp=(0d0,0d0)
  do i=1,NMAT
    do j=1,NMAT
      tmp = tmp + PhiMat(i,j,ls)*PHiMat(j,i,ls)
    enddo
  enddo
  tmp=dconjg(tmp)/dcmplx(dble(NMAT))
  radius=cdabs(tmp)
  !phase=atan2(dble(tmp),dble(tmp*(0d0,-1d0)))
  phase=atan2(dble(tmp*(0d0,-1d0)),dble(tmp))

  tmp_Acomp=tmp_Acomp + dcmplx(radius**ratio) * cdexp( (0d0,1d0)*dcmplx(phase*ratio) )
enddo

call MPI_REDUCE(tmp_Acomp,Acomp,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)
  
Acomp=Acomp/dcmplx(dble(global_num_sites))

end subroutine calc_phibar_compensator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! local phibar
subroutine calc_phibar_compensator_site(phibar,PhiMat)
use global_parameters
use parallel
implicit none

complex(kind(0d0)), intent(out) :: phibar(1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
double precision :: ratio,eular
double precision :: radius, phase
complex(kind(0d0)) :: tmp
integer :: ls,i,j

eular=global_num_sites-global_num_links+global_num_faces 
ratio=dble((NMAT*NMAT-1)*eular)/4d0 

do ls=1,num_necessary_sites
  tmp=(0d0,0d0)
  do i=1,NMAT
    do j=1,NMAT
      tmp = tmp + PhiMat(i,j,ls)*PHiMat(j,i,ls)
    enddo
  enddo
  tmp=dconjg(tmp)/dcmplx(dble(NMAT))
  radius=cdabs(tmp)
  !phase=atan2(dble(tmp),dble(tmp*(0d0,-1d0)))
  phase=atan2(dble(tmp*(0d0,-1d0)),dble(tmp))

  phibar(ls) = dcmplx(radius**ratio) * cdexp( (0d0,1d0)*dcmplx(phase*ratio) )
enddo

end subroutine calc_phibar_compensator_site

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Q( C_{\bar\phi} ) \Xi term 
!!  r/NS \sum_s ( B_s^{r-1} 2/N Tr(\bar\phi_s \eta_s) \Xi
!! where
!!  B_s=1/N Tr(\bar\phi^2)
!!  r = \chi dimG / 4
subroutine calc_QCphibar_Xi(&
    QC_Xi_site,QC_Xi_link,QC_Xi_face, &
    Geta_eta,Geta_lambda,Geta_chi,Umat,PhiMat)
implicit none

complex(kind(0d0)), intent(out) :: QC_Xi_site
complex(kind(0d0)), intent(out) :: QC_Xi_link
complex(kind(0d0)), intent(out) :: QC_Xi_face
complex(kind(0d0)), intent(in) ::  Geta_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_sites) 
complex(kind(0d0)), intent(in) :: Geta_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_links) 
complex(kind(0d0)), intent(in) :: Geta_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_faces) 
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)


complex(kind(0d0)) Xi_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) Xi_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) Xi_chi(1:NMAT,1:NMAT,1:num_necessary_faces)

complex(kind(0d0)) tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) DinvXi_site(1:NMAT,1:NMAT)
complex(kind(0d0)) DinvXi_link(1:NMAT,1:NMAT)
complex(kind(0d0)) DinvXi_face(1:NMAT,1:NMAT)
complex(kind(0d0)) trace
complex(kind(0d0)) tmp_QC_Xi_site
complex(kind(0d0)) tmp_QC_Xi_link
complex(kind(0d0)) tmp_QC_Xi_face
integer gs,ll,lf,ls,ls2,rank
integer i,j,k,l

complex(kind(0d0)) tr_phibar2
double precision :: ratio,eular
double precision :: radius, phase

eular=global_num_sites-global_num_links+global_num_faces 
ratio=dble((NMAT*NMAT-1)*eular)/4d0 


QC_Xi_site=(0d0,0d0)
QC_Xi_link=(0d0,0d0)
QC_Xi_face=(0d0,0d0)
tmp_QC_Xi_site=(0d0,0d0)
tmp_QC_Xi_link=(0d0,0d0)
tmp_QC_Xi_face=(0d0,0d0)
call make_XiVec(Xi_eta,Xi_lambda,Xi_chi,Umat,PhiMat)
!!!
do gs=1,global_num_sites
  ls2=local_site_of_global(gs)%label_
  rank=local_site_of_global(gs)%rank_

  !!!!!!!!!
  !! tr_phibar2 = (1/N Tr( \bar\phi^2 ))^{r-1} at ls2
  tr_phibar2=(0d0,0d0)
  if( MYRANK == rank ) then 
    trace=(0d0,0d0)
    do i=1,NMAT
      do j=1,NMAT
        trace = trace+PhiMat(j,i,ls2)*PhiMat(i,j,ls2)
      enddo
    enddo
    trace = dconjg(trace)/dcmplx(dble(NMAT))
  
    radius=cdabs(trace)
    !phase=atan2(dble(trace),dble(trace*(0d0,-1d0)))
    phase=atan2(dble(trace*(0d0,-1d0)),dble(trace))
    tr_phibar2 &
      = dcmplx(radius**(ratio-1d0)) * cdexp( (0d0,1d0)*dcmplx(phase*(ratio-1d0)) )
  endif

  !!!!!!!!!
  !! 2/N Tr(\bar\phi \eta) \Xi at ls2
  !! site
  tmpmat=(0d0,0d0)
  do ls=1,num_sites
    do k=1,NMAT
      do l=1,NMAT
        tmpmat=tmpmat+Geta_eta(:,:,l,k,gs,ls)*Xi_eta(k,l,ls)
      enddo
    enddo
  enddo
  DinvXi_site=(0d0,0d0)
  call MPI_REDUCE(tmpmat,DinvXi_site(:,:),NMAT*NMAT,MPI_DOUBLE_COMPLEX, &
    MPI_SUM,rank,MPI_COMM_WORLD,IERR)
  !! link
  tmpmat=(0d0,0d0)
  do ll=1,num_links
    do k=1,NMAT
      do l=1,NMAT
        tmpmat=tmpmat+Geta_lambda(:,:,l,k,gs,ll)*Xi_lambda(k,l,ll)
      enddo
    enddo
  enddo
  DinvXi_link=(0d0,0d0)
  call MPI_REDUCE(tmpmat,DinvXi_link(:,:),NMAT*NMAT,MPI_DOUBLE_COMPLEX, &
    MPI_SUM,rank,MPI_COMM_WORLD,IERR)
  !! face
  tmpmat=(0d0,0d0)
  do lf=1,num_faces
    do k=1,NMAT
      do l=1,NMAT
        tmpmat=tmpmat+Geta_chi(:,:,l,k,gs,lf)*Xi_chi(k,l,lf)
      enddo
    enddo
  enddo
  DinvXi_face=(0d0,0d0)
  call MPI_REDUCE(tmpmat,DinvXi_face(:,:),NMAT*NMAT,MPI_DOUBLE_COMPLEX, &
    MPI_SUM,rank,MPI_COMM_WORLD,IERR)
  
  !!!!!!!!!!!!!!!!!!
  !! trace = 2/N Tr( \bar\phi \eta ) \Xi at ls2
  !! site
  if( MYRANK == rank ) then
    trace=(0d0,0d0)
    do i=1,NMAT
      do j=1,NMAT
        trace = trace + DinvXi_site(i,j)*dconjg(PhiMat(i,j,ls2))
      enddo
    enddo
    trace = trace * dcmplx( 2d0/dble(Nmat) )
    tmp_QC_Xi_site = tmp_QC_Xi_site + tr_phibar2 * trace
  endif
  !! link
  if( MYRANK == rank ) then
    trace=(0d0,0d0)
    do i=1,NMAT
      do j=1,NMAT
        trace = trace + DinvXi_link(i,j)*dconjg(PhiMat(i,j,ls2))
      enddo
    enddo
    trace = trace * dcmplx( 2d0/dble(Nmat) )
    tmp_QC_Xi_link = tmp_QC_Xi_link + tr_phibar2 * trace
  endif
  !! face
  if( MYRANK == rank ) then
    trace=(0d0,0d0)
    do i=1,NMAT
      do j=1,NMAT
        trace = trace + DinvXi_face(i,j)*dconjg(PhiMat(i,j,ls2))
      enddo
    enddo
    trace = trace * dcmplx( 2d0/dble(Nmat) )
    tmp_QC_Xi_face = tmp_QC_Xi_face + tr_phibar2 * trace
  endif
enddo

!! site
call MPI_REDUCE(tmp_QC_Xi_site,QC_Xi_site,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)
QC_Xi_site = QC_Xi_site * dcmplx( ratio / dble(global_num_sites) )
!! link
call MPI_REDUCE(tmp_QC_Xi_link,QC_Xi_link,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)
QC_Xi_link = QC_Xi_link * dcmplx( ratio / dble(global_num_sites) )
!! face
call MPI_REDUCE(tmp_QC_Xi_face,QC_Xi_face,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)
QC_Xi_face = QC_Xi_face * dcmplx( ratio / dble(global_num_sites) )

end subroutine calc_QCphibar_Xi


