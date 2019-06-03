!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate D^{-1}_{*;X,i,j}
!!   X: Ftype, global_position
!!   i,j : matrix element
!subroutine calc_DinvVec(Dinv_eta,Dinv_lambda,Dinv_chi,Ftype,global_position,i,j,Phimat,Umat)
!use global_parameters
!use parallel
!implicit none
!
!complex(kind(0d0)), intent(out) :: Dinv_eta(1:NMAT,1:NMAT,1:num_sites)
!complex(kind(0d0)), intent(out) :: Dinv_lambda(1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)), intent(out) :: Dinv_chi(1:NMAT,1:NMAT,1:num_faces)
!character, intent(in) :: Ftype ! S:eta, L:lambda, F:chi
!integer, intent(in) :: global_position, i, j
!complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
!
!integer :: info
!complex(kind(0d0)) :: PFeta(1:NMAT,1:NMAT,1:num_necessary_sites)
!complex(kind(0d0)) :: PFlambda(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)) :: PFchi(1:NMAT,1:NMAT,1:num_necessary_faces)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! FOR TEST
!!complex(kind(0d0)) :: Dinv2_eta(1:NMAT,1:NMAT,1:num_sites)
!!complex(kind(0d0)) :: Dinv2_lambda(1:NMAT,1:NMAT,1:num_links)
!!complex(kind(0d0)) :: Dinv2_chi(1:NMAT,1:NMAT,1:num_faces)
!!integer :: s,l,f,a,b
!!double precision :: tmp
!!! END FOR TEST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!PFeta=(0d0,0d0)
!PFlambda=(0d0,0d0)
!PFchi=(0d0,0d0)
!
!if( Ftype=='S' ) then 
!  if( local_site_of_global(global_position)%rank_ == MYRANK ) then
!    PFeta(j,i,local_site_of_global(global_position)%label_) = (1d0,0d0)
!  endif
!  call syncronize_sites(PFeta)
!elseif( Ftype=='L' ) then
!  if( local_link_of_global(global_position)%rank_ == MYRANK ) then
!    PFlambda(j,i,local_link_of_global(global_position)%label_) = (1d0,0d0)
!  endif
!  call syncronize_links(PFlambda)
!elseif( Ftype=='F' ) then
!  if( local_face_of_global(global_position)%rank_ == MYRANK ) then
!    PFeta(j,i,local_face_of_global(global_position)%label_) = (1d0,0d0)
!  endif
!  call syncronize_faces(PFchi)
!endif
!
!!call calc_DinvF(Dinv_eta, Dinv_lambda, Dinv_chi, PFeta, PFlambda, PFchi, Umat, PhiMat, info)
!call calc_DinvF_direct(Dinv_eta, Dinv_lambda, Dinv_chi, PFeta, PFlambda, PFchi, Umat, PhiMat)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! FOR TEST
!!write(*,*) info
!!call calc_DinvF_direct(Dinv2_eta, Dinv2_lambda, Dinv2_chi, PFeta, PFlambda, PFchi, Umat, PhiMat)
!!
!!if( MYRANK==0 ) write(*,*) "#######  ", Ftype,global_position,i,j, "#######"
!!tmp=0d0
!!do s=1,num_sites
!!  do b=1,NMAT
!!    do a=1,NMAT
!!      tmp = tmp + &
!!       dble( (Dinv_eta(a,b,s)-Dinv2_eta(a,b,s)) &
!!             *dconjg(Dinv_eta(a,b,s)-Dinv2_eta(a,b,s)))
!!    enddo
!!  enddo
!!enddo
!!write(*,*) MYRANK, "site", tmp
!!
!!tmp=0d0
!!do l=1,num_links
!!  do b=1,NMAT
!!    do a=1,NMAT
!!      tmp = tmp + &
!!       dble( (Dinv_lambda(a,b,l)-Dinv2_lambda(a,b,l)) &
!!             *dconjg(Dinv_lambda(a,b,l)-Dinv2_lambda(a,b,l)))
!!    enddo
!!  enddo
!!enddo
!!write(*,*) MYRANK, "link", tmp
!!
!!tmp=0d0
!!do f=1,num_faces
!!  do b=1,NMAT
!!    do a=1,NMAT
!!      tmp = tmp + &
!!       dble( (Dinv_chi(a,b,f)-Dinv2_chi(a,b,f)) &
!!             *dconjg(Dinv_chi(a,b,f)-Dinv2_chi(a,b,f)))
!!    enddo
!!  enddo
!!enddo
!!write(*,*) MYRANK, "face", tmp
!!! END FOR TEST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!end subroutine calc_DinvVec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! calculate  rot[ Tr(\lambda \eta) ](gf)
!subroutine calc_rot_lambda_eta(trace, gf, PhiMat, Umat)
!use global_parameters
!use parallel
!implicit none
!
!complex(kind(0d0)), intent(out) :: trace
!integer, intent(in) :: gf ! global face
!complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
!
!complex(kind(0d0)) :: Dinv_eta(1:NMAT,1:NMAT,1:num_sites)
!complex(kind(0d0)) :: Dinv_lambda(1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)) :: Dinv_chi(1:NMAT,1:NMAT,1:num_faces)
!integer :: i,j
!integer :: kk,dir,gl
!integer :: gs,ll, rank
!complex(kind(0d0)) :: tmp
!
!
!trace=(0d0,0d0)
!tmp=(0d0,0d0)
!do kk=1, global_links_in_f(gf)%num_
!  gl=global_links_in_f(gf)%link_labels_(kk)
!  dir=global_links_in_f(gf)%link_dirs_(kk)
!  gs=global_link_org(gl)
!  rank=local_site_of_global(gs)%rank_
!  ll=local_link_of_global(gl)%label_
!  do i=1,NMAT
!    do j=1,NMAT
!      !! D^{-1}_{*;gs,i,j}
!      call calc_DinvVec(Dinv_eta,Dinv_lambda,Dinv_chi,'S',gs,i,j,Phimat,Umat)
!      if( MYRANK == rank ) then
!        tmp=tmp + dcmplx( dble(dir) * global_beta_f(gf) ) * Dinv_lambda(ll,j,i)
!      endif
!    enddo
!  enddo
!enddo
!
!call MPI_REDUCE(tmp,trace,1,MPI_DOUBLE_COMPLEX, &
!  MPI_SUM,0,MPI_COMM_WORLD,IERR)
!end subroutine calc_rot_lambda_eta
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! calculate div[ Tr(lambda \chi) ](gf)
!subroutine calc_div_lambda_chi(trace, gf, PhiMat, Umat)
!use global_parameters
!use parallel
!use global_subroutines, only : calc_prodUl_from_n1_to_n2_in_Uf
!use matrix_functions, only : matrix_3_product, matrix_product
!implicit none
!
!complex(kind(0d0)), intent(out) :: trace
!integer, intent(in) :: gf ! global link 
!complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
!
!complex(kind(0d0)) :: Dinv_eta(1:NMAT,1:NMAT,1:num_sites)
!complex(kind(0d0)) :: Dinv_lambda(1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)) :: Dinv_chi(1:NMAT,1:NMAT,1:num_faces)
!complex(kind(0d0)) :: Ucarry(1:NMAT,1:NMAT,1:num_sites), tmpmat(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: Utmp(1:NMAT,1:NMAT), Vtmp(1:NMAT,1:NMAT)
!integer :: i,j,k,l,kk
!integer :: local_f, Frank, Srank,Lrank
!integer :: gs, gl1, gl2
!integer :: ls,ll
!integer :: tag
!complex(kind(0d0)) :: tmp
!
!
!tmp=(0d0,0d0)
!trace=(0d0,0d0)
!
!! prepare Ucarry = U_{s(l)...f}
!Ucarry=(0d0,0d0)
!local_f = local_face_of_global(gf)%label_
!Frank = local_face_of_global(gf)%rank_
!do kk=1, global_sites_in_f(gf)%num_
!  gs=global_sites_in_f(gf)%label_(kk)
!  ! local data of the site
!  Srank=local_site_of_global(gs)%rank_
!  ls=local_site_of_global(gs)%label_
!  ! tag for sending/receiving Ucarry
!  tag=gs
!  if( MYRANK == Frank ) then
!    call calc_prodUl_from_n1_to_n2_in_Uf(tmpmat,local_f,1,kk-1,Umat)
!    if( MYRANK == Srank ) then
!      Ucarry(:,:,ls)=tmpmat
!    else
!      call MPI_SEND(tmpmat,NMAT*NMAT,MPI_DOUBLE_COMPLEX,Srank,tag,MPI_COMM_WORLD,IERR)
!    endif
!  else
!    call MPI_RECV(Ucarry(:,:,ls),NMAT*NMAT,MPI_DOUBLE_COMPLEX,Frank,tag,MPI_COMM_WORLD,ISTATUS,IERR)
!  endif
!enddo
!
!do l=1,NMAT
!  do k=1,NMAT
!    call calc_DinvVec(Dinv_eta,Dinv_lambda,Dinv_chi,'F',gf,k,l,Phimat,Umat)
!    !!!! links from s !!!!
!    do gs=1, global_sites_in_f(gf)%num_
!      ! local data of the site
!      Srank=local_site_of_global(gs)%rank_
!      ls=local_site_of_global(gs)%label_
!      if( MYRANK == Srank ) then
!        do kk=1, linktip_from_s(ls)%num_
!          ll=linktip_from_s(ls)%labels_(kk)
!          call matrix_3_product(tmpmat,&
!            Ucarry(:,:,ls),Dinv_lambda(:,:,ll),Ucarry(:,:,ls)&
!            ,'C','N','N')
!          tmp=tmp + tmpmat(l,k) / dcmplx(dble( num_faces_in_s(gs) ))
!        enddo
!      endif
!    enddo
!    !!!! links to s !!!!
!    ! prepare DinvVec( lambda_{l2}, ; chi_f,k,l )
!    do gs=1, global_sites_in_f(gf)%num_
!      ! local data of the site
!      Srank=local_site_of_global(gs)%rank_
!      ls=local_site_of_global(gs)%label_
!      do kk=1, global_linkorg_to_s(gs)%num_
!        gl2=global_linkorg_to_s(gs)%labels_(kk)
!        !!
!        Lrank=local_link_of_global(gl2)%rank_
!        ll=local_link_of_global(gl2)%label_
!        !!
!        tag=global_linkorg_to_s(gs)%labels_(kk)
!        if( MYRANK == Srank ) then
!          if( MYRANK == Lrank ) then
!            Vtmp = Dinv_lambda(:,:,ll)
!          else
!            call MPI_RECV(Vtmp(:,:),NMAT*NMAT,MPI_DOUBLE_COMPLEX,Lrank,tag,MPI_COMM_WORLD,ISTATUS,IERR)
!          endif
!        elseif( MYRANK == Lrank ) then
!          call MPI_SEND(Dinv_lambda(:,:,ll),NMAT*NMAT,MPI_DOUBLE_COMPLEX,Srank,tag,MPI_COMM_WORLD,IERR)
!        endif
!
!        if( MYRANK == Srank ) then
!          do ll=1, num_necessary_links
!            if( global_link_of_local(ll) == gl2 ) exit
!          enddo
!          call matrix_product(Utmp, Umat(:,:,ll),Ucarry(:,:,ls))
!          call matrix_3_product(tmpmat, Utmp,Vtmp,Utmp,'C','N','N')
!
!          tmp = tmp -  tmpmat(l,k) / dcmplx(dble( num_faces_in_s(gs) ))
!        endif
!      enddo
!    enddo
!  enddo
!enddo
!
!
!call MPI_REDUCE(tmp,trace,1,MPI_DOUBLE_COMPLEX, &
!  MPI_SUM,0,MPI_COMM_WORLD,IERR)
!end subroutine calc_div_lambda_chi
!    
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! calculate D^\mu J_\mu^{U(1)V}
!subroutine calc_divJ_U1V(divJ, gf, PhiMat, Umat)
!use global_parameters
!use parallel
!implicit none
!
!complex(kind(0d0)), intent(out) :: divJ
!integer, intent(in) :: gf ! global face 
!complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
!
!complex(kind(0d0)) :: trace1, trace2
!
!call calc_rot_lambda_eta(trace1, gf, PhiMat, Umat)
!call calc_div_lambda_chi(trace2, gf, PhiMat, Umat)
!
!divJ=(0d0,0d0)
!if( MYRANK==0 ) then
!  divJ = (0.5d0,0d0)*trace1 - trace2
!endif
!
!
!end subroutine calc_divJ_U1V
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate lambda.lambda(U.\bar{phi}.U^-1 + \bar{phi})
subroutine calc_SfL2(SfL2, PhiMat, Umat, Glambda_lambda)
use global_parameters
use parallel
use matrix_functions, only : matrix_3_product, trace_MM
use SUN_generators, only : trace_MTa, make_Mijkl_from_modes
implicit none

complex(kind(0d0)), intent(out) :: SfL2
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: Glambda_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_links)


complex(kind(0d0)) :: MMAT(1:NMAT,1:NMAT),tmpmat(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: DP(1:NMAT,1:NMAT,1:NMAT,1:NMAT)
!complex(kind(0d0)) :: modes(1:dimG,1:dimG)
complex(kind(0d0)) :: trace, tmpSfL2, tmp
integer :: gl,ll
integer :: a,b,i,j,k

SfL2=(0d0,0d0)
tmpSfL2=(0d0,0d0)
do ll=1,num_links
  gl=global_link_of_local(ll)
  !!!
  do j=1, NMAT
    do i=1, NMAT
      MMat(i,j)=dconjg( PhiMat(j,i,link_org(ll)) )
    enddo
  enddo
  call matrix_3_product(MMat,Umat(:,:,ll),PhiMat(:,:,link_tip(ll)),Umat(:,:,ll),&
    'N','C','C',(1d0,0d0),'ADD')
  MMat= MMat * alpha_l(ll)
  !!!!
  tmpmat=(0d0,0d0)
  do j=1,NMAT
    do i=1,NMAT
      tmp=(0d0,0d0)
      do k=1,NMAT
        tmp = tmp + Glambda_lambda(i,k,k,j,gl,ll) !DP(i,k,k,j)
      enddo
      trace = trace + tmp*MMAT(j,i)
    enddo
  enddo
  tmpSfL2=tmpSfL2+trace
enddo

call MPI_REDUCE(tmpSfL2,SfL2,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)

if( MYRANK == 0 ) then
  SfL2 = SfL2 * overall_factor
endif

end subroutine calc_SfL2

