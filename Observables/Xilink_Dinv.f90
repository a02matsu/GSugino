!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compute \Xi_l \lambda_l in the form that
!!   Xi_l Dinv_{l,I}
subroutine calc_Xilink_Dinv(XiDinv_eta, XiDinv_lambda, XiDinv_chi,&
    Geta_lambda,Glambda_lambda,Gchi_lambda,&
    PhiMat,Umat)
implicit none

!complex(kind(0d0)), intent(out) :: Xisite_QS
complex(kind(0d0)), intent(out) :: XiDinv_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(out) :: XiDinv_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(out) :: XiDinv_chi(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)), intent(in) :: Geta_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_necessary_links) 
complex(kind(0d0)), intent(in) :: Glambda_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_necessary_links) 
complex(kind(0d0)), intent(in) :: Gchi_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_necessary_links) 
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)


complex(kind(0d0)) :: Xi_lambda(1:NMAT,1:NMAT,1:num_necessary_links)

complex(kind(0d0)) :: tmp
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
integer :: ls,gs
integer :: ll,gl,l
integer :: lf,gf
integer :: i,j
integer :: rank


!! Xi_site = Xe_eta.eta
call make_XiVec_link(Xi_lambda,Umat,PhiMat)

!! 
tmpmat2=(0d0,0d0)
do gs=1,global_num_sites
  rank=local_site_of_global(gs)%rank_
  ls=local_site_of_global(gs)%label_
  tmpmat=(0d0,0d0)
  do l=1,num_links
    do j=1,NMAT
      do i=1,NMAT
        tmpmat=tmpmat - Xi_lambda(i,j,l)*Geta_lambda(:,:,j,i,gs,l)
      enddo
    enddo
  enddo
  call MPI_REDUCE(tmpmat,tmpmat2,NMAT*NMAT,MPI_DOUBLE_COMPLEX, &
    MPI_SUM,0,MPI_COMM_WORLD,IERR)
  if( MYRANK==rank ) XiDinv_eta(:,:,ls)=tmpmat2
enddo
!!
tmpmat2=(0d0,0d0)
do gl=1,global_num_links
  rank=local_link_of_global(gl)%rank_
  ll=local_link_of_global(gl)%label_
  tmpmat=(0d0,0d0)
  do l=1,num_links
    do j=1,NMAT
      do i=1,NMAT
        tmpmat=tmpmat - Xi_lambda(i,j,l)*Glambda_lambda(:,:,j,i,gl,l)
      enddo
    enddo
  enddo
  call MPI_REDUCE(tmpmat,tmpmat2,NMAT*NMAT,MPI_DOUBLE_COMPLEX, &
    MPI_SUM,0,MPI_COMM_WORLD,IERR)
  if( MYRANK==rank ) XiDinv_lambda(:,:,ll)=tmpmat2
enddo
!!
tmpmat2=(0d0,0d0)
do gf=1,global_num_faces
  rank=local_face_of_global(gf)%rank_
  lf=local_face_of_global(gf)%label_
  tmpmat=(0d0,0d0)
  do l=1,num_links
    do j=1,NMAT
      do i=1,NMAT
        tmpmat=tmpmat - Xi_lambda(i,j,l)*Gchi_lambda(:,:,j,i,gf,l)
      enddo
    enddo
  enddo
  call MPI_REDUCE(tmpmat,tmpmat2,NMAT*NMAT,MPI_DOUBLE_COMPLEX, &
    MPI_SUM,0,MPI_COMM_WORLD,IERR)
  if( MYRANK==rank ) XiDinv_chi(:,:,lf)=tmpmat2
enddo

end subroutine calc_Xilink_Dinv
