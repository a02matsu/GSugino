!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! fermion action link 1
!!    1/2g^2 \sum_l \alpha_l Tr( i\lambda_l D_l \eta )
subroutine calc_Sf_link1(Sflink,Geta_lambda,Umat,PhiMat)
implicit none

complex(kind(0d0)), intent(out) :: Sflink
complex(kind(0d0)), intent(in) ::  Geta_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_links) 
complex(kind(0d0)), intent(in) :: UMat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) trace,tmp
integer ll,ls,gs
integer i,j,k,l

Sflink=(0d0,0d0)
trace=(0d0,0d0)
do ll=1,num_links
  tmp=(0d0,0d0)
!!!!!!!
  ls=link_tip(ll)
  gs=global_site_of_local(ls)
  do l=1,NMAT
    do k=1,NMAT
      do j=1,NMAT
        do i=1,NMAT
          tmp=tmp&
            - Geta_lambda(i,j,k,l,gs,ll)*dconjg(Umat(k,j,ll))*Umat(l,i,ll)&
               * dconjg(U1Rfactor_link(ll)*U1R_ratio(ll))
        enddo
      enddo
    enddo
  enddo
!!!!!!!!
  ls=link_org(ll)
  gs=global_site_of_local(ls)
  do j=1,NMAT
    do i=1,NMAT
      tmp = tmp + Geta_lambda(j,i,i,j,gs,ll)
    enddo
  enddo
  trace = trace + tmp * (0d0,1d0) * dcmplx( alpha_l(ll) )
enddo

call MPI_REDUCE(trace,Sflink,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)

Sflink=Sflink * dcmplx( overall_factor ) 


end subroutine calc_Sf_link1

