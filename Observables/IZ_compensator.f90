subroutine calc_IZ_compensator(Acomp,Umat,PhiMat,Glambda_lambda)
use parallel
use global_parameters
use matrix_functions, only : matrix_power, trace_mm, matrix_3_product, matrix_inverse
implicit none
  
complex(kind(0d0)), intent(out) :: Acomp
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: Glambda_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_links) 
  
complex(kind(0d0)) :: trace, tmp
complex(kind(0d0)) :: phiinv(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
integer :: ratio
integer :: ls,lt,ll,gl
integer :: i,j,k,l

ratio = dimG * (global_num_sites - global_num_links + global_num_faces) / 2

tmp=(0d0,0d0)
do ll=1,num_links
  gl=global_link_of_local(ll)
  ls=link_org(ll)
  lt=link_tip(ll)

  !! phi^(-r-1)
  call matrix_power(phiinv, phimat(:,:,ls), ratio+1 )
  call matrix_inverse(phiinv)

  trace=(0d0,0d0)
  call matrix_3_product(tmpmat,Umat(:,:,ll),phimat(:,:,lt),Umat(:,:,ll),'N','N','C')
  call trace_mm(trace,tmpmat,phiinv)
  tmp = tmp + (0d0,1d0)*trace


  trace=(0d0,0d0)
  do k=1,NMAT
    do j=1,NMAT
      do i=1,NMAT
        trace = trace + Glambda_lambda(i,j,j,k,gl,ll)*phiinv(k,i)
      enddo
    enddo
  enddo
  tmp = tmp + (0d0,1d0)*trace
enddo

Acomp=(0d0,0d0)
call MPI_REDUCE(tmp,Acomp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,IERR)


end subroutine calc_IZ_compensator
  



