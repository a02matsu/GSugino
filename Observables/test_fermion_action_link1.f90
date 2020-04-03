!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! fermion action link 1
!!    1/2g^2 \sum_l \alpha_l Tr( i\lambda_l D_l \eta )
subroutine test_Sf_link1(Sflink,Umat,PhiMat,PFeta,PFlambda)
use matrix_functions, only : matrix_3_product, trace_mm
implicit none

complex(kind(0d0)), intent(out) :: Sflink
complex(kind(0d0)), intent(in) :: PFeta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: PFlambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: UMat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) trace,tmp
integer ll,ls,gs
integer i,j,k,l

Sflink=(0d0,0d0)
trace=(0d0,0d0)
do ll=1,num_links
  tmpmat=PFeta(:,:,link_org(ll))
  call matrix_3_product(tmpmat,Umat(:,:,ll),PFeta(:,:,link_tip(ll)),Umat(:,:,ll),'N','N','C',(1d0,0d0),'ADD')

  call trace_mm(tmp,PFlambda(:,:,ll),tmpmat)
  trace = trace + tmp * (0d0,1d0) * dcmplx( alpha_l(ll) )
enddo

call MPI_REDUCE(trace,Sflink,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)

Sflink=Sflink * dcmplx( overall_factor ) 

end subroutine test_Sf_link1

