!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate lambda.lambda(U.\bar{phi}.U^-1 + \bar{phi})
subroutine calc_Sf_link2(SfL2, PhiMat, Umat, Glambda_lambda)
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

end subroutine calc_Sf_link2

