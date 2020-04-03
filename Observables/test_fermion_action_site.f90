!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine test_Sf_site(Sfsite,PhiMat,PFeta)
use matrix_functions, only : matrix_commutator, trace_mm
implicit none

complex(kind(0d0)), intent(out) :: Sfsite
complex(kind(0d0)), intent(in) ::  PFeta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) trace,tmp
integer ls
integer i,j,k

Sfsite=(0d0,0d0)
trace=(0d0,0d0)
do ls=1,num_sites
  tmp=(0d0,0d0)
  call matrix_commutator(tmpmat,PhiMat(:,:,ls),PFeta(:,:,ls))
  call trace_mm(tmp,PFeta(:,:,ls),tmpmat)
  trace = trace + tmp * dcmplx( -0.25d0 * alpha_s(ls) )
enddo
call MPI_REDUCE(trace,Sfsite,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)

Sfsite=Sfsite * dcmplx( overall_factor )

end subroutine test_Sf_site

