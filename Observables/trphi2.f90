!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_trphi2(trphi2,PhiMat)
use parallel
implicit none

!double precision, intent(out) :: trphi2(1:num_sites)
double precision, intent(out) :: trphi2
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

double precision :: tmp
integer :: s,i,j

trphi2=0d0
tmp=0d0
do s=1,num_sites
  do j=1,NMAT
    do i=1,NMAT
      tmp=tmp+dble(PhiMat(i,j,s)*dconjg(PhiMat(i,j,s)))
    enddo
  enddo
enddo
tmp=tmp/(LatticeSpacing*LatticeSpacing*dble(NMAT))

call MPI_REDUCE(tmp,trphi2,1,MPI_DOUBLE_PRECISION, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)


end subroutine calc_trphi2
