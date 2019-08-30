!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_mass_reweight(mass_reweight,PhiMat)
use parallel
implicit none

double precision, intent(out) :: mass_reweight
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

double precision :: trphi2(1:num_sites)
double precision :: tmp
integer :: s,i,j

mass_reweight=0d0
trphi2=0d0
tmp=0d0
do s=1,num_sites
  do j=1,NMAT
    do i=1,NMAT
      trphi2(s)=trphi2(s)+dble(PhiMat(i,j,s)*dconjg(PhiMat(i,j,s)))
    enddo
  enddo
  tmp=tmp+dexp( mass_square_phi*0.5d0*overall_factor*( 1d0 - alpha_s(s) ) )
enddo

call MPI_REDUCE(tmp,mass_reweight,1,MPI_DOUBLE_PRECISION, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)

end subroutine calc_mass_reweight
