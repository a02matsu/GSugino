!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_trphi2(trphi2,PhiMat)
use parallel
implicit none

double precision, intent(out) :: trphi2(1:num_sites)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) :: tmp
integer :: s,i,j

trphi2=(0d0,0d0)
do s=1,num_sites
  do j=1,NMAT
    do i=1,NMAT
      trphi2(s)=trphi2(s)+dble(PhiMat(i,j,s)*dconjg(PhiMat(i,j,s)))
    enddo
  enddo
enddo
trphi2=trphi2/dble(NMAT)

end subroutine calc_trphi2

