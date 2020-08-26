!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_regularized_compensator(Acomp,PhiMat,reg)
use parallel
implicit none

complex(kind(0d0)), intent(out) :: Acomp
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
double precision, intent(in) :: reg
complex(kind(0d0)) tmp,A_tmp
double precision :: radius, phase, ratio
integer :: s,i,j,eular


eular=global_num_sites-global_num_links+global_num_faces 
ratio=dble(-(NMAT*NMAT-1)*eular)/4d0 
Acomp=(0d0,0d0)
A_tmp=(0d0,0d0)
do s=1,num_sites
  tmp=(0d0,0d0)
  do i=1,NMAT
    do j=1,NMAT
      if( i .ne. j ) then
        tmp=tmp+PhiMat(i,j,s)*PhiMat(i,j,s)
      else
        tmp=tmp+(PhiMat(i,i,s)+dcmplx(reg))*(PhiMat(i,i,s)+dcmplx(reg))
      endif
    enddo
  enddo
  tmp=(tmp/dcmplx(dble(NMAT)))
  radius=cdabs(tmp)
  phase=atan2(dble(tmp),dble(tmp*(0d0,-1d0)))

  A_tmp=A_tmp + dcmplx(radius**ratio) * cdexp( (0d0,1d0)*dcmplx(phase*ratio) )
enddo

call MPI_REDUCE(A_tmp,Acomp,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)
  
Acomp=Acomp/dcmplx(dble(global_num_sites))

end subroutine calc_regularized_compensator


