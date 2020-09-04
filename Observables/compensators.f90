!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_trace_compensator(Acomp,PhiMat)
use parallel
implicit none

complex(kind(0d0)), intent(out) :: Acomp
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
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
      tmp=tmp+PhiMat(i,j,s)*PhiMat(i,j,s)
    enddo
  enddo
  tmp=(tmp/dcmplx(dble(NMAT)))
  radius=cdabs(tmp)
  !phase=atan2(dble(tmp),dble(tmp*(0d0,-1d0)))
  phase=atan2(dble(tmp*(0d0,-1d0)),dble(tmp))

  A_tmp=A_tmp + dcmplx(radius**ratio) * cdexp( (0d0,1d0)*dcmplx(phase*ratio) )
enddo

call MPI_REDUCE(A_tmp,Acomp,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)
  
Acomp=Acomp/dcmplx(dble(global_num_sites))

end subroutine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_VM_compensator(Acomp,PhiMat)
use parallel
implicit none

complex(kind(0d0)), intent(out) :: Acomp
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: Atmp,cartan(1:NMAT-1),VMdet
integer :: s,i,j,eular

eular=global_num_sites-global_num_links+global_num_faces
Acomp=(0d0,0d0)
Atmp=(0d0,0d0)

do s=1, num_sites
  call calc_cartan(cartan,PhiMat(:,:,s))
  VMdet=(1d0,0d0)
  do I=1,NMAT-2
    do j=I+1,NMAT-1
      VMdet=VMdet*(cartan(i)-cartan(j))
    enddo
  enddo
  do i=1,NMAT-1
    VMdet=VMdet * cartan(i)**0.5d0
  enddo
  VMdet=VMdet**(-eular)

  Atmp=Atmp+VMdet
enddo

call MPI_REDUCE(Atmp,Acomp,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)
  
Acomp=Acomp/dcmplx(dble(global_num_sites))

end subroutine calc_VM_compensator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Cartan elements
!!   cartan(1:NMAT-1)
!!   MAT(1:NMAT,1:NMAT)
subroutine calc_cartan(cartan,MAT)
use SUN_generators, only : trace_MTa
implicit none

complex(kind(0d0)), intent(out) :: cartan(:)
complex(kind(0d0)), intent(in) :: MAT(:,:)
integer :: n,i

n=size(MAT,1)

do i=1,n-1
  call trace_MTa(cartan(i),MAT,i,n)
enddo

end subroutine calc_cartan




