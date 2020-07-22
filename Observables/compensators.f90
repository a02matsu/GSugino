!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! A = 1/NS sum_s ( 1/N Tr(\bar(\phi)^2) )^{dimG*\chi/4}
!subroutine calc_phibar_compensator(Acomp,PhiMat)
!use parallel
!use global_parameters
!use matrix_functions, only : matrix_product
!implicit none
!
!complex(kind(0d0)), intent(out) :: Acomp
!complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
!double precision :: ratio,eular
!double precision :: radius, phase
!
!complex(kind(0d0)) :: tmp_Acomp, tmp
!integer :: ls
!integer :: i,j
!
!eular=global_num_sites-global_num_links+global_num_faces 
!ratio=dble((NMAT*NMAT-1)*eular)/4d0 
!Acomp=(0d0,0d0)
!tmp_Acomp=(0d0,0d0)
!do ls=1,num_sites
!  tmp=(0d0,0d0)
!  do i=1,NMAT
!    do j=1,NMAT
!      tmp = tmp + dconjg(PhiMat(i,j,ls))*dconjg(PHiMat(j,i,ls))
!    enddo
!  enddo
!  tmp=(tmp/dcmplx(dble(NMAT)))
!  radius=cdabs(tmp)
!  phase=atan2(dble(tmp),dble(tmp*(0d0,-1d0)))
!
!  tmp_Acomp=tmp_Acomp + dcmplx(radius**ratio) * cdexp( (0d0,1d0)*dcmplx(phase*ratio) )
!enddo
!
!call MPI_REDUCE(tmp_Acomp,Acomp,1,MPI_DOUBLE_COMPLEX, &
!  MPI_SUM,0,MPI_COMM_WORLD,IERR)
!  
!Acomp=Acomp/dcmplx(dble(global_num_sites))
!
!end subroutine calc_phibar_compensator



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_face_compensator(Acomp,Umat,PhiMat,Geta_chi)
use parallel
use global_parameters
use matrix_functions, only : matrix_product
implicit none

complex(kind(0d0)), intent(out) :: Acomp
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: Geta_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_faces) 

complex(kind(0d0)) :: tmp_Acomp, tmp
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: phibar_p(1:NMAT,1:NMAT,0:dimG+2)
complex(kind(0d0)) :: Omega(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT)
integer :: lf,ls,gs,gf
integer :: i,j,k,l,p

Acomp=(0d0,0d0)
tmp_Acomp=(0d0,0d0)
do lf=1, num_faces
  ls=sites_in_f(lf)%label_(1)
  gf=global_face_of_local(lf)
  gs=global_sites_in_f(gf)%label_(1)

  !! phibar_p = \bar(\PhiMat)^p
  phibar_p(:,:,0)=(0d0,0d0)
  do i=1,NMAT
    phibar_p(i,i,0)=(1d0,0d0)
  enddo
  !!
  do j=1,NMAT
    do i=1,NMAT
      phibar_p(i,j,1)=dconjg(PhiMat(j,i,ls))
    enddo
  enddo
  !!1
  do k=2,dimG+2
    call matrix_product(phibar_p(:,:,k),phibar_p(:,:,k-1),PhiMat(:,:,ls),'N','C')
  enddo

  tmp=(0d0,0d0)
  do p=0,dimG+1
    do l=1,NMAT
      do k=1,NMAT
        do j=1,NMAT
          do i=1,NMAT
            tmp = tmp + phibar_p(i,j,dimG+1-p)*phibar_p(k,l,p)&
              *Geta_chi(j,k,l,i,gs,lf)
          enddo
        enddo
      enddo
    enddo
  enddo

  !! Omega
  call Make_face_variable(Uf,lf,UMAT)
  call Make_moment_map_adm(Omega,Uf)
  do j=1,NMAT
    do i=1,NMAT
      tmp = tmp + (0d0,0.5d0)*beta_f(lf)*phibar_p(i,j,dimG+2)*Omega(j,i)
    enddo
  enddo

  tmp_Acomp=tmp_Acomp + alpha_f(lf)*tmp
enddo

call MPI_REDUCE(tmp_Acomp,Acomp,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)
  
Acomp=Acomp/dcmplx(dble(global_num_faces))

end subroutine calc_face_compensator

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
  phase=atan2(dble(tmp),dble(tmp*(0d0,-1d0)))

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




