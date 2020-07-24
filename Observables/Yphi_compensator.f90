!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! A = 1/NF sum_f ( 1/N Tr(\phibar^{dimG*eular/2} ) 
subroutine calc_Yphi_compensator(Acomp,PhiMat,Umat)
use parallel
use global_parameters
use matrix_functions, only : matrix_product,matrix_inverse,matrix_power,trace_mm
implicit none

complex(kind(0d0)), intent(out) :: Acomp
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: UMat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: phiinv(1:NMAT,1:NMAT)
complex(kind(0d0)) :: phiinv_p(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Ymat(1:NMAT,1:NMAT)
integer :: ratio,eular
double precision :: radius, phase

complex(kind(0d0)) :: tmp_Acomp, tmp
integer :: lf,ls
integer :: i,j

eular=global_num_sites-global_num_links+global_num_faces 
ratio=(NMAT*NMAT-1)*eular/2
Acomp=(0d0,0d0)
tmp_Acomp=(0d0,0d0)
do lf=1,num_faces
  !! phi^{-dimG Eular/2}
  ls=sites_in_f(lf)%label_(1)
  call matrix_inverse(phiinv,Phimat(:,:,ls))
  call matrix_power(phiinv_p,phiinv,ratio)
  !! Omega
  call Make_face_variable(Uf,lf,UMAT)
  call Make_moment_map_adm(Ymat,Uf)
  Ymat = Ymat * (0d0,0.5d0)*beta_f(lf)*Ymat
  !! tr( phi^{-dimG r} Y )
  call trace_MM(tmp,phiinv_p, Ymat)
  !! 
  tmp_Acomp=tmp_Acomp + tmp/dcmplx(NMAT)
enddo

call MPI_REDUCE(tmp_Acomp,Acomp,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)
  
Acomp=Acomp/dcmplx(dble(global_num_faces))

end subroutine calc_Yphi_compensator


