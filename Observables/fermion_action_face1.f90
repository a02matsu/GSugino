!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_Sf_face1(Sf_face1,Gchi_chi,PhiMat)
implicit none

complex(kind(0d0)), intent(out) :: Sf_face1
complex(kind(0d0)), intent(in) ::  Gchi_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_faces) 
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) trace,tmp
integer lf,gf,ls
integer i,j,k

Sf_face1=(0d0,0d0)
trace=(0d0,0d0)
do lf=1,num_faces
  gf=global_face_of_local(lf)
  ls=sites_in_f(lf)%label_(1)
  tmp=(0d0,0d0)
  do k=1,NMAT
    do j=1,NMAT
      do i=1,NMAT
        tmp=tmp&
          +Gchi_chi(k,i,j,k,gf,lf)*Phimat(i,j,ls) &
          -Gchi_chi(i,k,k,j,gf,lf)*PhiMat(j,i,ls)
      enddo
    enddo
  enddo
  trace = trace + tmp * dcmplx( - alpha_f(lf) )
enddo
call MPI_REDUCE(trace,Sf_face1,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)

Sf_face1=Sf_face1 * dcmplx( overall_factor )


end subroutine calc_Sf_face1

