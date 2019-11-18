!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_Sf_site(Sfsite,Geta_eta,PhiMat)
implicit none

complex(kind(0d0)), intent(out) :: Sfsite
complex(kind(0d0)), intent(in) ::  Geta_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_sites) 
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) trace,tmp
integer ls,gs
integer i,j,k

Sfsite=(0d0,0d0)
trace=(0d0,0d0)
do ls=1,num_sites
  gs=global_site_of_local(ls)
  tmp=(0d0,0d0)
  do k=1,NMAT
    do j=1,NMAT
      do i=1,NMAT
        tmp=tmp&
          +Geta_eta(k,i,j,k,gs,ls)*Phimat(i,j,ls) &
          -Geta_eta(i,k,k,j,gs,ls)*PhiMat(j,i,ls)
      enddo
    enddo
  enddo
  trace = trace + tmp * dcmplx( -0.25d0 * alpha_s(ls) )
enddo
call MPI_REDUCE(trace,Sfsite,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)

Sfsite=Sfsite * dcmplx( overall_factor )


end subroutine calc_Sf_site

