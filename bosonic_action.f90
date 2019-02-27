!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! [\phi,\bar\phi]-part in Sb
subroutine calc_bosonic_action_site(Sb_S,PhiMat)
use parallel
implicit none

double precision, intent(out) :: Sb_S
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

double precision :: tmp

SB_S=0d0
tmp=0d0
call bosonic_action_site(tmp,PhiMat)

call MPI_REDUCE(tmp,Sb_S,1,MPI_DOUBLE_PRECISION, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)

end subroutine calc_bosonic_action_site

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! |D\Phi|^2-part in Sb
subroutine calc_bosonic_action_link(Sb_L,Umat,PhiMat)
use parallel
implicit none

double precision, intent(out) :: Sb_L
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

double precision :: tmp

SB_L=0d0
tmp=0d0
call bosonic_action_link(tmp,Umat,PhiMat)

call MPI_REDUCE(tmp,Sb_L,1,MPI_DOUBLE_PRECISION, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)
end subroutine calc_bosonic_action_link

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! \Omega^2-part in Sb
subroutine calc_bosonic_action_face(Sb_F,Umat)
use parallel
implicit none

double precision, intent(out) :: Sb_F
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)

double precision :: tmp

SB_F=0d0
tmp=0d0
call bosonic_action_face(tmp,Umat)

call MPI_REDUCE(tmp,Sb_F,1,MPI_DOUBLE_PRECISION, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)
end subroutine calc_bosonic_action_face

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_bosonic_action(Sb,UMAT,PhiMat)
!use hamiltonian
implicit none

double precision, intent(out) :: Sb
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

double precision :: SB_S,SB_L,SB_F,SB_local

SB_S=0d0
SB_L=0d0
SB_F=0d0
SB=0d0
call bosonic_action_site(SB_S,PhiMat)
call bosonic_action_link(SB_L,UMAT,PhiMat)
call bosonic_action_face(SB_F,UMAT)

Sb_local=SB_S+SB_L+SB_F

#ifdef PARALLEL
call MPI_REDUCE(SB_local,SB,1,MPI_DOUBLE_PRECISION, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)
#else
SB=SB_local
#endif

end subroutine calc_bosonic_action

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_TrX2(TrX2, PhiMat)
implicit none

double precision, intent(out) :: TrX2
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
double precision :: tmp

integer :: s,i,j

tmp=0d0
do s=1,num_sites
  do i=1,NMAT
    do j=1,NMAT
      tmp=tmp+ dble( PhiMat(i,j,s)*conjg( PhiMat(i,j,s) ) )
    enddo
  enddo
enddo
#ifdef PARALLEL
call MPI_REDUCE(tmp,TrX2,1,MPI_DOUBLE_PRECISION, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)
#else
TrX2=tmp
#endif

TrX2 = TrX2 / (2d0 * dble(global_num_sites) )

end subroutine calc_TrX2


