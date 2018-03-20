!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! set Remez data
subroutine set_Remez_data
implicit none

integer i

#ifdef PARALLEL
if (MYRANK==0) then
#endif

open(Remez4_FILE, file=Remez_1ovminus4, status='OLD',action='READ')
read(Remez4_FILE,*) N_Remez4
read(Remez4_FILE,*) Remez_min4
read(Remez4_FILE,*) Remez_max4
allocate( Remez_alpha4(0:N_Remez4) )
allocate( Remez_beta4(1:N_Remez4) )
do i=0,N_Remez4 
  read(Remez4_FILE,*) Remez_alpha4(i)
enddo
do i=1,N_Remez4 
  read(Remez4_FILE,*) Remez_beta4(i)
enddo
close(Remez4_FILE)

#ifdef PARALLEL
endif
call MPI_BCAST(N_Remez4,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(Remez_min4,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(Remez_max4,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
if( MYRANK .ne. 0) then
  allocate( Remez_alpha4(0:N_Remez4) )
  allocate( Remez_beta4(1:N_Remez4) )
endif
call MPI_BCAST(Remez_alpha4,N_Remez4+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(Remez_beta4,N_Remez4,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
#endif

!! rescaling
Remez_alpha4(0)=Remez_alpha4(0)*Remez_factor4**(-0.25d0)
do i=1,N_Remez4
  Remez_alpha4(i)=Remez_alpha4(i)*Remez_factor4**(0.75d0)
  Remez_beta4(i)=Remez_beta4(i)*Remez_factor4
enddo



!!
#ifdef PARALLEL
if (MYRANK==0) then
#endif
open(Remez8_FILE, file=Remez_1ov8, status='OLD',action='READ')
read(Remez8_FILE,*) N_Remez8
read(Remez8_FILE,*) Remez_min8
read(Remez8_FILE,*) Remez_max8
allocate( Remez_alpha8(0:N_Remez8) )
allocate( Remez_beta8(1:N_Remez8) )
do i=0,N_Remez8 
  read(Remez8_FILE,*) Remez_alpha8(i)
enddo
do i=1,N_Remez8 
  read(Remez8_FILE,*) Remez_beta8(i)
enddo
close(Remez8_FILE)

#ifdef PARALLEL
endif
call MPI_BCAST(N_Remez8,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(Remez_min8,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(Remez_max8,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
if( MYRANK .ne. 0) then
  allocate( Remez_alpha8(0:N_Remez8) )
  allocate( Remez_beta8(1:N_Remez8) )
endif
call MPI_BCAST(Remez_alpha8,N_Remez8+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(Remez_beta8,N_Remez8,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
#endif

!! rescaling
Remez_alpha8(0)=Remez_alpha8(0)*Remez_factor8**(0.125d0)
do i=1,N_Remez8
  Remez_alpha8(i)=Remez_alpha8(i)*Remez_factor8**(1.125d0)
  Remez_beta8(i)=Remez_beta8(i)*Remez_factor8
enddo
end subroutine set_Remez_data

