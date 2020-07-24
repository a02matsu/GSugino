!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Measurement 
!! % a.out [MEDFILE]
!! 
!! [MEDFILE] must be the name 
!!   "MEDCONF/medconfig_..."
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
use global_parameters
use initialization_calcobs
use simulation
use parallel
implicit none

integer, parameter :: N_MEDFILE=100


character(128) :: MEDFILE
character(128) :: DinvFILE
integer :: control,ite,iarg
integer :: gs,gl,ls,ll,rank,i,j,tag
complex(kind(0d0)), allocatable :: tmp(:,:)

INPUT_FILE_NAME="inputfile"
iarg=iargc()
if( iarg ==0 ) then
  if (MYRANK==0) write(*,*) "use as a.out [MEDFILE]"
  stop
endif
call getarg(1,MEDFILE)
call initialization 
allocate( tmp(1:NMAT,1:NMAT) )

if( MYRANK==0 ) then
  open(N_MEDFILE, file=MEDFILE, status='OLD',action='READ',form='unformatted')
  !write(*,*) Nmat, "   # Nc"
  !write(*,*) LatticeSpacing, "   # lattice spacing a"
endif
do    
  call read_config_from_medfile(Umat,PhiMat,ite,N_MEDFILE,control)
  if( control == 0 ) then 
    !if( MYRANK==0 ) then
      !write(*,*) ite, "   # iteration"
      !write(*,*) "# PhiMat(i,j)"
    !endif
    !!!
    do gs=1,global_num_sites
      ls=local_site_of_global(gs)%label_
      rank=local_site_of_global(gs)%rank_ 
      tag=gs
      if( MYRANK==0 ) then
        if( rank==0 ) then
          tmp=PhiMat(:,:,ls)
        else
          call MPI_RECV(tmp,NMAT*NMAT,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,ISTATUS,IERR)
        endif
      else
        if( MYRANK==rank ) then
          call MPI_SEND(PHiMat(:,:,ls),NMAT*NMAT,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,IERR)
        endif
      endif

      if( MYRANK==0 ) then
        do i=1,NMAT
          do j=1,NMAT
            write(*,'(I5,2X,A,2X,I5,2X,I5,2X,I5,2X,E23.15,2X,E23.15,2X)') &
              ite,"P",gs,i,j,dble(tmp(i,j)),dble( (0d0,-1d0)*tmp(i,j))
          enddo
        enddo
      endif
    enddo
    !!!!!!!!!!!!!!!
    !if( MYRANK==0 ) then
      !write(*,*) "# UMat(i,j)"
    !endif
    do gl=1,global_num_links
      ll=local_link_of_global(gl)%label_
      rank=local_link_of_global(gl)%rank_ 
      tag=global_num_sites+gl
      if( MYRANK==0 ) then
        if( rank==0 ) then
          tmp=UMat(:,:,ll)
        else
          call MPI_RECV(tmp,NMAT*NMAT,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,ISTATUS,IERR)
        endif
      else
        if( MYRANK==rank ) then
          call MPI_SEND(UMat(:,:,ll),NMAT*NMAT,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,IERR)
        endif
      endif

      if( MYRANK==0 ) then
        do i=1,NMAT
          do j=1,NMAT
            write(*,'(I5,2X,A,2X,I5,2X,I5,2X,I5,2X,E23.15,2X,E23.15,2X)') &
              ite,"U",gl,i,j,dble(tmp(i,j)),dble( (0d0,-1d0)*tmp(i,j))
          enddo
        enddo
      endif
    enddo
  else 
    stop
  endif
enddo
end program main


