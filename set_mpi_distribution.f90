subroutine set_mpi_distribution(local_site_list)
#ifdef PARALLEL
use parallel
#endif
implicit none

character(50), parameter :: DIST_FILE_NAME="mpi_distribution"
integer, parameter :: DIST_FILE=110
integer :: rank,tmp_num_sites!,tmp_num_links,tmp_num_faces
!integer,allocatable :: tmp_global_site_of_local(:)
integer :: s,l,f,i,part
type(SITE_DIST), intent(out) :: local_site_list(0:NPROCS-1)


if (MYRANK==0) then
  open(DIST_FILE, file=DIST_FILE_NAME, status='old',action='READ')
!! num_sub_SC
  read(DIST_FILE,'()')
  read(DIST_FILE,*) num_sub_SC
!! rank and local sites
  if( num_sub_SC == NPROCS) then
    do part=0,num_sub_SC-1
      read(DIST_FILE,'()') 
      read(DIST_FILE,'()') 
      read(DIST_FILE,*) rank, tmp_num_sites
      local_site_list(part)%rank_=rank 
      local_site_list(part)%num_sites_=tmp_num_sites
      allocate( local_site_list(part)%site_list_(1:tmp_num_sites) )
      read(DIST_FILE,'()') 
      read(DIST_FILE,*) (local_site_list(part)%site_list_(i),i=1,tmp_num_sites)
    enddo
  endif
  close(DIST_FILE)
endif

call MPI_BCAST(num_sub_SC,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
if (num_sub_SC .ne. NPROCS) then
  if(MYRANK==0) write(*,*) "number of core is mismatch."
  call stop_for_test
endif

end subroutine set_mpi_distribution
