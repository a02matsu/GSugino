subroutine set_simulation_parameters(local_site_list)
#ifdef PARALLEL
use parallel
#endif
implicit none

integer :: rank,tmp_num_sites!,tmp_num_links,tmp_num_faces
!integer,allocatable :: tmp_global_site_of_local(:)
integer :: s,l,f,i,part
type(SITE_DIST), intent(out) :: local_site_list(0:NPROCS-1)


#ifdef PARALLEL
if (MYRANK==0) then
  open(INPUT_FILE, file=INPUT_FILE_NAME, status='old',action='READ')
!! job_number
  read(INPUT_FILE,'()') 
  read(INPUT_FILE,*) job_number
!! Tau ; trajectory length
  read(INPUT_FILE,'()') 
  read(INPUT_FILE,*) Tau
!! Nfermion 
  read(INPUT_FILE,'()') 
  read(INPUT_FILE,*) Nfermion
!! FB_ratio ; force計算のfermion/boson比
  read(INPUT_FILE,'()') 
  read(INPUT_FILE,*) FB_ratio
!! num_ite
  read(INPUT_FILE,'()') 
  read(INPUT_FILE,*) num_ite
!! save_med_step
  read(INPUT_FILE,'()') 
  read(INPUT_FILE,*) save_med_step
!! save_config_step
  read(INPUT_FILE,'()') 
  read(INPUT_FILE,*) save_config_step
!! obs_step; step to compute observables
  read(INPUT_FILE,'()') 
  read(INPUT_FILE,*) obs_step
!! num_sub_SC
  read(INPUT_FILE,'()') 
  read(INPUT_FILE,*) num_sub_SC
!! rank and local sites
  if( num_sub_SC == NPROCS) then
    do part=0,num_sub_SC-1
      read(INPUT_FILE,'()') 
      read(INPUT_FILE,'()') 
      read(INPUT_FILE,*) rank, tmp_num_sites
      local_site_list(part)%rank_=rank 
      local_site_list(part)%num_sites_=tmp_num_sites
      allocate( local_site_list(part)%site_list_(1:tmp_num_sites) )
      read(INPUT_FILE,'()') 
      read(INPUT_FILE,*) (local_site_list(part)%site_list_(i),i=1,tmp_num_sites)
    enddo
  endif
  close(INPUT_FILE)
endif

call MPI_BCAST(job_number,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(Tau,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(Nfermion,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(FB_ratio,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(num_ite,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(save_med_step,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(save_config_step,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(obs_step,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(num_sub_SC,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
if (num_sub_SC .ne. NPROCS) then
  if(MYRANK==0) write(*,*) "number of core is mismatch."
  call stop_for_test
endif

Nboson = Nfermion * FB_ratio
Dtau_boson = Tau / dble(Nboson*Nfermion)
Dtau_fermion = Tau / dble(Nfermion)

Ntau = Nboson * Nfermion
Dtau_phi = Dtau_boson
Dtau_A = Dtau_boson
Dtau = Dtau_boson

#else
open(INPUT_FILE, file=INPUT_FILE_NAME, status='old',action='READ')
  read(INPUT_FILE,'()') 
  read(INPUT_FILE,*) job_number
  read(INPUT_FILE,'()') 
  read(INPUT_FILE,*) Tau
  read(INPUT_FILE,'()') 
  read(INPUT_FILE,*) Nfermion
  read(INPUT_FILE,'()') 
  read(INPUT_FILE,*) FB_ratio
  read(INPUT_FILE,'()') 
  read(INPUT_FILE,*) num_ite
  read(INPUT_FILE,'()') 
  read(INPUT_FILE,*) save_med_step
  read(INPUT_FILE,'()') 
  read(INPUT_FILE,*) save_config_step
  read(INPUT_FILE,'()') 
  read(INPUT_FILE,*) obs_step
  read(INPUT_FILE,'()') 
close(INPUT_FILE)
Nboson = Nfermion * FB_ratio
Dtau_boson = Tau / dble(Nboson*Nfermion)
Dtau_fermion = Tau / dble(Nfermion)

Ntau = Nboson * Nfermion
Dtau_phi = Dtau_boson
Dtau_A = Dtau_boson
Dtau = Dtau_boson
#endif
end subroutine set_simulation_parameters
