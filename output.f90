!module output
!use global_parameters
!implicit none

!integer, parameter :: num_obs=2
!character(10) :: obs_name(1:num_obs)
!data obs_name/ "Sb","TrX2" /
!double precision :: OBS(1:num_obs) ! 1) bosonic action SB

!contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_header(seed,UMAT,PhiMat)
#ifdef PARALLEL
use parallel
#endif
implicit none

integer, intent(in) :: seed 
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) :: min_eigen,max_eigen
integer :: output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compute the eigenvalues of D
if( eval_eigen/=0 ) then 
  call calc_smallset_and_largest_eigenvalues_of_D(min_eigen,max_eigen,UMAT,PhiMat)
endif

!!!!!!!!!!!!!!!!!!!!!!!!!
!! for standard output
output=6 
#ifdef PARALLEL
if( MYRANK == 0 ) then
#endif
call write_header_to(output,min_eigen,max_eigen)
call write_observable_list(output)

!!!!!!!!!!!!!!!!!!!!!!!!!
!! for output file
output=OUTPUT_FILE
call write_header_to(output,min_eigen,max_eigen)
!!!
!if( read_alpha == 0 ) then
!  write(output,*) "# alpha and beta for SC are set to 1.0"
!else
!  write(output,*) "# alpha and beta are set in", trim(ALPHA_BETA)
!endif
!write(output,*)
!!!
if( new_config == 0 ) then
  write(output,*) "# configs read from ", trim(Fconfigin)
  if( fix_seed == 0 ) then
    write(output,*) "# random seed is succeeded from the previous simulation"
  elseif( fix_seed == 1 ) then
    write(output,*) "# random seed is fixed to seed=",seed
  elseif( fix_seed == 2 ) then
    write(output,*) "# random seed is determined by the system time"
  endif
else
  write(output,*) "# new configs"
  if( fix_seed == 1 ) then
    write(output,*) "# random seed is fixed to seed=",seed
  else
    write(output,*) "# random seed is determined by the system time"
  endif
endif
!!!
call write_observable_list(output)

#ifdef PARALLEL
endif
#endif
end subroutine write_header

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_header_to(output,min_eigen,max_eigen)
implicit none

complex(kind(0d0)), intent(in) :: min_eigen,max_eigen
integer, intent(in) :: output

write(output,'(a,a)') "# simplicial complex : ",trim(SC_FILE_NAME)
write(output,'(A,F6.4)') "# lattice spacing (\lambda=1)= ",LatticeSpacing
write(output,'(A,I5)') "# NMAT= ",NMAT
write(output,'(A,E12.5)') "# mass_square_phi= ",mass_square_phi
write(output,'(A,E12.5)') "# phys_mass_square_phi= ",phys_mass_square_phi
write(output,'(A,E12.5)') "# mass_f= ",mass_f
write(output,'(A)') "######################"
write(output,'(A,I5)') "# job_number= ",job_number
write(output,'(a)') "#"
write(output,'(A,F10.8)') "# Tau= ",Tau
write(output,'(A,I5)') "# Nfermion= ",Nfermion
write(output,'(A,I5)') "# Nboson= ",Nboson
write(output,'(A,I5)') "# iterations= ", num_ite
write(output,'(a)') "#"
write(output,'(A,I5)') "# m_omega= ",m_omega
write(output,'(A,E12.5)') "# Remez_factor4= ",Remez_factor4
write(output,'(A,E12.5)') "# Remez_factor8= ",Remez_factor8
if( eval_eigen /= 0 ) then 
  write(output,'(A,E12.5)') "# minimal eigenvalue of DD^\dagger= ",dble(min_eigen*conjg(min_eigen))
  write(output,'(A,E12.5)') "# maximal eigenvalue of DD^\dagger= ",dble(max_eigen*conjg(max_eigen))
else
  write(output,'(A)') "# omitted the evaluation of the eigenvalues"
endif

end subroutine write_header_to


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_observable_list(output)
implicit none

integer, intent(in) :: output
integer :: i,fix_num

!obs_name(1)="Sb"

fix_num=4
write(output,'(a)',advance='no') "# 1) iteration, 2) delta Hamiltonian, 3) max||Uf-1||/max 4) CG ite, "
if( force_measurement == 1 ) then
  fix_num=6
  write(output,'(a)',advance='no') "5) FF/FB_phi, "
  write(output,'(a)',advance='no') "6) FF/FB_A, "
endif
do i=1,num_obs
  write(output,'(I3,a1,a,a1)',advance='no') i+fix_num,")",trim(obs_name(i) ),","
enddo
if( force_measurement /= 0 ) then 
  write(output,'(I3,a)') num_obs+fix_num,") acceptance rate"
else
  write(output,'(I3,a,a,a,a)') num_obs+fix_num,") acceptance rate, &
    FB_phi/FB_A, FF_phi/FB_phi, FF_A/FB_A)"
endif

end subroutine write_observable_list

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_observables(&
    PhiMat,UMAT,ite,accept,delta_Ham,total_ite,CGite,ratio)
!use observables
use SUN_generators, only : trace_MTa
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
double precision, intent(in) :: delta_Ham, ratio
integer, intent(in) :: ite, accept, total_ite, CGite
integer :: output,i,s,a

call calc_bosonic_action(OBS(1),UMAT,PhiMat)
call calc_TrX2(OBS(2),PhiMat)
!call calc_PCSC(OBS(3),OBS(4),Umat,PhiMat,1)

#ifdef PARALLEL
if( MYRANK == 0 ) then 
#endif
!! for standard output
call write_observables_to(6,ite,total_ite,accept,delta_Ham,CGite,ratio)
!! for output file
call write_observables_to(OUTPUT_FILE,ite,total_ite,accept,delta_Ham,CGite,ratio)
#ifdef PARALLEL
endif
#endif

end subroutine write_observables

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_observables_to(output,ite,total_ite,accept,delta_Ham,CGite,ratio)
implicit none

integer, intent(in) :: output,ite,accept,total_ite
double precision, intent(in) :: delta_Ham, ratio
integer, intent(in) :: CGite
integer i

write(output,'(I6,2X)',advance='no') ite
write(output,'(f12.5,2X)',advance='no') delta_Ham
write(output,'(f6.2,2X)',advance='no') ratio
write(output,'(I6,2X)',advance='no') CGite
if( force_measurement == 1 ) then
  write(output,'(f6.2,2X)',advance='no') fermionic_force_Phi/bosonic_force_Phi
  write(output,'(f6.2,2X)',advance='no') fermionic_force_A/bosonic_force_A
endif
do i=1,num_obs
  write(output,'(E12.5,2X)',advance='no') OBS(i)
enddo
if( force_measurement /= 0 ) then 
  write(output,"(F6.4)") dble(accept)/dble(ite-total_ite)
else
  write(output,"(F6.4)",advance='no') dble(accept)/dble(ite-total_ite)
  write(output,"(3f10.3)") &
    bosonic_force_Phi/bosonic_force_A, &
    fermionic_force_Phi/bosonic_force_Phi, &
    fermionic_force_A/bosonic_force_A 
endif


end subroutine write_observables_to

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_basic_info_to_medfile(MED_CONF_FILE)
implicit none

integer, intent(in) :: MED_CONF_FILE

write(MED_CONF_FILE) NMAT
write(MED_CONF_FILE) LatticeSpacing
write(MED_CONF_FILE) SC_FILE_NAME
!write(MED_CONF_FILE) FsubSC
!write(MED_CONF_FILE) ALPHA_BETA
write(MED_CONF_FILE) test_mode
write(MED_CONF_FILE) force_measurement
write(MED_CONF_FILE) new_config
write(MED_CONF_FILE) fix_seed
write(MED_CONF_FILE) reset_ite
!write(MED_CONF_FILE) read_alpha
write(MED_CONF_FILE) save_med_step
write(MED_CONF_FILE) save_config_step
write(MED_CONF_FILE) obs_step
!write(MED_CONF_FILE) seed
write(MED_CONF_FILE) m_omega
write(MED_CONF_FILE) phys_mass_square_phi
write(MED_CONF_FILE) mass_f
write(MED_CONF_FILE) Remez_factor4
write(MED_CONF_FILE) Remez_factor8
write(MED_CONF_FILE) epsilon
write(MED_CONF_FILE) CG_max
write(MED_CONF_FILE) num_ite
write(MED_CONF_FILE) Nfermion
write(MED_CONF_FILE) FB_ratio
write(MED_CONF_FILE) Tau
write(MED_CONF_FILE) Fconfigin
write(MED_CONF_FILE) Foutput
write(MED_CONF_FILE) Fconfigout
write(MED_CONF_FILE) Fmedconf

end subroutine write_basic_info_to_medfile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_config_to_medfile(ite,UMAT,PhiMat)
use parallel
implicit none

integer, intent(in) :: ite
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMAT(1:NMAT,1:NMAT,1:num_necessary_sites)

integer :: l,s,ll,ls,info,tag,rank
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)


!! write present iteration
write(MED_CONF_FILE) ite

!! write Umat
do l=1,global_num_links
  tag=l
  ll=local_link_of_global(l)%label_
  rank=local_link_of_global(l)%rank_

  if( MYRANK == 0 ) then
    if( rank == 0 ) then
      tmpmat=UMAT(:,:,ll)
    else
      call MPI_RECV(tmpmat,NMAT*NMAT,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,ISTATUS,IERR)
    endif
    write(MED_CONF_FILE) tmpmat
  elseif( MYRANK == rank ) then
    call MPI_SEND(Umat(:,:,ll),NMAT*NMAT,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,IERR)
  endif
enddo


!! write Phimat
do s=1,global_num_sites
  tag=global_num_links+s
  ls=local_site_of_global(s)%label_
  rank=local_site_of_global(s)%rank_

  if( MYRANK == 0 ) then
    if( rank == 0 ) then
      tmpmat=PhiMat(:,:,ls)
    else
      call MPI_RECV(tmpmat,NMAT*NMAT,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,ISTATUS,IERR)
    endif
    write(MED_CONF_FILE) tmpmat
  elseif( MYRANK == rank ) then
    call MPI_SEND(PhiMat(:,:,ls),NMAT*NMAT,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,IERR)
  endif
enddo

end subroutine write_config_to_medfile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_config_file(ite,UMAT,PhiMat,state,srepr)
#ifdef PARALLEL
use parallel
#endif
implicit none

integer, intent(in) :: ite
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
type(genrand_state), intent(inout) :: state
type(genrand_srepr), intent(inout) :: srepr
#ifdef PARALLEL
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
integer :: l,s,ll,ls,tag,rank,turn
type(genrand_srepr) :: tmp_srepr

if(MYRANK==0) then
  open(unit=OUT_CONF_FILE,status='replace',file=Fconfigout,action='write',form='unformatted')
  write(OUT_CONF_FILE) job_number
  write(OUT_CONF_FILE) ite-1
endif

do l=1,global_num_links
  tag=l
  ll=local_link_of_global(l)%label_
  rank=local_link_of_global(l)%rank_

  if( MYRANK == 0 ) then
    if( rank == 0 ) then
      tmpmat=UMAT(:,:,ll)
    else
      call MPI_RECV(tmpmat,NMAT*NMAT,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,ISTATUS,IERR)
    endif
    write(OUT_CONF_FILE) tmpmat
  elseif( MYRANK == rank ) then
    call MPI_SEND(Umat(:,:,ll),NMAT*NMAT,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,IERR)
  endif
enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do s=1,global_num_sites
  tag=global_num_links+s
  ls=local_site_of_global(s)%label_
  rank=local_site_of_global(s)%rank_

  if( MYRANK == 0 ) then
    if( rank == 0 ) then
      tmpmat=PhiMat(:,:,ls)
    else
      call MPI_RECV(tmpmat,NMAT*NMAT,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,ISTATUS,IERR)
    endif
    write(OUT_CONF_FILE) tmpmat
  elseif( MYRANK == rank ) then
    call MPI_SEND(PhiMat(:,:,ls),NMAT*NMAT,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,IERR)
  endif
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if( MYRANK == 0 ) then 
  call genrand_init( get=state )
  srepr=state ! mt95では"="がassignmentされている
  write(OUT_CONF_FILE) srepr%repr
  do rank=1,NPROCS-1
    call MPI_RECV(tmp_srepr,1,MPI_GENRAND_SREPR,rank,rank,MPI_COMM_WORLD,ISTATUS,IERR)
    write(OUT_CONF_FILE) tmp_srepr%repr
  enddo
else
  call genrand_init( get=state )
  srepr=state ! mt95では"="がassignmentされている
  call MPI_SEND(srepr,1,MPI_GENRAND_SREPR,0,MYRANK,MPI_COMM_WORLD,IERR)
endif 
if(MYRANK==0) then
  close( OUT_CONF_FILE)
  call system('cd CONFIG; FILE=$(ls inputconf* | tail -1); &
    LINK="latestconf"; if [ -e "$LINK" ]; then /bin/rm $LINK; &
    fi; ln -s $FILE $LINK; cd ..' )
endif

#else
open(unit=OUT_CONF_FILE,status='replace',file=Fconfigout,action='write',form='unformatted')
write(OUT_CONF_FILE) ite-1
write(OUT_CONF_FILE) UMAT
write(OUT_CONF_FILE) PhiMat
call genrand_init( get=state )
srepr=state ! mt95では"="がassignmentされている
write(OUT_CONF_FILE) srepr%repr
close( OUT_CONF_FILE)
call system('cd CONFIG; FILE=$(ls inputconf* | tail -1); &
  LINK="latestconf"; if [ -e "$LINK" ]; then /bin/rm $LINK; &
  fi; ln -s $FILE $LINK; cd ..' )
#endif 

end subroutine write_config_file





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_smallset_and_largest_eigenvalues_of_D(min_eigen,max_eigen,UMAT,PhiMat)
!use observables, only : calc_eigenvalues_Dirac
implicit none

complex(kind(0d0)), intent(out) :: min_eigen, max_eigen
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)),allocatable :: eigenvalues(:)
integer :: sizeD

#ifdef PARALLEL
sizeD=dimG*(global_num_sites+global_num_links+global_num_faces)
#else
sizeD=dimG*(num_sites+num_links+num_faces)
#endif
allocate(eigenvalues(1:sizeD))

call calc_eigenvalues_Dirac(eigenvalues,UMAT,PhiMat)

min_eigen=eigenvalues(1)
max_eigen=eigenvalues(sizeD)

end subroutine calc_smallset_and_largest_eigenvalues_of_D
!end module output
