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
  !call calc_smallest_and_largest_eigenvalues_of_DdagD(min_eigen,max_eigen,UMAT,PhiMat)
call max_eigen_DdagD(max_eigen,Umat,PhiMat)
call min_eigen_DdagD(min_eigen,Umat,PhiMat)
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
write(*,*) new_config
if( new_config == 0 ) then
  write(output,*) "# configs read from ", trim(Fconfigin)
  if( fix_seed == 0 ) then
    write(output,*) "# random seed is succeeded from the previous simulation"
  elseif( fix_seed == 1 ) then
    write(output,*) "# random seed is fixed to seed=",seed
  elseif( fix_seed == 2 ) then
    write(output,*) "# random seed is determined by the system time"
  endif
elseif( new_config == 1 ) then
  write(output,*) "# cold start: A=0, phi=0"
  if( fix_seed == 1 ) then
    write(output,*) "# random seed is fixed to seed=",seed
  else
    write(output,*) "# random seed is determined by the system time"
  endif
  write(*,*) Phimat, Umat
elseif( new_config == 2 ) then 
  write(output,*) "# cold start(A=0,phi=0) and all accept"
  if( fix_seed == 1 ) then
    write(output,*) "# random seed is fixed to seed=",seed
  else
    write(output,*) "# random seed is determined by the system time"
  endif
elseif( new_config == 3 ) then 
  write(output,*) "# configs read from ", trim(Fconfigin), "and all accept"
  if( fix_seed == 0 ) then
    write(output,*) "# random seed is succeeded from the previous simulation"
  elseif( fix_seed == 1 ) then
    write(output,*) "# random seed is fixed to seed=",seed
  elseif( fix_seed == 2 ) then
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
if( branch_mode==0 ) then 
  write(output,'(A,I5)') "# branch= ",branch_use
else
  write(output,'(A,I5,A,I5)') "# make branch from ",branch_root," to ",new_branch_label
endif
write(output,'(a)') "#"
write(output,'(A,F10.8)') "# Tau for A= ",Tau
write(output,'(A,F10.8)') "# Tau for Phi= ",Tau*ratio_DtauPhi_over_DtauA
write(output,'(A,F10.8)') "# DTau_phi/Dtau_A= ",ratio_DtauPhi_over_DtauA
write(output,'(A,I5)') "# Nfermion= ",Nfermion
write(output,'(A,I5)') "# Nboson= ",Nboson
write(output,'(A,I5)') "# iterations= ", num_ite
write(output,'(a)') "#"
write(output,'(A,I5)') "# m_omega= ",m_omega
write(output,'(A,E12.5,X,a,E12.5,a,E12.5)') "# Remez_factor4= ",Remez_factor4,&
  "range:",Remez_min4*Remez_factor4,"...",Remez_max4*Remez_factor4
write(output,'(A,E12.5,X,a,E12.5,a,E12.5)') "# Remez_factor8= ",Remez_factor8, &
  "range:",Remez_min8*Remez_factor8,"...",Remez_max8*Remez_factor8
if( eval_eigen /= 0 ) then 
  write(output,'(A,E12.5,a,E12.5)') "#          eigenvalue of DD^\dagger:",&
    cdabs(min_eigen),"...",cdabs(max_eigen)
  !dble(min_eigen*conjg(min_eigen)),"...",dble(max_eigen*conjg(max_eigen))
else
  write(output,'(A)') "# omitted the evaluation of the eigenvalues"
endif
write(output,'(a)') "#"

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
  fix_num=7
  write(output,'(a)',advance='no') "5) FB_A/FB_phi, "
  write(output,'(a)',advance='no') "6) FB_phi/FF_phi, "
  write(output,'(a)',advance='no') "7) FB_A/FF_A, "
endif
if( eigen_measurement == 1 ) then
  fix_num=9
  write(output,'(a)',advance='no') "8) min(DDdag), "
  write(output,'(a)',advance='no') "9) max(DDdag), "
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

complex(kind(0d0)) :: min_eigen,max_eigen

call calc_bosonic_action(OBS(1),UMAT,PhiMat)
!! OBS(1)= Sb / (NS+NL)*dimG/2
OBS(1)=OBS(1) * 2d0/dble( (NMAT*NMAT-1)*(global_num_sites+global_num_links) )
call calc_TrX2(OBS(2),PhiMat)
!call calc_PCSC(OBS(3),OBS(4),Umat,PhiMat,1)
if( eigen_measurement == 1 ) then
  call max_eigen_DdagD(max_eigen,Umat,PhiMat)
  call min_eigen_DdagD(min_eigen,Umat,PhiMat)
endif

#ifdef PARALLEL
if( MYRANK == 0 ) then 
#endif
!! for standard output
call write_observables_to(6,ite,total_ite,accept,delta_Ham,CGite,ratio,max_eigen,min_eigen)
!! for output file
call write_observables_to(OUTPUT_FILE,ite,total_ite,accept,delta_Ham,CGite,ratio,max_eigen,min_eigen)
#ifdef PARALLEL
endif
#endif

end subroutine write_observables

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_observables_to(output,ite,total_ite,accept,delta_Ham,CGite,ratio,maxDDdag,minDDdag)
implicit none

integer, intent(in) :: output,ite,accept,total_ite
double precision, intent(in) :: delta_Ham, ratio
integer, intent(in) :: CGite
integer i
complex(kind(0d0)) :: minDDdag, maxDDdag

write(output,'(I6,2X)',advance='no') ite
write(output,'(f12.5,2X)',advance='no') delta_Ham
write(output,'(f6.2,2X)',advance='no') ratio
write(output,'(I6,2X)',advance='no') CGite
if( force_measurement == 1 ) then
  write(output,'(f6.2,2X)',advance='no') bosonic_force_A/bosonic_force_Phi
  write(output,'(f6.2,2X)',advance='no') bosonic_force_Phi/fermionic_force_Phi
  write(output,'(f6.2,2X)',advance='no') bosonic_force_A/fermionic_force_A
endif
if( eigen_measurement == 1 ) then
  write(output,'(E12.5,2X)',advance='no') cdabs(minDDdag)
  write(output,'(E12.5,2X)',advance='no') cdabs(maxDDdag)
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
write(MED_CONF_FILE) test_mode
write(MED_CONF_FILE) force_measurement
write(MED_CONF_FILE) new_config
write(MED_CONF_FILE) fix_seed
write(MED_CONF_FILE) reset_ite
write(MED_CONF_FILE) save_med_step
write(MED_CONF_FILE) save_config_step
write(MED_CONF_FILE) obs_step
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
subroutine read_basic_info_from_medfile(N_MEDFILE)
implicit none

integer, intent(in) :: N_MEDFILE

  write(*,*) NMAT,Nfermion, FB_ratio
  read(N_MEDFILE) NMAT
  write(*,*) NMAT,Nfermion, FB_ratio
  read(N_MEDFILE) LatticeSpacing
  read(N_MEDFILE) SC_FILE_NAME
  read(N_MEDFILE) test_mode
  read(N_MEDFILE) force_measurement
  read(N_MEDFILE) new_config
  read(N_MEDFILE) fix_seed
  read(N_MEDFILE) reset_ite
  read(N_MEDFILE) save_med_step
  read(N_MEDFILE) save_config_step
  read(N_MEDFILE) obs_step
  read(N_MEDFILE) m_omega
  read(N_MEDFILE) phys_mass_square_phi
  read(N_MEDFILE) mass_f
  read(N_MEDFILE) Remez_factor4
  read(N_MEDFILE) Remez_factor8
  read(N_MEDFILE) epsilon
  read(N_MEDFILE) CG_max
  read(N_MEDFILE) num_ite
  read(N_MEDFILE) Nfermion
  read(N_MEDFILE) FB_ratio
  write(*,*) NMAT,Nfermion, FB_ratio
  read(N_MEDFILE) Tau
  read(N_MEDFILE) Fconfigin
  read(N_MEDFILE) Foutput
  read(N_MEDFILE) Fconfigout
  read(N_MEDFILE) Fmedconf

end subroutine read_basic_info_from_medfile



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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_config_from_medfile(UMAT,PhiMat,ite,N_MEDFILE,control)
use global_parameters
use parallel
implicit none

integer, intent(in) :: N_MEDFILE
integer, intent(out) :: ite
complex(kind(0d0)), intent(out) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(out) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_links)
integer, intent(out) :: control
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
integer :: l,ll,s,ls,tag,rank,i,ios

control=0
if( MYRANK == 0 ) then
  !read(N_MEDFILE, END=700) ite
  read(N_MEDFILE,iostat=ios) ite
  if( ios < 0 ) then
    control=1
    call MPI_BCAST(control, 1, MPI_INTEGER,0,MPI_COMM_WORLD,IERR) 
    return
  endif
endif

call MPI_BCAST(control, 1, MPI_INTEGER,0,MPI_COMM_WORLD,IERR)

if( control == 0) then 
  do l=1,global_num_links
    tag=l
    ll=local_link_of_global(l)%label_
    rank=local_link_of_global(l)%rank_
    if( MYRANK == 0 ) then
      read(N_MEDFILE) tmpmat
      !read(N_MEDFILE, END=700) tmpmat
      if( rank == 0 ) then
        Umat(:,:,ll) = tmpmat
      else
        call MPI_SEND(tmpmat,NMAT*NMAT,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,IERR)
      endif
    else
      if( rank == MYRANK ) then
        call MPI_RECV(Umat(:,:,ll),NMAT*NMAT,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,ISTATUS,IERR)
      endif
    endif
  enddo
  
  do s=1,global_num_sites
    tag=global_num_links+s
    ls=local_site_of_global(s)%label_
    rank=local_site_of_global(s)%rank_
    if( MYRANK == 0 ) then
      !read(N_MEDFILE, END=700) tmpmat
      read(N_MEDFILE) tmpmat
      if( rank == 0 ) then
        PhiMat(:,:,ls) = tmpmat
      else
        call MPI_SEND(tmpmat,NMAT*NMAT,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,IERR)
      endif
    else
      if( rank == MYRANK ) then
        call MPI_RECV(PhiMat(:,:,ls),NMAT*NMAT,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,ISTATUS,IERR)
      endif
    endif
  enddo
  !!! syncronize
  call syncronize_links(Umat)
  call syncronize_sites(PhiMat)
  return
else
  !write(*,*) MYRANK, control
  return
endif

!700 control=1
!write(*,*) MYRANK, control
!call MPI_BCAST(control, 1, MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
!return

end subroutine read_config_from_medfile


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
  call make_latestconf_link
  !call system('cd CONFIG; FILE=$(ls inputconf* | tail -1); &
    !LINK="latestconf"; if [ -e "$LINK" ]; then /bin/rm $LINK; &
    !fi; ln -s $FILE $LINK; cd ..' )
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

call make_latestconf_link
!call system('cd CONFIG; FILE=$(ls inputconf* | tail -1); &
  !LINK="latestconf"; if [ -e "$LINK" ]; then /bin/rm $LINK; &
  !fi; ln -s $FILE $LINK; cd ..' )
#endif 

end subroutine write_config_file


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine make_latestconf_link
implicit none
character(500) COMMAND
character(128) tmpc

if( branch_mode==0) then 
  if( branch_use==0) then 
    call system('cd CONFIG; FILE=$(ls inputconf* | tail -1); &
      LINK="latestconf"; if [ -e "$LINK" ]; then /bin/rm $LINK; &
      fi; ln -s $FILE $LINK; cd ..' )
  else
    write(tmpc,*) branch_use
    COMMAND='cd CONFIG'//trim(adjustl(tmpc)) &
      //'; FILE=$(ls inputconf* | tail -1); LINK="latestconf"; &
      if [ -e "$LINK" ]; then /bin/rm $LINK; fi; ln -s $FILE $LINK; cd ..' 
    call system(COMMAND)
  endif
endif
if( branch_mode==1 ) then
  write(tmpc,*) new_branch_label
  COMMAND='cd CONFIG'//trim(adjustl(tmpc))//'; FILE=$(ls inputconf* | tail -1); &
    LINK="latestconf"; if [ -e "$LINK" ]; then /bin/rm $LINK; &
    fi; ln -s $FILE $LINK; cd ..' 
  call system(COMMAND)
endif

end subroutine make_latestconf_link


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_smallest_and_largest_eigenvalues_of_DdagD(min_eigen,max_eigen,UMAT,PhiMat)
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

call calc_eigenvalues_DdagD(eigenvalues,UMAT,PhiMat)

min_eigen=eigenvalues(1)
max_eigen=eigenvalues(sizeD)

end subroutine calc_smallest_and_largest_eigenvalues_of_DdagD
!end module output
