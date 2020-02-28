!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine test_hamiltonian(UMAT,PhiMat,seed)
use SUN_generators, only : make_traceless_matrix_from_modes
use mt95
use Dirac_operator
implicit none

complex(kind(0d0)), intent(inout) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(inout) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
integer, intent(inout) :: seed

complex(kind(0d0)) :: P_PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: P_AMat(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: PF_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: PF_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: PF_chi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)) :: PhiMat_BAK(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: UMAT_BAK(1:NMAT,1:NMAT,1:num_necessary_links)
double precision :: Hold,Hnew
integer :: n
integer :: CGite,info,s,l,f,a
complex(kind(0d0)) :: tmp,tmp2
type(genrand_state) :: state
type(genrand_srepr) :: srepr

integer :: t_start, t_end, t_rate, t_max
double precision :: diff,ratio
integer gl,s1,s2,ll,rank

!complex(kind(0d0)) :: QS

integer i,j

#ifdef PARALLEL
if ( MYRANK == 0 ) then
#endif
  write(*,'(a)') " # test hamiltonian"
  write(*,'(a,I8)') " # NMAT=",NMAT
  write(*,'(a,I4)') " # fix_seed=",fix_seed
  !write(*,'(a,I2)') " # read_alpha=",read_alpha
  write(*,'(a,I3)') " # m_omega=",m_omega
  write(*,*) "# Tau=",Tau
  write(*,*) "# Nfermion=",Nfermion
  write(*,*) "# Nboson=",Nboson
  if( p_mass == 0 ) then
    write(*,*) "# mass_f=",mass_f
  else
    write(*,*) "# mass_f= 0"
  endif
  if( pf==1 ) then 
    write(*,*) "# only boson"
  else
    write(*,*) "# p1..p5=",p1,p2,p3,p4,p5
  endif
#ifdef PARALLEL
endif
#endif

!! backup
PhiMat_BAK=PhiMat
!Phi_BAK=Phi
UMAT_BAK=UMAT



call check_Dirac(Umat,Phimat)
call check_distance(info,ratio,Umat)
if(info .ne.0) then 
  write(*,*) "Umat is out of range.",ratio
  call stop_for_test
endif
!if( MYRANK==0 ) then
!do gl=1,global_num_links
!  write(*,*) gl, global_U1R_ratio(gl)
!enddo
!endif
!call stop_for_test

call check_QS(Umat,PhiMat)

call genrand_init( get=state )
n=1
do while ( n > 0 ) 
Ntau=Ntau*2
Dtau_phi=Dtau_phi/2d0
Dtau_A=Dtau_A/2d0

Nfermion=Nfermion*2
Dtau_boson=Dtau_boson/2d0
Dtau_fermion=Dtau_fermion/2d0

!! fix the random seed
!if( fix_seed == 1 ) then
  !seed=12345
!call genrand_init( put=state )
call genrand_init( put=seed )
!endif
!! set random momentum
call set_randomP_single(P_AMat,P_PhiMat)


! produce pseudo-fermion
if( pf==0 ) call make_pseudo_fermion(PF_eta,PF_lambda,PF_chi,UMAT,PhiMat)

call system_clock(t_start)
call Make_Hamiltonian(Hold,CGite,info,UMAT,PhiMat,PF_eta,PF_lambda,PF_chi,P_AMat,P_PhiMat)
!call check_local_vals(PhiMat,UMat)
!call stop_for_test

!if(MYRANK==0) write(*,*) "H=",Hold
!! molecular evolution
if( force_measurement == 1 ) then
  b_phi_count=0
  f_phi_count=0
  b_A_count=0
  f_A_count=0
  bosonic_force_Phi=0d0
  fermionic_force_Phi=0d0
  bosonic_force_A=0d0
  fermionic_force_A=0d0
endif
!call molecular_evolution_Omelyan(UMAT,PhiMat,PF_eta,PF_lambda,PF_chi,P_AMat,P_PhiMat,info)
call molecular_evolution_multistep(UMAT,PhiMat,PF_eta,PF_lambda,PF_chi,P_AMat,P_PhiMat,info)
if( force_measurement == 1 ) then
  bosonic_force_Phi= bosonic_force_Phi / dble(b_phi_count)
  fermionic_force_Phi= fermionic_force_Phi / dble(f_phi_count)
  bosonic_force_A= bosonic_force_A / dble(b_A_count)
  fermionic_force_A= fermionic_force_A / dble(f_A_count)
endif

call Make_Hamiltonian(Hnew,CGite,info,UMAT,PhiMat,PF_eta,PF_lambda,PF_chi,P_AMat,P_PhiMat)
!! write the simulation time
call system_clock(t_end, t_rate, t_max)
if ( t_end < t_start ) then
  diff = dble((t_max - t_start) + t_end + 1) / dble(t_rate)
else
  diff = dble(t_end - t_start) / dble(t_rate)
endif
!! metropolice
#ifdef PARALLEL
if( MYRANK == 0 ) then 
#endif 
write(*,'(I8,e20.4,5f15.4)') Nfermion, dabs(Hnew-Hold)/Hold,diff, &
  fermionic_force_Phi/bosonic_force_Phi, &
  fermionic_force_A/bosonic_force_A
#ifdef PARALLEL
endif
#endif 

!! return to the original values
PhiMat=PhiMat_bak
UMAT=UMAT_bak

!srepr=state ! mt95では"="がassignmentされている
!write(*,*) srepr%repr
n=n+1
!if( n>9 ) call stop_for_test
enddo


end subroutine test_hamiltonian


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to test Q-transformation of S
subroutine check_QS(Umat,PhiMat)
use matrix_functions, only : matrix_commutator, matrix_3_product
use Dirac_operator, only : Prod_Dirac
!use simulation, only : make_bosonic_force_nomass
#ifdef PARALLEL
use parallel
#endif
implicit none

complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)) :: Qeta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: Qlambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: Qchi(1:NMAT,1:NMAT,1:num_necessary_faces)

complex(kind(0d0)) :: DQeta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: DQlambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: DQchi(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: Omega(1:NMAT,1:NMAT,1:num_faces)

complex(kind(0d0)) :: Bforce_s(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: Bforce_l(1:NMAT,1:NMAT,1:num_links)

complex(kind(0d0)) :: ctmp
double precision :: tmp,QS
integer :: info,s,l,f,i,j,triger

!! fermion part
call make_Qfermion(Qeta,Qlambda,Qchi,Omega,Umat,PhiMat)
call Prod_Dirac(DQeta,DQlambda,DQchi,Qeta,Qlambda,Qchi,UMAT,Phimat)

!! boson part
call make_bosonic_force_nomass(Bforce_s,Bforce_l,Umat,PhiMat)
do s=1,num_sites
  Bforce_s(:,:,s)=Bforce_s(:,:,s)*dconjg(U1Rfactor_site(s))
enddo
do l=1,num_links
  Bforce_l(:,:,l)=Bforce_l(:,:,l)*U1Rfactor_site(link_org(l))
enddo


! Q^2 \Omega を care する
!do f=1,num_faces
!  call matrix_commutator(tmpmat,PhiMat(:,:,sites_in_f(f)%label_(1)),Omega(:,:,f))
!  DQchi(:,:,f)=DQchi(:,:,f)+(0d0,1d0)*dcmplx(beta_f(f))*tmpmat
!enddo

if(MYRANK==0) write(*,*) "# QS = 0 ?"
!write(*,*) DQchi
QS=0d0
tmp=0d0
!call make_bosonic_force_phi_site(Bforce_s,PhiMat)
do s=1,num_sites
  tmp=0d0
  do i=1,NMAT
    do j=1,NMAT
      ctmp=&
        - DQeta(i,j,s) &
        + dconjg(Bforce_s(j,i,s))!*U1Rfactor_site(s)
      tmp=tmp+dble( ctmp*dconjg(ctmp) )
    enddo
    !write(*,*) global_site_of_local(s), cdabs(ctmp), &
    !  dble(cdlog(site_U1Rfactor(s))/LatticeSpacing*(0d0,-1d0))
  enddo
enddo
call MPI_REDUCE(tmp,QS,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
if(MYRANK==0) write(*,*) "#   site:",QS

QS=0d0
tmp=0d0
do l=1,num_links
  do i=1,NMAT
    do j=1,NMAT
      !tmpmat(i,j) = &
      ctmp = &
        - DQlambda(i,j,l) & !*site_U1Rfactor(link_org(l))&
        + Bforce_l(i,j,l) !*site_U1Rfactor(link_org(l))
      tmp=tmp+dble( ctmp*dconjg(ctmp) )
    enddo
  enddo
  !write(*,*) "(L)",global_link_of_local(l), cdabs(ctmp)
enddo
call MPI_REDUCE(tmp,QS,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
if(MYRANK==0) write(*,*) "#   link:",QS

QS=0d0
tmp=0d0
do f=1,num_faces
  !tmp=0d0
  do i=1,NMAT
    do j=1,NMAT
      ctmp=-DQchi(j,i,f) !*site_U1Rfactor(sites_in_f(f)%label_(1))
      tmp=tmp+dble( ctmp*dconjg(ctmp) )
    enddo
  enddo
  !write(*,*) "(F)",global_face_of_local(f), cdabs(ctmp)
enddo
call MPI_REDUCE(tmp,QS,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
if(MYRANK==0) write(*,*) "#   face:",QS

end subroutine check_QS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to make Q\Psi
subroutine make_Qfermion(Qeta,Qlambda,Qchi,Omega,Umat,PhiMat)
use global_subroutines
use matrix_functions, only : matrix_commutator, matrix_3_product,matrix_power
use parallel
implicit none

complex(kind(0d0)), intent(out) :: Qeta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(out) :: Qlambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(out) :: Qchi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), intent(out) :: Omega(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
!complex(kind(0d0)), intent(in) :: PF_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)) :: comm(1:NMAT,1:NMAT)

integer :: s,l,f,i,j
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Ufm(1:NMAT,1:NMAT)

Qeta=(0d0,0d0)
do s=1,num_sites
  call matrix_commutator(Qeta(:,:,s),PhiMat(:,:,s),PhiMat(:,:,s),'N','C')
  Qeta(:,:,s)=Qeta(:,:,s)*U1Rfactor_site(s)
enddo

Qlambda=(0d0,0d0)
do l=1,num_links
  Qlambda(:,:,l)=(0d0,-1d0)*PhiMat(:,:,link_org(l))
  call matrix_3_product(Qlambda(:,:,l),&
    Umat(:,:,l),PhiMat(:,:,link_tip(l)),Umat(:,:,l),&
    'N','N','C',(0d0,1d0)*U1Rfactor_link(l)**2d0*U1R_ratio(l)**2d0,'ADD')
  Qlambda(:,:,l)=Qlambda(:,:,l)*U1Rfactor_site( link_org(l) )
enddo

Qchi=(0d0,0d0)
do f=1,num_faces
  call Make_face_variable(Uf,f,UMAT)
  if(m_omega == 0) then 
    call Make_moment_map0(Omega(:,:,f),Uf)
  elseif(m_omega == -1) then
    call Make_moment_map_adm(Omega(:,:,f),Uf)
  else
    call matrix_power(Ufm,Uf,m_omega)
    call Make_moment_map(Omega(:,:,f),Ufm)
  endif
  Qchi(:,:,f)=(0d0,-0.5d0)*dcmplx(beta_f(f))*Omega(:,:,f)
  Qchi(:,:,f)=Qchi(:,:,f)*U1Rfactor_site(sites_in_f(f)%label_(1))
enddo

#ifdef PARALLEL
call syncronize_sites(Qeta)
call syncronize_links(Qlambda)
call syncronize_faces(Qchi)
#endif

end subroutine make_Qfermion



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine HybridMonteCarlo(UMAT,PhiMat,seed,total_ite)
!use output
use SUN_generators, only : make_traceless_matrix_from_modes,trace_MTa
use matrix_functions, only : make_matrix_traceless
use mt95
implicit none

complex(kind(0d0)), intent(inout) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(inout) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
integer, intent(inout) :: seed
integer, intent(inout) :: total_ite

complex(kind(0d0)) :: P_AMat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: P_PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: PF_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: PF_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: PF_chi(1:NMAT,1:NMAT,1:num_necessary_faces)

complex(kind(0d0)) PhiMat_BAK(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) UMAT_BAK(1:NMAT,1:NMAT,1:num_necessary_links)
double precision Hold, Hnew, ratio, ave_dHam
integer :: ite
integer :: accept
type(genrand_state) :: state
type(genrand_srepr) :: srepr

integer :: t_start, t_end, t_rate, t_max
integer :: CGite1, CGite2, info1, info2, info
double precision :: diff

complex(kind(0d0)) :: min_eigen,max_eigen

integer s,a,l
complex(kind(0d0)) :: tmp

!! prepare intermediate file
#ifdef PARALLEL
if ( MYRANK == 0  ) then 
#endif
if( save_med_step /= 0 ) then
  open(unit=MED_CONF_FILE,status='replace',file=Fmedconf,action='write',form='unformatted')
endif
  !call write_basic_info_to_medfile(MED_CONF_FILE)
#ifdef PARALLEL
endif
#endif

!! prepare output file
#ifdef PARALLEL
if ( MYRANK == 0  ) then 
#endif
if( obs_step /= 0 ) then
  open(unit=OUTPUT_FILE,status='replace',file=Foutput,action='write')
endif
  write(*,*) "# start Monte Carlo simulation"
#ifdef PARALLEL
endif
#endif

ave_dHam=0d0
accept=0
call write_header(seed,UMAT,PhiMat)
!! measure the time used in this simulation
call system_clock(t_start)
do ite=total_ite+1,total_ite+num_ite
  !! set random momentuet
  call set_randomP(P_AMat,P_PhiMat)
  !! produce pseudo-fermion
  if(pf==0) call make_pseudo_fermion(PF_eta,PF_lambda,PF_chi,UMAT,PhiMat)
  !do l=1,num_links
    !call Make_traceless_matrix_from_modes(P_Amat(:,:,l), NMAT, dcmplx(P_A(:,l)))
  !enddo
  !! calculate Hamiltonian 
  call Make_Hamiltonian(Hold,CGite1,info1,UMAT,PhiMat,PF_eta,PF_lambda,PF_chi,P_AMat,P_PhiMat)

  !! backup
  PhiMat_BAK=PhiMat
  !Phi_BAK=Phi
  UMAT_BAK=UMAT
  !! molecular evolution
  !call mat_to_vec(PF,PF_eta,PF_lambda,PF_chi)
  if( force_measurement == 1 ) then
    b_phi_count=0
    f_phi_count=0
    b_A_count=0
    f_A_count=0
    bosonic_force_Phi=0d0
    fermionic_force_Phi=0d0
    bosonic_force_A=0d0
    fermionic_force_A=0d0
  endif
  !call molecular_evolution_Omelyan(UMAT,PhiMat,PF_eta,PF_lambda,PF_chi,P_AMat,P_PhiMat,info)
  info=0 ! check if CG successes during molecular evolution
  call molecular_evolution_multistep(UMAT,PhiMat,PF_eta,PF_lambda,PF_chi,P_AMat,P_PhiMat,info)

  if( force_measurement == 1 ) then
    bosonic_force_Phi= bosonic_force_Phi / dble(b_phi_count)
    fermionic_force_Phi= fermionic_force_Phi / dble(f_phi_count)
    bosonic_force_A= bosonic_force_A / dble(b_A_count)
    fermionic_force_A= fermionic_force_A / dble(f_A_count)
  endif
  !call vec_to_mat(PF_eta,PF_lambda,PF_chi,PF)
  do s=1,num_necessary_sites
    call make_matrix_traceless(PhiMat(:,:,s))
  enddo
  if( info == 1 .and. new_config <= 1 ) then
    PhiMat=PhiMat_BAK
    UMAT=UMAT_BAK
#ifdef PARALLEL
    if(MYRANK==0) then
#endif 
    write(*,*) "### CAUTION: CG iterations reaches to the maximal during molecular evolution."
#ifdef PARALLEL
    endif
#endif 
    
  else
    !call check_vacuum(UMAT)
    !! calculate Hamiltonian 
    call Make_Hamiltonian(Hnew,CGite2,info2,UMAT,PhiMat,PF_eta,PF_lambda,PF_chi,P_AMat,P_PhiMat)
    !! metropolice
#ifdef PARALLEL
    if(MYRANK==0) then
#endif 
    ave_dHam=ave_dHam+abs(Hnew-Hold)
#ifdef PARALLEL
    endif
#endif 
    call Metropolice_test(Hnew-Hold,PhiMat_Bak,UMAT_BAK,PhiMat,UMAT,accept)
    !!
    !! check distance of Uf from the origin
    call check_distance(info,ratio,UMAT)
    if ( info .ne. 0 ) then 
      write(*,*) "!!! CAUTION: plaquette variables are out of proper region !!!!"
    endif
  endif
  !! write out the configuration
  if ( save_med_step/=0 .and. mod(ite,save_med_step) == 0 ) then
#ifdef PARALLEL
     call write_config_to_medfile(ite,UMAT,PhiMat)
     !if( MYRANK == 0 ) then
     !  write(MED_CONF_FILE) ite
     !  write(MED_CONF_FILE) NPROCS
     !  write(MED_CONF_FILE) UMAT
     !  write(MED_CONF_FILE) PHIMAT
     !endif
#else
     write(MED_CONF_FILE) ite
     !write(MED_CONF_FILE) 1
     write(MED_CONF_FILE) UMAT
     write(MED_CONF_FILE) PHIMAT
#endif
   endif
   !! write out the observables 
   if ( obs_step/=0 .and. mod(ite,obs_step) == 0 ) then
     !if( eigen_measurement == 1 ) then
      !call max_eigen_DdagD(max_eigen,Umat,PhiMat)
      !call min_eigen_DdagD(min_eigen,Umat,PhiMat)
    !endif
     call write_observables(&
       PhiMat,UMAT,ite,accept,Hnew-Hold,total_ite,CGite1,ratio)
   endif
   !! write out config file
   if ( save_config_step/=0 .and. mod(ite,save_config_step) == 0 ) then
     call write_config_file(ite,UMAT,PhiMat,state,srepr)
   endif 
enddo

!! write the simulation time
#ifdef PARALLEL
if(MYRANK==0) then
#endif 
call system_clock(t_end, t_rate, t_max)
if ( t_end < t_start ) then
  diff = dble((t_max - t_start) + t_end + 1) / dble(t_rate)
else
  diff = dble(t_end - t_start) / dble(t_rate)
endif
ave_dHam = ave_dHam/dble(num_ite)
write(OUTPUT_FILE,*) "# Total time: ",diff,"[s]"
write(OUTPUT_FILE,*) "# Time/ite: ",diff/dble(num_ite),"[s]"
write(OUTPUT_FILE,*) "# Time/ite/Ntau: ",diff/dble(num_ite*Ntau),"[s]"
write(OUTPUT_FILE,*) "# ave_dHam=",ave_dHam
write(OUTPUT_FILE,*) "# acceptance=",dble(accept)/dble(num_ite)
write(*,*) "# Total time: ",diff,"[s]"
write(*,*) "# Step time: ",diff/dble(num_ite),"[s]"
write(*,*) "# Time/ite/Ntau: ",diff/dble(num_ite*Ntau),"[s]"
write(*,*) "# ave_dHam=",ave_dHam
write(*,*) "# acceptance=",dble(accept)/dble(num_ite)

#ifdef PARALLEL
endif
#endif 

call write_config_file(ite,UMAT,PhiMat,state,srepr)

#ifdef PARALLEL
if(MYRANK==0) then
#endif 
close( MED_CONF_FILE )
close( OUTPUT_FILE )
#ifdef PARALLEL
endif
#endif 

end subroutine HybridMonteCarlo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_randomP(P_AMat,P_PhiMat)
use SUN_generators, only : make_traceless_matrix_from_modes, make_sun_generators
#ifdef PARALLEL
use parallel
#endif
implicit none

complex(kind(0d0)), intent(inout) :: P_AMat(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(inout) :: P_PhiMat(1:NMAT,1:NMAT,1:num_sites)
double precision g_site(1:2*dimG,1:num_sites)
double precision g_rlink(1:dimG,1:num_links)
integer s,l,f,a
integer i,j
complex(kind(0d0)) :: TMAT(1:NMAT,1:NMAT,1:NMAT**2-1)
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT), trace
integer num
double precision :: rtmp


call make_SUN_generators(TMAT,NMAT)

call BoxMuller2(g_site,num_sites*dimG)
call BoxMuller2(g_rlink,num_links*dimG/2)

P_PHIMAT=(0d0,0d0)
do s=1,num_sites
  do a=1,dimG
    P_PHIMAT(:,:,s)=P_PhiMat(:,:,s) + &
      (dcmplx(g_site(2*a-1,s))+(0d0,1d0)*dcmplx(g_site(2*a,s)))*dcmplx(dsqrt(0.5d0))&
      *TMAT(:,:,a)
  enddo
enddo

! random number must be sufficiently small
P_AMAT=(0d0,0d0)
do l=1,num_links
  do a=1,dimG
    P_AMAT(:,:,l)=P_AMAT(:,:,l)+g_rlink(a,l)*TMAT(:,:,a)
  enddo
enddo

end subroutine set_randomP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_randomP_single(P_AMat,P_PhiMat)
use SUN_generators, only : make_traceless_matrix_from_modes, make_sun_generators
#ifdef PARALLEL
use parallel
#endif
implicit none

complex(kind(0d0)), intent(inout) :: P_AMat(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(inout) :: P_PhiMat(1:NMAT,1:NMAT,1:num_sites)
!double precision :: P_A(1:dimG)
!complex(kind(0d0)) :: P_Phi(1:dimG)
double precision, allocatable :: gauss(:)
double precision, allocatable :: g_rsite(:)
double precision, allocatable :: g_rlink(:,:)
integer s,l,f,a
integer i,j
complex(kind(0d0)) :: TMAT(1:NMAT,1:NMAT,1:NMAT**2-1)
complex(kind(0d0)) :: G_PhiMat(1:NMAT,1:NMAT,1:global_num_sites)
complex(kind(0d0)) :: G_AMat(1:NMAT,1:NMAT,1:global_num_links)
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT), trace
integer num,rank,ls,ll
double precision :: rtmp


call make_SUN_generators(TMAT,NMAT)

!! テスト用に、シングルコアと同じconfigurationを用意する。
#ifdef PARALLEL
!call genrand_real3(rsite)
if( MYRANK == 0 ) then 
#endif
  allocate(g_rsite(1:2*NMAT*NMAT*global_num_sites))
  allocate(g_rlink(1:dimG,1:global_num_links))

  call genrand_real3(g_rlink)
  call BoxMuller2(g_rsite,global_num_sites*NMAT*NMAT)

  num=0
  do s=1,global_num_sites
    do i=1,NMAT
      do j=1,NMAT
        num=num+1
        G_PHIMAT(i,j,s)=dcmplx(g_rsite(2*num-1))+(0d0,1d0)*dcmplx(g_rsite(2*num))
      enddo
    enddo
    trace=(0d0,0d0)
    do i=1,NMAT
      trace=trace+G_PHIMAT(i,i,s)
    enddo
    do i=1,NMAT
      G_PHIMAT(i,i,s)=G_PHIMAT(i,i,s)-trace/dcmplx(dble(NMAT))
    enddo
  enddo

  ! random number must be sufficiently small
  G_AMAT=(0d0,0d0)
  do l=1,global_num_links
    do a=1,dimG
      G_AMAT(:,:,l)=G_AMAT(:,:,l)+g_rlink(a,l)*TMAT(:,:,a)
    enddo
  enddo
#ifdef PARALLEL
endif
#endif

#ifdef PARALLEL
do s=1,global_num_sites
  rank=local_site_of_global(s)%rank_ 
  ls=local_site_of_global(s)%label_
  if( MYRANK == 0 ) then
    if( rank == 0 ) then 
      P_PhiMat(:,:,ls) = G_PhiMat(:,:,s)
    else
      call MPI_SEND(G_PhiMat(:,:,s),NMAT*NMAT,MPI_DOUBLE_COMPLEX,rank,s,MPI_COMM_WORLD,IERR)
    endif
  else
    if( MYRANK == rank ) then 
      call MPI_RECV(P_PhiMat(:,:,ls),NMAT*NMAT,MPI_DOUBLE_COMPLEX,0,s,MPI_COMM_WORLD,ISTATUS,IERR)
    endif
  endif
enddo


do l=1,global_num_links
  !if( MYRANK == 0 ) then
    !tmpmat = G_AMat(:,:,l)
  !endif
  rank=local_link_of_global(l)%rank_ 
  ll=local_link_of_global(l)%label_
  if( MYRANK == 0 ) then
    if( rank == 0 ) then 
      P_AMat(:,:,ll) = G_Amat(:,:,l)
    else
      call MPI_SEND(G_AMat(:,:,l),NMAT*NMAT,MPI_DOUBLE_COMPLEX,rank,global_num_sites+l,MPI_COMM_WORLD,IERR)
    endif
  else
    if( MYRANK == rank  ) then 
      call MPI_RECV(P_AMat(:,:,ll),NMAT*NMAT,MPI_DOUBLE_COMPLEX,0,global_num_sites+l,MPI_COMM_WORLD,ISTATUS,IERR)
    endif
  endif
enddo

#else
P_PhiMat=G_PhiMat
P_AMat=G_AMat
#endif
      
!rtmp=site_abs(P_PhiMat)
!if( MYRANK==0 ) write(*,*) "|P_PhiMat|^2=",rtmp
!rtmp=link_abs(P_AMat)
!if( MYRANK==0 ) write(*,*) "|P_AMat|^2=",rtmp
!
!call stop_for_test
!integer rank
!
!
!if( MYRANK == 0 ) then
!  call genrand_init( put=10000 )
!  allocate( gauss(1:dimG*global_num_sites) )
!  call BoxMuller(gauss,(dimG)*global_num_sites)
!endif
!
!do s=1,global_num_sites
!  rank = local_site_of_global(s)%rank_
!
!  if( MYRANK == 0) then 
!    do a=1,dimG
!      i=dimG*(s-1)+a
!      P_Phi(a)=dcmplx(dsqrt(0.5d0)*gauss(2*i-1)) &
!          + im_unit*dcmplx(dsqrt(0.5d0)*gauss(2*i))
!    enddo
!    if( rank == 0 ) then
!      call make_traceless_matrix_from_modes(P_PhiMat(:,:,local_site_of_global(s)%label_),NMAT,P_Phi)
!    else
!      call MPI_SEND(P_Phi,dimG,MPI_DOUBLE_COMPLEX,rank,s,MPI_COMM_WORLD,IERR)
!    endif
!  elseif( rank == MYRANK ) then
!    call MPI_RECV(P_Phi,dimG,MPI_DOUBLE_COMPLEX,0,s,MPI_COMM_WORLD,ISTATUS,IERR)
!    call make_traceless_matrix_from_modes(P_PhiMat(:,:,local_site_of_global(s)%label_),NMAT,P_Phi)
!  endif
!enddo
!
!if( MYRANK == 0 ) then
!  call BoxMuller( gauss,(dimG*global_num_links+1)/2 )
!endif
!
!do l=1,global_num_links
!  rank = local_link_of_global(l)%rank_
!
!  if( MYRANK == 0 ) then
!    do a=1,dimG
!      i=dimG*(l-1)+a
!      P_A(a)=gauss(i)
!    enddo
!    if( rank == 0 ) then
!      call make_traceless_matrix_from_modes(P_AMat(:,:,l),NMAT,dcmplx( P_A(:) ))
!    else
!      call MPI_SEND(P_A,dimG,MPI_DOUBLE_PRECISION,rank,global_num_sites+l,MPI_COMM_WORLD,IERR)
!    endif
!  elseif( rank == MYRANK ) then
!    call MPI_RECV(P_A,dimG,MPI_DOUBLE_PRECISION,0,global_num_sites+l,MPI_COMM_WORLD,ISTATUS,IERR)
!    call make_traceless_matrix_from_modes(P_AMat(:,:,local_link_of_global(l)%label_),NMAT,dcmplx( P_A(:) ))
!  endif
!enddo
!
!
!#else
!!! set random momentum
!allocate(g_rsite(1:2*NMAT*NMAT*global_num_sites))
!allocate(g_rlink(1:dimG*global_num_links))
!call BoxMuller(gauss,(dimG)*num_sites)
!do s=1,num_sites
!  do a=1,dimG
!    i=dimG*(s-1)+a
!    P_Phi(a)=dcmplx(dsqrt(0.5d0)*gauss(2*i-1)) &
!        + im_unit*dcmplx(dsqrt(0.5d0)*gauss(2*i))
!  enddo
!  call make_traceless_matrix_from_modes(P_PhiMat(:,:,s),NMAT,P_Phi)
!enddo
!
!call BoxMuller( gauss,(dimG*num_links+1)/2 )
!do l=1,num_links
!  do a=1,dimG
!    i=dimG*(l-1)+a
!    P_A(a)=gauss(i)
!  enddo
!  call make_traceless_matrix_from_modes(P_AMat(:,:,l),NMAT,dcmplx( P_A(:) ))
!enddo
!#endif


end subroutine set_randomP_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Molecular Evolution by multistep Leap Frog ingegrator
subroutine molecular_evolution_multistep(UMAT,PhiMat,PF_eta,PF_lambda,PF_chi,P_AMat,P_PhiMat,info)
#ifdef PARALLEL
use parallel   
#endif         
implicit none  
               
complex(kind(0d0)), intent(inout) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(inout) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: PF_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: PF_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PF_chi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), intent(inout) :: P_PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(inout) :: P_AMat(1:NMAT,1:NMAT,1:num_links)
integer, intent(inout) :: info
integer :: i,j,step,local_info
double precision :: rtmp


info=0
!!!!!!!!!!!!!!!!!!
!!!! Update Phi
!! first step
!write(*,*) "test"
call update_momentum_boson(P_PhiMat,P_AMat,PhiMat,UMAT,Dtau_boson*0.5d0)
if(pf==0) call update_momentum_fermion(P_PhiMat,P_AMat,PhiMat,UMAT,PF_eta,PF_lambda,PF_chi,local_info,Dtau_fermion*0.5d0)
if(local_info==1) then
  info=1
  return
endif

!! main step
step=0
do i=1,Nfermion
  do j=1,Nboson
    call update_PhiMat(PhiMat,P_phiMat,Dtau_boson)
    call update_Umat(UMat,P_AMat,Dtau_boson)
    step=step+1
    if( step .ne. Nfermion*Nboson ) then
      call update_momentum_boson(P_PhiMat,P_AMat,PhiMat,UMAT,Dtau_boson)
    endif
  enddo
  if( step .ne. Nfermion*Nboson ) then
    if(pf==0) call update_momentum_fermion(P_PhiMat,P_AMat,PhiMat,UMAT,PF_eta,PF_lambda,PF_chi,local_info,Dtau_fermion)
    if(local_info==1) then
      info=1
      return
    endif
  endif
enddo

!! final step
call update_momentum_boson(P_PhiMat,P_AMat,PhiMat,UMAT,Dtau_boson*0.5d0)
if(pf==0) call update_momentum_fermion(P_PhiMat,P_AMat,PhiMat,UMAT,PF_eta,PF_lambda,PF_chi,local_info,Dtau_fermion*0.5d0)
if(local_info==1) then
  info=1
  return
endif

end subroutine molecular_evolution_multistep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Molecular Evolution by Omelyan integrator
!subroutine molecular_evolution_Omelyan(UMAT,PhiMat,PF_eta,PF_lambda,PF_chi,P_AMat,P_PhiMat,info)
!implicit none
!
!double precision, parameter :: lambda=0.1931833275d0
!complex(kind(0d0)), intent(inout) :: UMAT(1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)), intent(inout) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
!complex(kind(0d0)), intent(in) :: PF_eta(1:NMAT,1:NMAT,1:num_sites)
!complex(kind(0d0)), intent(in) :: PF_lambda(1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)), intent(in) :: PF_chi(1:NMAT,1:NMAT,1:num_faces)
!complex(kind(0d0)), intent(inout) :: P_PhiMat(1:NMAT,1:NMAT,1:num_sites)
!complex(kind(0d0)), intent(inout) :: P_AMat(1:NMAT,1:NMAT,1:num_links)
!integer, intent(inout) :: info
!
!
!integer :: f_skip ! fermion forceを計算するかどうか 0/1
!integer :: s,l,f,a
!integer :: i
!integer :: j,k,ii
!
!complex(kind(0d0)) :: minimal, maximal,tmp
!
!f_skip=0
!do i=1, Ntau
!! momentum
!  call update_momentum(P_PhiMat,P_AMat,PhiMat,UMAT,PF_eta,PF_lambda,PF_chi,info,Dtau_phi*lambda,Dtau_A*lambda,f_skip)
!  if ( info == 1 ) return
!! variables
!  call update_PhiMat(PhiMat,P_phiMat,Dtau_phi*0.5d0)
!  call update_UMAT(UMAT,P_AMat,0.5d0*Dtau_A)
!! momentum
!  call update_momentum(P_PhiMat,P_AMat,PhiMat,UMAT,PF_eta,PF_lambda,PF_chi,info,Dtau_phi*(1d0-2d0*lambda),Dtau_A*(1d0-2d0*lambda),f_skip)
!  if( info == 1 ) return
!! variables
!  call update_PhiMat(PhiMat,P_phimat,Dtau_phi*0.5d0)
!  call update_UMAT(UMAT,P_AMat,0.5d0*Dtau_A)
!! momentum
!  call update_momentum(P_PhiMat,P_AMat,PhiMat,UMAT,PF_eta,PF_lambda,PF_chi,info,Dtau_phi*lambda,Dtau_A*lambda,f_skip)
!  if ( info == 1 ) return
!enddo
!
!end subroutine molecular_evolution_Omelyan


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! update UMAT
subroutine update_UMAT(UMAT,P_AMat,Dtau)
use matrix_functions, only : MATRIX_EXP,make_matrix_traceless
use SUN_generators, only : make_traceless_matrix_from_modes
#ifdef PARALLEL
use parallel
#endif
implicit none

complex(kind(0d0)),intent(inout) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
!double precision,intent(in) :: P_A(1:dimG,1:num_links)
complex(kind(0d0)), intent(in) :: P_AMat(1:NMAT,1:NMAT,1:num_links)
double precision, intent(in) :: Dtau
complex(kind(0d0)) :: Plmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: dU(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
integer :: l

!do l=1,num_links
  !call Make_traceless_matrix_from_modes(P_Amat(:,:,l), NMAT, dcmplx(P_A(:,l)))
!enddo
do l=1,num_links
  Plmat=(0d0,1d0)*Dtau*P_Amat(:,:,l)
  !call MATRIX_EXP(NMAT,Plmat,dU)
  !call traceless_projection(Plmat)
  call make_matrix_traceless(Plmat)
  call MATRIX_EXP(dU,Plmat)
  tmpmat=UMAT(:,:,l)
  call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
    dU, NMAT, &
    tmpmat, NMAT, &
    (0d0,0d0), UMAT(:,:,l), NMAT)
enddo

call syncronize_links(UMat)

end subroutine update_UMAT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! update PhiMat
subroutine update_PhiMat(PhiMat,P_phimat,Dtau)
use SUN_generators, only : make_traceless_matrix_from_modes
#ifdef PARALLEL
use parallel
#endif
implicit none

complex(kind(0d0)),intent(inout) :: PhiMAT(1:NMAT,1:NMAT,1:num_necessary_sites)
!complex(kind(0d0)),intent(inout) :: Phi(1:dimG,1:num_links)
complex(kind(0d0)),intent(in) :: P_PhiMat(1:NMAT,1:NMAT,1:num_sites)
double precision, intent(in) :: Dtau
integer :: s,i,j

do s=1,num_sites
  !do a=1,dimG
  do j=1,NMAT
    do i=1,NMAT
      PhiMat(i,j,s)=PhiMat(i,j,s)+Dtau*dconjg( P_PhiMat(j,i,s) )
      !Phi(a,s)=Phi(a,s) + Dtau * dconjg( P_phi(a,s) )
    enddo
  enddo
enddo
!do s=1,num_sites
  !call make_traceless_matrix_from_modes(PhiMat(:,:,s),NMAT,Phi(:,s))
!enddo

call syncronize_sites(PhiMat)

end subroutine update_PhiMat



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! update momentum
!subroutine update_momentum(P_PhiMat,P_AMat,PhiMat,UMAT,PF_eta,PF_lambda,PF_chi,info,deltaPPhi,deltaA,f_skip)
!use SUN_generators, only : make_traceless_matrix_from_modes, trace_mta
!implicit none
!
!double precision :: P_A(1:dimG,1:num_links)
!complex(kind(0d0)), intent(inout) :: P_AMat(1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)), intent(inout) :: P_PhiMat(1:NMAT,1:NMAT,1:num_sites)
!complex(kind(0d0)), intent(inout) :: UMAT(1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)), intent(inout) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)
!!complex(kind(0d0)) :: PF(1:sizeD)
!complex(kind(0d0)), intent(in) :: PF_eta(1:NMAT,1:NMAT,1:num_sites)
!complex(kind(0d0)), intent(in) :: PF_lambda(1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)), intent(in) :: PF_chi(1:NMAT,1:NMAT,1:num_faces)
!complex(kind(0d0)) :: tmp
!double precision, intent(in) :: deltaPPhi,deltaA
!integer, intent(out) :: info
!integer, intent(in) :: f_skip ! fermion forceを計算するかどうか 0/1
!
!
!complex(kind(0d0)) :: dSdPhi(1:NMAT,1:NMAT,1:num_sites)
!complex(kind(0d0)) :: dSdA(1:NMAT,1:NMAT,1:num_links)
!
!call Make_force(dSdPhi,dSdA,UMAT,PhiMat,PF_eta,PF_lambda,PF_chi,info,f_skip)
!
!P_PhiMat(:,:,:)=P_PhiMat(:,:,:) - dSdPhi(:,:,:) * deltaPPhi
!P_AMat(:,:,:) = P_AMat(:,:,:) - dSdA(:,:,:) * dcmplx(deltaA)
!
!end subroutine update_momentum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! update momentum of A by boson
subroutine update_momentum_boson(P_PhiMat,P_AMat,PhiMat,UMAT,delta_b)
#ifdef PARALLEL
use parallel
#endif
implicit none

complex(kind(0d0)), intent(inout) :: P_PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(inout) :: P_AMat(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(inout) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(inout) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
double precision, intent(in) :: delta_b
complex(kind(0d0)) :: tmp

complex(kind(0d0)) :: dSdA(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: dSdPhi(1:NMAT,1:NMAT,1:num_sites)
integer :: s,l

call Make_bosonic_force(dSdPhi,dSdA,UMAT,PhiMat)

P_PhiMat(:,:,:)=P_PhiMat(:,:,:) - dSdPhi(:,:,:) * delta_b
P_AMat(:,:,:) = P_AMat(:,:,:) - dSdA(:,:,:) * dcmplx(delta_b)

end subroutine update_momentum_boson

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! update momentum by fermion
subroutine update_momentum_fermion(P_Phimat,P_AMat,PhiMat,UMAT,PF_eta,PF_lambda,PF_chi,info,delta_f)
#ifdef PARALLEL
use parallel
#endif
implicit none

complex(kind(0d0)), intent(inout) :: P_PhiMat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(inout) :: P_AMat(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(inout) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(inout) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
!complex(kind(0d0)) :: PF(1:sizeD)
complex(kind(0d0)), intent(in) :: PF_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: PF_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PF_chi(1:NMAT,1:NMAT,1:num_necessary_faces)
double precision, intent(in) :: delta_f
integer, intent(out) :: info
complex(kind(0d0)) :: tmp


complex(kind(0d0)) :: dSdA(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: dSdPhi(1:NMAT,1:NMAT,1:num_sites)

info=0
call Make_fermionic_force(dSdPhi,dSdA,UMAT,PhiMat,PF_eta,PF_lambda,PF_chi,info)
if( info == 1) return

P_PhiMat(:,:,:)=P_PhiMat(:,:,:) - dSdPhi(:,:,:) * delta_f
P_AMat(:,:,:) = P_AMat(:,:,:) - dSdA(:,:,:) * dcmplx(delta_f)

end subroutine update_momentum_fermion



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! make pseudo fermion
subroutine make_pseudo_fermion(PF_eta,PF_lambda,PF_chi,UMAT,PhiMat)
!use matrix_functions,  only : BoxMuller2
use SUN_generators, only : trace_MTa
use Dirac_operator
use rational_algorithm
#ifdef PARALLEL
use parallel
!use global_subroutines, only : syncronize_sites, syncronize_links, syncronize_faces
#endif
!use global_subroutines, only : mat_to_vec, vec_to_mat
implicit none


complex(kind(0d0)), intent(out) :: PF_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(out) :: PF_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(out) :: PF_chi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
!complex(kind(0d0)) :: Phi(1:dimG,1:num_sites)

!complex(kind(0d0)) :: PF(1:sizeD)
complex(kind(0d0)) :: G_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: G_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: G_chi(1:NMAT,1:NMAT,1:num_necessary_faces)
!#ifdef PARALLEL
!complex(kind(0d0)) :: tmp_G_eta(1:NMAT,1:NMAT,1:num_sites)
!complex(kind(0d0)) :: tmp_G_lambda(1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)) :: tmp_G_chi(1:NMAT,1:NMAT,1:num_faces)
!#endif
complex(kind(0d0)), allocatable :: gauss(:)
double precision, allocatable :: gauss2(:)
double precision :: rtmp
!complex(kind(0d0)) :: gauss(1:sizeD)
!double precision :: gauss2(1:2*sizeD),rtmp
integer :: s,l,f,k,i,j,info,CGite,trace
complex(kind(0d0)) :: tmp, tmp2

if( test_mode == 0 ) then
  allocate( gauss(1:(NMAT*NMAT-1)*(num_sites+num_links+num_faces) ) )
  allocate( gauss2(1:2*(NMAT*NMAT-1)*(num_sites+num_links+num_faces) ) )
  
  call BoxMuller2(gauss2, (NMAT*NMAT-1)*(num_sites+num_links+num_faces) )
  
  do i=1,(NMAT*NMAT-1)*(num_sites+num_links+num_faces)
    gauss(i)=(gauss2(2*i-1) + gauss2(2*i)*(0d0,1d0)) * dcmplx(dsqrt(0.5d0))
  enddo
  
  call vec_to_mat(G_eta(:,:,1:num_sites),G_lambda(:,:,1:num_links),G_chi(:,:,1:num_faces),gauss)

#ifdef PARALLEL
  call syncronize_sites(G_eta)
  call syncronize_links(G_lambda)
  call syncronize_faces(G_chi)
#endif
else
  !!!!!!!!!!!!!!!
  !!for test
  allocate( gauss(1:(NMAT*NMAT-1)*(global_num_sites+global_num_links+global_num_faces) ) )
  allocate( gauss2(1:2*(NMAT*NMAT-1)*(global_num_sites+global_num_links+global_num_faces) ) )
  if( MYRANK == 0 ) then 
    call genrand_init( put=10000 )
    call BoxMuller2(gauss2, (NMAT*NMAT-1)*(global_num_sites+global_num_links+global_num_faces) )
    do i=1,(NMAT*NMAT-1)*(global_num_sites+global_num_links+global_num_faces)
      gauss(i)=(gauss2(2*i-1) + gauss2(2*i)*(0d0,1d0)) * dcmplx(dsqrt(0.5d0))
    enddo
  endif

  call MPI_BCAST(gauss,(NMAT*NMAT-1)*(global_num_sites+global_num_links+global_num_faces),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)

  call globalvec_to_localmat(G_eta,G_lambda,G_chi,gauss)
  !!!!!!!!!!!!!!!!!!!!
  
endif

!!!!!!!!!!!!!!!
!! for test
!if( test_mode /= 0 ) then 
!  rtmp=site_abs(PhiMat(:,:,1:num_sites))
!  if( MYRANK == 0 ) then
!    write(*,'(a,E25.18)') "# |PhiMat|^2=",rtmp
!  endif
!  rtmp=link_abs(UMat(:,:,1:num_links))
!  if( MYRANK == 0 ) then
!    write(*,'(a,E25.18)') "# |UMat|^2=",rtmp
!  endif
!  rtmp=site_abs(G_eta(:,:,1:num_sites))
!  if( MYRANK == 0 ) then
!    write(*,'(a,E25.18)') "# |G_eta|^2=",rtmp
!  endif
!  rtmp=link_abs(G_lambda(:,:,1:num_links))
!  if( MYRANK == 0 ) then
!    write(*,'(a,E25.18)') "# |G_lambda|^2=",rtmp
!  endif
!  rtmp=face_abs(G_chi(:,:,1:num_faces))
!  if( MYRANK == 0 ) then
!    write(*,'(a,E25.18)') "# |G_chi|^2=",rtmp
!  endif
!endif
!
!call stop_for_test
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!do i=1,sizeD
!write(*,*) MYRANK, gauss(i)
!enddo
info=0
call calc_matrix_rational_power(&
  PF_eta(:,:,1:num_sites), PF_lambda(:,:,1:num_links), PF_chi(:,:,1:num_faces), &
  G_eta(:,:,1:num_necessary_sites), G_lambda(:,:,1:num_necessary_links), G_chi(:,:,1:num_necessary_faces), &
  epsilon, CG_max, info, CGite, &
  Remez_alpha8, Remez_beta8, UMAT,PhiMat, prod_DdagD)
!write(*,*) MYRANK,info
#ifdef PARALLEL
call syncronize_sites(PF_eta)
call syncronize_links(PF_lambda)
call syncronize_faces(PF_chi)
#endif


!if( test_mode /= 0 ) then 
!  rtmp=site_abs(PF_eta(:,:,1:num_sites))
!  if( MYRANK == 0 ) then
!    write(*,'(a,E25.18)') "# |PF_eta|^2=",rtmp
!  endif
!  rtmp=link_abs(PF_lambda(:,:,1:num_links))
!  if( MYRANK == 0 ) then
!    write(*,'(a,E25.18)') "# |PF_lambda|^2=",rtmp
!  endif
!  rtmp=face_abs(PF_chi(:,:,1:num_faces))
!  if( MYRANK == 0 ) then
!    write(*,'(a,E25.18)') "# |PF_chi|^2=",rtmp
!  endif
!endif

end subroutine make_pseudo_fermion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Metropolice_test(delta_Ham,PhiMat_Bak,UMAT_BAK,PhiMat,UMAT,accept)
#ifdef PARALLEL
use parallel
#endif
implicit none

double precision, intent(in) :: delta_Ham
!complex(kind(0d0)), intent(in) :: Phi_BAK(1:dimG,1:num_sites)
complex(kind(0d0)), intent(in) :: PhiMat_BAK(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: UMAT_BAK(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)), intent(inout) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)), intent(inout) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(inout) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
integer, intent(inout) :: accept
integer :: local_info

double precision :: random_num

#ifdef PARALLEL
if (MYRANK == 0) then 
#endif 

local_info=1
if( delta_Ham <= 0d0 .or. new_config >= 2 ) then 
  local_info=0
  accept=accept+1
else
  call genrand_real3(random_num)
  if( exp( -delta_Ham ) > random_num ) then
    local_info=0
    accept=accept+1
  endif
endif
#ifdef PARALLEL
endif
call MPI_BCAST(local_info,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
#endif 
if( local_info == 1 ) then 
  PhiMat=PhiMat_BAK
  UMAT=UMAT_BAK
endif
return
end subroutine Metropolice_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine check_vacuum(UMAT)
use matrix_functions, only : matrix_norm, make_unit_matrix
!use global_subroutines
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
integer l,i,j,f
double precision norm
complex(kind(0d0)) UNITM(1:NMAT,1:NMAT)
complex(kind(0d0)) tmp(1:NMAT,1:NMAT)
complex(kind(0d0)) Uf(1:NMAT,1:NMAT)
double precision min_dist
double precision dist(1:NMAT-1)

!write(*,*) "===== check distance from 1 ==========="
if( NMAT <= 4 ) then 
  min_dist=2d0*dsqrt(2d0)
else
  min_dist=dsin(PI/dble(NMAT))*2d0*sqrt(2d0)
endif


call make_unit_matrix(UNITM)
do f=1,num_faces
  call Make_face_variable(Uf,f,UMAT)
  tmp=UNITM-Uf
  !call matrix_norm(tmp,NMAT,norm)
  call matrix_norm(norm,tmp)
  write(*,*) f,min_dist-norm
enddo
end subroutine check_vacuum

