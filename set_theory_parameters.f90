!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11!!!!
!! SUBROUTINES for initialization
subroutine set_theory_parameters(seed)
implicit none

integer, intent(inout) :: seed
! Open parameter file
#ifdef PARALLEL
if (MYRANK==0) then
#endif
open(PAR_FILE, file=PAR_FILE_NAME, status='old',action='READ')
  !! NMAT
  !  data of the simplicial complex
  !read(PAR_FILE,*) FsubSC
  !read(PAR_FILE,*) ALPHA_BETA
  !read(PAR_FILE,*) read_alpha
!! NMAT
  read(PAR_FILE,'()') 
  read(PAR_FILE,*) NMAT
!! mass_phi2 : physical mass square of \Phi
  read(PAR_FILE,'()') 
  read(PAR_FILE,*) phys_mass_square_phi
!! pf, p1, p2, p3, p4, p5, p_mass : set 1 for omitting each fermion action(pf means all femrions)
  read(PAR_FILE,'()') 
  read(PAR_FILE,*) pf, p1, p2, p3, p4, p5, p_mass
!! mass_f; fermion mass
  read(PAR_FILE,'()') 
  read(PAR_FILE,*) mass_f
!!  data of the simplicial complex
  read(PAR_FILE,'()') 
  read(PAR_FILE,*) SC_FILE_NAME
!! test_mode ; 0:Simulation mode 1:Test mode, 
  read(PAR_FILE,'()') 
  read(PAR_FILE,*) test_mode
!! branch_mode !! 0:normal, 1:make branch 
  read(PAR_FILE,'()') 
  read(PAR_FILE,*) branch_mode
!! branch_root ! make branch from this config. default:0
  read(PAR_FILE,'()') 
  read(PAR_FILE,*) branch_root
!! new_branch_label !! label of branches to make
  read(PAR_FILE,'()') 
  read(PAR_FILE,*) new_branch_label
!! force_measurement ; 1:measure forces  
  read(PAR_FILE,'()') 
  read(PAR_FILE,*) force_measurement
!! new_config ; 0:Old Config 1:New Config
  read(PAR_FILE,'()') 
  read(PAR_FILE,*) new_config
!!  cold_hot : 0: cold start, 1: hot start
  read(PAR_FILE,'()') 
  read(PAR_FILE,*) cold_hot
!! omit_metroplis : 0: with metropolice test, 1: all accept
  read(PAR_FILE,'()') 
  read(PAR_FILE,*) omit_metropolis
!! fix_seed ; 0:previous value 1:fix 2:system time 
  read(PAR_FILE,'()') 
  read(PAR_FILE,*) fix_seed
!! reset_ite ; 0:continue total ite, 1:reset total ite
  read(PAR_FILE,'()') 
  read(PAR_FILE,*) reset_ite
!! seed (works when fix_seed=1)
  read(PAR_FILE,'()') 
  read(PAR_FILE,*) seed
!! m_omega ;integer to avoid vacuum degeneracy (mm>=NMAT/4)
  read(PAR_FILE,'()') 
  read(PAR_FILE,*) m_omega
!!  Remez data for x^{-1/4}
  read(PAR_FILE,'()') 
  read(PAR_FILE,*) Remez_1ovminus4
!!  Remez data for x^{1/8}
  read(PAR_FILE,'()') 
  read(PAR_FILE,*) Remez_1ov8
!! Remez_factor4 ; range of remez approx( min*factor .. max*factor) 
  read(PAR_FILE,'()') 
  read(PAR_FILE,*) Remez_factor4
!! Remez_factor8 ; range of remez approx( min*factor .. max*factor) 
  read(PAR_FILE,'()') 
  read(PAR_FILE,*) Remez_factor8
!! epsilon ; range to stop CG solver
  read(PAR_FILE,'()') 
  read(PAR_FILE,*) epsilon
!! CG_max ; maximum number of CG iteration
  read(PAR_FILE,'()') 
  read(PAR_FILE,*) CG_max
!!! Nfermion 
!  read(PAR_FILE,'()') 
!  read(PAR_FILE,*) Nfermion
!!! FB_ratio ; force計算のfermion/boson比
!  read(PAR_FILE,'()') 
!  read(PAR_FILE,*) FB_ratio
close(PAR_FILE)

if( p_mass /= 0 ) mass_f=0d0

#ifdef PARALLEL
endif
! send all parameters to all the other rank
!  read(PAR_FILE,*) NMAT
call MPI_BCAST(NMAT,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) pf, p1, p2, p3, p4, p5
call MPI_BCAST(pf,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(p1,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(p2,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(p3,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(p4,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(p5,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(p_mass,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) mass_f
call MPI_BCAST(mass_f,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) SC_FILE_NAME
call MPI_BCAST(SC_FILE_NAME,128,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) test_mode
call MPI_BCAST(test_mode,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) branch_mode
call MPI_BCAST(branch_mode,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) branch_root
call MPI_BCAST(branch_root,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) new_branch_label
call MPI_BCAST(new_branch_label,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) force_measurement
call MPI_BCAST(force_measurement,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) new_config
call MPI_BCAST(new_config,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) cold_hot
call MPI_BCAST(cold_hot,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) omit_metropolis
call MPI_BCAST(omit_metropolis,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) fix_seed
call MPI_BCAST(fix_seed,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) reset_ite
call MPI_BCAST(reset_ite,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) seed
call MPI_BCAST(seed,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) m_omega
call MPI_BCAST(m_omega,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) phys_mass_square_phi
call MPI_BCAST(phys_mass_square_phi,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) Remez_1ovminus4
call MPI_BCAST(Remez_1ovminus4,128,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) Remez_1ov8
call MPI_BCAST(Remez_1ov8,128,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) Remez_factor4
call MPI_BCAST(Remez_factor4,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) Remez_factor8
call MPI_BCAST(Remez_factor8,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) epsilon
call MPI_BCAST(epsilon,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !read(PAR_FILE,*) CG_max
call MPI_BCAST(CG_max,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
!  !dimG=NMAT*NMAT-1
!call MPI_BCAST(dimG,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
!  !maximal_dist = 2d0*dsqrt(2d0)
!call MPI_BCAST(maximal_dist,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!  !e_max
!call MPI_BCAST(e_max,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
#endif

!!!!!!!!!!!!!!
dimG=NMAT*NMAT-1
!!!!!!!!!!!!!!
if (NMAT<=4) then
  maximal_dist = 2d0*dsqrt(2d0)
else
  maximal_dist = 2d0*sqrt(dble(NMAT))*dsin(3.1415926535898d0/dble(NMAT))
endif
!!!!!!!!!!!!!!
if( NMAT <= 4 ) then 
  e_max = dcmplx( 2d0 * dsqrt(2d0) )
else
  e_max = dcmplx( 2d0 * dsqrt(dble(NMAT)) * dsin( 3.1415926835898d0/dble(NMAT) ) )
endif
e_max=e_max*(0.8,0d0)

end subroutine set_theory_parameters


