!! UMAT, Phi, 
program write_down_config
use SUN_generators, only : make_traceless_matrix_from_modes
use global_parameters
use initialization
use simplicial_complex
use hamiltonian
use Dirac_operator, only : make_Dirac, make_DdagD
use matrix_functions, only : PfaffianLog, matrix_eigenvalues
implicit none

integer :: ios, num_data
integer :: num, iargc

integer,parameter :: CONFIG_IN=10 ! input config from medfie (binary)
integer,parameter :: CONFIG_OUT_UMAT=11 ! UMAT
integer,parameter :: CONFIG_OUT_Phi=14 ! Phi 
integer,parameter :: EIGEN_OUT=12 ! UMATの固有値
integer,parameter :: RESULT_OUT=13 ! bosonic actionの成分など
integer,parameter :: CONFIG_OUT_UF=15
character(128) :: config_in_filename
character(128),parameter :: config_out_UMAT_filename="config_UMAT_out.dat"
character(128),parameter :: config_out_UF_filename="config_Uf_out.dat"
character(128),parameter :: config_out_Phi_filename="config_PHI_out.dat"
character(128),parameter :: eigen_out_filename="Uf_eigen.dat"
character(128),parameter :: result_out_filename="result_out.dat"

!! variables
complex(kind(0d0)), allocatable :: UMAT(:,:,:) ! unitary link variables
complex(kind(0d0)), allocatable :: Uf(:,:) ! face variable
complex(kind(0d0)), allocatable :: Uf_eigen(:) ! face variable
complex(kind(0d0)), allocatable :: PHI(:,:) ! complex scalar at sites (su(N))
complex(kind(0d0)), allocatable :: PHI_MAT(:,:,:) ! complex scalar at sites (matrix)
integer :: ite

real(8) :: SB_M, SB_S, SB_L, SB_F
complex(kind(0d0)), allocatable :: Dirac(:,:) 
complex(kind(0d0)), allocatable :: DDdag(:,:)
complex(kind(0d0)), allocatable :: DDdag2(:,:)
complex(kind(0d0)), allocatable :: Dirac_inv(:,:) 
complex(kind(0d0)), allocatable :: eigenvalues(:)
complex(kind(0d0)), allocatable :: eigenvalues2(:)
complex(kind(0d0)), allocatable :: Ufall(:,:,:) !(1:NMAT,1:NMAT), Ufm(1:NMAT,1:NMAT)


real(8) :: logPF
complex(kind(0d0)) :: phasePF, PF



integer i,j,k
integer s,t,l,l1,l2,f,f1,f2
integer a,b

num = iargc()
call getarg(1,config_in_filename)

open(unit=CONFIG_IN,file=config_in_filename,action='read',form='unformatted')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! initialization
read(CONFIG_IN) NMAT
read(CONFIG_IN) LatticeSpacing
read(CONFIG_IN) SC_FILE_NAME
read(CONFIG_IN) ALPHA_BETA
read(CONFIG_IN) test_mode
read(CONFIG_IN) new_config
read(CONFIG_IN) fix_seed
read(CONFIG_IN) read_alpha
read(CONFIG_IN) save_med_step
read(CONFIG_IN) obs_step
read(CONFIG_IN) m_omega
read(CONFIG_IN) mass_square_phi
read(CONFIG_IN) mass_f
read(CONFIG_IN) Remez_factor4
read(CONFIG_IN) Remez_factor8
read(CONFIG_IN) epsilon
read(CONFIG_IN) CG_max
read(CONFIG_IN) num_ite
read(CONFIG_IN) Ntau
read(CONFIG_IN) Dtau
read(CONFIG_IN) R_phi
read(CONFIG_IN) R_A
read(CONFIG_IN) Fconfigin
read(CONFIG_IN) Foutput
read(CONFIG_IN) Fconfigout
read(CONFIG_IN) Fmedconf
dimG=NMAT*NMAT-1
Dtau_phi = R_phi * Dtau
Dtau_A = R_A * Dtau
one_ov_2g2N=1d0/(2d0*LatticeSpacing*LatticeSpacing)
overall_factor=dble(NMAT)/(2d0*LatticeSpacing*LatticeSpacing)

call set_NZF
call set_sc
call set_alpha_beta
allocate( UMAT(1:NMAT,1:NMAT,1:num_links) )
allocate( Uf(1:NMAT,1:NMAT) )
allocate( Uf_eigen(1:NMAT) )
allocate( PHI(1:dimG, 1:num_sites) )
allocate( PHI_MAT(1:NMAT,1:NMAT, 1:num_sites) )
allocate( Dirac(1:sizeD,1:sizeD) )
allocate( Dirac_inv(1:sizeD,1:sizeD) )
allocate( DDdag(1:sizeD,1:sizeD) )
allocate( DDdag2(1:sizeD,1:sizeD) )
allocate( eigenvalues(1:sizeD) )
allocate( eigenvalues2(1:sizeD) )

allocate( Ufall(1:NMAT,1:NMAT,1:num_faces) )


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 最初のconfigurationを読む。
read(CONFIG_IN,iostat=ios) ite
read(CONFIG_IN,iostat=ios) UMAT
read(CONFIG_IN,iostat=ios) PHI


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UMATとPhiを書き出し

do s=1,num_sites
  call make_traceless_matrix_from_modes(phi_mat(:,:,s),NMAT,Phi(:,s))
enddo

open(unit=CONFIG_OUT_PHi,status='replace',file=config_out_Phi_filename,action='write')
do s=1,num_sites
  do i=1,NMAT
    do j=1,NMAT
      write(CONFIG_OUT_Phi,'(I3,I3,I3,2(E15.5))') s,i,j, &
        dble(Phi_MAT(i,j,s)),dble((0d0,-1d0)*Phi_MAT(i,j,s)) 
    enddo
  enddo
enddo
close(CONFIG_OUT_PHI)

open(unit=CONFIG_OUT_UMAT,status='replace',file=config_out_UMAT_filename,action='write')
do l=1,num_links
  do i=1,NMAT
    do j=1,NMAT
      write(CONFIG_OUT_UMAT,'(3(I3),2(E15.5))') l,i,j, &
        dble(UMAT(i,j,l)), dble( (0d0,-1d0)*UMAT(i,j,l) ) 
    enddo
  enddo
enddo
close(CONFIG_OUT_UMAT)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Ufを書き出し
do f=1,num_faces
  call Make_face_variable(Ufall(:,:,f),f,UMAT)
enddo

open(unit=CONFIG_OUT_UF,status='replace',file=config_out_UF_filename,action='write')
do f=1,num_faces
  do i=1,NMAT
    do j=1,NMAT
      write(CONFIG_OUT_UF,'(3(I3),2(E15.5))') f,i,j,Ufall(i,j,f)
    enddo
  enddo
enddo
close(CONFIG_OUT_UF)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! bosonic actionとDiracのPfaffianを書き出し

call bosonic_action_mass(SB_M,Phi)
call bosonic_action_site(SB_S,Phi)
call bosonic_action_link(SB_L,UMAT,Phi)
call bosonic_action_face(SB_F,UMAT)



call make_Dirac(Dirac,UMAT,Phi)
call make_DdagD(DDdag2,UMAT,Phi)
DDdag=(0d0,0d0)
do i=1,sizeD
  do j=1,sizeD
    do k=1,sizeD
      DDdag(i,j)=DDdag(i,j)+conjg(Dirac(k,i))*(Dirac(k,j))
    enddo
  enddo
enddo
!write(*,*) DDdag-DDdag2
call Matrix_Eigenvalues(eigenvalues2,DDdag)
call Matrix_Eigenvalues(eigenvalues,Dirac)
call make_Dirac(Dirac,UMAT,Phi)
call PfaffianLog(DIrac,logPF,phasePF)
Pf=cmplx(exp(LogPF))*PhasePF

open(unit=RESULT_OUT,status='replace',file=result_out_filename,action='write')
write(RESULT_OUT,'(a,E15.5)') "SB_mass=", SB_M
write(RESULT_OUT,'(a,E15.5)') "SB_site=", SB_S
write(RESULT_OUT,'(a,E15.5)') "SB_link=", SB_L
write(RESULT_OUT,'(a,E15.5)') "SB_face=", SB_F
write(RESULT_OUT,'(a,2(E15.5))') "Pfaffian=", dble(Pf), dble((0d0,-1d0)*Pf)
write(RESULT_OUT,*) "# eigenvalues of Dirac"
do i=1,sizeD
  write(RESULT_OUT,'(3(E15.5))') dble(eigenvalues(i)),&
    dble((0d0,-1d0)*eigenvalues(i)), &
    sqrt(dble(eigenvalues(i)*conjg(eigenvalues(i))))
enddo
write(RESULT_OUT,*) "# eigenvalues of D.D^\dagger"
do i=1,sizeD
  write(RESULT_OUT,'(3(E15.5))') dble(eigenvalues2(i)),&
    dble((0d0,-1d0)*eigenvalues2(i)), &
    sqrt(dble(eigenvalues2(i)*conjg(eigenvalues2(i))))
enddo
write(RESULT_OUT,*) "# components of Dirac:site-site"
do s=1,num_sites
  do a=1,dimG
    do t=1,num_sites
      do b=1,dimG
        write(RESULT_OUT,'(a,2I3,a,a,2I3,a,2(E15.5))') "(",s,a,")","(",t,b,")",&
          dble(Dirac(dimG*(s-1)+a,dimG*(t-1)+b)),&
          dble((0d0,-1d0)*Dirac(dimG*(s-1)+a,dimG*(t-1)+b))
      enddo
    enddo
  enddo
enddo
write(RESULT_OUT,*) "# components of Dirac:face-face"
do f1=1,num_faces
  do a=1,dimG
    do f2=1,num_faces
      do b=1,dimG
        write(RESULT_OUT,'(a,2I3,a,a,2I3,a,2(E15.5))') "(",f1,a,")","(",f2,b,")",&
          dble(Dirac(dimG*(num_sites+num_links+f1-1)+a,&
                         dimG*(num_sites+num_links+f2-1)+b)),&
          dble((0d0,-1d0)*&
                   Dirac(dimG*(num_sites+num_links+f1-1)+a,&
                         dimG*(num_sites+num_links+f2-1)+b))
      enddo
    enddo
  enddo
enddo
write(RESULT_OUT,*) "# components of Dirac:link-link"
do l1=1,num_links
  do a=1,dimG
    do l2=1,num_links
      do b=1,dimG
        write(RESULT_OUT,'(a,2I3,a,a,2I3,a,2(E15.5))') "(",l1,a,")","(",l2,b,")",&
          dble(Dirac(dimG*(num_sites+l1-1)+a,&
                         dimG*(num_sites+l2-1)+b)),&
          dble((0d0,-1d0)*&
                   Dirac(dimG*(num_sites+l1-1)+a,&
                         dimG*(num_sites+l2-1)+b))
      enddo
    enddo
  enddo
enddo
write(RESULT_OUT,*) "# components of Dirac:site-link"
do s=1,num_sites
  do a=1,dimG
    do l2=1,num_links
      do b=1,dimG
        write(RESULT_OUT,'(a,2I3,a,a,2I3,a,2(E15.5))') "(",s,a,")","(",l2,b,")",&
          dble(Dirac(dimG*(s-1)+a,&
                           dimG*(num_sites+l2-1)+b)),&
          dble((0d0,-1d0)*Dirac(dimG*(s-1)+a,&
                           dimG*(num_sites+l2-1)+b))
      enddo
    enddo
  enddo
enddo

write(RESULT_OUT,*) "# components of Dirac:link-face"
do l1=1,num_links
  do a=1,dimG
    do f2=1,num_faces
      do b=1,dimG
        write(RESULT_OUT,'(a,2I3,a,a,2I3,a,2(E15.5))') "(",l1,a,")","(",f2,b,")",&
          dble(Dirac(dimG*(num_sites+l1-1)+a,&
                         dimG*(num_sites+num_links+f2-1)+b)),&
          dble((0d0,-1d0)*&
                   Dirac(dimG*(num_sites+l1-1)+a,&
                         dimG*(num_sites+num_links+f2-1)+b))
      enddo
    enddo
  enddo
enddo

close(Result_Out)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! U_fの固有値を書き出し
open(unit=EIGEN_OUT,status='replace',file=eigen_out_filename,action='write')
do f=1,num_faces
  call Make_face_variable(Uf(:,:),f,UMAT)
  call Matrix_eigenvalues(Uf_eigen, Uf)
  do i=1,NMAT
    write(EIGEN_OUT,'(2(E15.5))') dble(Uf_eigen(i)), dble((0d0,-1d0)*Uf_eigen(i))
  enddo
enddo
!! U_fの固有値は最後まで書き出す
do while (ios == 0)
  num_data=num_data+1
  read(CONFIG_IN,iostat=ios) ite
  read(CONFIG_IN,iostat=ios) UMAT
  read(CONFIG_IN,iostat=ios) PHI

  if( ios /= 0 ) exit
  do f=1,num_faces
    call Make_face_variable(Uf(:,:),f,UMAT)
    call Matrix_eigenvalues(Uf_eigen, Uf)
    do i=1,NMAT
      write(EIGEN_OUT,'(2(E15.5))') dble(Uf_eigen(i)), dble((0d0,-1d0)*Uf_eigen(i))
    enddo
  enddo
enddo
close(EIGEN_OUT)

close(CONFIG_IN)


end program write_down_config

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine calc_Pfaffian2(Pf,Dirac)
!implicit none
!
!complex(kind(0d0)), intent(in) :: Dirac(:,:) ! anti-symmetric matrix
!complex(kind(0d0)), intent(out) :: Pf ! pfaffian
!
!integer :: MATSIZE
!
!integer, allocatable :: IWORK(:)
!integer  LWORK, LDA, NB
!complex(kind(0d0)), allocatable :: WORK(:)
!real(8), allocatable :: RWORK(:)
!complex(kind(0d0)) Pf, phase
!real(8) LogPf, RePfPhase, ImPfPhase
!
!MATSIZE=size(Dirac,1)
!
!NB=NMAT*NMAT
!LWORK = MATSIZE*NB
!allocate ( IWORK(1:MATSIZE) )
!allocate ( RWORK(1:MATSIZE-1) )
!allocate ( WORK(LWORK) )
!LDA=MATSIZE
!call ZSKPFA('U','P',MATSIZE,Dirac,LDA,LogPf,phase,&
!             IWORK,WORK,LWORK,RWORK,INFO)
!Pf= cmplx(exp(LogPf)) * phase
!
!end subroutine calc_Pfaffian
