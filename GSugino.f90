program main
use mt95
use global_parameters
use SUN_generators, only : make_traceless_matrix_from_modes
use simulation
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! variables
complex(kind(0d0)), allocatable :: UMAT(:,:,:) ! unitary link variables
complex(kind(0d0)), allocatable :: PHIMAT(:,:,:) ! complex scalar at sites

integer :: total_ite ! total iteration
integer :: seed, time
type(genrand_state) :: state_mt95 ! state of mt95
double precision :: tmp

integer :: s

integer num_para, iargc
num_para = iargc()
if(num_para > 0) then
 call getarg(1,PAR_FILE_NAME)
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! set the parameters in the module "global_parameters"
! read input parameter from the parameter file
    call set_parameters(seed)
! set the structure constant of SU(NMAT)
    call set_NZF
! set the data of the simplicial complex
!  "sizeD" is set here
    call set_sc 
    call set_alpha_beta
! remove comment if test of the module is needed
!    call test_module_simplicial_complex(sc)
! initialize the size of the variables
allocate( UMAT(1:NMAT,1:NMAT,1:num_links) )
allocate( PHIMAT(1:NMAT,1:NMAT, 1:num_sites) )

! set the seed of the random number generators
! The default value is set in the parameter file.
! When fix_seed=2, seed is set by the system time.
  if (fix_seed==2 .or. fix_seed==0) then
    call system_clock(count=time)
    seed=time
  endif
! at this stage, seed is fixed to some value

! set the variables depending on simulation_mode and test_mode
  !if (test_mode==1 .or. new_config==1) then
  if (new_config==1) then
    call set_random_config(UMAT,PHIMAT) 
    total_ite=0
  else
    call read_config(total_ite,UMAT,PHIMat,state_mt95) 
  endif
  if( fix_seed == 0 ) then
    call genrand_init( put=state_mt95 )
  else
    call genrand_init( put=seed )
  endif
  if( reset_ite == 1 ) then
    total_ite=0
  endif
! set Remez data
    call set_Remez_data
!! check unitaryty
    !call check_unitary
!! check the distance from 1_N 
    !call check_distance
!! check mode expansion
    !call check_Ta_expansion
    !stop
!! check anti-symmetricity of Dirac and Hermiticity of D\dag D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!do s=1,num_sites
!call make_traceless_matrix_from_modes(PhiMat(:,:,s),NMAT,Phi(:,s))
!enddo

  !if (writedown_mode==1) then
    !call writedown_config_action_and_fores(UMAT,PhiMat,seed)
    !stop
  !endif

  if (test_mode==1) then
    !call check_Dirac(UMAT,Phi)
    call test_hamiltonian(UMAT,PhiMat)
  else
    call HybridMonteCarlo(UMAT,PhiMat,seed,total_ite)
  endif



end program main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! routines for Initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! read the preevious configuration
subroutine read_config(total_ite,UMAT,PhiMat,state_mt95)
use mt95
use global_parameters
implicit none

integer, intent(inout) :: total_ite
type(genrand_state), intent(inout) :: state_mt95
complex(kind(0d0)), intent(inout) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(inout) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)

type(genrand_srepr) :: char_mt95

open(IN_CONF_FILE, file=Fconfigin, status='OLD',action='READ',form='unformatted')
read(IN_CONF_FILE) total_ite
read(IN_CONF_FILE) UMAT
read(IN_CONF_FILE) PHIMat
read(IN_CONF_FILE) char_mt95
close(IN_CONF_FILE)
state_mt95=char_mt95

end subroutine read_config


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! set random configuration
!! 
!! For SU(N) case, the initial UMAT must be nearer than
!! all the Z_N centers, Omega_n=diag( exp( 2\pi i n / N ) ) from 1_N. 
!! We see
!!   min_n( || 1 - Omega_n || ) = 2 sin( \pi/N ). 
!! On the other hand, 
!!   || 1 - U || = 4/N Tr( sin^2( \theat T / 2 ) ) \sim < \theta^2 > 
!! Thus the random number must satisfy 
!!   <\theta^2> < 2\pi/N
!! 
subroutine set_random_config(UMAT,PhiMat)
use SUN_generators, only : Make_SUN_generators
use matrix_functions, only : matrix_exp
use simulation, only : BoxMuller2
use mt95
use global_parameters
implicit none

complex(kind(0d0)), intent(inout) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(inout) :: PhiMat(1:NMAT,1:NMAT,1:num_sites)

!double precision, parameter :: PI=dacos(-1d0)
complex(kind(0d0)) :: TMAT(1:NMAT,1:NMAT,1:NMAT**2-1)
double precision :: rsite(1:2*dimG*num_sites) ! for PHI
double precision :: rlink(1:dimG,1:num_links) ! for UMAT
complex(kind(0d0)) :: AMAT(1:NMAT,1:NMAT,1:num_links)
integer :: s,l,a,f,i,j,num

call make_SUN_generators(TMAT,NMAT)

!call genrand_real3(rsite)
call BoxMuller2(rsite,num_sites*dimG)
call genrand_real3(rlink)

rsite=rsite * 0.01d0 !/ mass_square_phi
do s=1,num_sites
  num=0
  do i=1,NMAT
    do j=1,NMAT
      num=num+1
      if ( i.ne.NMAT .or. j.ne.NMAT ) then
        PHIMAT(i,j,s)=dcmplx(rsite(2*num-1))+(0d0,1d0)*dcmplx(rsite(2*num))
      endif
    enddo
  enddo
  PhiMat(NMAT,NMAT,s)=(0d0,0d0)
  do i=1,NMAT-1
    PhiMat(NMAT,NMAT,s)=PhiMat(NMAT,NMAT,s)-PhiMat(i,i,s)
  enddo
enddo


! random number must be sufficiently small
if( m_omega == 0 ) then 
  rlink=rlink * ( 1d0/dble(NMAT*NMAT) )
else
  rlink=rlink * ( 1d0/dble(NMAT*NMAT*m_omega) )
endif
AMAT=(0d0,0d0)
do l=1,num_links
  do a=1,dimG
    AMAT(:,:,l)=AMAT(:,:,l)+rlink(a,l)*TMAT(:,:,a)
  enddo
enddo


do l=1,num_links
!call matrix_exp(NMAT,(0d0,1d0)*AMAT(:,:,l),UMAT(:,:,l))
call matrix_exp(UMAT(:,:,l),(0d0,1d0)*AMAT(:,:,l))
enddo

end subroutine set_random_config



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! routines for check
subroutine test_module_simplicial_complex(sc)
use simplicial_complex
implicit none
type(SmpCom) sc
integer numsites,numlinks,numfaces
integer origin,tip
integer link_label, link_dir
integer,allocatable :: link_labels(:),tips(:),origins(:),faces(:),sites(:)
integer s,l,f
integer i,j,k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call get_numsite_sc(sc,numsites)
write(*,*) "test:get_numsite_sc: numsite=",numsites
write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call get_numlink_sc(sc,numlinks)
write(*,*) "test:get_numlink_sc: numlink=",numlinks
write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call get_numface_sc(sc,numfaces)
write(*,*) "test:get_numface_sc: numface=",numfaces
write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*) "test:get_link_sc"
do l=1,numlinks
call get_link_sc(sc,l,origin,tip)
  write(*,*) l,"'th link=(",origin,tip,")"
enddo
write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*) "test:get_linklabel_sc"
do i=1,numsites
do j=1,numsites
  if(i .ne. j) then
  call get_linklabel_sc(sc,i,j,link_label,link_dir)
  write(*,*) "(",i,j,"): label=",link_label,",dir=",link_dir
endif
enddo
enddo
write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*) "test:get_links_from_s_sc"
do s=1,numsites
  call get_links_from_s_sc(sc,s,link_labels,tips)
  write(*,*) "links from",s,"going to"
  do i=1,size(link_labels)
    write(*,*) "  ",tips(i),"is label=",link_labels(i)
  enddo
enddo
write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*) "test:get_links_to_s_sc"
do s=1,numsites
  call get_links_to_s_sc(sc,s,link_labels,origins)
  write(*,*) "links to",s,"coming from"
  do i=1,size(link_labels)
    write(*,*) "  ",origins(i),"is label=",link_labels(i)
  enddo
enddo
write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*) "test:get_face_components_sc"
do f=1,numfaces
  call get_face_components_sc(sc,f,sites)
  write(*,*) "face",f,"is made of", sites
enddo
write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*) "test:get_faces_in_site_sc"
do s=1,numsites
  call get_faces_in_site_sc(sc,s,faces)
  write(*,*) "site",s,"is included in the faces",faces 
enddo
write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*) "test:get_faces_in_link_sc"
do l=1,numlinks
  call get_faces_in_link_sc(sc,l,faces)
  write(*,*) "link",l,"is included in the faces",faces
enddo
write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

end subroutine test_module_simplicial_complex


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! check unitaryty
subroutine  check_unitary(UMAT)
use global_parameters, only : NMAT,num_links
use matrix_functions, only : matrix_norm
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
integer l,i,j,k
double precision norm
complex(kind(0d0)) tmp(1:NMAT,1:NMAT)

write(*,*) "===== check unitarity of UMAT ========="
do l=1,num_links
  tmp=(0d0,0d0)
  do i=1,NMAT
    do j=1,NMAT
        do k=1,NMAT
            tmp(i,j)=tmp(i,j)+UMAT(i,k,l)*dconjg(UMAT(j,k,l))
        enddo
    enddo
  enddo
  do i=1,NMAT
    tmp(i,i)=tmp(i,i)-(1d0,0d0)
  enddo
  call matrix_norm(norm,tmp)
    write(*,*) "link:",l,":||U_l.U_l^\dagger - 1||=",norm
enddo

end subroutine check_unitary



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! test  Make_traceless_matrix_from_modes
subroutine  check_Ta_expansion
use mt95
!use global_parameters
use simulation
use SUN_generators, only : make_traceless_matrix_from_modes, trace_MTa
implicit none

double precision MAT_re(1:NMAT,1:NMAT), MAT_im(1:NMAT,1:NMAT)
double precision MODES_re(1:dimG),MODES_im(1:dimG)
complex(kind(0d0)) MAT(1:NMAT,1:NMAT),tmp
complex(kind(0d0)) MAT2(1:NMAT,1:NMAT)
complex(kind(0d0)) MODES(1:dimG),MODES2(1:dimG)
integer i,j,k,a

call genrand_real3(MAT_re)
call genrand_real3(MAT_im)

MAT=MAT_re + (0d0,1d0)*MAT_im
tmp=(0d0,0d0)
do i=1,NMAT-1
    tmp=tmp+MAT(i,i)
enddo
MAT(NMAT,NMAT)=-tmp

do a=1,dimG
    call Trace_MTa(MODES(a),MAT,a,NMAT)
    write(*,*) MODES(a)
enddo

call Make_traceless_matrix_from_modes(MAT2,NMAT,MODES)

write(*,*) "### TEST1 ###"
do i=1,NMAT
do j=1,NMAT
write(*,*) i,j, MAT(i,j)-MAT2(i,j)
enddo
enddo


call genrand_real3(MODES_re)
call genrand_real3(MODES_im)
MODES=dcmplx(MODES_re) + (0d0,1d0)*dcmplx(MODES_im)
call Make_traceless_matrix_from_modes(MAT,NMAT,MODES)
do a=1,dimG
    call Trace_MTa(MODES2(a),MAT,a,NMAT)
enddo

write(*,*) "### TEST2 ###"
do a=1,dimG
write(*,*) MODES(a) - MODES2(a)
enddo


end subroutine check_Ta_expansion


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! check anti-symmetricity of D
!subroutine check_Dirac(UMAT,Phi)
!use Dirac_operator
!implicit none 
!
!complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)), intent(in) :: Phi(1:NMAT,1:NMAT,1:num_sites)
!
!complex(kind(0d0)) :: PF(1:sizeD)
!complex(kind(0d0)) :: Dirac(1:sizeD,1:sizeD), Det
!double precision :: rtmp,dtmp
!integer :: seed
!integer :: i,j
!integer :: s,l,a,b
!
!
!!seed=12345
!!call genrand_init( put=seed )
!! produce pseudo-fermion
!!call make_pseudo_fermion(PF,UMAT,Phi)
!!! calculate Hamiltonian 
!
!call make_Dirac(Dirac,UMAT,Phi)
!
!rtmp=0d0
!dtmp=0d0
!do i=1,sizeD-1
!  rtmp=rtmp+dble(Dirac(i,i)*dconjg(Dirac(i,i)))
!  do j=i+1,sizeD
!    rtmp=rtmp+dble( (Dirac(i,j)+Dirac(j,i)) * dconjg(Dirac(i,j)+Dirac(j,i)) ) 
!  enddo
!enddo
!write(*,*) "# D's diagonal elements?", dtmp+dble(Dirac(sizeD,sizeD)*dconjg(Dirac(sizeD,sizeD)))
!write(*,*) "# Is D anti-symmetric?", rtmp
!
!call make_DdagD(Dirac,UMAT,Phi)
!rtmp=0d0
!do i=1,sizeD
!  do j=i,sizeD
!    rtmp=rtmp&
!        +dble ( &
!          ( Dirac(i,j) - dconjg( Dirac(j,i) ) ) &
!          *dconjg( Dirac(i,j) - dconjg( Dirac(j,i) ) ) )
!  enddo
!enddo
!write(*,*) "# Is D^\dag D hermitian?", rtmp
!
!end subroutine check_Dirac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! check 







