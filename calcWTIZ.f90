!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! code for computing WT identity with IZ compensator
!
program calcwt
use global_parameters
use initialization
use simplicial_complex
use observables
use Dirac_operator, only : make_Dirac
use matrix_functions, only : matrix_inverse, matrix_product
use SUN_generators, only : make_traceless_matrix_from_modes
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! observables 
integer, parameter :: num_obs=12
double precision :: obs(1:num_obs)

character(100) :: obs_name(1:num_obs)
data obs_name/ &
"Re:WT", & ! obs(1)
"Im:WT", & ! obs(2)
"Sb*|A_IZ|", &  ! obs(3)
"-dim(G)*num*|A_IZ|", &  ! obs(4)
"Re(pcsc_mass)*|A_IZ|", &  ! obs(5)
"Re(pcsc_correction*e^{-iarg(A)})", &  ! obs(6)
"Re(correction*e^{-iarg(A)})", &  ! obs(7)
"Im:pcsc_mass * |A_IZ|", &  ! obs(8)
"Im(pcsc_correction*e^{-iarg(A)})", &  ! obs(9)
"Im(correction*e^{-iarg(A)})", &  ! obs(10)
"|A_IZ|", & ! obs(11)
"arg(A_IZ)"/  ! obs(12)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


integer :: ios, num_data, i, j, s, l, f
integer :: num, iargc

integer,parameter :: MED_FILE=10
integer :: OUT_FILE
character(128) :: medfile,outputfile,char_testmode,char_numcon
integer :: testmode, numcon

!! variables
complex(kind(0d0)), allocatable :: UMAT(:,:,:) ! unitary link variables
complex(kind(0d0)), allocatable :: PHI(:,:) ! complex scalar at sites
integer :: ite
complex(kind(0d0)), allocatable ::PhiMat(:,:) ! Phi in a matrix form

complex(kind(0d0)), allocatable :: Dirac(:,:) 
complex(kind(0d0)), allocatable :: Dirac_inv(:,:) 

complex(kind(0d0)), allocatable :: Bfermi(:) ! 1/Nc Tr(lambda lambda(UPhiU+Phi)_l
complex(kind(0d0)), allocatable :: Bbose(:) ! 1/Nc Tr(2 Phi U Phi U^inv)_l
complex(kind(0d0)), allocatable :: CC(:,:) ! 1/Nc (UPhiU+Phi)_(l,a)

!! for anomaly-phase quench
complex(kind(0d0)) :: A_IZ, pcsc_mass,  corr, WT, pcsc_corr
real(8) :: Sb

!! general observable
complex(kind(0d00)) :: observable

num = iargc()
call getarg(1,medfile)

OUT_FILE=6

open(unit=MED_FILE,file=medfile,action='read',form='unformatted')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! initialization
read(MED_FILE) NMAT
read(MED_FILE) LatticeSpacing
read(MED_FILE) SC_FILE_NAME
read(MED_FILE) ALPHA_BETA
read(MED_FILE) test_mode
read(MED_FILE) new_config
read(MED_FILE) fix_seed
read(MED_FILE) read_alpha
read(MED_FILE) config_step
read(MED_FILE) obs_step
read(MED_FILE) m_omega
read(MED_FILE) mass_square_phi
read(MED_FILE) mass_f
read(MED_FILE) Remez_factor4
read(MED_FILE) Remez_factor8
read(MED_FILE) epsilon
read(MED_FILE) CG_max
read(MED_FILE) num_ite
read(MED_FILE) Ntau
read(MED_FILE) Dtau
read(MED_FILE) R_phi
read(MED_FILE) R_A
read(MED_FILE) Fconfigin
read(MED_FILE) Foutput
read(MED_FILE) Fconfigout
read(MED_FILE) Fmedconf
dimG=NMAT*NMAT-1
Dtau_phi = R_phi * Dtau
Dtau_A = R_A * Dtau
!one_ov_2g2N=1d0/(2d0*LatticeSpacing*LatticeSpacing)
overall_factor=dble(NMAT)/(2d0*LatticeSpacing*LatticeSpacing)

call set_NZF
call set_sc
call set_alpha_beta
allocate( UMAT(1:NMAT,1:NMAT,1:num_links) )
allocate( PHI(1:dimG, 1:num_sites) )
allocate( PHIMAT(1:NMAT, 1:NMAT) )
allocate( Dirac(1:sizeD,1:sizeD) )
allocate( Dirac_inv(1:sizeD,1:sizeD) )
allocate( Bbose(1:num_links) )
allocate( Bfermi(1:num_links) )
allocate( CC(1:num_links,1:dimG) )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!do i=1,num_obs
   !write(OUT_FILE,'(I3,a1,a)',advance='no') i+1,")",trim(obs_name(i))
 !enddo
 !write(OUT_FILE,*) 


num_data=0
ios=0
do while (ios == 0)
  num_data=num_data+1
  read(MED_FILE,iostat=ios) ite
  read(MED_FILE,iostat=ios) UMAT
  read(MED_FILE,iostat=ios) PHI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! prepare dirac inverse
  call make_Dirac(Dirac_inv,UMAT,Phi)
  call matrix_inverse(Dirac_inv)      ! compute inverse

  call calc_B(Bbose,Bfermi,CC,Phi,UMAT,Dirac_inv)

  call calc_IZcompensator(A_IZ,Bbose,Bfermi)
  call calc_bosonic_action(Sb,UMAT,Phi)
  call calc_correction(corr,Bbose,Bfermi)
  call calc_pcsc_mass(pcsc_mass,Phi,UMAT,Dirac_inv)
  call calc_pcsc_correction(pcsc_corr,Phi,UMAT,Dirac_inv,Bbose,CC)


  WT=cmplx(dble(-dimG*(num_sites+num_links))/dble(2) * abs(A_IZ))  &
    + cmplx(Sb*abs(A_IZ)) &
    + pcsc_mass*cmplx(abs(A_IZ)) &
    + pcsc_corr*exp( (0d0,-1d0)*arg(A_IZ) ) &
    + corr*exp( (0d0,-1d0)*arg(A_IZ) )

!"Re:A_IZ", &  ! obs(1)
!"Re:Sb*A_IZ", &  ! obs(2)
!"Re:mass contribution", &  ! obs(3)
!"Re:correction", &  ! obs(4)
  obs(1)=dble(WT)
  obs(2)=dble((0d0,-1d0)*WT)
  !!!!
  obs(3)=Sb*abs(A_IZ)
  obs(4)=dble(-dimG*(num_sites+num_links))/dble(2)*abs(A_IZ) 
  obs(5)=dble(pcsc_mass)*abs(A_IZ)
  obs(6)=dble(pcsc_corr*exp( (0d0,-1d0)*arg(A_IZ)))
  obs(7)=dble(corr*exp( (0d0,-1d0)*arg(A_IZ))) 
  !!!!
  obs(8)=dble((0d0,-1d0)*pcsc_corr)*abs(A_IZ)
  obs(9)=dble((0d0,-1d0)*pcsc_corr*exp( (0d0,-1d0)*arg(A_IZ)))
  obs(10)=dble((0d0,-1d0)*corr*exp( (0d0,-1d0)*arg(A_IZ))) 
  !!!!
  obs(11)=abs(A_IZ)
  obs(12)=arg(A_IZ)


  write(OUT_FILE,'(I6,2X)',advance='no') ite
  do i=1,num_obs
    write(OUT_FILE,'(E12.5,2X)',advance='no') OBS(i)
  enddo
  write(OUT_FILE,*) 

enddo

end program calcwt


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! construct 
!!    CC(1:num_links,1:dimG) ! (Ul Phi Ul^inv + Phi)_(l,c)
!!    BB(1:num_links) ! 1/NMAT Tr( lambda lambda CC)
!!    DD(1:num_links) ! 1/NMAT Tr( 2 Phi U Phi U^inv )
subroutine calc_B(Bbose,Bfermi,CC,Phi,UMAT,Dirac_inv)
use global_parameters
use global_subroutines
use SUN_generators, only : make_traceless_matrix_from_modes, trace_MTa
use matrix_functions, only : matrix_product
implicit none

complex(kind(0d0)), intent(out) :: Bbose(1:num_links)
complex(kind(0d0)), intent(out) :: Bfermi(1:num_links)
complex(kind(0d0)), intent(out) :: CC(1:num_links,1:dimG)
complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: Dirac_inv(1:sizeD,1:sizeD)
complex(kind(0d0)) :: phi_tip(1:NMAT,1:NMAT), phi_org(1:NMAT,1:NMAT)
complex(kind(0d0)) :: UphiUdag(1:NMAT,1:NMAT), basemat(1:NMAT,1:NMAT)
!complex(kind(0d0)) :: CC(1:dimG), tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)

!! data of simplicial complex
integer :: l,a,b,c,n,i,j
double precision :: f_abc

!! construct CC and prepare BB
Bbose=(0d0,0d0)
Bfermi=(0d0,0d0)
do l=1,num_links
  call make_traceless_matrix_from_modes(phi_org,NMAT,Phi(:,link_org(l)))
  call make_traceless_matrix_from_modes(phi_tip,NMAT,Phi(:,link_tip(l)))

  ! U Phi U^inv
  call matrix_product(tmpmat,UMAT(:,:,l),phi_tip,'N','N')
  call matrix_product(UphiUdag,tmpmat,UMAT(:,:,l),'N','C')

  !! Bfermi and CC
  do a=1,dimG
    call trace_MTa(CC(l,a), UphiUdag+phi_org, a, NMAT)
    CC(l,a)=CC(l,a)/cmplx(dble(NMAT))
  enddo
  do n=1,NZF
    a=NZF_index(1,n)
    b=NZF_index(2,n)
    c=NZF_index(3,n)
    f_abc=cmplx(NZF_value(n))

    Bfermi(l)=Bfermi(l)&
      +(0d0,0.5d0)*f_abc*Dirac_inv( link_index(a,l), link_index(b,l) )*CC(l,c)
  enddo
  !Bfermi(l)=Bfermi(l)/cmplx(dble(NMAT))

  !! Bbose
  do i=1,NMAT
    do j=1,NMAT
      Bbose(l)=Bbose(l)+(2d0,0d0)*phi_org(i,j)*UphiUdag(j,i)
    enddo
  enddo
  Bbose(l)=Bbose(l)/cmplx(dble(NMAT))

enddo

end subroutine calc_B

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate IZ compensator
subroutine calc_IZcompensator(comp,Bbose,Bfermi)
use global_parameters
use global_subroutines
implicit none

complex(kind(0d0)), intent(in) :: Bbose(1:num_links), Bfermi(1:num_links)
complex(kind(0d0)), intent(out) :: comp
double precision :: r_num
complex(kind(0d0)) :: f_abc
double precision :: g_combination
complex(kind(0d0)) :: c_power

integer :: l,k

r_num = dble( dimG*(num_sites-num_links+num_faces) ) / 4d0 

comp=(0d0,0d0)
do l=1,num_links
  do k=0, dimG/2
    comp = comp & 
      + g_combination( -r_num, k ) & 
        * c_power( Bbose(l), -r_num - dble(k) ) &
        * c_power( Bfermi(l), dble(k) ) 
  enddo
enddo
comp=comp/cmplx(dble(num_links))

end subroutine calc_IZcompensator


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate correction 
subroutine calc_correction(corr,Bbose,Bfermi)
use global_parameters
use global_subroutines
implicit none

complex(kind(0d0)), intent(out) :: corr
complex(kind(0d0)), intent(in) :: Bbose(1:num_links)
complex(kind(0d0)), intent(in) :: Bfermi(1:num_links)

integer :: l,a,b,c,n,k
real(8) :: r_num
complex(kind(0d0)) :: f_abc
double precision :: g_combination
complex(kind(0d0)) :: c_power

r_num = dble( -dimG*(num_sites-num_links+num_faces) ) / 4d0 

corr=(0d0,0d0)
do l=1,num_links
  do k=0, dimG/2-1 
    corr=corr &
      + g_combination( r_num-1, k) &
        * c_power( Bfermi(l), r_num - dble(k+1) ) &
        * c_power( Bbose(l), dble(k+1) )
  enddo
enddo
corr=corr * cmplx(  r_num  / dble(num_links) )

end subroutine calc_correction


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! complex power: power = c^p
complex(kind(0d0)) function c_power( c, p )
use global_subroutines
implicit none

complex(kind(0d0)), intent(in) :: c
real(8), intent(in) :: p

real(8) :: abs_c, arg_c
complex(kind(0d0)) :: f_abc


abs_c=abs(c)
arg_c=arg(c)

c_power = cmplx(abs_c**( p )) * exp( (0d0,1d0) * cmplx( arg_c * p ) )

end function c_power


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! general combination rCk
double precision function g_combination(r,k)
implicit none

double precision, intent(in) :: r
integer, intent(in) :: k

integer :: i, s
double precision :: tmp

s=1
tmp=0d0
do i=1,k 
  if( r-i+1 < 0d0 ) then
    s=s*(-1) 
  endif
  tmp=tmp+log(dble(abs(r-i+1))) - log(dble(i))
enddo
g_combination=exp( tmp )*dble(s)

end function g_combination

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! contruction part of <[mass-terms] A_IZ> 
subroutine calc_pcsc_correction(pcsc_corr,Phi,UMAT,Dirac_inv,Bbose,CC)
use global_parameters
use global_subroutines
use matrix_functions, only : matrix_inverse, matrix_commutator, matrix_power
use SUN_generators, only : make_traceless_matrix_from_modes, trace_mta
implicit none

complex(kind(0d0)), intent(in) :: Phi(1:dimG,1:num_sites)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: Dirac_inv(1:sizeD,1:sizeD)
complex(kind(0d0)), intent(in) :: Bbose(1:num_links)
complex(kind(0d0)), intent(in) :: CC(1:num_links,1:dimG)
complex(kind(0d0)), intent(out) :: pcsc_corr

complex(kind(0d0)) :: AIZ1(1:num_links,1:dimG,1:dimG)

complex(kind(0d0)) :: phi_mat(1:NMAT,1:NMAT),tmp(1:NMAT,1:NMAT)
complex(kind(0d0)) :: comm(1:num_sites,1:dimG)

complex(kind(0d0)) :: bPhi(1:dimG,1:num_sites)
complex(kind(0d0)) :: dphi_mat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: dphi_a(1:num_links,1:dimG)

complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT), Ufm(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Omega(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Omega_a(1:num_faces,1:dimG)

complex(kind(0d0)) :: c_power

complex(kind(0d0)) :: f_abc
integer :: l,a,b,c,d,n,f,l1,l2,s,t

AIZ1=(0d0,0d0)
do l=1,num_links
  do n=1,NZF
    a=NZF_index(1,n)
    b=NZF_index(2,n)
    c=NZF_index(3,n)
    f_abc=cmplx(NZF_value(n))

    AIZ1(l,a,b)=AIZ1(l,a,b) &
      + (0d0,-0.125d0)*cmplx( dble( dimG*(num_sites-num_links+num_faces) )&
        / dble(num_links) ) &
       * f_abc &
       * c_power(Bbose(l), &
                 dble(-dimG*(num_sites-num_links+num_faces))/4d0-1d0)  &
       * CC(l,c)!1/NMATが含まれていることに注意！
  enddo
enddo

bPhi=conjg(Phi)
do l=1,num_links
  call make_diff_phi(dphi_mat,l,UMAT,bPhi)
  do a=1,dimG
    call trace_MTa(dphi_a(l,a),dphi_mat,a,NMAT)
  enddo
enddo

do s=1,num_sites
  call make_traceless_matrix_from_modes(phi_mat, NMAT, Phi(:,s))
  call matrix_commutator(tmp, phi_mat, phi_mat, 'N', 'C')
  do a=1,dimG
    call trace_MTa(comm(s,a), tmp, a, NMAT)
  enddo
enddo

do f=1,num_faces
  call Make_face_variable(Uf,f,UMAT)
  call matrix_power(Ufm,Uf,m_omega)
  call make_moment_map(Omega,Ufm)
  do a=1,dimG
    call trace_MTa(Omega_a(f,a),Omega,a,NMAT)
  enddo
enddo



pcsc_corr=(0d0,0d0)
!Site
do s=1,num_sites
  do t=1,num_sites
    do l=1,num_links
      do a=1,dimG
      do b=1,dimG
        do c=1,dimG
        do d=1,dimG
          pcsc_corr=pcsc_corr&
           + (0.25d0,0d0)*cmplx(alpha_s(s))*comm(s,a)*Phi(t,b)&
             *AIZ1(l,c,d)&
             *Dirac_inv(site_index(b,t),link_index(c,l)) &
             *Dirac_inv(site_index(a,s),link_index(d,l))
        enddo
        enddo
      enddo
      enddo
    enddo
  enddo
enddo

!Link
do s=1,num_sites
  do l1=1,num_links
    do l2=1,num_links
      do a=1,dimG
      do b=1,dimG
        do c=1,dimG
        do d=1,dimG
          pcsc_corr=pcsc_corr&
           + (0d0,-1d0)*cmplx(alpha_l(l1))*dphi_a(l1,a)*Phi(s,b)&
             *AIZ1(l2,c,d)&
             *Dirac_inv(site_index(b,s),link_index(c,l2)) &
             *Dirac_inv(link_index(a,l1),link_index(d,l2))
        enddo
        enddo
      enddo
      enddo
    enddo
  enddo
enddo

!Face
do s=1,num_sites
  do f=1,num_faces
    do l=1,num_links
      do a=1,dimG
      do b=1,dimG
        do c=1,dimG
        do d=1,dimG
          pcsc_corr=pcsc_corr&
           + (0d0,-0.5d0)*cmplx(alpha_f(f)*beta_f(f))&
             *Omega_a(f,a)*Phi(s,b)&
             *AIZ1(l,c,d)&
             *Dirac_inv(site_index(b,s),link_index(c,l)) &
             *Dirac_inv(face_index(a,f),link_index(d,l))
        enddo
        enddo
      enddo
      enddo
    enddo
  enddo
enddo

pcsc_corr=pcsc_corr*cmplx( 2d0 * mass_square_phi * overall_factor**2 )


end subroutine calc_pcsc_correction





