module global_writeDirac
implicit none

character(128), allocatable :: MEDFILE(:)
character(128), allocatable :: DiracFILE(:)

integer, parameter :: N_MEDFILE=100
integer, parameter :: N_DiracFILE=101
integer, parameter :: N_tmpfile=102
character(128), parameter :: tmpfile='abcdefg'
character(128) :: COMMAND

contains
end module global_writeDirac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
use global_parameters
use global_writeDirac
use simulation
use SUN_generators, only : make_SUN_generators
use global_subroutines, only : syncronize_sites, syncronize_links, syncronize_faces
use Dirac_operator, only : Prod_Dirac
use parallel
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! variables
complex(kind(0d0)), allocatable :: UMAT(:,:,:) ! unitary link variables
complex(kind(0d0)), allocatable :: PHIMAT(:,:,:) ! complex scalar at sites
complex(kind(0d0)), allocatable :: PFeta(:,:,:) 
complex(kind(0d0)), allocatable :: PFlambda(:,:,:) 
complex(kind(0d0)), allocatable :: PFchi(:,:,:) 
complex(kind(0d0)), allocatable :: DPFeta(:,:,:) 
complex(kind(0d0)), allocatable :: DPFlambda(:,:,:) 
complex(kind(0d0)), allocatable :: DPFchi(:,:,:) 
complex(kind(0d0)), allocatable :: T(:,:,:) ! SU(N) generators

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! FOR PARALLEL
type(SITE_DIST), allocatable,save :: local_site_list(:) !(0:NPROCS-1)

! for type genrand_srepr
!integer :: MPI_GENRAND_SREPR, MPI_LOCAL_LABEL, MPI_GENRAND_STATE
integer IBLOCK1(1),IDISP1(1),ITYPE1(1)
integer IBLOCK2(1),IDISP2(1),ITYPE2(1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: seed
integer :: gs,gl,gf,s,l,f
integer :: i,j,k,narg
integer :: a,rank
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: iarg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: control, info, ite
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
character(128) :: FMT1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
iarg=iargc()
allocate( MEDFILE(1:iarg) )
allocate( DiracFILE(1:iarg) )
if( iarg ==0 ) then
  if (MYRANK==0) write(*,*) "use as a.out [MEDFILE-1] [MEDFILE-2]..."
  stop
endif
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call MPI_INIT(IERR)
call MPI_COMM_SIZE(MPI_COMM_WORLD,NPROCS,IERR)
call MPI_COMM_RANK(MPI_COMM_WORLD,MYRANK,IERR)

! setup MPI_GENRAND_SREPR
IBLOCK1(1)=1
IDISP1(1)=0
ITYPE1(1)=MPI_CHARACTER
call MPI_TYPE_STRUCT(1,IBLOCK1,IDISP1,ITYPE1,MPI_GENRAND_SREPR,IERR)
call MPI_TYPE_COMMIT(MPI_GENRAND_SREPR,IERR)

! setup MPI_LOCAL_LABEL
IBLOCK2(1)=2
IDISP2(1)=0
ITYPE2(1)=MPI_INTEGER
call MPI_TYPE_STRUCT(1,IBLOCK2,IDISP2,ITYPE2,MPI_LOCAL_LABEL,IERR)
call MPI_TYPE_COMMIT(MPI_LOCAL_LABEL,IERR)

allocate(local_site_list(0:NPROCS-1) )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call set_theory_parameters(seed)
! set the structure constant of SU(NMAT)
call set_NZF
! set global data of the simplicial complex
call set_global_simplicial_complex 
! set global parameters
overall_factor=dble(NMAT)/(2d0*LatticeSpacing*LatticeSpacing)
mass_square_phi = phys_mass_square_phi * (LatticeSpacing*latticeSpacing)
! set simulation data
!   ここで、inputfileを読み込み、各コアに必要な情報を設定して、
!   Rank=0にたまっているalphaを振り分ける。
!   この段階で、
!      num_sites,num_links,num_faces
!      local_{site,link,face}_of_globalの(1:num_{sites,links,faces})が確定する。
call set_simulation_parameters(local_site_list)

! set local data
call set_local_data(local_site_list)

allocate( UMAT(1:NMAT,1:NMAT,1:num_necessary_links) )
allocate( PHIMAT(1:NMAT,1:NMAT,1:num_necessary_sites) )
allocate( PFeta(1:NMAT,1:NMAT,1:num_necessary_sites) )
allocate( PFlambda(1:NMAT,1:NMAT,1:num_necessary_links) )
allocate( PFchi(1:NMAT,1:NMAT,1:num_necessary_faces) )
allocate( DPFeta(1:NMAT,1:NMAT,1:num_sites) )
allocate( DPFlambda(1:NMAT,1:NMAT,1:num_links) )
allocate( DPFchi(1:NMAT,1:NMAT,1:num_faces) )
allocate( T(1:NMAT,1:NMAT,1:dimG) )

! SU(N) generators
call make_SUN_generators(T,NMAT)

if( MYRANK==0 ) then
do narg=1,iarg
  call getarg(narg,MEDFILE(narg))
  COMMAND = 'echo "' // trim(adjustl(MEDFILE(narg))) // '" | sed "s/medconfig/Dirac/" > ' // trim(adjustl(tmpfile) )
  call system(COMMAND)
  !!!
  if( MYRANK == 0 ) then
    open(N_tmpfile, file=tmpfile, status='OLD',action='READ')
    READ(N_tmpfile,'(a)') DiracFILE(narg)
    close(N_tmpfile)
  endif
enddo
endif

do narg=1,iarg
  if( MYRANK==0 ) then
    open(N_MEDFILE, file=MEDFILE(narg), status='OLD',action='READ',form='unformatted')
    open(N_DiracFILE, file=DiracFILE(narg), status='REPLACE')
  endif
  !!!!
  do
    call read_config_from_medfile(Umat,PhiMat,ite,N_MEDFILE,control)
    if( control == 0 ) then 
      if( MYRANK==0 )  then
        write(N_DiracFILE,'(a,2X,10I)') 'I', ite !END 
      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! make Dirac matrix
      J=0
      do gs=1, global_num_sites
        s=local_site_of_global(gs)%label_
        rank=local_site_of_global(gs)%rank_

        PFlambda=(0d0,0d0)
        PFchi=(0d0,0d0)
        do a=1, dimG
          J=J+1
          PFeta=(0d0,0d0)
          if( MYRANK == rank ) then
            PFeta(:,:,s) = T(:,:,a)
          endif
          call syncronize_sites(PFeta)

          !call calc_DinvF(DPFeta, DPFlambda, DPFchi, &
            !PFeta, PFlambda, PFchi, Umat, PhiMat, info)

          call Prod_Dirac(&
            DPFeta, DPFlambda, DPFchi, &
            PFeta,PFlambda,PFchi, UMAT,PhiMat)

          !!!
          call write_PF(N_DiracFILE,DPFeta,DPFlambda,DPFchi,J)
        enddo
      enddo
      !!!!!!!!!!!!!!
      do gl=1, global_num_links
        l=local_link_of_global(gl)%label_
        rank=local_link_of_global(gl)%rank_

        PFeta=(0d0,0d0)
        PFchi=(0d0,0d0)
        do a=1, dimG
          J=J+1
          PFlambda=(0d0,0d0)
          if( MYRANK == rank ) then
            PFlambda(:,:,l) = T(:,:,a)
          endif
          call syncronize_links(PFlambda)

          !call calc_DinvF(DPFeta, DPFlambda, DPFchi, &
            !PFeta, PFlambda, PFchi, Umat, PhiMat, info)
          call Prod_Dirac(&
            DPFeta, DPFlambda, DPFchi, &
            PFeta,PFlambda,PFchi, UMAT,PhiMat)
          !!!
          call write_PF(N_DiracFILE,DPFeta,DPFlambda,DPFchi,J)
        enddo
      enddo
      !!!!!!!!!!!!!!
      do gf=1, global_num_faces
        f=local_face_of_global(gf)%label_
        rank=local_face_of_global(gf)%rank_

        PFeta=(0d0,0d0)
        PFlambda=(0d0,0d0)
        do a=1, dimG
          J=J+1
          PFchi=(0d0,0d0)
          if( MYRANK == rank ) then
            PFchi(:,:,f) = T(:,:,a)
          endif
          call syncronize_faces(PFchi)

          !call calc_DinvF(DPFeta, DPFlambda, DPFchi, &
            !PFeta, PFlambda, PFchi, Umat, PhiMat, info)
          call Prod_Dirac(&
            DPFeta, DPFlambda, DPFchi, &
            PFeta,PFlambda,PFchi, UMAT,PhiMat)
          !!!
          call write_PF(N_DiracFILE,DPFeta,DPFlambda,DPFchi,J)
        enddo
      enddo
      !!!!!!!!!!!!!!
      if( MYRANK==0 )  then
        write(N_DiracFILE,'(a)') 'E' !END 
      endif
    else
      exit
    endif
  enddo
  if( MYRANK == 0 ) then
    close(N_MEDFILE)
  endif
enddo
end program main

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! write eta, lambda, chi to N_FILE (which is opened already)
subroutine write_PF(N_FILE,eta,lambda,chi,J)
use global_parameters
!use matrix_functions, only : matrix_norm
use parallel
implicit none

integer, intent(in) :: N_FILE
complex(kind(0d0)), intent(in) :: eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(in) :: lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: chi(1:NMAT,1:NMAT,1:num_faces)
integer, intent(in) :: J !! J of D_{IJ}
integer :: I !! I of D_{IJ}

integer :: gs,gl,gf
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
double precision :: norm

I=0
do gs=1,global_num_sites
  call sendPFeta(tmpmat,gs,eta)
  !call matrix_norm(norm, tmpmat)
  !if( norm > epsilon ) then 
    call write_modes(N_FILE,tmpmat,I,J)
  !else
    !I=I+dimG
  !endif
enddo
!!!
do gl=1,global_num_links
  call sendPFlambda(tmpmat,gl,lambda)
  !if( norm > epsilon ) then 
    call write_modes(N_FILE,tmpmat,I,J)
  !else
    !I=I+dimG
  !endif
enddo
!!!
do gf=1,global_num_faces
  call sendPFchi(tmpmat,gf,chi)
  !if( norm > epsilon ) then 
    call write_modes(N_FILE,tmpmat,I,J)
  !else
    !I=I+dimG
  !endif
enddo

end subroutine write_PF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! write matrix in mode to N_FILE
subroutine write_modes(N_FILE, mat,I,J)
use global_parameters
use SUN_generators, only : trace_MTa
use parallel
implicit none


complex(kind(0d0)), intent(in) :: mat(1:NMAT,1:NMAT)
integer, intent(in) :: N_FILE
integer, intent(in) :: J !! J of D_{IJ}
integer, intent(inout) :: I !! I of D_{IJ}

integer a
complex(kind(0d0)) :: tmp

if( MYRANK==0 ) then
  do a=1,dimG
    I=I+1
    call trace_MTa(tmp,mat,a,NMAT)
    if( cdabs(tmp) > epsilon ) then 
      write(N_FILE,'(a,2X,I10,2X,I10,2X,E15.8,2X,E15.8,2X)') &
        'D',I,J,dble(tmp), dble((0d0,-1d0)*tmp)
    endif
  enddo
endif

end subroutine write_modes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! write eta, lambda, chi to N_FILE (which is opened already)
!subroutine PFtoRank0(mat,gp,FTYPE,PsF)
subroutine sendPFeta(mat,gp,PFeta)
use global_parameters
use parallel
implicit none

complex(kind(0d0)), intent(out) :: mat(1:NMAT,1:NMAT)
integer, intent(in) :: gp
complex(kind(0d0)), intent(in) :: PFeta(1:NMAT,1:NMAT,1:num_sites)

integer :: tag, lp,rank

  tag=gp
  lp=local_site_of_global(gp)%label_
  rank=local_site_of_global(gp)%rank_

if( MYRANK == 0 ) then 
  if(rank==0) then
    mat=PFeta(:,:,lp)
  else
    call MPI_RECV(mat,NMAT*NMAT,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,ISTATUS,IERR)
  endif
else
  if( MYRANK == rank ) then
    call MPI_SEND(PFeta(:,:,lp),NMAT*NMAT,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,IERR)
  endif
endif
end subroutine sendPFeta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
subroutine sendPFlambda(mat,gp,PFlambda)
use global_parameters
use parallel
implicit none

complex(kind(0d0)), intent(out) :: mat(1:NMAT,1:NMAT)
integer, intent(in) :: gp
complex(kind(0d0)), intent(in) :: PFlambda(1:NMAT,1:NMAT,1:num_links)

integer :: tag, lp,rank

  tag=global_num_sites+gp
  lp=local_link_of_global(gp)%label_
  rank=local_link_of_global(gp)%rank_

if( MYRANK == 0 ) then 
  if(rank==0) then
    mat=PFlambda(:,:,lp)
  else
    call MPI_RECV(mat,NMAT*NMAT,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,ISTATUS,IERR)
  endif
else
  if( MYRANK == rank ) then
    call MPI_SEND(PFlambda(:,:,lp),NMAT*NMAT,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,IERR)
  endif
endif
end subroutine sendPFlambda


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
subroutine sendPFchi(mat,gp,PFchi)
use global_parameters
use parallel
implicit none

complex(kind(0d0)), intent(out) :: mat(1:NMAT,1:NMAT)
integer, intent(in) :: gp
complex(kind(0d0)), intent(in) :: PFchi(1:NMAT,1:NMAT,1:num_faces)

integer :: tag, lp,rank

  tag=global_num_sites+global_num_links+gp
  lp=local_face_of_global(gp)%label_
  rank=local_face_of_global(gp)%rank_

if( MYRANK == 0 ) then 
  if(rank==0) then
    mat=PFchi(:,:,lp)
  else
    call MPI_RECV(mat,NMAT*NMAT,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,ISTATUS,IERR)
  endif
else
  if( MYRANK == rank ) then
    call MPI_SEND(PFchi(:,:,lp),NMAT*NMAT,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,IERR)
  endif
endif
end subroutine sendPFchi








