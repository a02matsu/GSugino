!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Test code to use library PZGETRI in PBLAS
! which compute the inverse matrix using the LU factorization
! computed by PZGETRF.
! We now use it to compute the inverse of Dirac matrix
! So Matsuura 2019/3/12
!
! ここでは16x16行列Aを4コアで並列化して計算することにする
! 1. Initialize the process grid
! 2. Distribute the matrix on the process grid
! 3. Call ScaLAPACK routine
! 4. Release the process grid 
!
!! compile
! % mpiifort -mkl=cluster inverse_matrix_pblas.f90
program calcDinv
use mpi
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!
! 読み込み／書き込みファイル
character(128) :: INFILE !="MEDCONF/Dirac_029810+000300.dat"
character(128) :: OUTFILE !="MEDCONF/Dinv_029810+000300.dat"
integer :: N_INFILE=100
integer :: N_OUTFILE=101

character :: TRIG
integer :: C_TRIG
integer :: N_TRIG
integer :: ite


!!!!!!!!!!!!!!!!!!!!!!!!!!
! simulation data
integer :: NMAT !=3
integer :: num_sites !=32
integer :: num_links !=56
integer :: num_faces !=26
character(50) :: C_NMAT !=3
character(50) :: C_num_sites !=32
character(50) :: C_num_links !=56
character(50) :: C_num_faces !=26

integer :: sizeD

integer :: s,l,f
integer :: ss,ll,ff

double precision :: real, imag
complex(kind(0d0)) :: ele

integer :: control

!!!!!!!!!!!!!!!!!!!!!!!!!!
! process gridの初期化
integer :: NPROW != 2 ! 行のプロセス数
integer :: NPCOL != 2 ! 列のプロセス数
integer :: MB  != 8 ! 行のブロック数 
integer :: NB  != 8 ! 列のブロック数
integer :: TEST
character(50) :: C_NPROW != 2 ! 行のプロセス数
character(50) :: C_NPCOL != 2 ! 列のプロセス数
character(50) :: C_MB != 8 ! 行のブロック数 
character(50) :: C_NB != 8 ! 列のブロック数
character(50) :: C_TEST


!!!!!!!!!!!!!!!!!!!!!!!!!!
integer ICTXT ! context (BLACS で通信する process grid のラベル)
integer INFO, IAM, NPROCS ! BLACS_PINFO用 
integer, parameter :: RSRC = 0, CSRC = 0 ! 最初のブロックをどのプロセスに配置するか。
!! for MPI
integer :: MYRANK, IERR, ISTATUS, tag, RANK

integer :: MXLLD ! 各gridに配置される部分行列Aの中での最大の大きさ
integer MYROW, MYCOL ! 各プロセスのラベル
integer L_ROW, L_COL !それぞれの process grid にあるlocal matrixの行と列のサイズ
integer NUMROC ! L_ROW と L_COL を計算する関数
integer BLACS_PNUM ! process numberを返す函数   
integer LLD

!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Matrix
integer, parameter :: control_out=1 
complex(kind(0d0)), allocatable :: Dirac(:,:) ! Diracを配置する配列
complex(kind(0d0)), allocatable :: Dinv(:,:) ! D^{-1}を配置する配列
complex(kind(0d0)), allocatable :: prod(:,:) ! D.D^{-1}を配置する配列
integer DESC_A(9) ! 行列Aの descriptor
! descriptorとは、global arrayと対応するプロセスやメモリとの
! mapping情報を格納する変数 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! for output
integer IPW ! なぞ。PDLAPRNT 用
integer :: NOUT = 6 ! 標準出力
complex(kind(0d0)), allocatable :: work(:), work1(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! for PZGETRI
integer :: IA
integer :: JA
integer, allocatable :: IPIV(:)
integer :: LWORK
integer, allocatable :: IWORK(:)
integer :: LIWORK

integer :: LOCr, LOCc
integer :: LCM
integer :: NQ ! =NUMROC(sizeD+mod(JA-1,NB), NB, MYCOL, CSRC, NPCOL)
integer :: MQ ! =NUMROC(sizeD+mod(JA-1,NB), MB, MYCOL, CSRC, NPCOL)
integer :: MP ! =NUMROC(sizeD+mod(JA-1,NB), MB, MYROW, RSRC, NPROW)
integer :: LIWORK1, LIWORK2
integer :: ILCM
integer :: ICEIL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Others
integer Pr, Pc, X, Y ! Row process and local position corresponding to global index i,j
integer i,j
integer IROW, ICOL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: iarg

! 0) initialization
iarg=iargc()
if( iarg <= 9 ) stop
call getarg(1,INFILE)
call getarg(2,OUTFILE)
call getarg(3,C_NMAT)
call getarg(4,C_num_sites)
call getarg(5,C_num_links)
call getarg(6,C_num_faces)
call getarg(7,C_NPROW)    
call getarg(8,C_NPCOL)
call getarg(9,C_MB)    
call getarg(10,C_NB)
if( iarg >= 11 ) then
  call getarg(11,C_TEST)    
  read(C_TEST,*) TEST
else
  TEST=0
endif

read(C_NMAT,*) NMAT
read(C_num_sites,*) num_sites
read(C_num_links,*) num_links
read(C_num_faces,*) num_faces
read(C_NPROW,*) NPROW
read(C_NPCOL,*) NPCOL
read(C_MB,*) MB
read(C_NB,*) NB

sizeD = (NMAT*NMAT-1)*(num_sites+num_links+num_faces)
NPROCS=NPROW*NPCOL

! 1) Create Process Grid
call SL_INIT( ICTXT, NPROW, NPCOL )

call MPI_COMM_SIZE(MPI_COMM_WORLD, NPROCS, info)
call MPI_COMM_RANK(MPI_COMM_WORLD, MYRANK, IERR)

!call BLACS_PINFO( IAM, NPROCS )
call BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )

L_ROW = NUMROC( sizeD, MB, MYROW, RSRC, NPROW ) ! number of row
L_COL = NUMROC( sizeD, NB, MYCOL, CSRC, NPCOL ) ! number of col
MXLLD = max(L_ROW, L_COL)
LLD=max(1,L_ROW)

!write(*,*) L_ROW, L_COL
! 2) Distribute the matrix on the process grid
call DESCINIT( DESC_A, sizeD, sizeD, MB, NB, RSRC, CSRC, ICTXT, MXLLD, INFO )

! 3) allocate memory
allocate( Dinv(1:L_ROW, 1:L_COL) )
if( TEST==1 ) then
  allocate( Dirac(1:L_ROW, 1:L_COL) )
  allocate( prod(1:L_ROW, 1:L_COL) )
  allocate( work1(1:sizeD) )
endif

! 4) open the Dirac file and output file
if( MYROW==0 .and. MYCOL==0 ) then
  open(N_INFILE, file=INFILE, status='OLD',action='READ')
  open(N_OUTFILE, file=OUTFILE, status='REPLACE')

  control=0
endif

! 5) setup the matrix A
!! for PZGETRF
IA=1
JA=1
!allocate( IPIV( L_ROW+MB ) )
allocate( IPIV( MB+L_ROW ) )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! for PZGETRI
 LWORK=NUMROC(sizeD+mod(IA-1,MB), MB, MYROW, RSRC, NPCOL)*NB
 allocate( WORK(LWORK) )
 MP=NUMROC(sizeD, MB, MYROW, RSRC, NPROW)
 MQ=NUMROC(sizeD, MB, MYCOL, CSRC, NPCOL)
 NQ=NUMROC(sizeD, NB, MYCOL, CSRC, NPCOL)
 LCM=ILCM(NPROW,NPCOL)
 !NQ=NUMROC(sizeD+mod(JA-1,NB), NB, MYCOL, CSRC, NPCOL)
 !MQ=NUMROC(sizeD+mod(JA-1,NB), MB, MYCOL, CSRC, NPCOL)
 !MP=NUMROC(sizeD+mod(JA-1,NB), MB, MYROW, RSRC, NPROW)
 if( NPROW==NPCOL ) then
   LIWORK=NQ+NB
 else
   LIWORK1=MB*ICEIL(ICEIL(MP,MB), LCM/NPROW)
   LIWORK2=NB*ICEIL(ICEIL(MQ,NB), LCM/NPCOL)
   LIWORK=NQ+max(LIWORK1,LIWORK2,NB)
 endif
 allocate( IWORK(LIWORK) )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do 
  if( MYROW==0 .and. MYCOL==0 ) then
    read(N_INFILE,*, END=1111) trig, ite
    write(N_OUTFILE,'(I,2X)',advance='no') ite
  endif
  call MPI_BCAST(control,1,MPI_INTEGER,BLACS_PNUM(0,0),MPI_COMM_WORLD,IERR)
  call MPI_BCAST(ite,1,MPI_INTEGER,BLACS_PNUM(0,0),MPI_COMM_WORLD,IERR)
    !call IGEBS2D(ICTXT,'ALL','i-ring',1,1,control,1)
  !else
    !call IGEBR2D(ICTXT,'ALL','i-ring',1,1,control,1,0,0)
  !endif
  if( control == 1 ) stop

  Dinv=(0d0,0d0)
  if( TEST==1 ) Dirac=(0d0,0d0)
  do
    if( MYROW==0 .and. MYCOL==0 ) then
      read(N_INFILE,'(a1)',advance='no') trig
      if( trig == 'D' ) then
        C_trig = 0
      elseif( trig == 'E' ) then
        read(N_INFILE,*)
        C_trig = 1
      else
        C_trig = 2
      endif
!      call IGEBS2D(ICTXT,'ALL','i-ring',1,1,C_trig,1)
!    else
!      call IGEBR2D(ICTXT,'ALL','i-ring',1,1,C_trig,1,0,0)
    endif
    call MPI_BCAST(C_trig,1,MPI_INTEGER,BLACS_PNUM(0,0),MPI_COMM_WORLD,IERR)
    if( C_trig == 0 ) then
      N_trig=0
      if( MYROW==0 .and. MYCOL==0 ) then
        read(N_INFILE,*) I,J,real,imag
        ele = dcmplx(real) + (0d0,1d0)*dcmplx(imag)
!        call IGEBS2D(ICTXT,'ALL','i-ring',1,1,I,1)
!        call IGEBS2D(ICTXT,'ALL','i-ring',1,1,J,1)
!        call ZGEBS2D(ICTXT,'ALL','i-ring',1,1,ele,1)
!      else
!        call IGEBR2D(ICTXT,'ALL','i-ring',1,1,I,1,0,0)
!        call IGEBR2D(ICTXT,'ALL','i-ring',1,1,J,1,0,0)
!        call ZGEBR2D(ICTXT,'ALL','i-ring',1,1,ele,1,0,0)
      endif
      call MPI_BCAST(I,1,MPI_INTEGER,BLACS_PNUM(0,0),MPI_COMM_WORLD,IERR)
      call MPI_BCAST(J,1,MPI_INTEGER,BLACS_PNUM(0,0),MPI_COMM_WORLD,IERR)
      call MPI_BCAST(ele,1,MPI_DOUBLE_COMPLEX,BLACS_PNUM(0,0),MPI_COMM_WORLD,IERR)
      call PZELSET( Dinv, I, J, DESC_A, ele )
    elseif( C_trig==1 ) then
      exit
    elseif( C_trig==2 ) then
      write(*,*) "something happened"
      stop
    endif
  enddo
  !if( .not. allocated(work1) ) allocate(work1(1:sizeD))
  !call PZLAPRNT(sizeD, sizeD, Dinv, 1, 1, DESC_A, 0, 0, 'Dirac', NOUT, work1 )
  !deallocate(work1)
  !stop
  !if( MYROW==0 .and. MYCOL==0 ) write(*,*) "================="


! 5) compute matrix inverse
  if( TEST==1 ) then 
    Dirac=Dinv
  endif
  call PZGETRF(sizeD,sizeD,Dinv,IA,JA,DESC_A,IPIV,INFO)
  !!
  call PZGETRI(sizeD,Dinv,IA,JA,DESC_A, IPIV, WORK, LWORK, IWORK, LIWORK, INFO)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! TEST 
if( TEST==1 ) then
  call PZGEMM( 'N', 'N', sizeD, sizeD, sizeD, (1d0,0d0), &
       Dirac, 1, 1, DESC_A, &
       Dinv, 1, 1, DESC_A, &
       (0d0,0d0), prod, 1, 1, DESC_A )
  call PZLAPRNT(sizeD, sizeD, prod, 1, 1, DESC_A, 0, 0, 'D.Dinv', NOUT, work1 )
  if( MYROW==0 .and. MYCOL==0 ) then
    write(*,*) "==========================================="
  endif
endif
!stop
!! END TEST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! 6) write out Dinv
!  deallocate(work)
!  allocate(work(1:sizeD))
!  call PZLAPRNT(sizeD, sizeD, Dinv, 1, 1, DESC_A, 0, 0, 'D.Dinv', NOUT, work )

  !if( .not. allocated(work1) ) allocate(work1(1:sizeD))
  !call PZLAPRNT(sizeD, sizeD, Dinv, 1, 1, DESC_A, 0, 0, 'Dinv', NOUT, work1 )
  !deallocate(work1)
  do i=1,sizeD
    do j=1,sizeD
      call PZELGET('','i-ring',ele,Dinv,i,j,DESC_A)
      !CALL INFOG2L( I, J, DESC_A, NPROW, NPCOL, MYROW, MYCOL, IA, JA, IROW, ICOL )
      IROW=mod((I-1)/MB,NPROW)
      ICOL=mod((J-1)/NB,NPCOL)
      if( IROW/=0 .or. ICOL/=0 ) then
        if( MYROW==IROW .and. MYCOL==ICOL ) then
          call MPI_SEND(ele,1,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,IERR)
          !call ZGEBS2D(ICTXT,'ALL','i-ring',1,1,ele,1)
          !call ZGESD2D( ICTXT,1,1,ele,1,0,0)
        elseif( MYROW==0 .and. MYCOL==0 ) then
          call MPI_RECV(ele,1,MPI_DOUBLE_COMPLEX,&
            BLACS_PNUM(ICTXT,IROW,ICOL),tag,MPI_COMM_WORLD,ISTATUS,IERR)
          !call ZGERV2D( ICTXT,1,1,ele,1,IROW,ICOL )
          !call ZGEBR2D(ICTXT,'ALL','i-ring',1,1,ele,1,IROW,ICOL)
        endif
      endif
      call MPI_BARRIER(MPI_COMM_WORLD,IERR)
      !call BLACS_BARRIER( ICTXT, 'A' )
      if( MYROW==0 .and. MYCOL==0 ) then 
        write(N_OUTFILE,'(E15.8,2X,E15.8,2X)',advance='no') ele
        !write(N_OUTFILE,*) i,j,ele
      endif
    enddo
  enddo
  if( MYROW==0 .and. MYCOL==0 ) then
    write(N_OUTFILE,*)
    !write(*,*) "====="
  endif
  !stop
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
1111 control = 1
call IGEBS2D(ICTXT,'ALL','i-ring',1,1,control,1)
stop

end program calcDinv

