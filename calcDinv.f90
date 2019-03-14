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
integer :: MB != 8 ! 行のブロック数 
integer :: NB != 8 ! 列のブロック数
character(50) :: C_NPROW != 2 ! 行のプロセス数
character(50) :: C_NPCOL != 2 ! 列のプロセス数
character(50) :: C_MB != 8 ! 行のブロック数 
character(50) :: C_NB != 8 ! 列のブロック数

!!!!!!!!!!!!!!!!!!!!!!!!!!
integer ICTXT ! context (BLACS で通信する process grid のラベル)
integer INFO, IAM, NPROCS ! BLACS_PINFO用 
integer, parameter :: RSRC = 0, CSRC = 0 ! 最初のブロックをどのプロセスに配置するか。

integer :: MXLLD ! 各gridに配置される部分行列Aの中での最大の大きさ
integer MYROW, MYCOL ! 各プロセスのラベル
integer L_ROW, L_COL !それぞれの process grid にあるlocal matrixの行と列のサイズ
integer NUMROC ! L_ROW と L_COL を計算する関数
integer LLD

!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Matrix
integer, parameter :: control_out=1 
!complex(kind(0d0)), allocatable :: Dirac(:,:) ! Diracを配置する配列
complex(kind(0d0)), allocatable :: Dinv(:,:) ! D^{-1}を配置する配列
!complex(kind(0d0)), allocatable :: prod(:,:) ! D.D^{-1}を配置する配列
integer DESC_A(9) ! 行列Aの descriptor
! descriptorとは、global arrayと対応するプロセスやメモリとの
! mapping情報を格納する変数 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! for output
integer IPW ! なぞ。PDLAPRNT 用
integer :: NOUT = 6 ! 標準出力
complex(kind(0d0)), allocatable :: work(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! for PZGETRI
integer :: IA
integer :: JA
integer, allocatable :: IPIV(:)
integer :: LWORK
integer, allocatable :: IWORK(:)
integer :: LIWORK

integer :: LOCr, LOCc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Others
integer Pr, Pc, X, Y ! Row process and local position corresponding to global index i,j
integer i,j
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: iarg

! 0) initialization
iarg=iargc()
if( iarg .ne. 10 ) stop
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

read(C_NMAT,*) NMAT
read(C_num_sites,*) num_sites
read(C_num_links,*) num_links
read(C_num_faces,*) num_faces
read(C_NPROW,*) NPROW
read(C_NPCOL,*) NPCOL
read(C_MB,*) MB
read(C_NB,*) NB

sizeD = (NMAT*NMAT-1)*(num_sites+num_links+num_faces)


! 1) Create Process Grid
call BLACS_PINFO( IAM, NPROCS )
call SL_INIT( ICTXT, NPROW, NPCOL )
call BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )

L_ROW = NUMROC( sizeD, MB, MYROW, RSRC, NPROW ) ! number of row
L_COL = NUMROC( sizeD, NB, MYCOL, CSRC, NPCOL ) ! number of col
MXLLD = max(L_ROW, L_COL)
LLD=max(1,L_ROW)

! 2) Distribute the matrix on the process grid
call DESCINIT( DESC_A, sizeD, sizeD, MB, NB, RSRC, CSRC, ICTXT, MXLLD, INFO )

! 3) allocate memory
allocate( Dinv(1:L_ROW, 1:L_COL) )
!allocate( Dirac(1:L_ROW, 1:L_COL) )
!allocate( prod(1:L_ROW, 1:L_COL) )

! 4) open the Dirac file and output file
Dinv=(0d0,0d0)
if( MYROW==0 .and. MYCOL==0 ) then
  open(N_INFILE, file=INFILE, status='OLD',action='READ')
  open(N_OUTFILE, file=OUTFILE, status='REPLACE')

  control=0
endif

! 5) setup the matrix A
!! for PZGETRF
IA=1
JA=1
allocate( IPIV( L_ROW+MB ) )
!! for PZGETRI
 LWORK=NUMROC(sizeD+mod(IA-1,MB), MB, MYROW, RSRC, NPROW)*NB
 allocate( WORK(LWORK) )
 LIWORK=NUMROC(sizeD+mod(JA-1,NB), NB, MYCOL, CSRC, NPROW)+NB
 allocate( IWORK(LIWORK) )
do 
  if( MYROW==0 .and. MYCOL==0 ) then
    read(N_INFILE,*, END=1111) trig, ite
    write(N_OUTFILE,'(I,2X)',advance='no') ite
    call IGEBS2D(ICTXT,'ALL','i-ring',1,1,control,1)
  else
    call IGEBR2D(ICTXT,'ALL','i-ring',1,1,control,1,0,0)
    if( control == 1 ) stop
  endif

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
      call IGEBS2D(ICTXT,'ALL','i-ring',1,1,C_trig,1)
    else
      call IGEBR2D(ICTXT,'ALL','i-ring',1,1,C_trig,1,0,0)
    endif
    if( C_trig == 0 ) then
      N_trig=0
      if( MYROW==0 .and. MYCOL==0 ) then
        read(N_INFILE,*) I,J,real,imag
        call global_to_local(Pr,Pc,X,Y,I,J,MB,NB,NPROW,NPCOL)
        ele = dcmplx(real) + (0d0,1d0)*dcmplx(imag)

        call IGEBS2D(ICTXT,'ALL','i-ring',1,1,Pr,1)
        call IGEBS2D(ICTXT,'ALL','i-ring',1,1,Pc,1)
        !!!!
        if( Pr/=0 .or. Pc/=0 ) then 
          call IGESD2D(ICTXT,1,1,X,1,Pr,Pc)
          call IGESD2D(ICTXT,1,1,Y,1,Pr,Pc)
          call ZGESD2D(ICTXT,1,1,ele,1,Pr,Pc)
        else
          Dinv(X,Y) = ele
        endif
      else
        call IGEBR2D(ICTXT,'ALL','i-ring',1,1,Pr,1,0,0)
        call IGEBR2D(ICTXT,'ALL','i-ring',1,1,Pc,1,0,0)
        !!
        if( MYROW==Pr .and. MYCOL==Pc ) then
          call IGERV2D(ICTXT,1,1,X,1,0,0)
          call IGERV2D(ICTXT,1,1,Y,1,0,0)
          call ZGERV2D(ICTXT,1,1,ele,1,0,0)
          Dinv(X,Y) = ele
        endif
      endif
    elseif( C_trig==1 ) then
      exit
    elseif( C_trig==2 ) then
      write(*,*) "something happened"
      stop
    endif
  enddo
  !allocate(work(1:sizeD))
  !call PZLAPRNT(sizeD, sizeD, Dinv, 1, 1, DESC_A, 0, 0, 'Dirac', NOUT, work )
  !deallocate(work)
  !stop


! 5) compute matrix inverse
  !Dirac=Dinv
  call PZGETRF(sizeD,sizeD,Dinv,IA,JA,DESC_A,IPIV,INFO)
  !!
  call PZGETRI(sizeD,Dinv,IA,JA,DESC_A, IPIV, WORK, LWORK, IWORK, LIWORK, INFO)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! TEST 
!call PZGEMM( 'N', 'N', sizeD, sizeD, sizeD, (1d0,0d0), &
!     Dirac, 1, 1, DESC_A, &
!     Dinv, 1, 1, DESC_A, &
!     (0d0,0d0), prod, 1, 1, DESC_A )
!deallocate(work)
!allocate(work(1:sizeD))
!call PZLAPRNT(sizeD, sizeD, prod, 1, 1, DESC_A, 0, 0, 'D.Dinv', NOUT, work )
!stop
!! END TEST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! 6) write out Dinv
!  deallocate(work)
!  allocate(work(1:sizeD))
!  call PZLAPRNT(sizeD, sizeD, Dinv, 1, 1, DESC_A, 0, 0, 'D.Dinv', NOUT, work )

  do i=1,sizeD
    do j=1,sizeD
      call GLOBAL_TO_LOCAL(Pr,Pc,X,Y,I,J,MB,NB,NPROW,NPCOL)
      if( MYROW==0 .and. MYCOL==0 ) then
        if( Pr/=0 .or. Pc/=0 ) then 
          call ZGERV2D(ICTXT,1,1,ele,1,Pr,Pc)
        else
          ele=Dinv(X,Y)
        endif
        write(N_OUTFILE,'(E15.8,2X,E15.8,2X)',advance='no') ele
      else
        if( MYROW==Pr .and. MYCOL==Pc ) then
          ele=Dinv(X,Y)
          call ZGESD2D(ICTXT,1,1,ele,1,0,0)  
        endif
      endif
    enddo
  enddo
  if( MYROW==0 .and. MYCOL==0 ) then
    write(N_OUTFILE,*)
  endif
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
1111 control = 1
call IGEBS2D(ICTXT,'ALL','i-ring',1,1,control,1)
stop

end program calcDinv



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute local position from a global position
SUBROUTINE GLOBAL_TO_LOCAL(Pr,Pc,X,Y,I,J,MB,NB,NPROW,NPCOL)
implicit none

integer, INTENT(OUT) :: Pr,Pc,X,Y
integer, INTENT(IN) :: I,J,MB,NB,NPROW,NPCOL

integer R1, R2, tmp

Pr = mod( (I-1)/MB, NPROW )
Pc = mod( (J-1)/NB, NPCOL )

tmp = NPROW*MB
R1 = (I-1)/tmp

tmp = NPCOL*NB
R2 = (J-1)/tmp


X = R1 * MB + mod(I-1,MB) + 1
Y = R2 * NB + mod(J-1,NB) + 1

END SUBROUTINE GLOBAL_TO_LOCAL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute global position (I,J) from a local position (Pr,Pc,X,Y)
SUBROUTINE LOCAL_TO_GLOBAL(I,J,Pr,Pc,X,Y,MB,NB,NPROW,NPCOL)
implicit none

integer, INTENT(IN) :: Pr,Pc,X,Y,MB,NB,NPROW,NPCOL
integer, INTENT(OUT) :: I,J

integer r1,r2

r1 = (X-1)/MB
r2 = (Y-1)/NB

I = MB * NPROW * r1 + MB * Pr + mod(X-1,MB) + 1
J = NB * NPCOL * r2 + NB * Pc + mod(Y-1,NB) + 1

END SUBROUTINE LOCAL_TO_GLOBAL


!***********************************************************
!***********************************************************
! Box-Muller method for generating Gaussian random number
SUBROUTINE BoxMuller(p,q,seed)  

  implicit none 
  
  integer seed
  double precision p,q,r,s,Pi

  Pi=2d0*DASIN(1d0)
  call RandomNum(seed, r)
  call RandomNum(seed, s)

  p=dsqrt(-2d0*dlog(r))*DSIN(2d0*Pi*s)
  q=dsqrt(-2d0*dlog(r))*DCOS(2d0*Pi*s)

  return

END SUBROUTINE BoxMuller

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Random number generater
SUBROUTINE RandomNum(seed, r)
implicit none

integer seed
doubleprecision r
INTEGER, PARAMETER :: mask = 2147483647, a = 48828125

seed = IAND( a*seed, mask )
r = DBLE(seed)/DBLE(mask)
END SUBROUTINE RandomNum





