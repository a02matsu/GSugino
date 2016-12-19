!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! module on making quantities on face variables
module face_variables
use global_parameters
!use simplicial_complex
implicit none
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to make plaquette variable of face f
subroutine Make_face_variable(f,Uf)
use simplicial_complex, only : get_links_in_face_sc
implicit none

integer, intent(in) :: f
complex(kind(0d0)), intent(inout) :: Uf(1:NMAT,1:NMAT)
integer :: FaceSize
integer, allocatable :: sites(:),link_labels(:),link_dirs(:)
character(1) :: char1
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
integer :: l,i,j

call get_links_in_face_sc(sc,f,FaceSize,sites,link_labels,link_dirs)
!write(*,*) "### face ", f
!write(*,*) "sites:", sites
!write(*,*) "links:", link_labels
!write(*,*) "directions:", link_dirs

if(link_dirs(1)==1) then
  Uf=UMAT(:,:,link_labels(1))
elseif(link_dirs(1)==-1) then
  do i=1,NMAT
  do j=1,NMAT
    UF(i,j)=dconjg(UMAT(j,i,link_labels(1)))
  enddo
  enddo
else
  write(*,*) "There is a blanck link in the face",f
  stop 1
endif
!
do l=2,FaceSize
  tmpmat=Uf
  if(link_dirs(l)==1) then 
    char1='N'
  elseif(link_dirs(l)==-1) then
    char1='C'
  else
    write(*,*) "There is a blanck link in the face",f
    stop
  endif
  !! Uf=tmpmat*Uf or tmpmat*(Uf)^\dagger
  call ZGEMM('N',char1,NMAT,NMAT,NMAT,(1d0,0d0), &
    tmpmat,NMAT, &
    UMAT(:,:,link_labels(l)),NMAT, &
    (0d0,0d0), Uf, NMAT)
enddo
end subroutine Make_face_variable



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to make moment map \Omega(Uf)
subroutine moment_map(f,Omega)
use matrix_functions, only : matrix_power, matrix_inverse
implicit none

integer, intent(in) :: f
complex(kind(0d0)), intent(inout) :: Omega(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Ufm(1:NMAT,1:NMAT),Uf(1:NMAT,1:NMAT)
complex(kind(0d0)) :: SMAT(1:NMAT,1:NMAT)
complex(kind(0d0)) :: CMAT(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
integer :: i,j

call Make_face_variable(f,Uf)
call matrix_power(Ufm,Uf,mm)

do i=1,NMAT
  do j=1,NMAT
    SMAT(i,j)=-im_unit*( Ufm(i,j) - dconjg( Ufm(j,i) ) )
    CMAT(i,j)=( Ufm(i,j) + dconjg( Ufm(j,i) ) )
  enddo
enddo

! CMAT --> CMAT^{-1}
call matrix_inverse(CMAT)

call ZGEMM('N','N',NMAT,NMAT,NMAT,dcmplx(1d0/dble(mm)), &
  SMAT,NMAT,&
  CMAT,NMAT,&
  (0d0,0d0),tmpmat,NMAT)

do i=1,NMAT
  do j=1,NMAT
    Omega(i,j)=tmpmat(i,j)+dconjg(tmpmat(j,i))
  enddo
enddo

end subroutine moment_map

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to compute \frac{d}{d A_l^a} U_f^m
subroutine Make_diff_Ufm(diff_Ufm, f, l)
!use matrix_functions, only : matrix_inverse
implicit none

complex(kind(0d0)), intent(inout) :: diff_Ufm(1:NMAT,1:NMAT)
integer, intent(in) :: f, l
integer :: FaceSize
integer, allocatable :: sites(:),link_labels(:),link_dirs(:)
integer :: nl, info
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)) :: XMAT(1:NMAT,1:NMAT), YMAT(1:NMAT,1:NMAT), tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Ufla(1:NMAT,1:NMAT)
integer :: i,j

! make face variable
call Make_face_varible(f,Uf)
! get information of the face
call get_links_in_face_sc(sc,f,FaceSize,sites,link_labels,link_dirs)
! find the place of link l in the face l
info=0
do nl=1,FaceSize
  if( link_labels(i) == l ) then 
    info=1
    exit
  endif
enddo
if ( info == 0 ) then 
  write(*,*) "link",l,"is not included in the face",f
  stop
endif

!! when Ul is the first link
if( nl == 1 ) then 
!!!
! X=1, Y=U_l^{-1} U_f ( dir=1 )
  if( link_dirs(nl) == 1 ) then
    call TtimesM(Ufla, Uf, a, NMAT)
    Ufla = (0d0,1d0)*Ufla
!!!
! X=1, Y=U_l U_f      ( dir=-1 )
  elseif( link_dirs(nl) == -1 ) then
    do i=1,NMAT
      do j=1,NMAT
        ! X \equiv U_l^{-1}
        XMAT(i,j)=dconjg( UMAT(j,i,l) )
      enddo
    enddo
    ! Y=U_l.U_f
    call ZGEMM('N','N',NMAT,NMAT,NMAT,(1d0,0d0), &
      UMAT(:,:,l),NMAT, &
      Uf, NMAT, &
      (0d0,0d0), YMAT, NMAT)
    ! Uf_l^a = -i Ul^{-1}.T_a.Y
    call MTN(Ufla,XMAT,YMAT,a,NMAT)
    Ufla = (0d0,-1d0) * Ufla
  endif
! when Ul is the last link
elseif( nl == FaceSize ) then
  if( link_dir(nl) == 1) then
!  X=U_f U_l^{-1}, Y=1 (dir=1)
    call ZGEMM('N','C',NMAT,NMAT,NMAT,(1d0,0d0), &
      Uf,NMAT, &
      UMAT(:,:,l), NMAT, &
      (0d0,0d0), XMAT, NMAT)
    call MTN(Ufla, XMAT, UMAT(:,:,l), NMAT)
    Ufla = (0d0,1d0) * Ufla
  elseif( link_dir(nl) == -1 ) then
!  X=U_f U_l,      Y=1 (dir=-1)
    call MtimesT(Ulfa,Uf,a,NMAT)
    Ulfa=Ulfa*(0d0,-1d0)
  endif
!!!
else
! X=U_1 ... U_{nl-1}
  if(link_dirs(1)==1) then
    XMAT=UMAT(:,:,link_labels(1))
  elseif(link_dirs(1)==-1) then
    do i=1,NMAT
    do j=1,NMAT
      XMAT(i,j)=dconjg(UMAT(j,i,link_labels(1)))
    enddo
    enddo
  endif
  !
  do l=2,nl-1
    tmpmat=XMAT
    if(link_dirs(l)==1) then 
      char1='N'
    elseif(link_dirs(l)==-1) then
      char1='C'
    endif
    !! Uf=tmpmat*Uf or tmpmat*(Uf)^\dagger
    call ZGEMM('N',char1,NMAT,NMAT,NMAT,(1d0,0d0), &
      tmpmat,NMAT, &
      UMAT(:,:,link_labels(l)),NMAT, &
      (0d0,0d0), XMAT, NMAT)
  enddo
! Y=U_{nl+1}...U_{FaceSize}
  if(link_dirs(nl+1)==1) then
    YMAT=UMAT(:,:,link_labels(nl+1))
  elseif(link_dirs(nl+1)==-1) then
    do i=1,NMAT
    do j=1,NMAT
      YMAT(i,j)=dconjg(UMAT(j,i,link_labels(nl+1)))
    enddo
    enddo
  endif
  !
  do l=nl+1,FaceSize
    tmpmat=YMAT
    if(link_dirs(l)==1) then 
      char1='N'
    elseif(link_dirs(l)==-1) then
      char1='C'
    endif
    !! Uf=tmpmat*Uf or tmpmat*(Uf)^\dagger
    call ZGEMM('N',char1,NMAT,NMAT,NMAT,(1d0,0d0), &
      tmpmat,NMAT, &
      UMAT(:,:,link_labels(l)),NMAT, &
      (0d0,0d0), YMAT, NMAT)
  enddo

  if( link_dirs(nl) == 1 ) then
    call ZGEMM('N','N',NMAT,NMAT,NMAT,(0d0,1d0), &
      UMAT(:,:,l),NMAT, &
      YMAT,NMAT, &
      (0d0,0d0), tmpmat, NMAT)
    call MTN(Ulfa,XMAT,tmpmat,a,NMAT)
  elseif( link_dirs(nl) == -1 ) then
    call ZGEMM('N','C',NMAT,NMAT,NMAT,(0d0,-1d0), &
      XMAT,NMAT, &
      UMAT(:,:,l),NMAT, &
      (0d0,0d0), tmpmat, NMAT)
    call MTN(Ulfa,tmpmat,YMAT,a,NMAT)
  endif







end subroutine Make_diff_Ufm


end module face_variables
