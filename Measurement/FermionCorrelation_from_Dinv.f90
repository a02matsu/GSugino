!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
subroutine make_fermion_correlation_from_Dinv(&
    Geta_eta, Glambda_eta, Gchi_eta, &
    Geta_lambda, Glambda_lambda, Gchi_lambda, &
    Geta_chi, Glambda_chi, Gchi_chi, &
    Dinv)
use global_parameters
use SUN_generators, only : make_SUN_generators
use parallel
implicit none

complex(kind(0d0)), intent(out) :: Geta_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_sites) 
complex(kind(0d0)), intent(out) :: Glambda_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_sites) 
complex(kind(0d0)), intent(out) :: Gchi_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_sites) 
!!!
complex(kind(0d0)), intent(out) :: Geta_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_links) 
complex(kind(0d0)), intent(out) :: Glambda_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_links) 
complex(kind(0d0)), intent(out) :: Gchi_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_links) 
!!!
complex(kind(0d0)), intent(out) :: Geta_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_faces) 
complex(kind(0d0)), intent(out) :: Glambda_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_faces) 
complex(kind(0d0)), intent(out) :: Gchi_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_faces) 

complex(kind(0d0)), intent(in) :: Dinv(:,:)
integer :: numF

complex(kind(0d0)), allocatable :: mode_Geta_eta(:,:,:,:)
complex(kind(0d0)), allocatable :: mode_Geta_lambda(:,:,:,:)
complex(kind(0d0)), allocatable :: mode_Geta_chi(:,:,:,:)
complex(kind(0d0)), allocatable :: mode_Glambda_eta(:,:,:,:)
complex(kind(0d0)), allocatable :: mode_Glambda_lambda(:,:,:,:)
complex(kind(0d0)), allocatable :: mode_Glambda_chi(:,:,:,:)
complex(kind(0d0)), allocatable :: mode_Gchi_eta(:,:,:,:)
complex(kind(0d0)), allocatable :: mode_Gchi_lambda(:,:,:,:)
complex(kind(0d0)), allocatable :: mode_Gchi_chi(:,:,:,:)
complex(kind(0d0)), allocatable :: Dinv_check(:,:)

complex(kind(0d0)) :: T(1:NMAT,1:NMAT,1:dimG) 
complex(kind(0d0)) :: TT(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:dimG,1:dimG) 

complex(kind(0d0)) :: tmp
integer :: gs,gl,gf
integer :: gs2,gl2,gf2,DL1,DL2
integer :: ls,ll,lf
integer :: ls2,ll2,lf2,rank2
integer :: a,b,i,j,k,l
integer :: Dlabel
integer :: rank,tag

numF=(global_num_sites+global_num_links+global_num_faces)*(NMAT*NMAT-1)
allocate( mode_Geta_eta(1:dimG,1:dimG,1:global_num_sites,1:num_sites) )
allocate( mode_Geta_lambda(1:dimG,1:dimG,1:global_num_sites,1:num_links) )
allocate( mode_Geta_chi(1:dimG,1:dimG,1:global_num_sites,1:num_faces) )
allocate( mode_Glambda_eta(1:dimG,1:dimG,1:global_num_links,1:num_sites) )
allocate( mode_Glambda_lambda(1:dimG,1:dimG,1:global_num_links,1:num_links) )
allocate( mode_Glambda_chi(1:dimG,1:dimG,1:global_num_links,1:num_faces) )
allocate( mode_Gchi_eta(1:dimG,1:dimG,1:global_num_faces,1:num_sites) )
allocate( mode_Gchi_lambda(1:dimG,1:dimG,1:global_num_faces,1:num_links) )
allocate( mode_Gchi_chi(1:dimG,1:dimG,1:global_num_faces,1:num_faces) )
!allocate( mode_mat_eta(1:numF,1:NMAT,1:NMAT,1:num_sites) )
!allocate( mode_mat_lambda(1:numF,1:NMAT,1:NMAT,1:num_links) )
!allocate( mode_mat_chi(1:numF,1:NMAT,1:NMAT,1:num_faces) )
allocate( Dinv_check(1:numF,1:numF) )

if( size(Dinv,1) .ne. numF) then
  write(*,*) "size of Dinv is incorrect."
  call stop_for_test
endif

Geta_eta=(0d0,0d0)
Glambda_eta=(0d0,0d0)
Gchi_eta=(0d0,0d0)
Geta_lambda=(0d0,0d0)
Glambda_lambda=(0d0,0d0)
Gchi_lambda=(0d0,0d0)
Geta_chi=(0d0,0d0)
Glambda_chi=(0d0,0d0)
Gchi_chi=(0d0,0d0)

! SU(N) generators
call make_SUN_generators(T,NMAT)
do b=1,dimG
  do a=1,dimG
    do l=1,NMAT
      do k=1,NMAT
        do j=1,NMAT
          do i=1,NMAT
            TT(i,j,k,l,a,b)=T(i,j,a)*T(k,l,b)
          enddo
        enddo
      enddo
    enddo
  enddo
enddo


!! Distribute Dinv to MPI_COMM_WORLD
DL2=0
tag=0
do gs=1,global_num_sites
  ls=local_site_of_global(gs)%label_
  rank=local_site_of_global(gs)%rank_ 
  do a=1,dimG
    DL2=DL2+1
    !!!!!!!!!!!
    DL1=0
    do gs2=1,global_num_sites
      do b=1,dimG
        DL1=DL1+1
        tag=tag+1
        !call set_mode(mode_Geta_eta(b,a,gs2,ls),Dinv(DL1,DL2),rank,tag)
        call set_mode(tmp,Dinv(DL1,DL2),rank,tag)
        if( MYRANK==rank ) mode_Geta_eta(b,a,gs2,ls)=tmp
      enddo
    enddo
    do gl2=1,global_num_links
      do b=1,dimG
        DL1=DL1+1
        tag=tag+1
        !call set_mode(mode_Glambda_eta(b,a,gl2,ls),Dinv(DL1,DL2),rank,tag)
        call set_mode(tmp,Dinv(DL1,DL2),rank,tag)
        if( MYRANK==rank ) mode_Glambda_eta(b,a,gl2,ls)=tmp
      enddo
    enddo
    do gf2=1,global_num_faces
      do b=1,dimG
        DL1=DL1+1
        tag=tag+1
        !call set_mode(mode_Gchi_eta(b,a,gf2,ls),Dinv(DL1,DL2),rank,tag)
        call set_mode(tmp,Dinv(DL1,DL2),rank,tag)
        if( MYRANK==rank ) mode_Gchi_eta(b,a,gf2,ls)=tmp
      enddo
    enddo
  enddo
enddo
do gl=1,global_num_links
  ll=local_link_of_global(gl)%label_
  rank=local_link_of_global(gl)%rank_ 
  do a=1,dimG
    DL2=DL2+1
    !!!!!!!!!!!
    DL1=0
    do gs2=1,global_num_sites
      do b=1,dimG
        DL1=DL1+1
        tag=tag+1
        !call set_mode(mode_Geta_lambda(b,a,gs2,ll),Dinv(DL1,DL2),rank,tag)
        call set_mode(tmp,Dinv(DL1,DL2),rank,tag)
        if( MYRANK==rank ) mode_Geta_lambda(b,a,gs2,ll)=tmp
      enddo
    enddo
    do gl2=1,global_num_links
      do b=1,dimG
        DL1=DL1+1
        tag=tag+1
        !call set_mode(mode_Glambda_lambda(b,a,gl2,ll),Dinv(DL1,DL2),rank,tag)
        call set_mode(tmp,Dinv(DL1,DL2),rank,tag)
        if(MYRANK==rank) mode_Glambda_lambda(b,a,gl2,ll)=tmp
      enddo
    enddo
    do gf2=1,global_num_faces
      do b=1,dimG
        DL1=DL1+1
        tag=tag+1
        !call set_mode(mode_Gchi_lambda(b,a,gf2,ll),Dinv(DL1,DL2),rank,tag)
        call set_mode(tmp,Dinv(DL1,DL2),rank,tag)
        if( MYRANK==rank ) mode_Gchi_lambda(b,a,gf2,ll)=tmp
      enddo
    enddo
  enddo
enddo
do gf=1,global_num_faces
  lf=local_face_of_global(gf)%label_
  rank=local_face_of_global(gf)%rank_ 
  do a=1,dimG
    DL2=DL2+1
    !!!!!!!!!!!
    DL1=0
    do gs2=1,global_num_sites
      do b=1,dimG
        DL1=DL1+1
        tag=tag+1
        !call set_mode(mode_Geta_chi(b,a,gs2,lf),Dinv(DL1,DL2),rank,tag)
        call set_mode(tmp,Dinv(DL1,DL2),rank,tag)
        if( MYRANK==rank ) mode_Geta_chi(b,a,gs2,lf)=tmp
      enddo
    enddo
    do gl2=1,global_num_links
      do b=1,dimG
        DL1=DL1+1
        tag=tag+1
        !call set_mode(mode_Glambda_chi(b,a,gl2,lf),Dinv(DL1,DL2),rank,tag)
        call set_mode(tmp,Dinv(DL1,DL2),rank,tag)
        if( MYRANK==rank ) mode_Glambda_chi(b,a,gl2,lf)=tmp
      enddo
    enddo
    do gf2=1,global_num_faces
      do b=1,dimG
        DL1=DL1+1
        tag=tag+1
        !call set_mode(mode_Gchi_chi(b,a,gf2,lf),Dinv(DL1,DL2),rank,tag)
        call set_mode(tmp,Dinv(DL1,DL2),rank,tag)
        if( MYRANK==rank ) mode_Gchi_chi(b,a,gf2,lf)=tmp
      enddo
    enddo
  enddo
enddo


!! convert modes to matrices
do ls=1,num_sites
  do gs=1,global_num_sites
    do a=1,dimG
      do b=1,dimG
        Geta_eta(:,:,:,:,gs,ls) = Geta_eta(:,:,:,:,gs,ls) &
          +TT(:,:,:,:,a,b)*mode_Geta_eta(a,b,gs,ls)
      enddo
    enddo
  enddo
  do gl=1,global_num_links
    do a=1,dimG
      do b=1,dimG
        Glambda_eta(:,:,:,:,gl,ls)=Glambda_eta(:,:,:,:,gl,ls)&
          +TT(:,:,:,:,a,b)*mode_Glambda_eta(a,b,gl,ls)
      enddo
    enddo
  enddo
  do gf=1,global_num_faces
    do a=1,dimG
      do b=1,dimG
        Gchi_eta(:,:,:,:,gf,ls)=Gchi_eta(:,:,:,:,gf,ls) &
          +TT(:,:,:,:,a,b)*mode_Gchi_eta(a,b,gf,ls)
      enddo
    enddo
  enddo
enddo
!!
do ll=1,num_links
  do gs=1,global_num_sites
    do a=1,dimG
      do b=1,dimG
        Geta_lambda(:,:,:,:,gs,ll)= Geta_lambda(:,:,:,:,gs,ll) &
          +TT(:,:,:,:,a,b)*mode_Geta_lambda(a,b,gs,ll)
      enddo
    enddo
  enddo
  do gl=1,global_num_links
    do a=1,dimG
      do b=1,dimG
        Glambda_lambda(:,:,:,:,gl,ll)= Glambda_lambda(:,:,:,:,gl,ll) &
          +TT(:,:,:,:,a,b)*mode_Glambda_lambda(a,b,gl,ll)
      enddo
    enddo
  enddo
  do gf=1,global_num_faces
    do a=1,dimG
      do b=1,dimG
        Gchi_lambda(:,:,:,:,gf,ll)= Gchi_lambda(:,:,:,:,gf,ll) &
          +TT(:,:,:,:,a,b)*mode_Gchi_lambda(a,b,gf,ll)
      enddo
    enddo
  enddo
enddo
!!
do lf=1,num_faces
  do gs=1,global_num_sites
    do a=1,dimG
      do b=1,dimG
        Geta_chi(:,:,:,:,gs,lf)= Geta_chi(:,:,:,:,gs,lf) &
          +TT(:,:,:,:,a,b)*mode_Geta_chi(a,b,gs,lf)
      enddo
    enddo
  enddo
  do gl=1,global_num_links
    do a=1,dimG
      do b=1,dimG
        Glambda_chi(:,:,:,:,gl,lf)= Glambda_chi(:,:,:,:,gl,lf) &
          +TT(:,:,:,:,a,b)*mode_Glambda_chi(a,b,gl,lf)
      enddo
    enddo
  enddo
  do gf=1,global_num_faces
    do a=1,dimG
      do b=1,dimG
        Gchi_chi(:,:,:,:,gf,lf)= Gchi_chi(:,:,:,:,gf,lf) &
          +TT(:,:,:,:,a,b)*mode_Gchi_chi(a,b,gf,lf)
      enddo
    enddo
  enddo
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! check if mode_*** is the same with Dinv !!!!!
!call check_mode_Dinv(&
!    Dinv, &
!    mode_Geta_eta, mode_Glambda_eta, mode_Gchi_eta, &
!    mode_Geta_lambda, mode_Glambda_lambda, mode_Gchi_lambda, &
!    mode_Geta_chi, mode_Glambda_chi, mode_Gchi_chi)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! check if Gany_any is the same with Dinv !!!!!
!call correlation_to_Dinv(&
!    Dinv, &
!    Geta_eta, Glambda_eta, Gchi_eta, &
!    Geta_lambda, Glambda_lambda, Gchi_lambda, &
!    Geta_chi, Glambda_chi, Gchi_chi,T)


contains
  subroutine set_mode(mode,element,rank,tag)
  implicit none

  complex(kind(0d0)), intent(out) :: mode
  complex(kind(0d0)), intent(in) :: element
  integer, intent(in) :: rank,tag

  if( rank==0 ) then
    if( MYRANK==0 ) then
      mode=Dinv(DL1,DL2)
    endif
  else
    if( MYRANK==0 ) then
      call MPI_SEND(element, &
        1,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,IERR)
    elseif( MYRANK==rank ) then
      call MPI_RECV(mode, &
        1,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,ISTATUS,IERR)
    endif
  endif
  end subroutine set_mode

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine pull_mode(mode,element,rank,tag)
  implicit none

  complex(kind(0d0)), intent(out) :: mode
  complex(kind(0d0)), intent(in) :: element
  integer, intent(in) :: rank,tag

  if( rank==0 ) then
    if( MYRANK==0 ) then
      mode=element
    endif
  else
    if( MYRANK==rank ) then
      call MPI_SEND(element, &
        1,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,IERR)
    elseif( MYRANK==0 ) then
      call MPI_RECV(mode, &
        1,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,ISTATUS,IERR)
    endif
  endif
  end subroutine pull_mode

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine check_mode_Dinv(&
    Dinv, &
    mode_Geta_eta, mode_Glambda_eta, mode_Gchi_eta, &
    mode_Geta_lambda, mode_Glambda_lambda, mode_Gchi_lambda, &
    mode_Geta_chi, mode_Glambda_chi, mode_Gchi_chi)
  implicit none

  complex(kind(0d0)), intent(in) :: Dinv(:,:)
  complex(kind(0d0)), intent(in) :: mode_Geta_eta(1:dimG,1:dimG,1:global_num_sites,1:num_sites) 
  complex(kind(0d0)), intent(in) :: mode_Glambda_eta(1:dimG,1:dimG,1:global_num_links,1:num_sites) 
  complex(kind(0d0)), intent(in) :: mode_Gchi_eta(1:dimG,1:dimG,1:global_num_faces,1:num_sites) 
  !!!
  complex(kind(0d0)), intent(in) :: mode_Geta_lambda(1:dimG,1:dimG,1:global_num_sites,1:num_links) 
  complex(kind(0d0)), intent(in) :: mode_Glambda_lambda(1:dimG,1:dimG,1:global_num_links,1:num_links) 
  complex(kind(0d0)), intent(in) :: mode_Gchi_lambda(1:dimG,1:dimG,1:global_num_faces,1:num_links) 
  !!!
  complex(kind(0d0)), intent(in) :: mode_Geta_chi(1:dimG,1:dimG,1:global_num_sites,1:num_faces) 
  complex(kind(0d0)), intent(in) :: mode_Glambda_chi(1:dimG,1:dimG,1:global_num_links,1:num_faces) 
  complex(kind(0d0)), intent(in) :: mode_Gchi_chi(1:dimG,1:dimG,1:global_num_faces,1:num_faces) 
  
  integer :: DL1,DL2,tag,rank
  complex(kind(0d0)) :: tmp

  if(MYRANK==0) then
    write(*,*) "==== check Dinv-mode ===="
  endif

  DL2=0
  tag=0
  do gs=1,global_num_sites
    ls=local_site_of_global(gs)%label_
    rank=local_site_of_global(gs)%rank_ 
    do a=1,dimG
      DL2=DL2+1
      !!!!!!!!!!!
      DL1=0
      do gs2=1,global_num_sites
        do b=1,dimG
          DL1=DL1+1
          tag=tag+1
          call pull_mode(tmp,mode_Geta_eta(b,a,gs2,ls),rank,tag)
          if(MYRANK==0) then
            if( cdabs(tmp-Dinv(DL1,DL2))>1d-3 ) then 
              write(*,*) DL1,DL2,tmp,Dinv(DL1,DL2)
            endif
          endif
        enddo
      enddo
      do gl2=1,global_num_links
        do b=1,dimG
          DL1=DL1+1
          tag=tag+1
          call pull_mode(tmp,mode_Glambda_eta(b,a,gl2,ls),rank,tag)
          if(MYRANK==0) then
            if( cdabs(tmp-Dinv(DL1,DL2))>1d-3 ) then 
              write(*,*) DL1,DL2,tmp,Dinv(DL1,DL2)
            endif
          endif
        enddo
      enddo
      do gf2=1,global_num_faces
        do b=1,dimG
          DL1=DL1+1
          tag=tag+1
          call pull_mode(tmp,mode_Gchi_eta(b,a,gf2,ls),rank,tag)
          if(MYRANK==0) then
            if( cdabs(tmp-Dinv(DL1,DL2))>1d-3 ) then 
              write(*,*) DL1,DL2,tmp,Dinv(DL1,DL2)
            endif
          endif
        enddo
      enddo
    enddo
  enddo
  do gl=1,global_num_links
    ll=local_link_of_global(gl)%label_
    rank=local_link_of_global(gl)%rank_ 
    do a=1,dimG
      DL2=DL2+1
      !!!!!!!!!!!
      DL1=0
      do gs2=1,global_num_sites
        do b=1,dimG
          DL1=DL1+1
          tag=tag+1
          call pull_mode(tmp,mode_Geta_lambda(b,a,gs2,ll),rank,tag)
          if(MYRANK==0) then
            if( cdabs(tmp-Dinv(DL1,DL2))>1d-3 ) then 
              write(*,*) DL1,DL2,tmp,Dinv(DL1,DL2)
            endif
          endif
        enddo
      enddo
      do gl2=1,global_num_links
        do b=1,dimG
          DL1=DL1+1
          tag=tag+1
          call pull_mode(tmp,mode_Glambda_lambda(b,a,gl2,ll),rank,tag)
          if(MYRANK==0) then
            if( cdabs(tmp-Dinv(DL1,DL2))>1d-3 ) then 
              write(*,*) DL1,DL2,tmp,Dinv(DL1,DL2)
            endif
          endif
        enddo
      enddo
      do gf2=1,global_num_faces
        do b=1,dimG
          DL1=DL1+1
          tag=tag+1
          call pull_mode(tmp,mode_Gchi_lambda(b,a,gf2,ll),rank,tag)
          if(MYRANK==0) then
            if( cdabs(tmp-Dinv(DL1,DL2))>1d-3 ) then 
              write(*,*) DL1,DL2,tmp,Dinv(DL1,DL2)
            endif
          endif
        enddo
      enddo
    enddo
  enddo
  do gf=1,global_num_faces
    lf=local_face_of_global(gf)%label_
    rank=local_face_of_global(gf)%rank_ 
    do a=1,dimG
      DL2=DL2+1
      !!!!!!!!!!!
      DL1=0
      do gs2=1,global_num_sites
        do b=1,dimG
          DL1=DL1+1
          tag=tag+1
          call pull_mode(tmp,mode_Geta_chi(b,a,gs2,lf),rank,tag)
          if(MYRANK==0) then
            if( cdabs(tmp-Dinv(DL1,DL2))>1d-3 ) then 
              write(*,*) DL1,DL2,tmp,Dinv(DL1,DL2)
            endif
          endif
        enddo
      enddo
      do gl2=1,global_num_links
        do b=1,dimG
          DL1=DL1+1
          tag=tag+1
          call pull_mode(tmp,mode_Glambda_chi(b,a,gl2,lf),rank,tag)
          if(MYRANK==0) then
            if( cdabs(tmp-Dinv(DL1,DL2))>1d-3 ) then 
              write(*,*) DL1,DL2,tmp,Dinv(DL1,DL2)
            endif
          endif
        enddo
      enddo
      do gf2=1,global_num_faces
        do b=1,dimG
          DL1=DL1+1
          tag=tag+1
          call pull_mode(tmp,mode_Gchi_chi(b,a,gf2,lf),rank,tag)
          if(MYRANK==0) then
            if( cdabs(tmp-Dinv(DL1,DL2))>1d-3 ) then 
              write(*,*) DL1,DL2,tmp,Dinv(DL1,DL2)
            endif
          endif
        enddo
      enddo
    enddo
  enddo

  if(MYRANK==0) then
    write(*,*) "==== END check Dinv-mode ===="
  endif
  end subroutine check_mode_Dinv


end subroutine make_fermion_correlation_from_Dinv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
subroutine make_Dinv_matrix(&
    Geta_eta, Glambda_eta, Gchi_eta, &
    Geta_lambda, Glambda_lambda, Gchi_lambda, &
    Geta_chi, Glambda_chi, Gchi_chi)
use global_parameters
use SUN_generators, only : make_SUN_generators
use parallel
implicit none

complex(kind(0d0)), intent(in) :: Geta_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_sites) 
complex(kind(0d0)), intent(in) :: Glambda_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_sites) 
complex(kind(0d0)), intent(in) :: Gchi_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_sites) 
!!!
complex(kind(0d0)), intent(in) :: Geta_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_links) 
complex(kind(0d0)), intent(in) :: Glambda_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_links) 
complex(kind(0d0)), intent(in) :: Gchi_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_links) 
!!!
complex(kind(0d0)), intent(in) :: Geta_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_faces) 
complex(kind(0d0)), intent(in) :: Glambda_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_faces) 
complex(kind(0d0)), intent(in) :: Gchi_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_faces) 

complex(kind(0d0)), allocatable :: Dinv(:,:)
integer :: numFmat

integer :: gs1,gl1,gf1,gs2,gl2,gf2
integer :: ls2,ll2,lf2
integer :: i,j,k,l
integer :: DL1,DL2
integer :: rank,tag


numFmat= (global_num_sites+global_num_links+global_num_faces)*(NMAT*NMAT)
if(MYRANK==0) allocate(Dinv(1:numFmat,1:numFmat))

tag=0
DL1=0
do gs1=1,global_num_sites; do i=1,NMAT; do j=1,NMAT
  DL1=DL1+1
  DL2=0
  do gs2=1,global_num_sites
    ls2=local_site_of_global(gs2)%label_
    rank=local_site_of_global(gs2)%rank_ 
    do k=1,NMAT; do l=1,NMAT
      DL2=DL2+1
      if( rank==0 ) then 
        if( MYRANK==0 ) then
          Dinv(DL1,DL2)=Geta_eta(i,j,k,l,gs1,ls2)
        endif
      else
        tag=tag+1
        if( MYRANK==rank ) then
          call MPI_SEND(Geta_eta(i,j,k,l,gs1,ls2), &
            1,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,IERR)
        elseif( MYRANK==0 ) then
          call MPI_RECV(Dinv(DL1,DL2), &
            1,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,ISTATUS,IERR)
        endif
      endif
    enddo;enddo;
  enddo
      
  !!
  do gl2=1,global_num_links
    ll2=local_link_of_global(gl2)%label_
    rank=local_link_of_global(gl2)%rank_ 
    do k=1,NMAT; do l=1,NMAT
      DL2=DL2+1
      if( rank==0 ) then 
        if( MYRANK==0 ) then
          Dinv(DL1,DL2)=Geta_lambda(i,j,k,l,gs1,ll2)
        endif
      else
        tag=tag+1
        if( MYRANK==rank ) then
          call MPI_SEND(Geta_lambda(i,j,k,l,gs1,ll2), &
            1,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,IERR)
        elseif( MYRANK==0 ) then
          call MPI_RECV(Dinv(DL1,DL2), &
            1,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,ISTATUS,IERR)
        endif
      endif
    enddo;enddo;
  enddo
  !!
  do gf2=1,global_num_faces
    lf2=local_face_of_global(gf2)%label_
    rank=local_face_of_global(gf2)%rank_ 
    do k=1,NMAT; do l=1,NMAT
      DL2=DL2+1
      if( rank==0 ) then 
        if( MYRANK==0 ) then
          Dinv(DL1,DL2)=Geta_chi(i,j,k,l,gs1,lf2)
        endif
      else
        tag=tag+1
        if( MYRANK==rank ) then
          call MPI_SEND(Geta_chi(i,j,k,l,gs1,lf2), &
            1,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,IERR)
        elseif( MYRANK==0 ) then
          call MPI_RECV(Dinv(DL1,DL2), &
            1,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,ISTATUS,IERR)
        endif
      endif
    enddo;enddo;
  enddo
enddo;enddo;enddo
!!!!!!!!!!!!
do gl1=1,global_num_links; do i=1,NMAT; do j=1,NMAT
  DL1=DL1+1
  DL2=0
  do gs2=1,global_num_sites
    ls2=local_site_of_global(gs2)%label_
    rank=local_site_of_global(gs2)%rank_ 
    do k=1,NMAT; do l=1,NMAT
      DL2=DL2+1
      if( rank==0 ) then 
        if( MYRANK==0 ) then
          Dinv(DL1,DL2)=Glambda_eta(i,j,k,l,gl1,ls2)
        endif
      else
        tag=tag+1
        if( MYRANK==rank ) then
          call MPI_SEND(Glambda_eta(i,j,k,l,gl1,ls2), &
            1,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,IERR)
        elseif( MYRANK==0 ) then
          call MPI_RECV(Dinv(DL1,DL2), &
            1,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,ISTATUS,IERR)
        endif
      endif
    enddo;enddo;
  enddo
  !!
  do gl2=1,global_num_links
    ll2=local_link_of_global(gl2)%label_
    rank=local_link_of_global(gl2)%rank_ 
    do k=1,NMAT; do l=1,NMAT
      DL2=DL2+1
      if( rank==0 ) then 
        if( MYRANK==0 ) then
          Dinv(DL1,DL2)=Glambda_lambda(i,j,k,l,gl1,ll2)
        endif
      else
        tag=tag+1
        if( MYRANK==rank ) then
          call MPI_SEND(Glambda_lambda(i,j,k,l,gl1,ll2), &
            1,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,IERR)
        elseif( MYRANK==0 ) then
          call MPI_RECV(Dinv(DL1,DL2), &
            1,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,ISTATUS,IERR)
        endif
      endif
    enddo;enddo;
  enddo
  !!
  do gf2=1,global_num_faces
    lf2=local_face_of_global(gf2)%label_
    rank=local_face_of_global(gf2)%rank_ 
    do k=1,NMAT; do l=1,NMAT
      DL2=DL2+1
      if( rank==0 ) then 
        if( MYRANK==0 ) then
          Dinv(DL1,DL2)=Glambda_chi(i,j,k,l,gl1,lf2)
        endif
      else
        tag=tag+1
        if( MYRANK==rank ) then
          call MPI_SEND(Glambda_chi(i,j,k,l,gl1,lf2), &
            1,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,IERR)
        elseif( MYRANK==0 ) then
          call MPI_RECV(Dinv(DL1,DL2), &
            1,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,ISTATUS,IERR)
        endif
      endif
    enddo;enddo;
  enddo
enddo;enddo;enddo
!!!!!!!!!!!!
do gf1=1,global_num_faces; do i=1,NMAT; do j=1,NMAT
  DL1=DL1+1
  DL2=0
  do gs2=1,global_num_sites
    ls2=local_site_of_global(gs2)%label_
    rank=local_site_of_global(gs2)%rank_ 
    do k=1,NMAT; do l=1,NMAT
      DL2=DL2+1
      if( rank==0 ) then 
        if( MYRANK==0 ) then
          Dinv(DL1,DL2)=Gchi_eta(i,j,k,l,gf1,ls2)
        endif
      else
        tag=tag+1
        if( MYRANK==rank ) then
          call MPI_SEND(Gchi_eta(i,j,k,l,gf1,ls2), &
            1,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,IERR)
        elseif( MYRANK==0 ) then
          call MPI_RECV(Dinv(DL1,DL2), &
            1,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,ISTATUS,IERR)
        endif
      endif
    enddo;enddo;
  enddo
  !!
  do gl2=1,global_num_links
    ll2=local_link_of_global(gl2)%label_
    rank=local_link_of_global(gl2)%rank_ 
    do k=1,NMAT; do l=1,NMAT
      DL2=DL2+1
      if( rank==0 ) then 
        if( MYRANK==0 ) then
          Dinv(DL1,DL2)=Gchi_lambda(i,j,k,l,gf1,ll2)
        endif
      else
        tag=tag+1
        if( MYRANK==rank ) then
          call MPI_SEND(Gchi_lambda(i,j,k,l,gf1,ll2), &
            1,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,IERR)
        elseif( MYRANK==0 ) then
          call MPI_RECV(Dinv(DL1,DL2), &
            1,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,ISTATUS,IERR)
        endif
      endif
    enddo;enddo;
  enddo
  !!
  do gf2=1,global_num_faces
    lf2=local_face_of_global(gf2)%label_
    rank=local_face_of_global(gf2)%rank_ 
    do k=1,NMAT; do l=1,NMAT
      DL2=DL2+1
      if( rank==0 ) then 
        if( MYRANK==0 ) then
          Dinv(DL1,DL2)=Gchi_chi(i,j,k,l,gf1,lf2)
        endif
      else
        tag=tag+1
        if( MYRANK==rank ) then
          call MPI_SEND(Gchi_chi(i,j,k,l,gf1,lf2), &
            1,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,IERR)
        elseif( MYRANK==0 ) then
          call MPI_RECV(Dinv(DL1,DL2), &
            1,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,ISTATUS,IERR)
        endif
      endif
    enddo;enddo;
  enddo
enddo;enddo;enddo
  
if( MYRANK==0 ) then
  write(*,*) "======================="
  do i=1,numFmat
    do j=1,numFmat
      if( cdabs(Dinv(i,j)+Dinv(j,i))>1d-3) write(*,*) i,j,Dinv(i,j), Dinv(j,i)
    enddo
  enddo
  write(*,*) "======================="
endif


end subroutine make_Dinv_matrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
subroutine correlation_to_Dinv(&
    Dinv, &
    Geta_eta, Glambda_eta, Gchi_eta, &
    Geta_lambda, Glambda_lambda, Gchi_lambda, &
    Geta_chi, Glambda_chi, Gchi_chi,T)
use global_parameters
use SUN_generators, only : make_SUN_generators
use parallel
implicit none

complex(kind(0d0)), intent(in) :: Dinv(:,:)
complex(kind(0d0)), intent(in) :: Geta_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_sites) 
complex(kind(0d0)), intent(in) :: Glambda_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_sites) 
complex(kind(0d0)), intent(in) :: Gchi_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_sites) 
!!!
complex(kind(0d0)), intent(in) :: Geta_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_links) 
complex(kind(0d0)), intent(in) :: Glambda_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_links) 
complex(kind(0d0)), intent(in) :: Gchi_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_links) 
!!!
complex(kind(0d0)), intent(in) :: Geta_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_faces) 
complex(kind(0d0)), intent(in) :: Glambda_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_faces) 
complex(kind(0d0)), intent(in) :: Gchi_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_faces) 
complex(kind(0d0)), intent(in) :: T(1:NMAT,1:NMAT,1:dimG) 

integer :: numF

complex(kind(0d0)), allocatable :: mode_eta(:,:,:)
complex(kind(0d0)), allocatable :: mode_lambda(:,:,:)
complex(kind(0d0)), allocatable :: mode_chi(:,:,:)
complex(kind(0d0)), allocatable :: Dinv_check(:,:)




integer :: gs,gl,gf,gs1,gl1,gf1
integer :: ls,ll,lf
integer :: a,b,i,j,k,l
integer :: Dlabel
integer :: rank,tag

numF=(global_num_sites+global_num_links+global_num_faces)*(NMAT*NMAT-1)
allocate(mode_eta(1:numF,1:dimG,1:num_sites))
allocate(mode_lambda(1:numF,1:dimG,1:num_links))
allocate(mode_chi(1:numF,1:dimG,1:num_faces))
allocate(Dinv_check(1:numF,1:numF))
mode_eta=(0d0,0d0)
mode_lambda=(0d0,0d0)
mode_chi=(0d0,0d0)
Dinv_check=(0d0,0d0)

do ls=1,num_sites; do a=1,dimG
  Dlabel=0
  do gs=1,global_num_sites; do b=1,dimG
    Dlabel=Dlabel+1 
    do l=1,NMAT; do k=1,NMAT; do j=1,NMAT; do i=1,NMAT
      mode_eta(Dlabel,a,ls)=mode_eta(Dlabel,a,ls)&
        + Geta_eta(i,j,k,l,gs,ls)*T(j,i,b)*T(l,k,a)
    enddo; enddo; enddo; enddo
    !write(*,*) MYRANK,Dlabel, (global_site_of_local(ls)-1)+a,  mode_eta(Dlabel,a,ls)
  enddo; enddo
  !!
  do gl=1,global_num_links; do b=1,dimG
    Dlabel=Dlabel+1 
    do l=1,NMAT; do k=1,NMAT; do j=1,NMAT; do i=1,NMAT
      mode_eta(Dlabel,a,ls)=mode_eta(Dlabel,a,ls)&
        + Glambda_eta(i,j,k,l,gl,ls)*T(j,i,b)*T(l,k,a)
    enddo; enddo; enddo; enddo
  enddo; enddo
  !!
  do gf=1,global_num_faces; do b=1,dimG
    Dlabel=Dlabel+1 
    do l=1,NMAT; do k=1,NMAT; do j=1,NMAT; do i=1,NMAT
      mode_eta(Dlabel,a,ls)=mode_eta(Dlabel,a,ls)&
        + Gchi_eta(i,j,k,l,gf,ls)*T(j,i,b)*T(l,k,a)
    enddo; enddo; enddo; enddo
  enddo; enddo
enddo; enddo
!!!!!!!!!!!!!!!!!
do ll=1,num_links; do a=1,dimG
  Dlabel=0
  do gs=1,global_num_sites; do b=1,dimG
    Dlabel=Dlabel+1 
    do l=1,NMAT; do k=1,NMAT; do j=1,NMAT; do i=1,NMAT
      mode_lambda(Dlabel,a,ll)=mode_lambda(Dlabel,a,ll)&
        + Geta_lambda(i,j,k,l,gs,ll)*T(j,i,b)*T(l,k,a)
    enddo; enddo; enddo; enddo
  enddo; enddo
  !!
  do gl=1,global_num_links; do b=1,dimG
    Dlabel=Dlabel+1 
    do l=1,NMAT; do k=1,NMAT; do j=1,NMAT; do i=1,NMAT
      mode_lambda(Dlabel,a,ll)=mode_lambda(Dlabel,a,ll)&
        + Glambda_lambda(i,j,k,l,gl,ll)*T(j,i,b)*T(l,k,a)
    enddo; enddo; enddo; enddo
  enddo; enddo
  !!
  do gf=1,global_num_faces; do b=1,dimG
    Dlabel=Dlabel+1 
    do l=1,NMAT; do k=1,NMAT; do j=1,NMAT; do i=1,NMAT
      mode_lambda(Dlabel,a,ll)=mode_lambda(Dlabel,a,ll)&
        + Gchi_lambda(i,j,k,l,gf,ll)*T(j,i,b)*T(l,k,a)
    enddo; enddo; enddo; enddo
  enddo; enddo
enddo; enddo
!write(*,*) "======", myrank
!write(*,*) mode_lambda
!write(*,*) "======"
!!!!!!!!!!!!!!!!!
do lf=1,num_faces; do a=1,dimG
  Dlabel=0
  do gs=1,global_num_sites; do b=1,dimG
    Dlabel=Dlabel+1 
    do l=1,NMAT; do k=1,NMAT; do j=1,NMAT; do i=1,NMAT
      mode_chi(Dlabel,a,lf)=mode_chi(Dlabel,a,lf)&
        + Geta_chi(i,j,k,l,gs,lf)*T(j,i,b)*T(l,k,a)
    enddo; enddo; enddo; enddo
  enddo; enddo
  !!
  do gl=1,global_num_links; do b=1,dimG
    Dlabel=Dlabel+1 
    do l=1,NMAT; do k=1,NMAT; do j=1,NMAT; do i=1,NMAT
      mode_chi(Dlabel,a,lf)=mode_chi(Dlabel,a,lf)&
        + Glambda_chi(i,j,k,l,gl,lf)*T(j,i,b)*T(l,k,a)
    enddo; enddo; enddo; enddo
  enddo; enddo
  !!
  do gf=1,global_num_faces; do b=1,dimG
    Dlabel=Dlabel+1 
    do l=1,NMAT; do k=1,NMAT; do j=1,NMAT; do i=1,NMAT
      mode_chi(Dlabel,a,lf)=mode_chi(Dlabel,a,lf)&
        + Gchi_chi(i,j,k,l,gf,lf)*T(j,i,b)*T(l,k,a)
    enddo; enddo; enddo; enddo
  enddo; enddo
enddo; enddo
!!!!!!!!!!!!!!!!!

Dlabel=1
tag=0
do gs=1,global_num_sites
  ls=local_site_of_global(gs)%label_
  rank=local_site_of_global(gs)%rank_ 
  if( rank==0 ) then
    if( MYRANK==0 ) then
      Dinv_check(:,Dlabel:Dlabel+dimG-1)=mode_eta(:,:,ls)
    endif
  else
    tag=tag+1
    if( MYRANK==rank ) then
      call MPI_SEND(mode_eta(:,:,ls), &
        numF*dimG,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,IERR)
    elseif( MYRANK==0 ) then
      call MPI_RECV(Dinv_check(:,Dlabel:Dlabel+dimG-1), &
        numF*dimG,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,ISTATUS,IERR)
    endif
  endif
  Dlabel=Dlabel+dimG
enddo
!!
do gl=1,global_num_links
  ll=local_link_of_global(gl)%label_
  rank=local_link_of_global(gl)%rank_ 
  if( rank==0 ) then
    if( MYRANK==0 ) then
      Dinv_check(:,Dlabel:Dlabel+dimG-1)=mode_lambda(:,:,ll)
    endif
  else
    tag=tag+1
    if( MYRANK==rank ) then
      call MPI_SEND(mode_lambda(:,:,ll), &
        numF*dimG,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,IERR)
    elseif( MYRANK==0 ) then
      call MPI_RECV(Dinv_check(:,Dlabel:Dlabel+dimG-1), &
        numF*dimG,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,ISTATUS,IERR)
    endif
  endif
  Dlabel=Dlabel+dimG
enddo
!!
do gf=1,global_num_faces
  lf=local_face_of_global(gf)%label_
  rank=local_face_of_global(gf)%rank_ 
  if( rank==0 ) then
    if( MYRANK==0 ) then
      Dinv_check(:,Dlabel:Dlabel+dimG-1)=mode_chi(:,:,lf)
    endif
  else
    tag=tag+1
    if( MYRANK==rank ) then
      call MPI_SEND(mode_chi(:,:,lf), &
        numF*dimG,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,IERR)
    elseif( MYRANK==0 ) then
      call MPI_RECV(Dinv_check(:,Dlabel:Dlabel+dimG-1), &
        numF*dimG,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,ISTATUS,IERR)
    endif
  endif
  Dlabel=Dlabel+dimG
enddo

if( MYRANK==0 ) then
  write(*,*) "=== check if correlation is Dinv  ==="
  do i=1,numF
    do j=1,numF
      if( cdabs(Dinv(i,j)-Dinv_check(i,j)) > 1d-5 ) then 
        write(*,*) i,j,Dinv(i,j), Dinv_check(i,j)
      endif
    enddo
  enddo
  write(*,*) "=== END check if correlation is Dinv  ==="
endif


end subroutine correlation_to_Dinv




