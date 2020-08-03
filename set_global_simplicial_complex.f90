!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_global_simplicial_complex
use simplicial_complex
implicit none
!!!!!!!!!!!!!!!!!!
integer l,origin,tip
integer s,f
integer :: fsize
integer, allocatable :: fsites(:), faces_l(:), sites_f(:)
integer :: FaceSize
integer, allocatable :: sites(:)
integer k,j,i
character(128) tmp
double precision :: alpha, tmp_U1Rmass
! open SC_FILE 
if(MYRANK==0) then
  open(SC_FILE, file=SC_FILE_NAME, status='old',action='READ')
  read(SC_FILE,'()') ! skip 1 line
  read(SC_FILE,*) global_num_sites
  read(SC_FILE,*) global_num_links
  read(SC_FILE,*) global_num_faces
  read(SC_FILE,*) LatticeSpacing
  read(SC_FILE,'()') ! skip 1 line
  !write(*,*) num_sites,num_links,num_faces,LatticeSpacing
endif
call MPI_BCAST(global_num_sites,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(global_num_links,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(global_num_faces,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(LatticeSpacing,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!ひとまずglobalなalpha,betaを設定してしまう。
  allocate( global_alpha_s(1:global_num_sites) )
  allocate( global_alpha_l(1:global_num_links) )
  allocate( global_alpha_f(1:global_num_faces) )
  allocate( global_beta_f(1:global_num_faces) )
! U(1)_R factor on global links
  allocate( global_U1Rfactor_link(1:global_num_links) )
  allocate( global_U1R_ratio(1:global_num_links) )
! U(1)_R factor on global sites
  allocate( global_U1Rfactor_site(1:global_num_sites) )
! initialize the simplicial complex
call init_smpcom(SC,global_num_sites,global_num_links,global_num_faces)


!allocate( local_site_of_global(1:global_num_sites) )
!allocate( local_link_of_global(1:global_num_links) )
!allocate( local_face_of_global(1:global_num_faces) )

!!!!!!!!!!!!!!!!!!!!!!!!!
!! Set alpha_s
if( MYRANK==0 ) then
  do k=1,global_num_sites
    read(SC_FILE,*) s,alpha!,tmp_U1Rmass
    global_alpha_s(s)=alpha
  enddo
  read(SC_FILE,'()') ! skip 1 line
endif

!!!!!!!!!!!!!!!!!!!!!!!!!
!! Set alpha_l and U(1)_R mass
do k=1,global_num_links
  if( MYRANK==0 ) then 
    read(SC_FILE,*) l,origin,tip,alpha,tmp_U1Rmass
    global_alpha_l(l)=alpha
    !global_U1Rfactor_link(l)=global_U1Rfactor_site(tip)/global_U1Rfactor_site(origin)
    global_U1Rfactor_link(l)=cdexp( (0d0,2d0)*dacos(-1d0)*tmp_U1Rmass )
  endif
  call MPI_BCAST(l,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(origin,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(tip,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call put_link_sc(SC,l,origin,tip)
enddo
call MPI_BCAST(global_U1Rfactor_link,global_num_links,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
if( MYRANK==0 ) then 
  read(SC_FILE,'()') ! skip 1 line
endif

!!!!!!!!!!!!!!!!!!!!!!!!!
!! Set alpha_f
do k=1,global_num_faces
  if(MYRANK==0) then
    read(SC_FILE,"(I6,X,I6,X)",advance="no") f,fsize
    !write(*,*) f,fsize
  endif
  call MPI_BCAST(f,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(fsize,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)

  allocate( fsites(1:fsize ) )

! set sites constructing i'th face
  if (MYRANK==0) then
    do j=1,fsize
      read(SC_FILE,'(I6,X)',advance='no') fsites(j) 
    enddo
    !write(*,*) fsites
    read(SC_FILE,*) global_alpha_f(f),global_beta_f(f)
  endif
  call MPI_BCAST(fsites,fsize,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
!write(*,*) MYRANK, f,fsites
  call put_face_sc(SC,f,fsites)
!write(*,*) "2) ", MYRANK, f,fsites

  deallocate( fsites )
enddo
if (MYRANK==0) then
  close(SC_FILE)
endif


!write(*,*) "test"
!ここまでで、global_alpha,betaがRANK 0に設定された。

! set links
allocate(global_link_org(1:global_num_links))
allocate(global_link_tip(1:global_num_links))
do l=1,global_num_links
  call get_link_sc(sc,l,global_link_org(l),global_link_tip(l))
enddo

! tips of links from s
allocate(global_linktip_from_s(1:global_num_sites))
do s=1,global_num_sites
  call get_links_from_s_sc(sc,s,&
    global_linktip_from_s(s)%labels_,&
    global_linktip_from_s(s)%sites_)
  global_linktip_from_s(s)%num_=size(global_linktip_from_s(s)%labels_)
enddo

! origins of links to s
allocate(global_linkorg_to_s(1:global_num_sites))
do s=1,global_num_sites
  call get_links_to_s_sc(sc,s,&
    global_linkorg_to_s(s)%labels_,&
    global_linkorg_to_s(s)%sites_)
  global_linkorg_to_s(s)%num_=size(global_linkorg_to_s(s)%labels_)
enddo

! links included in a face f
allocate(global_links_in_f(1:global_num_faces))
do f=1,global_num_faces
  call get_links_in_face_sc(sc,f,&
    global_links_in_f(f)%num_,&
    sites,&
    global_links_in_f(f)%link_labels_,&
    global_links_in_f(f)%link_dirs_)
enddo

!write(*,*) "test4"
! faces included in a link l
allocate(global_face_in_l(1:global_num_links))
do l=1,global_num_links
  call get_faces_in_link_sc(sc,l,faces_l)
  global_face_in_l(l)%num_=size(faces_l)
  allocate(global_face_in_l(l)%label_(1:size(faces_l)))
  global_face_in_l(l)%label_=faces_l
enddo
  
! sites included in a face f
allocate(global_sites_in_f(1:global_num_faces))
do f=1,global_num_faces
  call get_face_components_sc(sc,f,sites_f)
  global_sites_in_f(f)%num_=size(sites_f)
  allocate(global_sites_in_f(f)%label_(1:size(sites_f)))
  global_sites_in_f(f)%label_=sites_f
enddo

! faces touching site s
allocate(global_face_in_s(1:global_num_sites))
do s=1,global_num_sites
  global_face_in_s(s)%num_=0
enddo
do f=1,global_num_faces
  do i=1,global_sites_in_f(f)%num_
    s=global_sites_in_f(f)%label_(i)
    global_face_in_s(s)%num_= global_face_in_s(s)%num_ + 1
  enddo
enddo
do s=1,global_num_sites
  allocate(global_face_in_s(s)%label_(1:global_face_in_s(s)%num_))
  global_face_in_s(s)%num_=0
enddo
do f=1,global_num_faces
  do i=1,global_sites_in_f(f)%num_
    s=global_sites_in_f(f)%label_(i)
    global_face_in_s(s)%num_= global_face_in_s(s)%num_ + 1
    global_face_in_s(s)%label_( global_face_in_s(s)%num_ ) = f
  enddo
enddo

!write(*,*) "test3"
!! set global_U1Rfactor_site from global_U1Rfactor_link
call set_global_U1Rfactor_site
!write(*,*) "test1"
!! set global_U1R_ratio
call set_global_U1R_ratio
!write(*,*) "test2"

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! make global_U1Rfactor_site from global_U1Rfacor_link in all ranks
subroutine set_global_U1Rfactor_site
implicit none

!integer :: U1Rfactor_site(1:global_num_sites)
integer :: groupS(1:global_num_sites)
integer :: groupT(1:global_num_sites)
integer :: groupU(1:global_num_sites)

integer :: s,t,l,ls
integer :: i,j,k
integer :: numS,numT,numU
integer :: info

double precision :: tmp

global_U1Rfactor_site=(1d0,0d0)
!!!!!!!!!!!
groupS=0
numS=0
!!!!!!!!!!!
groupT=0
numT=1
groupT(1)=1
!!!!!!!!!!!
groupU=0
numU=0
do while (numS < global_num_sites) 
  do i=1,numT
    s=groupT(i)
    !! sから出るlink
    do j=1,global_linktip_from_s(s)%num_
      t=global_linktip_from_s(s)%sites_(j)
      l=global_linktip_from_s(s)%labels_(j)

      info=1
      do k=1,numS
        if( t==groupS(k) ) then
          info=0
          exit
        endif
      enddo
      do k=1,numT
        if( t==groupT(k) ) then
          info=0
          exit
        endif
      enddo
      do k=1,numU
        if( t==groupU(k) ) then
          info=0
          exit
        endif
      enddo
      if( info==1 ) then
        numU=numU+1
        global_U1Rfactor_site(t)=global_U1Rfactor_site(s)*global_U1Rfactor_link(l)
        groupU(numU)=t
      endif
    enddo
    !! sに向かうlink
    do j=1,global_linkorg_to_s(s)%num_
      t=global_linkorg_to_s(s)%sites_(j)
      l=global_linkorg_to_s(s)%labels_(j)

      info=1
      do k=1,numS
        if( t==groupS(k) ) then
          info=0
          exit
        endif
      enddo
      do k=1,numT
        if( t==groupT(k) ) then
          info=0
          exit
        endif
      enddo
      do k=1,numU
        if( t==groupU(k) ) then
          info=0
          exit
        endif
      enddo
      if( info==1 ) then
        numU=numU+1
        global_U1Rfactor_site(t)=global_U1Rfactor_site(s)*dconjg(global_U1Rfactor_link(l))
        groupU(numU)=t
      endif
    enddo



  enddo
  !! update groupS
  groupS(numS+1:numS+numT)=groupT(1:numT)
  numS=numS+numT
  !!
  !write(*,*) numS
  !! update groupT
  groupT=groupU
  numT=numU
  !! initialize groupU
  groupU=0
  numU=0
  !! 
enddo
end subroutine set_global_U1Rfactor_site

subroutine set_global_U1R_ratio
implicit none

integer :: l,s,t

do l=1,global_num_links
  s=global_link_org(l)
  t=global_link_tip(l)
  global_U1R_ratio(l)=&
    global_U1Rfactor_site(t)/global_U1Rfactor_site(s)/global_U1Rfactor_link(l)
enddo

end subroutine set_global_U1R_ratio



end subroutine set_global_simplicial_complex


