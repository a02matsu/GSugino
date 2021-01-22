subroutine set_local_data(local_site_list)
#ifdef PARALLEL
use parallel
#endif
implicit none
type(SITE_DIST), intent(in) :: local_site_list(0:NPROCS-1)
integer :: s,l,f,i,j

#ifdef PARALLEL
  !! global {site,link,face} から local {site,link,face}へのmapを作成
  !!   global_num_{sites,links,faces} 
  !!   num_sites, num_links, num_faces
!  write(*,*) MYRANK, "test1"
  call set_global_to_local(local_site_list)
  !! 各ノードが計算に必要なsiteの情報と、site情報の送り先、受信元を設定する
  !!   num_necessary_sites
  !!   global_site_of_local(s)
  !!   send_sites(s)
  !!   recv_sites(s)
!  write(*,*) MYRANK, "test2"
  call set_local_sites
  !! 各ノードが計算に必要なlinkの情報と、link情報の送り先、受信元を設定する
  !!   num_necessary_links
  !!   global_link_of_local(l)
  !!   send_links(l)
  !!   recv_links(l)
!  write(*,*) MYRANK, "test3"
  call set_local_links
  !! 各ノードが計算に必要なfaceの情報を設定する
  !!   num_necessary_faces
  !!   global_face_of_local(l)
!  write(*,*) MYRANK, "test4"
  call set_local_faces
  !! localなlinkのorgとtipを設定
!  write(*,*) MYRANK, "test5"
  call set_local_link_org_and_tip
  !! localなsiteから出入りするlocal linkのラベルを設定
!  write(*,*) MYRANK, "test6"
  call set_local_link_fromto_site
  !! localなfaceに含まれるlinkのlocalラベルを設定
!  write(*,*) MYRANK, "test7"
  call set_local_links_in_f
  !! localなlinkを共有するfaceのlocalラベルを設定
!  write(*,*) MYRANK, "test8"
  call set_local_face_in_l
!  do l=1,num_links
!    do i=1,face_in_l(l)%num_
!      f=face_in_l(l)%label_(i)
!      do j=1,links_in_f(f)%num_
!        if( links_in_f(f)%link_labels_(j)==l ) then
!          write(*,*) MYRANK, global_face_of_local(f), global_link_of_local(l), links_in_f(f)%link_dirs_(j)
!        endif
!      enddo
!    enddo
!  enddo
!  call stop_for_test
  !! localなfaceに含まれるsiteのlocalラベルを設定
  !! ただし、使うのは sites_in_f(f)%label_(1)だけなので、
  !! はじめから size=1 にしておく。
!  write(*,*) MYRANK, "test9"
  call set_local_sites_in_f

!  write(*,*) MYRANK, "test10"
  call set_num_local_faces_in_s

  !! alpha と beta を割り振る
  !write(*,*) MYRANK, "test11"
  call set_local_alpha_beta
  !! U(1)_R mass を割り振る
  !write(*,*) MYRANK, "test12"
  call set_local_U1Rmass
  !call set_U1Rfactor_on_sites
#else
num_sites=global_num_sites
num_links=global_num_links
num_faces=global_num_faces
!!
allocate(alpha_s(1:num_sites))
allocate(alpha_l(1:num_links))
allocate(alpha_f(1:num_faces))
allocate(beta_f(1:num_faces))
alpha_s=globa_alpha_s
alpha_l=globa_alpha_l
alpha_f=globa_alpha_f
beta_f=globa_beta_f

! set links
allocate(link_org(1:num_links))
allocate(link_tip(1:num_links))
do l=1,num_links
  call get_link_sc(sc,l,link_org(l),link_tip(l))
enddo

! tips of links from s
allocate(linktip_from_s(1:num_sites))
do s=1,num_sites
  call get_links_from_s_sc(sc,s,&
    linktip_from_s(s)%labels_,&
    linktip_from_s(s)%sites_)
  linktip_from_s(s)%num_=size(linktip_from_s(s)%labels_)
enddo

! origins of links to s
allocate(linkorg_to_s(1:num_sites))
do s=1,num_sites
  call get_links_to_s_sc(sc,s,&
    linkorg_to_s(s)%labels_,&
    linkorg_to_s(s)%sites_)
  linkorg_to_s(s)%num_=size(linkorg_to_s(s)%labels_)
enddo

! links included in a face f
allocate(links_in_f(1:num_faces))
do f=1,num_faces
  call get_links_in_face_sc(sc,f,&
    links_in_f(f)%num_,&
    sites,&
    links_in_f(f)%link_labels_,&
    links_in_f(f)%link_dirs_)
enddo

! faces included in a link l
allocate(face_in_l(1:num_links))
do l=1,num_links
  call get_faces_in_link_sc(sc,l,faces_l)
  face_in_l(l)%num_=size(faces_l)
  allocate(face_in_l(l)%label_(1:size(faces_l)))
  face_in_l(l)%label_=faces_l
enddo
  
! sites included in a face f
allocate(sites_in_f(1:num_faces))
do f=1,num_faces
  call get_face_components_sc(sc,f,sites_f)
  sites_in_f(f)%num_=size(sites_f)
  allocate(sites_in_f(f)%label_(1:size(sites_f)))
  sites_in_f(f)%label_=sites_f
enddo
#endif

end subroutine set_local_data

#ifdef PARALLEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! global {site,link,face} から local {site,link,face}へのmapを作成
subroutine set_global_to_local(local_site_list)
use parallel
implicit none

type(SITE_DIST), intent(in) :: local_site_list(0:NPROCS-1)
integer, allocatable :: tmp_global_site_of_local(:)
integer :: tmp_num_sites, tmp_num_links, tmp_num_faces
integer :: s,l,f,rank,i

allocate(local_site_of_global(1:global_num_sites))
allocate(local_link_of_global(1:global_num_links))
allocate(local_face_of_global(1:global_num_faces))

!! local_site_of_global(s)
if(MYRANK==0)  then
  do rank=0,NPROCS-1
    !write(*,*) rank, local_site_list(rank)%site_list_
    tmp_num_sites=local_site_list(rank)%num_sites_
    if(rank==0) then 
      num_sites=tmp_num_sites
    else
      call MPI_SEND(tmp_num_sites,1,MPI_INTEGER,rank,1,MPI_COMM_WORLD,IERR)
    endif
    do i=1,tmp_num_sites
      s=local_site_list(rank)%site_list_(i)
      local_site_of_global(s)%rank_=rank
      local_site_of_global(s)%label_=i
    enddo
  enddo
  !do s=1,global_num_sites
    !write(*,*) s, local_site_of_global(s)
  !enddo
else
  call MPI_RECV(num_sites,1,MPI_INTEGER,0,1,MPI_COMM_WORLD,ISTATUS,IERR)
endif
call MPI_BCAST(local_site_of_global,global_num_sites,MPI_LOCAL_LABEL,0,MPI_COMM_WORLD,IERR)

!! local_link_of_global(l)
if(MYRANK==0) then
  do rank=0,NPROCS-1
    tmp_num_links=0
    do l=1,global_num_links
      do i=1,local_site_list(rank)%num_sites_
      s=local_site_list(rank)%site_list_(i)
        if( s == global_link_org(l) ) then
          tmp_num_links=tmp_num_links+1
          local_link_of_global(l)%rank_=rank
          local_link_of_global(l)%label_=tmp_num_links
        endif
      enddo
    enddo
    if(rank==0) then
      num_links=tmp_num_links
    else
      call MPI_SEND(tmp_num_links,1,MPI_INTEGER,rank,2,MPI_COMM_WORLD,IERR)
    endif
  enddo
else
  call MPI_RECV(num_links,1,MPI_INTEGER,0,2,MPI_COMM_WORLD,ISTATUS,IERR)
endif
call MPI_BCAST(local_link_of_global,global_num_links,MPI_LOCAL_LABEL,0,MPI_COMM_WORLD,IERR)
!if(MYRANK==0) then
!  do l=1,global_num_links
!    write(*,*) l,local_link_of_global(l)
!  enddo
!endif



!! local_face_of_global(f)
if(MYRANK==0) then
  do rank=0,NPROCS-1
    tmp_num_faces=0
    do f=1,global_num_faces
      do i=1,local_site_list(rank)%num_sites_
        s=local_site_list(rank)%site_list_(i)
        if( s == global_sites_in_f(f)%label_(1) ) then
          tmp_num_faces=tmp_num_faces+1
          local_face_of_global(f)%rank_=rank
          local_face_of_global(f)%label_=tmp_num_faces
        endif
      enddo
    enddo
    !write(*,*) MYRANK,rank,tmp_num_faces
    if(rank==0) then
      num_faces=tmp_num_faces
    else
      call MPI_SEND(tmp_num_faces,1,MPI_INTEGER,rank,3,MPI_COMM_WORLD,IERR)
    endif
  enddo
else
  call MPI_RECV(num_faces,1,MPI_INTEGER,0,3,MPI_COMM_WORLD,ISTATUS,IERR)
endif
call MPI_BCAST(local_face_of_global,global_num_faces,MPI_LOCAL_LABEL,0,MPI_COMM_WORLD,IERR)



end subroutine set_global_to_local


!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! 計算に必要なsiteの情報を設定
subroutine set_local_sites
use parallel
implicit none

integer s,l,f,i,nsend,nrecv,info,k,kk
integer :: tmp_global_site_of_local(1:global_num_sites)
!integer :: tmp_send_sites(1:num_sites)
!integer :: tmp_recv_sites(1:num_sites)
type(LOCAL_LABEL) :: tmp_send_sites(1:global_num_sites)
type(LOCAL_LABEL) :: tmp_recv_sites(1:global_num_sites)

!! 計算を担当するsiteを tmp_global_site_of_local の最初の方に設定
num_necessary_sites=0
tmp_global_site_of_local=0
do s=1,global_num_sites
  if( local_site_of_global(s)%rank_ == MYRANK ) then
    num_necessary_sites = num_necessary_sites + 1
    tmp_global_site_of_local(local_site_of_global(s)%label_)=s
    !tmp_global_site_of_local(num_necessary_sites)=s
  endif
enddo

!write(*,*) MYRANK, tmp_global_site_of_local
!! この段階で、num_necesarry_sites=num_sites
!! linkの始点と終点が別のnodeに属していたら、
!! tipの位置を、orgを担当するnodeに送る必要がある
!! (1) 担当するlinkの始点と終点
!! (2) 担当するサイトを始点とするlinkの終点
!! (3) 担当するサイトを終点とするlinkの始点
!! (4) 担当するfaceに属する全site
!! 重複があり得るので、逐一チェックしながら進める
nsend=0
nrecv=0
!! (1) 担当するlinkの終点を始点を持つrankに送る
do l=1,global_num_links
  !write(*,*) MYRANK, l, global_link_org(l), global_link_tip(l), &
             !local_site_of_global(global_link_org(l))%rank_, &
             !local_site_of_global(global_link_tip(l))%rank_
  if( local_site_of_global(global_link_org(l))%rank_ &
        /= local_site_of_global(global_link_tip(l))%rank_ ) then
    !!!!!!!!!!!
    if( local_site_of_global(global_link_tip(l))%rank_ == MYRANK ) then 
      !! 重複チェック
      info=0
      do i=1,nsend
        if( tmp_send_sites(i)%rank_ == local_site_of_global(global_link_org(l))%rank_ &
          .and. &
            tmp_send_sites(i)%label_ == global_link_tip(l) ) then
          info=1
          exit
        endif
      enddo
      !write(*,*) MYRANK, l, local_site_of_global(global_link_org(l))%rank_, &
            !global_link_tip(l), info
      if (info==0) then 
        nsend=nsend+1
        tmp_send_sites(nsend)%rank_ = local_site_of_global(global_link_org(l))%rank_ 
        tmp_send_sites(nsend)%label_ = global_link_tip(l) 
        !write(*,*) MYRANK, nsend, tmp_send_sites(nsend)%rank_, tmp_send_sites(nsend)%label_ 
      endif
    !!!!!!!!!!!
    elseif( local_site_of_global(global_link_org(l))%rank_ == MYRANK ) then
      !! 重複チェック
      info=0
      do i=1,nrecv
        if( tmp_recv_sites(i)%rank_ == local_site_of_global(global_link_tip(l))%rank_ &
          .and. &
            tmp_recv_sites(i)%label_ == global_link_tip(l) ) then
          info=1
          exit
        endif
      enddo
      if (info==0) then 
        nrecv=nrecv+1
        tmp_recv_sites(nrecv)%rank_ = local_site_of_global(global_link_tip(l))%rank_ 
        tmp_recv_sites(nrecv)%label_ = global_link_tip(l) 
        tmp_global_site_of_local(num_sites+nrecv)=global_link_tip(l)
      endif
    endif
  endif
enddo
!write(*,*) MYRANK, tmp_send_sites
!! (2) 担当するサイトを終点とするlinkの始点
do l=1,global_num_links
  if( local_site_of_global(global_link_org(l))%rank_ &
        /= local_site_of_global(global_link_tip(l))%rank_ ) then
    !!!!!!!!!!!
    if( local_site_of_global(global_link_org(l))%rank_ == MYRANK ) then
      !! 重複チェック
      info=0
      do i=1,nsend
        !write(*,*) MYRANK, l, tmp_send_sites(i)
        if( tmp_send_sites(i)%rank_ == local_site_of_global(global_link_tip(l))%rank_ &
          .and. &
            tmp_send_sites(i)%label_ == global_link_org(l) ) then
          info=1
          exit
        endif
      enddo
      !write(*,*) MYRANK, l, local_site_of_global(global_link_tip(l))%rank_, &
            !global_link_org(l), info
      if (info==0) then 
        nsend=nsend+1
        tmp_send_sites(nsend)%rank_ = local_site_of_global(global_link_tip(l))%rank_ 
        tmp_send_sites(nsend)%label_ = global_link_org(l) 
        !write(*,*) MYRANK, nsend, tmp_send_sites(nsend)%rank_, tmp_send_sites(nsend)%label_ 
      endif
    !!!!!!!!!!!
    elseif( local_site_of_global(global_link_tip(l))%rank_ == MYRANK ) then
      !! 重複チェック
      info=0
      do i=1,nrecv
        if( tmp_recv_sites(i)%rank_ == local_site_of_global(global_link_org(l))%rank_ &
          .and. &
            tmp_recv_sites(i)%label_ == global_link_org(l) ) then
          info=1
          exit
        endif
      enddo
      if (info==0) then 
        nrecv=nrecv+1
        tmp_recv_sites(nrecv)%rank_ = local_site_of_global(global_link_org(l))%rank_ 
        tmp_recv_sites(nrecv)%label_ = global_link_org(l) 
        tmp_global_site_of_local(num_sites+nrecv)=global_link_org(l)
      endif
    endif
  endif
enddo

!do s=1,global_num_sites
!  do k=1,global_linkorg_to_s(s)%num_
!    if( local_site_of_global(global_linkorg_to_s(s)%sites_(k))%rank_ &
!        /= local_site_of_global(s)%rank_ ) then
!      if(local_site_of_global(global_linkorg_to_s(s)%sites_(k))%rank_ == MYRANK) then
!        !! 重複チェック
!        info=0
!        do i=1,nsend
!          if( tmp_send_sites(i)%rank_ &
!              == local_site_of_global(s)%rank_ &
!            .and. &
!              tmp_send_sites(i)%label_  &
!              == global_linkorg_to_s(s)%sites_(k) &
!              ) then
!            info=1
!            exit
!          endif
!        enddo
!        if (info==0) then 
!          !write(*,*) "(S)site",linkorg_to_s(s)%sites_(k),"from",MYRANK,"to",local_site_of_global(s)%rank_
!          nsend=nsend+1
!          tmp_send_sites(nsend)%rank_ = local_site_of_global(s)%rank_ 
!          tmp_send_sites(nsend)%label_ = global_linkorg_to_s(s)%sites_(k)
!        endif
!  !!!!!!!!!!!
!      elseif( local_site_of_global(s)%rank_ == MYRANK ) then
!        !! 重複チェック
!        info=0
!        do i=1,nrecv
!          if( tmp_recv_sites(i)%rank_ &
!            == local_site_of_global(global_linkorg_to_s(s)%sites_(k))%rank_ &
!            .and. &
!              tmp_recv_sites(i)%label_ == global_linkorg_to_s(s)%sites_(k) ) then
!            info=1
!            exit
!          endif
!        enddo
!        if (info==0) then 
!          !write(*,*) "(R)site",linkorg_to_s(s)%sites_(k),"from",local_site_of_global(linkorg_to_s(s)%sites_(k))%rank_,"to",MYRANK
!          nrecv=nrecv+1
!          tmp_recv_sites(nrecv)%rank_  &
!            = local_site_of_global(global_linkorg_to_s(s)%sites_(k))%rank_ 
!          tmp_recv_sites(nrecv)%label_ = global_linkorg_to_s(s)%sites_(k)
!          tmp_global_site_of_local(num_sites+nrecv) = global_linkorg_to_s(s)%sites_(k)
!        endif
!      endif
!    endif
!  enddo
!enddo

!! (3) 担当するサイトを始点とするlinkの終点
!! これは(1)で尽きている
!do s=1,global_num_sites
!!write(*,*) s,global_linktip_from_s(s)%num_
!  do k=1,global_linktip_from_s(s)%num_
!    if( local_site_of_global(global_linktip_from_s(s)%sites_(k))%rank_  &
!        /= local_site_of_global(s)%rank_ ) then
!      if(local_site_of_global(global_linktip_from_s(s)%sites_(k))%rank_ == MYRANK) then
!        !! 重複チェック
!        info=0
!        do i=1,nsend
!          if( tmp_send_sites(i)%rank_ &
!              == local_site_of_global(s)%rank_ &
!            .and. &
!              tmp_send_sites(i)%label_  &
!              == global_linktip_from_s(s)%sites_(k)) then
!            info=1
!            exit
!          endif
!        enddo
!        write(*,*) info
!        if (info==0) then 
!          nsend=nsend+1
!          tmp_send_sites(nsend)%rank_ = local_site_of_global(s)%rank_ 
!          tmp_send_sites(nsend)%label_ = global_linktip_from_s(s)%sites_(k)
!        endif
!  !!!!!!!!!!!
!      elseif( local_site_of_global(s)%rank_ == MYRANK ) then
!        !! 重複チェック
!        info=0
!        do i=1,nrecv
!          if( tmp_recv_sites(i)%rank_ &
!            == local_site_of_global(global_linktip_from_s(s)%sites_(k))%rank_ &
!            .and. &
!              tmp_recv_sites(i)%label_ == global_linktip_from_s(s)%sites_(k) ) then
!            info=1
!            exit
!          endif
!        enddo
!        write(*,*) info
!        if (info==0) then 
!          nrecv=nrecv+1
!          tmp_recv_sites(nrecv)%rank_  &
!            = local_site_of_global(global_linktip_from_s(s)%sites_(k))%rank_ 
!          tmp_recv_sites(nrecv)%label_ = global_linktip_from_s(s)%sites_(k)
!          tmp_global_site_of_local(num_sites+nrecv)=global_linktip_from_s(s)%sites_(k)
!        endif
!      endif
!    endif
!  enddo
!enddo

!! (4) 担当するfaceを作る全てのsite
do f=1,global_num_faces
  do k=1, global_sites_in_f(f)%num_
    s=global_sites_in_f(f)%label_(k)
    if( local_face_of_global(f)%rank_ /= local_site_of_global(s)%rank_ ) then
      if( MYRANK == local_site_of_global(s)%rank_ ) then
        !! 重複チェック
        info=0
        do i=1,nsend
          if( tmp_send_sites(i)%rank_ &
              == local_face_of_global(f)%rank_ &
            .and. &
              tmp_send_sites(i)%label_ == s ) then
            info=1
            exit
          endif
        enddo
        if (info==0) then 
          nsend=nsend+1
          tmp_send_sites(nsend)%rank_ = local_face_of_global(f)%rank_ 
          tmp_send_sites(nsend)%label_ = s
        endif
  !!!!!!!!!!!
      elseif( MYRANK == local_face_of_global(f)%rank_ ) then
        !! 重複チェック
        info=0
        do i=1,nrecv
          if( tmp_recv_sites(i)%rank_ &
            == local_site_of_global(s)%rank_ &
            .and. &
              tmp_recv_sites(i)%label_ == s ) then
            info=1
            exit
          endif
        enddo
        if (info==0) then 
          nrecv=nrecv+1
          tmp_recv_sites(nrecv)%rank_  &
            = local_site_of_global(s)%rank_
          tmp_recv_sites(nrecv)%label_ = s
          tmp_global_site_of_local(num_sites+nrecv)=s
        endif
      endif
    endif
  enddo
enddo

!! (5) 担当するlinkを共有するface の代表点
do l=1,global_num_links
  do kk=1,global_face_in_l(l)%num_
    f=global_face_in_l(l)%label_(kk)

    do k=1, 1 !global_sites_in_f(f)%num_
      s=global_sites_in_f(f)%label_(k)
      if( local_link_of_global(l)%rank_ /= local_site_of_global(s)%rank_ ) then
        if( MYRANK == local_site_of_global(s)%rank_ ) then
          !! 重複チェック
          info=0
          do i=1,nsend
            if( tmp_send_sites(i)%rank_ &
                == local_link_of_global(l)%rank_ &
              .and. &
                tmp_send_sites(i)%label_ == s ) then
              info=1
              exit
            endif
          enddo
          if (info==0) then 
            nsend=nsend+1
            tmp_send_sites(nsend)%rank_ = local_link_of_global(l)%rank_ 
            tmp_send_sites(nsend)%label_ = s
          endif
  !!!!!!!!!!!
        elseif( MYRANK == local_link_of_global(l)%rank_ ) then
          !! 重複チェック
          info=0
          do i=1,nrecv
            if( tmp_recv_sites(i)%rank_ &
              == local_site_of_global(s)%rank_ &
              .and. &
                tmp_recv_sites(i)%label_ == s ) then
              info=1
              exit
            endif
          enddo
          if (info==0) then 
            nrecv=nrecv+1
            tmp_recv_sites(nrecv)%rank_  &
              = local_site_of_global(s)%rank_
            tmp_recv_sites(nrecv)%label_ = s
            tmp_global_site_of_local(num_sites+nrecv)=s
          endif
        endif
      endif
    enddo
  enddo
enddo

!! 重複度を除いて配列を定義
num_send_sites=nsend
num_recv_sites=nrecv
allocate( send_sites(num_send_sites) )
allocate( recv_sites(num_recv_sites) )
num_necessary_sites = num_sites + num_recv_sites
allocate( global_site_of_local(1:num_necessary_sites) )

!! 代入
do s=1,num_send_sites
  send_sites(s)%rank_ = tmp_send_sites(s)%rank_
  send_sites(s)%label_ = local_site_of_global(tmp_send_sites(s)%label_)%label_
  !write(*,*) MYRANK, send_sites(s)%rank_, send_sites(s)%label_, tmp_send_sites(s)%label_
enddo
do s=1,num_recv_sites
  recv_sites(s)%rank_ = tmp_recv_sites(s)%rank_
  recv_sites(s)%label_ = num_sites+s
enddo
do s=1,num_necessary_sites
  global_site_of_local(s) = tmp_global_site_of_local(s)
enddo

end subroutine set_local_sites

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! 計算に必要なlinkの情報を設定
subroutine set_local_links
use parallel
implicit none

integer s,l,f,i,j,k,ll,ll_label
integer :: tmp_global_link_of_local(1:global_num_links)
type(LOCAL_LABEL) :: tmp_send_links(1:global_num_links)
type(LOCAL_LABEL) :: tmp_recv_links(1:global_num_links)
integer :: nsend,nrecv,info


do l=1,global_num_links
  tmp_global_link_of_local(l)=0
  tmp_send_links(l)%rank_=-1
  tmp_send_links(l)%label_=-1
  tmp_recv_links(l)%rank_=-1
  tmp_recv_links(l)%label_=-1
enddo

!! 計算を担当するlinkを tmp_global_link_of_local の最初の方に設定
num_necessary_links=0
tmp_global_link_of_local=0
do l=1,global_num_links
  if( local_link_of_global(l)%rank_ == MYRANK ) then
    num_necessary_links = num_necessary_links + 1
    tmp_global_link_of_local(local_link_of_global(l)%label_)=l
    !tmp_global_link_of_local(num_necessary_links)=l
  endif
enddo
!! この段階で、num_necesarry_links=num_links
!! 計算に必要なのは、
!!   (1) 担当するsiteをtipとするlink
!!   (2) 担当するfaceを構成するlink
!!   (3) 担当するlinkを共有するfaceに含まれる全link
!!   (4) 担当するfaceに含まれるsiteに入る全link
!!   (5) 担当するfaceに含まれるsiteから出る全link
!! 重複があり得るので、逐一チェックしながら進める

!!   (1) 担当するsiteをtipとするlink
nsend=0
nrecv=0
do l=1,global_num_links
  if( local_site_of_global(global_link_tip(l))%rank_ /= &
      local_site_of_global(global_link_org(l))%rank_ ) then
    if( local_site_of_global(global_link_org(l))%rank_ == MYRANK ) then
      !! 重複チェック
      info=0
      do i=1,nsend
        if( tmp_send_links(i)%rank_ == local_site_of_global(global_link_tip(l))%rank_ &
          .and. &
            tmp_send_links(i)%label_ == l ) then
          info=1
          exit
        endif
      enddo
      if( info == 0 ) then
        nsend=nsend+1
        tmp_send_links(nsend)%rank_ = local_site_of_global(global_link_tip(l))%rank_
        tmp_send_links(nsend)%label_ = l
      endif
    !!!!!!!!!!!
    elseif( local_site_of_global(global_link_tip(l))%rank_ == MYRANK ) then
      !! 重複チェック
      info=0
      do i=1,nrecv
        if( tmp_recv_links(i)%rank_ == local_site_of_global(global_link_org(l))%rank_ &
          .and. &
            tmp_recv_links(i)%label_ == l ) then
          info=1
          exit
        endif
      enddo
      if (info==0) then
        nrecv=nrecv+1
        tmp_recv_links(nrecv)%rank_ = local_site_of_global(global_link_org(l))%rank_ 
        tmp_recv_links(nrecv)%label_ = l
        tmp_global_link_of_local(num_links+nrecv)=l
      endif
    endif
  endif
enddo

!!   (2) 担当するfaceを構成するlink
do f=1,global_num_faces
  do j=1,global_links_in_f(f)%num_
    l=global_links_in_f(f)%link_labels_(j)
    if( local_face_of_global(f)%rank_ /= local_link_of_global(l)%rank_ ) then
      if( local_link_of_global(l)%rank_ == MYRANK ) then
        !! 重複チェック
        info=0
        do i=1,nsend
          if( tmp_send_links(i)%rank_ == local_face_of_global(f)%rank_ &
            .and. &
              tmp_send_links(i)%label_ == l ) then
            info=1
            exit
          endif
        enddo
        if( info == 0 ) then
          nsend=nsend+1
          tmp_send_links(nsend)%rank_ = local_face_of_global(f)%rank_ 
          tmp_send_links(nsend)%label_ = l
        endif
      !!!!!!!!!!!
      elseif( local_face_of_global(f)%rank_ == MYRANK ) then
        !! 重複チェック
        info=0
        do i=1,nrecv
          if( tmp_recv_links(i)%rank_ == local_link_of_global(l)%rank_ &
            .and. &
              tmp_recv_links(i)%label_ == l ) then
            info=1
            exit
          endif
        enddo
        if (info==0) then 
          nrecv=nrecv+1
          tmp_recv_links(nrecv)%rank_ = local_link_of_global(l)%rank_
          tmp_recv_links(nrecv)%label_ = l
          tmp_global_link_of_local(num_links+nrecv)=l
        endif
      endif
    endif
  enddo
enddo

!!   (3) 担当するlinkを共有するfaceに含まれる全link
do l=1,global_num_links
  do k=1,global_face_in_l(l)%num_
    f=global_face_in_l(l)%label_(k)
    do ll_label=1,global_links_in_f(f)%num_
      ll=global_links_in_f(f)%link_labels_(ll_label)
      if( local_link_of_global(l)%rank_ /= local_link_of_global(ll)%rank_ ) then 
        if( local_link_of_global(ll)%rank_ == MYRANK ) then
          !! 重複チェック
          info=0
          do i=1,nsend
            if( tmp_send_links(i)%rank_ == local_link_of_global(l)%rank_ &
                .and. &
                tmp_send_links(i)%label_ == ll ) then 
              info=1
              exit
            endif
          enddo
          if (info == 0 ) then
            nsend=nsend+1
            tmp_send_links(nsend)%rank_ = local_link_of_global(l)%rank_ 
            tmp_send_links(nsend)%label_ = ll
            !write(*,*) MYRANK,tmp_send_links(nsend)%rank_,tmp_send_links(nsend)%label_
          endif
        !!!!!!!!!!!
        elseif( local_link_of_global(l)%rank_ == MYRANK ) then
          !! 重複チェック
          info=0
          do i=1,nrecv
            if( tmp_recv_links(i)%rank_ == local_link_of_global(ll)%rank_ &
                .and. &
                tmp_recv_links(i)%label_ == ll ) then 
              info=1
              exit
            endif
          enddo
          if (info == 0 ) then
            nrecv=nrecv+1
            tmp_recv_links(nrecv)%rank_ = local_link_of_global(ll)%rank_ 
            tmp_recv_links(nrecv)%label_ = ll
            tmp_global_link_of_local(num_links+nrecv)=ll
          endif
        endif
      endif
    enddo
  enddo
enddo

!!   (4) necessary siteに入る全link
do i=1,num_necessary_sites
  s=global_site_of_local(i)
  do j=1,global_linkorg_to_s(s)%num_
    ll=global_linkorg_to_s(s)%labels_(j)
    if( local_face_of_global(f)%rank_ /= local_link_of_global(ll)%rank_ ) then 
      if( local_link_of_global(ll)%rank_ == MYRANK ) then
        !! 重複チェック
        info=0
        do i=1,nsend
          if( tmp_send_links(i)%rank_ == local_face_of_global(f)%rank_ &
              .and. &
              tmp_send_links(i)%label_ == ll ) then 
            info=1
            exit
          endif
        enddo
        if (info == 0 ) then
          nsend=nsend+1
          tmp_send_links(nsend)%rank_ = local_face_of_global(f)%rank_ 
          tmp_send_links(nsend)%label_ = ll
          !write(*,*) MYRANK,tmp_send_links(nsend)%rank_,tmp_send_links(nsend)%label_
        endif
      !!!!!!!!!!!
      elseif( local_face_of_global(f)%rank_ == MYRANK ) then
        !! 重複チェック
        info=0
        do i=1,nrecv
          if( tmp_recv_links(i)%rank_ == local_link_of_global(ll)%rank_ &
              .and. &
              tmp_recv_links(i)%label_ == ll ) then 
            info=1
            exit
          endif
        enddo
        if (info == 0 ) then
          nrecv=nrecv+1
          tmp_recv_links(nrecv)%rank_ = local_link_of_global(ll)%rank_ 
          tmp_recv_links(nrecv)%label_ = ll
          tmp_global_link_of_local(num_links+nrecv)=ll
        endif
      endif
    endif
  enddo
enddo

!!   (5) necessary siteから出る全link
do i=1,num_necessary_sites
  s=global_site_of_local(i)
  do j=1,global_linktip_from_s(s)%num_
    ll=global_linktip_from_s(s)%labels_(j)
    if( local_face_of_global(f)%rank_ /= local_link_of_global(ll)%rank_ ) then 
      if( local_link_of_global(ll)%rank_ == MYRANK ) then
        !! 重複チェック
        info=0
        do i=1,nsend
          if( tmp_send_links(i)%rank_ == local_face_of_global(f)%rank_ &
              .and. &
              tmp_send_links(i)%label_ == ll ) then 
            info=1
            exit
          endif
        enddo
        if (info == 0 ) then
          nsend=nsend+1
          tmp_send_links(nsend)%rank_ = local_face_of_global(f)%rank_ 
          tmp_send_links(nsend)%label_ = ll
          !write(*,*) MYRANK,tmp_send_links(nsend)%rank_,tmp_send_links(nsend)%label_
        endif
      !!!!!!!!!!!
      elseif( local_face_of_global(f)%rank_ == MYRANK ) then
        !! 重複チェック
        info=0
        do i=1,nrecv
          if( tmp_recv_links(i)%rank_ == local_link_of_global(ll)%rank_ &
              .and. &
              tmp_recv_links(i)%label_ == ll ) then 
            info=1
            exit
          endif
        enddo
        if (info == 0 ) then
          nrecv=nrecv+1
          tmp_recv_links(nrecv)%rank_ = local_link_of_global(ll)%rank_ 
          tmp_recv_links(nrecv)%label_ = ll
          tmp_global_link_of_local(num_links+nrecv)=ll
        endif
      endif
    endif
  enddo
enddo


!!!   (4) 担当するfaceに含まれるsiteに入る全link
!do f=1,global_num_faces
!  do k=1,global_sites_in_f(f)%num_
!    s=global_sites_in_f(f)%label_(k)
!    do j=1,global_linkorg_to_s(s)%num_
!      ll=global_linkorg_to_s(s)%labels_(j)
!      if( local_face_of_global(f)%rank_ /= local_link_of_global(ll)%rank_ ) then 
!        if( local_link_of_global(ll)%rank_ == MYRANK ) then
!          !! 重複チェック
!          info=0
!          do i=1,nsend
!            if( tmp_send_links(i)%rank_ == local_face_of_global(f)%rank_ &
!                .and. &
!                tmp_send_links(i)%label_ == ll ) then 
!              info=1
!              exit
!            endif
!          enddo
!          if (info == 0 ) then
!            nsend=nsend+1
!            tmp_send_links(nsend)%rank_ = local_face_of_global(f)%rank_ 
!            tmp_send_links(nsend)%label_ = ll
!            !write(*,*) MYRANK,tmp_send_links(nsend)%rank_,tmp_send_links(nsend)%label_
!          endif
!        !!!!!!!!!!!
!        elseif( local_face_of_global(f)%rank_ == MYRANK ) then
!          !! 重複チェック
!          info=0
!          do i=1,nrecv
!            if( tmp_recv_links(i)%rank_ == local_link_of_global(ll)%rank_ &
!                .and. &
!                tmp_recv_links(i)%label_ == ll ) then 
!              info=1
!              exit
!            endif
!          enddo
!          if (info == 0 ) then
!            nrecv=nrecv+1
!            tmp_recv_links(nrecv)%rank_ = local_link_of_global(ll)%rank_ 
!            tmp_recv_links(nrecv)%label_ = ll
!            tmp_global_link_of_local(num_links+nrecv)=ll
!          endif
!        endif
!      endif
!    enddo
!  enddo
!enddo

!!   (5) 担当するfaceに含まれるsiteから出る全link
!do f=1,global_num_faces
!  do k=1,global_sites_in_f(f)%num_
!    s=global_sites_in_f(f)%label_(k)
!    do j=1,global_linktip_from_s(s)%num_
!      ll=global_linktip_from_s(s)%labels_(j)
!      if( local_face_of_global(f)%rank_ /= local_link_of_global(ll)%rank_ ) then 
!        if( local_link_of_global(ll)%rank_ == MYRANK ) then
!          !! 重複チェック
!          info=0
!          do i=1,nsend
!            if( tmp_send_links(i)%rank_ == local_face_of_global(f)%rank_ &
!                .and. &
!                tmp_send_links(i)%label_ == ll ) then 
!              info=1
!              exit
!            endif
!          enddo
!          if (info == 0 ) then
!            nsend=nsend+1
!            tmp_send_links(nsend)%rank_ = local_face_of_global(f)%rank_ 
!            tmp_send_links(nsend)%label_ = ll
!            !write(*,*) MYRANK,tmp_send_links(nsend)%rank_,tmp_send_links(nsend)%label_
!          endif
!        !!!!!!!!!!!
!        elseif( local_face_of_global(f)%rank_ == MYRANK ) then
!          !! 重複チェック
!          info=0
!          do i=1,nrecv
!            if( tmp_recv_links(i)%rank_ == local_link_of_global(ll)%rank_ &
!                .and. &
!                tmp_recv_links(i)%label_ == ll ) then 
!              info=1
!              exit
!            endif
!          enddo
!          if (info == 0 ) then
!            nrecv=nrecv+1
!            tmp_recv_links(nrecv)%rank_ = local_link_of_global(ll)%rank_ 
!            tmp_recv_links(nrecv)%label_ = ll
!            tmp_global_link_of_local(num_links+nrecv)=ll
!          endif
!        endif
!      endif
!    enddo
!  enddo
!enddo


!! 重複度を除いて配列を定義
num_send_links=nsend
num_recv_links=nrecv
allocate( send_links(num_send_links) )
allocate( recv_links(num_recv_links) )
num_necessary_links = num_links + num_recv_links
allocate( global_link_of_local(1:num_necessary_links) )

!! 代入
do l=1,num_send_links
  send_links(l)%rank_ = tmp_send_links(l)%rank_
  send_links(l)%label_ = local_link_of_global(tmp_send_links(l)%label_)%label_
  !write(*,*) MYRANK,l,tmp_send_links(l)%rank_, tmp_send_links(l)%label_,send_links(l)%label_
enddo
do l=1,num_recv_links
  recv_links(l)%rank_ = tmp_recv_links(l)%rank_
  recv_links(l)%label_ = num_links+l
enddo
do l=1,num_necessary_links
  global_link_of_local(l) = tmp_global_link_of_local(l)
enddo

end subroutine set_local_links


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! 計算に必要なfaceの情報を設定
subroutine set_local_faces
use parallel
implicit none

integer f,l,i,k
type(LOCAL_LABEL) :: tmp_send_faces(1:global_num_faces)
type(LOCAL_LABEL) :: tmp_recv_faces(1:global_num_faces)
integer :: nsend,nrecv,info
integer :: tmp_global_face_of_local(1:global_num_faces)

!! 計算を担当するfaceを tmp_global_face_of_local の最初の方に設定
num_necessary_faces=0
do f=1,global_num_faces
  if( local_face_of_global(f)%rank_ == MYRANK ) then
    num_necessary_faces = num_necessary_faces + 1
    tmp_global_face_of_local(local_face_of_global(f)%label_)=f
    !tmp_global_face_of_local(num_necessary_faces)=f
  endif
enddo

nsend=0
nrecv=0
!! 担当するlinkを共有するfaceが全て必要
do l=1,global_num_links
  do k=1,global_face_in_l(l)%num_
    f=global_face_in_l(l)%label_(k)
    if( local_link_of_global(l)%rank_ &
        /= local_face_of_global(f)%rank_ ) then
      if( local_face_of_global(f)%rank_ == MYRANK ) then
        !重複チェック
        info=0
        do i=1,nsend
          if( tmp_send_faces(i)%rank_ == local_link_of_global(l)%rank_ &
              .and. &
              tmp_send_faces(i)%label_ == f ) then 
            info=1
            exit
          endif
        enddo
        if (info == 0 ) then
          nsend=nsend+1
          tmp_send_faces(nsend)%rank_ = local_link_of_global(l)%rank_ 
          tmp_send_faces(nsend)%label_ = f
        endif
      !!!!!!!!!!!
      elseif( local_link_of_global(l)%rank_ == MYRANK ) then
        !! 重複チェック
        info=0
        do i=1,nrecv
          if( tmp_recv_faces(i)%rank_ == local_face_of_global(f)%rank_ &
              .and. &
              tmp_recv_faces(i)%label_ == f ) then 
            info=1
            exit
          endif
        enddo
        if (info == 0 ) then
          nrecv=nrecv+1
          tmp_recv_faces(nrecv)%rank_ = local_face_of_global(f)%rank_ 
          tmp_recv_faces(nrecv)%label_ = f
          tmp_global_face_of_local(num_faces+nrecv)=f
        endif
      endif
    endif
  enddo
enddo


!! 重複度を除いて配列を定義
num_send_faces=nsend
num_recv_faces=nrecv
allocate( send_faces(num_send_faces) )
allocate( recv_faces(num_recv_faces) )
num_necessary_faces = num_faces + num_recv_faces
allocate( global_face_of_local(1:num_necessary_faces) )
!write(*,*) MYRANK,"size",size(global_face_of_local),global_face_of_local(:)

!! 代入
do i=1,num_send_faces
  send_faces(i)%rank_ = tmp_send_faces(i)%rank_
  send_faces(i)%label_ = local_face_of_global(tmp_send_faces(i)%label_)%label_
enddo
do i=1,num_recv_faces
  recv_faces(i)%rank_ = tmp_recv_faces(i)%rank_
  recv_faces(i)%label_ = num_faces+i
enddo
do f=1,num_necessary_faces
  global_face_of_local(f) = tmp_global_face_of_local(f)
  !write(*,*) MYRANK,num_necessary_faces,f,global_face_of_local(f)
  !write(*,*) MYRANK,f,global_face_of_local(f)
enddo
!write(*,*) MYRANK,"##############"
!write(*,*) "check_tmp",MYRANK,size(tmp_global_face_of_local),tmp_global_face_of_local(:)
!write(*,*) "check3",MYRANK,global_face_of_local(:)

end subroutine set_local_faces


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! それぞれのnode上でのlink_tipとlink_orgを設定
subroutine set_local_link_org_and_tip
use parallel
implicit none

integer info_o, info_t
integer l,s

allocate( link_org(1:num_necessary_links) )
allocate( link_tip(1:num_necessary_links) )

do l=1,num_necessary_links
  info_o=0
  info_t=0
  do s=1,num_necessary_sites
    if( global_site_of_local(s) == global_link_org(global_link_of_local(l)) ) then
      link_org(l) = s
      info_o=1
    elseif( global_site_of_local(s) == global_link_tip(global_link_of_local(l)) ) then
      link_tip(l) = s
      info_t=1
    endif
    if( info_o * info_t == 1 ) exit
  enddo
  if( info_o==0 ) link_org(l)=0
  if( info_t==0 ) link_tip(l)=0
  !if( info_o * info_t == 0 ) then 
  !  write(*,*) "something is wrong in local link", l, "in RANK",MYRANK
  !  write(*,*) "  global_link:",global_link_of_local(l),&
  !    "org:",(local_link_org(l)),&
  !    "tip:",(local_link_tip(l))
  !  call stop_for_test
  !endif
enddo


end subroutine set_local_link_org_and_tip



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! それぞれのnode上で、siteに出入りするlinkのラベルを設定
subroutine set_local_link_fromto_site
use parallel
implicit none

integer s,ss,l,num,info,i

!allocate(linktip_from_s(1:num_sites))
!allocate(linkorg_to_s(1:num_sites))
allocate(linktip_from_s(1:num_necessary_sites)) !HERE
allocate(linkorg_to_s(1:num_necessary_sites)) !HERE
!do s=1,num_sites !HERE
do s=1,num_necessary_sites
  !! local_linktip
  num=global_linktip_from_s(global_site_of_local(s))%num_
  linktip_from_s(s)%num_=num
  allocate(linktip_from_s(s)%labels_(1:num) )
  allocate(linktip_from_s(s)%sites_(1:num) )!使わないので設定はしない
  do i=1,num
    info=0
    do l=1,num_necessary_links
      if( global_linktip_from_s(global_site_of_local(s))%labels_(i) &
          == global_link_of_local(l) ) then
        linktip_from_s(s)%labels_(i)=l
        info=1
        exit
      endif
    enddo
    if( info==0 ) then 
      write(*,*) MYRANK,"Something happpend in set_local_linktip_from_s"
      call stop_for_test
    endif
  enddo
  !! local_linkorg
  num=global_linkorg_to_s(global_site_of_local(s))%num_
  linkorg_to_s(s)%num_= num
  allocate(linkorg_to_s(s)%labels_(1:num) )
  allocate(linkorg_to_s(s)%sites_(1:num) ) !使わないので設定はしない
  do i=1,num
    do l=1,num_necessary_links
      info=0
      if( global_linkorg_to_s(global_site_of_local(s))%labels_(i) &
          == global_link_of_local(l) ) then
        linkorg_to_s(s)%labels_(i)=l
        info=1
        exit
      endif
    enddo
    if( info==0 ) then 
      write(*,*) "Something happpend in set_local_linkorg_to_s in link", l, "in RANK ",MYRANK
      call stop_for_test
    endif
  enddo
enddo

end subroutine set_local_link_fromto_site

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! それぞれのnode上で、faceを構成するlinkのラベルを設定
subroutine  set_local_links_in_f
use parallel
implicit none

integer f,l,ll, num, g_link, dir

allocate( links_in_f(1:num_necessary_faces) )
do f=1,num_necessary_faces
  !write(*,*) "test2",MYRANK,num_necessary_faces,f,global_face_of_local(f)
  num = global_links_in_f(global_face_of_local(f))%num_
  links_in_f(f)%num_ = num
  allocate( links_in_f(f)%link_labels_(1:num) )
  allocate( links_in_f(f)%link_dirs_(1:num))

  do l=1,num
    g_link = global_links_in_f(global_face_of_local(f))%link_labels_(l)
    dir = global_links_in_f(global_face_of_local(f))%link_dirs_(l)
    do ll=1,num_necessary_links
      if( global_link_of_local(ll) == g_link ) then
        links_in_f(f)%link_labels_(l) = ll
        links_in_f(f)%link_dirs_(l) = dir
      endif
    enddo
  enddo
enddo

end subroutine  set_local_links_in_f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! それぞれのnode上で、linkに共有されるfaceのラベルを設定
subroutine  set_local_face_in_l
use parallel
implicit none

integer l,num ,i,f,g_face,info

allocate( face_in_l(1:num_links) )
do l=1,num_links
  num = global_face_in_l(global_link_of_local(l))%num_
  face_in_l(l)%num_ = num
  allocate( face_in_l(l)%label_(1:num) )

  do i=1,num
    info=0
    g_face = global_face_in_l(global_link_of_local(l))%label_(i)
    do f=1,num_necessary_faces
      if( global_face_of_local(f) == g_face ) then
        face_in_l(l)%label_(i) = f
        info=1
      endif
    enddo
    if( info == 0 ) then 
      write(*,*) "something happened in setting local_face_in_l."
      write(*,*) MYRANK, f,i
      call stop_for_test
    endif
  enddo
enddo


end subroutine  set_local_face_in_l


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! それぞれのnode上で担当するfaceのlocalな代表点を設定
subroutine set_local_sites_in_f
use parallel
implicit none

integer f,s,global_s,info
integer gf,i
integer Nsites
integer num

allocate( sites_in_f(1:num_necessary_faces) )
do f=1,num_necessary_faces
  gf=global_face_of_local(f)

  Nsites = global_sites_in_f(gf)%num_
  
  sites_in_f(f)%num_ = Nsites
  allocate( sites_in_f(f)%label_(1:Nsites) )

  !! 担当するfaceは全てのサイトを、それ以外は代表点のみ
  if ( f <= num_faces ) then 
    num = Nsites
  else 
    num =1
  endif 

  do i=1, num
    global_s=global_sites_in_f(global_face_of_local(f))%label_(i)
    info=0
    do s=1,num_necessary_sites
      if( global_s == global_site_of_local(s) ) then
        sites_in_f(f)%label_(i) = s
        info=1
        exit
      endif
    enddo
  enddo
  if( info == 0 ) then
    write(*,*) MYRANK, gf, "something happened in set_local_sites_in_f"
    call stop_for_test
  endif
enddo

end subroutine set_local_sites_in_f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_num_local_faces_in_s
!use global_subroutines, only : syncronize_siteval_int
use parallel
implicit none
integer s,global_s

allocate( num_faces_in_s(1:num_necessary_sites) )
do s=1,num_necessary_sites
  global_s = global_site_of_local(s)
  num_faces_in_s(s)=global_face_in_s(global_s)%num_
enddo

!call syncronize_siteval_int(num_faces_in_s)


end subroutine set_num_local_faces_in_s


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! set alpha and beta
subroutine set_local_alpha_beta
use parallel
implicit none
integer s,l,f
integer :: local, rank,  tag

allocate ( alpha_s(1:num_necessary_sites) )
allocate ( alpha_l(1:num_necessary_links) )
allocate ( alpha_f(1:num_necessary_faces) )
allocate ( beta_f(1:num_necessary_faces) )
 
do s=1,global_num_sites
  local=local_site_of_global(s)%label_
  rank=local_site_of_global(s)%rank_
  tag=s
  if( MYRANK == 0 ) then
    if( rank == 0 ) then
      alpha_s(local)=global_alpha_s(s)
    else
      call MPI_SEND(global_alpha_s(s),1,MPI_DOUBLE_PRECISION,rank,tag,MPI_COMM_WORLD,IERR)
    endif
  elseif( rank == MYRANK ) then
    call MPI_RECV(alpha_s(local),1,MPI_DOUBLE_PRECISION, 0,tag,MPI_COMM_WORLD,ISTATUS,IERR)
  endif
enddo

do l=1,global_num_links
  local=local_link_of_global(l)%label_
  rank=local_link_of_global(l)%rank_
  tag=l+global_num_sites
  if( MYRANK == 0 ) then
    if( rank == 0 ) then
      alpha_l(local)=global_alpha_l(l)
    else
      call MPI_SEND(global_alpha_l(l),1,MPI_DOUBLE_PRECISION,rank,tag,MPI_COMM_WORLD,IERR)
    endif
  elseif( rank == MYRANK ) then
    call MPI_RECV(alpha_l(local),1,MPI_DOUBLE_PRECISION, 0,tag,MPI_COMM_WORLD,ISTATUS,IERR)
  endif
enddo

do f=1,global_num_faces
  local=local_face_of_global(f)%label_
  rank=local_face_of_global(f)%rank_
  tag=f+global_num_sites+global_num_links
  if( MYRANK == 0 ) then
    if( rank == 0 ) then
      alpha_f(local)=global_alpha_f(f)
    else
      call MPI_SEND(global_alpha_f(f),1,MPI_DOUBLE_PRECISION,rank,tag,MPI_COMM_WORLD,IERR)
    endif
  elseif( rank == MYRANK ) then
    call MPI_RECV(alpha_f(local),1,MPI_DOUBLE_PRECISION, 0,tag,MPI_COMM_WORLD,ISTATUS,IERR)
  endif
enddo

do f=1,global_num_faces
  local=local_face_of_global(f)%label_
  rank=local_face_of_global(f)%rank_
  tag=f+global_num_sites+global_num_links+global_num_faces
  if( MYRANK == 0 ) then
    if( rank == 0 ) then
      beta_f(local)=global_beta_f(f)
    else
      call MPI_SEND(global_beta_f(f),1,MPI_DOUBLE_PRECISION,rank,tag,MPI_COMM_WORLD,IERR)
    endif
  elseif( rank == MYRANK ) then
    call MPI_RECV(beta_f(local),1,MPI_DOUBLE_PRECISION, 0,tag,MPI_COMM_WORLD,ISTATUS,IERR)
  endif
enddo

call syncronize_ab(alpha_s,'S')
call syncronize_ab(alpha_l,'L')
call syncronize_ab(alpha_f,'F')
call syncronize_ab(beta_f,'F')


end subroutine set_local_alpha_beta


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! set U(1)_R mass 
subroutine set_local_U1Rmass
use parallel
implicit none
integer s,l
!integer :: local, rank,  tag
!double precision :: tmp
!double precision :: tmpr(1:num_necessary_links)
!double precision :: tmpi(1:num_necessary_links)

allocate ( U1Rfactor_link(1:num_necessary_links) )
allocate ( U1Rfactor_site(1:num_necessary_sites) )
allocate ( U1R_ratio(1:num_necessary_links) )

do s=1,num_necessary_sites
  U1Rfactor_site(s)=global_U1Rfactor_site( global_site_of_local(s) )
  !write(*,*) global_site_of_local(s), U1Rfactor_site(s)
enddo

do l=1,num_necessary_links
  U1Rfactor_link(l)=global_U1Rfactor_link( global_link_of_local(l) )
  U1R_ratio(l)=global_U1R_ratio( global_link_of_local(l) )
enddo

 
!tmpr=0d0
!tmpi=0d0
!do l=1,global_num_links
!  local=local_link_of_global(l)%label_
!  rank=local_link_of_global(l)%rank_
!  tag=l+global_num_sites
!  if( MYRANK == 0 ) then
!    if( rank == 0 ) then
!      tmp=global_U1Rmass_phys(l)*LatticeSpacing
!      tmpr(local)=dcos(tmp)
!      tmpi(local)=dsin(tmp)
!      !U1Rfactor(local)=dcmplx(dcos(tmp))+(0d0,-1d0)*dcmplx(dsin(tmp))
!    else
!      call MPI_SEND(global_U1Rmass_phys(l),1,MPI_DOUBLE_PRECISION,rank,tag,MPI_COMM_WORLD,IERR)
!    endif
!  elseif( rank == MYRANK ) then
!    call MPI_RECV(tmp,1,MPI_DOUBLE_PRECISION, 0,tag,MPI_COMM_WORLD,ISTATUS,IERR)
!    tmp=tmp*LatticeSpacing
!    tmpr(local)=dcos(tmp)
!    tmpi(local)=dsin(tmp)
!    !U1Rfactor(local)=dcmplx(dcos(tmp))+(0d0,-1d0)*dcmplx(dsin(tmp))
!  endif
!enddo
!call syncronize_ab(tmpr,'L')
!call syncronize_ab(tmpi,'L')
!
!U1Rfactor=dcmplx( tmpr ) + (0d0,1d0)*tmpi


end subroutine set_local_U1Rmass


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! alpha, betaを通信するsubroutine
subroutine syncronize_ab(alpha,C)
use parallel
double precision, intent(inout) :: alpha(:)
character, intent(in) :: C

integer :: num_send, num_recv
integer :: s_send
integer :: s_recv
integer :: local, rank, tag
integer, allocatable :: ISEND(:), IRECV(:) ! for MPI_WAIT 

!!!!!!!!

if( C=='S' ) then 
  num_send=num_send_sites
  num_recv=num_recv_sites
elseif( C=='L' ) then 
  num_send=num_send_links
  num_recv=num_recv_links
elseif( C=='F' ) then 
  num_send=num_send_faces
  num_recv=num_recv_faces
endif


allocate(ISEND(1:num_send))
allocate(IRECV(1:num_recv))
do s_send=1,num_send
  if( C=='S' ) then
    local=send_sites(s_send)%label_
    rank=send_sites(s_send)%rank_
    tag=10000*rank + global_site_of_local(local)
  elseif( C=='L' ) then
    local=send_links(s_send)%label_
    rank=send_links(s_send)%rank_
    tag=10000*rank + global_link_of_local(local)
  elseif( C=='F' ) then
    local=send_faces(s_send)%label_
    rank=send_faces(s_send)%rank_
    tag=10000*rank + global_face_of_local(local)
  endif

  call MPI_ISEND(alpha(local),1,MPI_DOUBLE_PRECISION,rank,tag,MPI_COMM_WORLD,ISEND(s_send),IERR)
enddo
do s_recv=1,num_recv
  if( C=='S' ) then
    local=recv_sites(s_recv)%label_
    rank=recv_sites(s_recv)%rank_
    tag=10000*MYRANK + global_site_of_local(local)
  elseif( C=='L' ) then
    local=recv_links(s_recv)%label_
    rank=recv_links(s_recv)%rank_
    tag=10000*MYRANK + global_link_of_local(local)
  elseif( C=='F' ) then
    local=recv_faces(s_recv)%label_
    rank=recv_faces(s_recv)%rank_
    tag=10000*MYRANK + global_face_of_local(local)
  endif
  call MPI_IRECV(alpha(local),1,MPI_DOUBLE_PRECISION,rank,tag,MPI_COMM_WORLD,IRECV(s_recv),IERR)
enddo

do s_send=1,num_send
  call MPI_WAIT(ISEND(s_send),ISTATUS,IERR)
enddo
do s_recv=1,num_recv
  call MPI_WAIT(IRECV(s_recv),ISTATUS,IERR)
enddo

deallocate(ISEND, IRECV)
end subroutine syncronize_ab

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine set_U1Rfactor_on_sites
!implicit none
!
!!integer :: U1Rfactor_site(1:global_num_sites)
!integer :: groupS(1:global_num_sites)
!integer :: groupT(1:global_num_sites)
!integer :: groupU(1:global_num_sites)
!
!integer :: s,t,l,ls
!integer :: i,j,k
!integer :: numS,numT,numU
!integer :: info
!
!double precision :: tmp
!
!allocate( global_site_U1Rfactor(1:global_num_sites) )
!allocate( site_U1Rfactor(1:num_necessary_sites) )
!global_site_U1Rfactor=(1d0,0d0)
!site_U1Rfactor=(1d0,0d0)
!
!!! set global_site_U1Rfactor in all ranks
!!!!!!!!!!!!
!groupS=0
!numS=0
!!!!!!!!!!!!
!groupT=0
!numT=1
!groupT(1)=1
!!!!!!!!!!!!
!groupU=0
!numU=0
!do while (numS < global_num_sites) 
!  do i=1,numT
!    s=groupT(i)
!    do j=1,global_linktip_from_s(s)%num_
!      t=global_linktip_from_s(s)%sites_(j)
!      l=global_linktip_from_s(s)%labels_(j)
!
!      info=1
!      do k=1,numS
!        if( t==groupS(k) ) then
!          info=0
!          exit
!        endif
!      enddo
!      do k=1,numT
!        if( t==groupT(k) ) then
!          info=0
!          exit
!        endif
!      enddo
!      do k=1,numU
!        if( t==groupU(k) ) then
!          info=0
!          exit
!        endif
!      enddo
!      if( info==1 ) then
!        numU=numU+1
!        global_site_U1Rfactor(t)=global_site_U1Rfactor(s)*global_U1Rfactor(l)
!        groupU(numU)=t
!      endif
!    enddo
!  enddo
!  !! update groupS
!  groupS(numS+1:numS+numT)=groupT(1:numT)
!  numS=numS+numT
!  !!
!  !write(*,*) numS
!  !! update groupT
!  groupT=groupU
!  numT=numU
!  !! initialize groupU
!  groupU=0
!  numU=0
!  !! 
!enddo
!
!do ls=1,num_necessary_sites
!  site_U1Rfactor(ls)=global_site_U1Rfactor(global_site_of_local(ls))
!enddo
!
!!do s=1,num_necessary_sites
!!  tmp=dble((0d0,-1d0)*cdlog(site_U1Rfactor(s)))/LatticeSpacing
!!  if( tmp < 0 ) then 
!!    write(*,*) global_site_of_local(s), site_U1Rfactor(s), tmp+dacos(-1d0)
!!  else
!!    write(*,*) global_site_of_local(s), site_U1Rfactor(s), tmp
!!  endif
!!  !write(*,*) global_site_of_local(s), site_U1Rfactor(s), &
!!    !dble((0d0,-1d0)*cdlog(site_U1Rfactor(s)))/LatticeSpacing
!!enddo
!
!
!end subroutine set_U1Rfactor_on_sites


#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! local なラベルを出力する
subroutine check_local_sc2
use parallel
implicit none

integer :: rank, turn
integer :: local,global

turn=0
if ( turn .ne. MYRANK ) then
  call MPI_RECV(turn,1,MPI_INTEGER,MYRANK-1,MYRANK,MPI_COMM_WORLD,ISTATUS,IERR)
endif
write(*,*) "#####",MYRANK,"#####"

write(*,*) "#!!!! local data to global data !!!!"
write(*,*) "### site ###"
do local=1,num_necessary_sites
  write(*,*) "local site",local,"->","global",global_site_of_local(local)
enddo
!write(*,*) "### link ###"
!do local=1,num_necessary_links
!  write(*,*) "local link",local,"->","global",global_link_of_local(local)
!enddo
!write(*,*) "### face ###"
!do local=1,num_necessary_faces
!  write(*,*) "local face",local,"->","global",global_face_of_local(local)
!enddo


write(*,*) "#!!!! SEND and RECV !!!!"
write(*,*) "### send site ###"
do local=1,size(send_sites,1)
  write(*,*) "send site",global_site_of_local(send_sites(local)%label_),"(local",send_sites(local)%label_,") to RANK",send_sites(local)%rank_
enddo
!write(*,*) "### send link ###"
!do local=1,size(send_links,1)
!  write(*,*) "send local link",send_links(local)%label_,"to RANK",send_links(local)%rank_
!enddo
!write(*,*) "### send face ###"
!do local=1,size(send_faces,1)
!  write(*,*) "send local face",send_faces(local)%label_,"to RANK",send_faces(local)%rank_
!enddo
write(*,*) "### recv site ###"
do local=1,size(recv_sites,1)
  write(*,*) "recv site",global_site_of_local(recv_sites(local)%label_),"(local",recv_sites(local)%label_,") from RANK",recv_sites(local)%rank_
enddo
!write(*,*) "### recv link ###"
!do local=1,size(recv_links,1)
!  write(*,*) "recv local link",recv_links(local)%label_,"from RANK",recv_links(local)%rank_
!enddo
!write(*,*) "### recv face ###"
!do local=1,size(recv_faces,1)
!  write(*,*) "recv local face",recv_faces(local)%label_,"from RANK",recv_faces(local)%rank_
!enddo

!write(*,*) "#!!!! SITE-LINK  !!!!"
!do local = 1,num_links
!  write(*,*) "local link",local,"is",link_org(local),link_tip(local)
!enddo
!
!write(*,*) "###"
!do local=1,num_sites
!  write(*,*) "local links from",local,";",linktip_from_s(local)%labels_
!enddo
!write(*,*) "###"
!do local=1,num_sites
!  write(*,*) "local links to",local,";",linkorg_to_s(local)%labels_
!enddo


write(*,*)
write(*,*)
turn=MYRANK+1
if( MYRANK .ne. NPROCS-1 ) then
  call MPI_SEND(turn,1,MPI_INTEGER,MYRANK+1,MYRANK+1,MPI_COMM_WORLD,IERR)
endif

end subroutine check_local_sc2
