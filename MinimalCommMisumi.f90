module global
use simplicial_complex
implicit none


type SITE_LINKDAT
  integer :: num_
  integer, allocatable :: labels_(:)
  integer, allocatable :: sites_(:)
end type SITE_LINKDAT

type FACEDAT
  integer :: num_ ! 要素の数
  integer, allocatable :: label_(:) ! faceのラベル label_(1:num_)
end type FACEDAT

type LINKDAT
  integer :: num_ ! linkの数
  integer, allocatable :: link_labels_(:) ! linkのラベル
  integer, allocatable :: link_dirs_(:) ! linkの方向
end type LINKDAT

type SITEDAT
  integer :: num_ ! 要素の数
  integer, allocatable :: label_(:) ! siteのラベル label_(1:num_)
end type SITEDAT

type SITE_DIST
  integer :: rank_
  integer :: num_sites_
  integer, allocatable :: site_list_(:)
end type SITE_DIST

character(20) :: SC_FILE_NAME
integer, parameter :: SC_FILE=10

!!!!!!!!!!!!!!!!!!! global !!!!!!!!!!!!!!!!!!!!!
integer, save :: NPROCS
type(SmpCom), save :: SC
integer global_num_sites, global_num_links, global_num_faces
integer, allocatable, save :: global_link_org(:) ! l番目のリンクのorigin
integer, allocatable, save :: global_link_tip(:) ! l番目のリンクのtop
type(SITE_LINKDAT), allocatable, save :: global_linktip_from_s(:) ! sから出るリンクのtip
type(SITE_LINKDAT), allocatable, save :: global_linkorg_to_s(:) ! sに入るリンクのorigin
type(LINKDAT), allocatable, save :: global_links_in_f(:) ! fに含まれるリンクのデータ
type(FACEDAT), allocatable, save :: global_face_in_l(:) ! lに含まれるface
type(SITEDAT), allocatable, save :: global_sites_in_f(:) ! face f に含まれるサイト

!!!!!!!!!!!!!!!!!!! local !!!!!!!!!!!!!!!!!!!!!
type LOCAL_LABEL
  integer :: rank_ ! nodes
  integer :: label_ ! label
end type LOCAL_LABEL

integer,allocatable,save :: num_sites(:),num_links(:),num_faces(:)

type(LOCAL_LABEL), allocatable, save :: local_site_of_global(:) ! map from global site to local site
type(LOCAL_LABEL), allocatable, save :: local_link_of_global(:) ! map from global link to local link
type(LOCAL_LABEL), allocatable, save :: local_face_of_global(:) ! map from global face to local face


integer, save :: num_send_sites, num_recv_sites
integer, save :: num_send_links, num_recv_links
integer, save :: num_send_faces, num_recv_faces

type(LOCAL_LABEL), allocatable, save :: send_sites(:) 
type(LOCAL_LABEL), allocatable, save :: send_links(:) 
type(LOCAL_LABEL), allocatable, save :: send_faces(:) 

type(LOCAL_LABEL), allocatable, save :: recv_sites(:) 
type(LOCAL_LABEL), allocatable, save :: recv_links(:) 
type(LOCAL_LABEL), allocatable, save :: recv_faces(:) 


!!!!!!!!!!!!!!!!!1
!type YoungTableau
  !integer :: num_row
  !integer,allocatable :: num_box(:)
!end type YoungTableau


end module global


!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Young図を扱うモジュール
module YoungTableau
implicit none

type array
  !integer, allocatable :: len
  integer, allocatable :: arr(:)
end type array

!! Young図。枠だけ
type YT
  integer :: num_row
  integer,allocatable :: num_box(:)
end type YT

!! Young Diagram。枠に数字が入ったもの
type YD
  integer :: num_row
  type(array), allocatable :: row(:)
end type YD

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 構造体arrayを長さ len に設定し、
!! 中身を0に初期化する。 
subroutine init_array(arr,len)
implicit none

type(array), intent(out) :: arr
integer, intent(in) :: len

!arr%len=len
if( allocated(arr%arr) ) deallocate( arr%arr )
allocate( arr%arr(1:len) )

!! 0で初期化
arr%arr=0

end subroutine init_array


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 構造体 YT の列数を num に設定する
subroutine init_YT(T,num)
implicit none

type(YT), intent(inout) :: T
integer, intent(in) :: num

T%num_row=num
if( allocated( T%num_box ) ) deallocate( T%num_box )
allocate( T%num_box(1:num) )
end subroutine init_YT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ヤング図（構造体 YD）の枠を、ヤング盤（構造体YT） T の形に合わせる
subroutine init_YD(D,T)
implicit none

type(YD), intent(out) :: D
type(YT), intent(in) :: T
integer :: r 

D%num_row=T%num_row
allocate( D%row(1:T%num_row) )

do r=1,T%num_row
  call init_array(D%row(r),T%num_box(r))
enddo

end subroutine init_YD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 構造体 YT をコピーする
subroutine copy_YT(Tout,Tin)
implicit none

type(YT), intent(in) :: Tin
type(YT), intent(out) :: Tout
integer i

call init_YT(Tout,Tin%num_row)

do i=1,Tin%num_row
  Tout%num_box(i)=Tin%num_box(i)
enddo

end subroutine copy_YT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 構造体 YT のlistに新しいものを追加する
subroutine add_YTlist(YTlist,YTlist_add)
implicit none

type(YT), allocatable, intent(inout) :: YTlist(:)
type(YT), intent(in) :: YTlist_add(:)

integer :: num,num_add
type(YT), allocatable :: tmp_YTlist(:)
integer :: i 

if( allocated(YTlist) ) then 
  num=size(YTlist)
else
  write(*,*) "allocate YTlist"
  stop
endif
num_add=size(YTlist_add)

allocate( tmp_YTlist(1:num) )
do i=1,num
  call copy_YT( tmp_YTlist(i), YTlist(i) )
enddo

deallocate( YTlist )
allocate( YTlist(1:num+num_add) )

do i=1,num
  call copy_YT( YTlist(i), tmp_YTlist(i) )
enddo
do i=1,num_add
  call copy_YT( YTlist(num+i), YTlist_add(i) )
enddo


end subroutine add_YTlist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 構造体 YT1とYT2を比較し、同じなら0, 違えば1を返す関数
integer function compare_YT(YT1, YT2) 
implicit none

type(YT), intent(in) :: YT1, YT2
integer i

compare_YT=0
if( YT1%num_row /= YT2%num_row ) then
  compare_YT=1
  return
endif
do i=1,YT1%num_row
  if(YT1%num_box(i) /= YT2%num_box(i) ) then
    compare_YT=1
    return
  endif
enddo

end  function compare_YT





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_possible_partitions(YT_list,num_total,num_div)
implicit none

type(YT), intent(inout),allocatable :: YT_list(:)
integer, intent(in) :: num_total,num_div
!type(YT), allocatable :: tmp_YT_list(:)!, YT_add(:)
!type(YT) :: tmptab
type(YT) :: tmptab(1)

integer :: num_tableau
integer :: i,j,k,l,info,same
integer :: num_scan0,num_scan1

if( num_total < num_div ) then 
  write(*,*) "there is no possible partition"
  stop
endif

! first step
num_tableau=1
if( allocated(YT_list) ) deallocate(YT_list)
allocate( YT_list(1:1))
call init_YT(YT_list(1),num_div)
YT_list(1)%num_box(1)=num_total-num_div+1
do k=2,num_div
  YT_list(1)%num_box(k)=1
enddo

! next step
info=0
num_scan0=1
num_scan1=1
do while(info==0)
  !! count the add list
  info=1
  do k=num_scan0, num_scan1
    do j=1,num_div-1
      if( YT_list(k)%num_box(j)-1 >= YT_list(k)%num_box(j+1)+1 ) then 
        call copy_YT(tmptab(1), YT_list(k))
        tmptab(1)%num_box(j)= tmptab(1)%num_box(j)-1
        tmptab(1)%num_box(j+1)= tmptab(1)%num_box(j+1)+1
        !!!!!
        same=1
        do l=1,size(YT_list)
          same=same * compare_YT( tmptab(1), YT_list(l) )
         if( same==0 ) exit
        enddo
        if( same == 1 ) then 
          info=0
          call add_YTlist( YT_list, tmptab )
        endif
      endif
      if( YT_list(k)%num_box(j)-1 >= YT_list(k)%num_box(j+1) ) then 
        do i=j+2,num_div
          if( YT_list(k)%num_box(i)+1 <= YT_list(k)%num_box(i-1) ) then
            call copy_YT(tmptab(1), YT_list(k))
            tmptab(1)%num_box(j)= tmptab(1)%num_box(j)-1
            tmptab(1)%num_box(i)= tmptab(1)%num_box(i)+1
            !!!!!
            same=1
            do l=1,size(YT_list)
              same=same * compare_YT( tmptab(1), YT_list(l) )
              if( same==0 ) exit
            enddo
            if( same == 1 ) then
              info=0
              call add_YTlist( YT_list, tmptab )
            endif
          endif
        enddo
      endif
    enddo
  enddo
  num_scan0=num_scan1+1
  num_scan1=size(YT_list)
enddo



           
end subroutine get_possible_partitions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! N-boxのヤング盤に、自然数1～Nを、
!! 各列が単調増加するように配置させる。
subroutine get_diagrams(YD_list,T)
implicit none

type(YD), allocatable,intent(out) :: YD_list(:)
type(YT), intent(in) :: T

integer :: num_YD
integer :: num_allbox
type(YD), allocatable:: tmp_YD_list(:)

integer :: r,i,j,k,l
integer :: counter

num_allbox=0
do r=1,T%num_row
  num_allbox=num_allbox+T%num_box(r)
enddo
!write(*,*) num_allbox

!! まず1を入れる
num_YD=T%num_row
allocate(YD_list(1:num_YD))
do r=1,T%num_row
  call init_YD(YD_list(r),T) !この段階で、rowは全てゼロ
  YD_list(r)%row(r)%arr(1)=1
enddo
allocate( tmp_YD_list(1:num_YD) )

do j=2,num_allbox
  deallocate( tmp_YD_list )
  allocate( tmp_YD_list(1:num_YD) )
  tmp_YD_list=YD_list
  deallocate( YD_list )
  !! count the number of the next YD_list
  num_YD=0
  do i=1,size(tmp_YD_list)
    do r=1,T%num_row
      do l=1,T%num_box(r)
        if( tmp_YD_list(i)%row(r)%arr(l) == 0 ) then
          num_YD=num_YD+1
          go to 10
        endif
      enddo
      10 continue
    enddo
  enddo
  !! make the diagrams
  allocate( YD_list(1:num_YD) )
  counter=0
  do i=1,size(tmp_YD_list)
    do r=1,T%num_row
      do l=1,T%num_box(r)
        if( tmp_YD_list(i)%row(r)%arr(l) == 0 ) then
          counter=counter+1
          YD_list(counter)=tmp_YD_list(i)
          YD_list(counter)%row(r)%arr(l)=j
          go to 20
        endif
      enddo
      20 continue
    enddo
  enddo
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!!!!!! write down all the Young diagrams
!do i=1,num_YD
!  write(*,*) "### diagram",i,"######"
!  do r=1,T%num_row
!    write(*,*) YD_list(i)%row(r)%arr
!  enddo
!enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
end subroutine get_diagrams

end module YoungTableau




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Program to determine such discretization which minimize 
!! the number of communication in the Misumi-discretization
!! ＜設計＞
!! ・単体的複体の情報と分割数をインプット
!! ・可能な分割の方法を総当りし、通信の数が最小になるような分割を決定
!! ・結果をinputfileの形式に合わせて標準出力する。
program main
use global
use simplicial_complex
use YoungTableau
implicit none


type(YT), allocatable :: P_list(:)
type(YD) :: D1,D2
type(YD), allocatable :: YD_list(:)
type(YD), allocatable :: YD_MIN(:)

integer :: num_part
integer s,l,f,i,j,k

integer :: tmp_num, diff

integer :: NUM_COMM
integer, allocatable :: NUM_COMM_MIN(:),num_max_calc(:)
integer :: NUMDIV=20
SC_FILE_NAME="S2MisumiM20N20R1.dat"
!integer :: NUMDIV=3
!SC_FILE_NAME="S2Misumi.dat"

call set_global_sc

!call get_possible_partitions(P_list,global_num_sites,NUMDIV)
!num_part=size(P_list)

num_part=1
allocate(P_list(1:1))
call init_YT(P_list(1),NUMDIV)
tmp_num = global_num_sites/NUMDIV
diff = global_num_sites - tmp_num*NUMDIV
do i=1,NUMDIV
  P_list(1)%num_box(i)=tmp_num
enddo
if( diff > 0 ) then
  do i=2,NUMDIV
    P_list(1)%num_box(i)=P_list(1)%num_box(i)+1
    diff=diff-1
    if( diff==0 ) exit
  enddo
endif


!write(*,*) num_part
allocate(YD_MIN(1:num_part))
allocate(NUM_COMM_MIN(1:num_part),num_max_calc(1:num_part))

do k=1,num_part
  !write(*,*) P_list(k)%num_box(NUMDIV), P_list(k)%num_box(1)
  write(*,*) P_list(k)%num_box
enddo

NUM_COMM_MIN=0
do k=1,num_part
  if( -P_list(k)%num_box(NUMDIV)+P_list(k)%num_box(1) < 2 ) then 
    write(*,*) "#start"
    call get_diagrams(YD_list,P_list(k))
    write(*,*) "#(YD)=",size(YD_list)
    call init_YD(YD_MIN(k),P_list(k))
    do j=1,size(YD_list)
      write(*,*) "## j=",j
      call count_num_comm(NUM_COMM,num_max_calc(k),YD_list(j))
      if(j==1) then 
        NUM_COMM_MIN(k)=NUM_COMM
        YD_MIN(k)=YD_list(1)
      endif
      if( NUM_COMM_MIN(k) > NUM_COMM ) then
        NUM_COMM_MIN(k)=NUM_COMM
        !call COPY_YD(YD_MIN(k),YD_list(j))
        YD_MIN(k)=YD_list(j)
      endif
   enddo
 endif
enddo

do k=1,num_part
  if( -P_list(k)%num_box(NUMDIV)+P_list(k)%num_box(1) < 2 ) then 
    write(*,*) "#########################################"
    write(*,*) "#########################################"
    write(*,*) "# #(calculation)=",num_max_calc(k)
    write(*,*) "# #(communication)=",num_comm_min(k)
    write(*,*) "##"
    call writedown_inputfile(YD_MIN(k))
  endif
enddo

end program main

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine writedown_inputfile(Y)
use global
use YoungTableau
implicit none

type(YD), intent(in) :: Y

integer :: rank

write(*,*) "############################"
do rank=1,Y%num_row
  write(*,'(a,I3,a)') "## Rank ",rank,", number of sites,"
  write(*,'(I3,a,I3)') rank," ",size(Y%row(rank)%arr(:))
  write(*,'(a)') "## list of sites"
  write(*,*) Y%row(rank)%arr(:)
enddo



end subroutine writedown_inputfile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_global_sc
use global
implicit none

integer fsize
integer, allocatable :: fsites(:), faces_l(:), sites_f(:), sites(:)
double precision LatticeSpacing,alpha,beta
integer origin,tip
integer s,l,f,j,k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Read SC_FILE_NAME and set SC
!! Initialize
open(SC_FILE, file=SC_FILE_NAME, status='old',action='READ')
read(SC_FILE,'()') ! skip 1 line
read(SC_FILE,*) global_num_sites
read(SC_FILE,*) global_num_links
read(SC_FILE,*) global_num_faces
read(SC_FILE,*) LatticeSpacing
read(SC_FILE,'()') ! skip 1 line
call init_smpcom(SC,global_num_sites,global_num_links,global_num_faces)

!!!!!!!!!!!!!!!!!!!!!!!!!
!! Set alpha_s(DUMMY)
do k=1,global_num_sites
  read(SC_FILE,*) s,alpha
enddo
read(SC_FILE,'()') ! skip 1 line

!!!!!!!!!!!!!!!!!!!!!!!!!
!! Set global links
do k=1,global_num_links
  read(SC_FILE,*) l,origin,tip,alpha
  call put_link_sc(SC,l,origin,tip)
enddo
read(SC_FILE,'()') ! skip 1 line

!!!!!!!!!!!!!!!!!!!!!!!!!
!! Set global faces
do k=1,global_num_faces
  read(SC_FILE,"(I6,X,I6,X)",advance="no") f,fsize
  allocate( fsites(1:fsize ) )

  do j=1,fsize
    read(SC_FILE,'(I6,X)',advance='no') fsites(j) 
  enddo
  read(SC_FILE,*) alpha,beta
  call put_face_sc(SC,f,fsites)

  deallocate( fsites )
enddo
close(SC_FILE)

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine set_global_sc



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine count_num_comm(num_comm,num_max_calc,local_site_list)
use global
use simplicial_complex
use YoungTableau
implicit none

integer, intent(out) :: num_comm,num_max_calc
type(YD), intent(in) :: local_site_list

integer :: num_comm_sites, num_comm_links, num_comm_faces

num_comm=0
num_comm_faces=0
!! global {site,link,face} から local {site,link,face}へのmapを作成
!!   global_num_{sites,links,faces} 
!!   num_sites, num_links, num_faces
call set_global_to_local(local_site_list,num_max_calc)
!! 各ノードが計算に必要なsiteの情報と、site情報の送り先、受信元を設定する
!!   num_necessary_sites
!!   global_site_of_local(s)
!!   send_sites(s)
!!   recv_sites(s)
call set_local_sites(num_comm_sites)
!! 各ノードが計算に必要なlinkの情報と、link情報の送り先、受信元を設定する
!!   num_necessary_links
!!   global_link_of_local(l)
!!   send_links(l)
!!   recv_links(l)
call set_local_links(num_comm_links)
!! 各ノードが計算に必要なfaceの情報を設定する
!!   num_necessary_faces
!!   global_face_of_local(l)
call set_local_faces(num_comm_faces)

num_comm=num_comm_sites+num_comm_links+num_comm_faces


end subroutine count_num_comm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! global {site,link,face} から local {site,link,face}へのmapを作成
subroutine set_global_to_local(local_site_list,num_max_calc)
use global
use YoungTableau
implicit none

type(YD), intent(in) :: local_site_list
integer, intent(out) :: num_max_calc
integer :: tmp_num_sites, tmp_num_links, tmp_num_faces
!integer,allocatable :: num_calc_site(:),num_calc_link(:),num_calc_face(:)
integer :: s,l,f,rank,i,tmp

NPROCS=local_site_list%num_row
!allocate(num_calc_site(0:PROCS-1))
!allocate(num_calc_link(0:PROCS-1))
!allocate(num_calc_face(0:PROCS-1))
num_max_calc=0
!num_calc_site=0
!num_calc_link=0
!num_calc_face=0


!write(*,*) size(local_site_of_global) 
!write(*,*) global_num_sites
if(.not. allocated(local_site_of_global)) allocate(local_site_of_global(1:NPROCS*global_num_sites))
if(.not. allocated(local_link_of_global)) allocate(local_link_of_global(1:NPROCS*global_num_links))
if(.not. allocated(local_face_of_global)) allocate(local_face_of_global(1:NPROCS*global_num_faces))

if(.not. allocated(num_sites)) allocate(num_sites(0:NPROCS-1))
if(.not. allocated(num_links)) allocate(num_links(0:NPROCS-1))
if(.not. allocated(num_faces)) allocate(num_faces(0:NPROCS-1))

!!! initialization
do i=1,NPROCS*global_num_sites
  local_site_of_global(i)%rank_=0
  local_site_of_global(i)%label_=0
enddo
do i=1,NPROCS*global_num_links
  local_link_of_global(i)%rank_=0
  local_link_of_global(i)%label_=0
enddo
do i=1,NPROCS*global_num_faces
  local_face_of_global(i)%rank_=0
  local_face_of_global(i)%label_=0
enddo
num_sites=0
num_links=0
num_faces=0


!! local_site_of_global(s)
do rank=0,NPROCS-1
  num_sites(rank)=size(local_site_list%row(rank+1)%arr)
  do i=1,num_sites(rank)
    s=local_site_list%row(rank+1)%arr(i)
    local_site_of_global(s)%rank_=rank
    local_site_of_global(s)%label_=i
  enddo
enddo

!! local_link_of_global(l)
do rank=0,NPROCS-1
  num_links(rank)=0
  do l=1,global_num_links
    do i=1,num_sites(rank)
      s=local_site_list%row(rank+1)%arr(i)
      if( s == global_link_org(l) ) then
        num_links(rank)=num_links(rank)+1
        local_link_of_global(l)%rank_=rank
        local_link_of_global(l)%label_=num_links(rank)
      endif
    enddo
  enddo
enddo


!! local_face_of_global(f)
do rank=0,NPROCS-1
  num_faces(rank)=0
  do f=1,global_num_faces
    do i=1,num_sites(rank)
      s=local_site_list%row(rank+1)%arr(i)
      if( s == global_sites_in_f(f)%label_(1) ) then
        num_faces(rank)=num_faces(rank)+1
        local_face_of_global(f)%rank_=rank
        local_face_of_global(f)%label_=num_faces(rank)
      endif
    enddo
  enddo
enddo

num_max_calc=num_sites(0)+num_links(0)+num_faces(0)
do rank=1,NPROCS-1
  tmp=num_sites(rank)+num_links(rank)+num_faces(rank)
  if(tmp>num_max_calc) num_max_calc=tmp
enddo



end subroutine set_global_to_local

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! 計算に必要なsiteの情報を設定
subroutine set_local_sites(num_comm_sites)
use global
use YoungTableau
implicit none

integer, intent(out) :: num_comm_sites

integer :: tmp_num_comm_sites(0:NPROCS-1)
integer s,l,i,info,k,rank
type(LOCAL_LABEL) :: tmp_send_sites(1:NPROCS*global_num_sites,0:NPROCS-1)

!! この段階で、num_necesarry_sites=num_sites
!! linkの始点と終点が別のnodeに属していたら、
!! tipの位置を、orgを担当するnodeに送る必要がある
!! (1) 担当するlinkの始点と終点
!! (2) 担当するサイトを始点とするlinkの終点
!! (3) 担当するサイトを終点とするlinkの始点
!! 重複があり得るので、逐一チェックしながら進める

tmp_num_comm_sites=0
do rank=0,NPROCS-1
  do i=1,global_num_sites
    tmp_send_sites(i,rank)%rank_=0
    tmp_send_sites(i,rank)%label_=0
  enddo
enddo
!! (1) 担当するlinkの始点と終点
do rank=0,NPROCS-1
  do l=1,global_num_links
    if( local_site_of_global(global_link_org(l))%rank_ &
          /= local_site_of_global(global_link_tip(l))%rank_ ) then
    !!!!!!!!!!!
      if( local_site_of_global(global_link_tip(l))%rank_ == rank ) then
        !! 重複チェック
        info=0
        do i=1,tmp_num_comm_sites(rank)
          if( tmp_send_sites(i,rank)%rank_ == local_site_of_global(global_link_org(l))%rank_ &
            .and. &
              tmp_send_sites(i,rank)%label_ == global_link_tip(l) ) then
            info=1
            exit
          endif
        enddo
        if (info==0) then 
          tmp_num_comm_sites(rank)=tmp_num_comm_sites(rank)+1
          tmp_send_sites(tmp_num_comm_sites(rank),rank)%rank_ = local_site_of_global(global_link_org(l))%rank_ 
          tmp_send_sites(tmp_num_comm_sites(rank),rank)%label_ = global_link_tip(l) 
        endif
      endif
    endif
  enddo
enddo
!! (2) 担当するサイトを終点とするlinkの始点
do rank=0,NPROCS-1
  do s=1,global_num_sites
    do k=1,global_linkorg_to_s(s)%num_
      if( local_site_of_global(global_linkorg_to_s(s)%sites_(k))%rank_ &
          /= local_site_of_global(s)%rank_ ) then
        if(local_site_of_global(global_linkorg_to_s(s)%sites_(k))%rank_ == rank) then
          !! 重複チェック
          info=0
          do i=1,tmp_num_comm_sites(rank)
            if( tmp_send_sites(i,rank)%rank_ &
                == local_site_of_global(s)%rank_ &
              .and. &
                tmp_send_sites(i,rank)%label_  &
                == global_linkorg_to_s(s)%sites_(k) &
                ) then
              info=1
              exit
            endif
          enddo
          if (info==0) then 
            tmp_num_comm_sites(rank)=tmp_num_comm_sites(rank)+1
            tmp_send_sites(tmp_num_comm_sites(rank),rank)%rank_ = local_site_of_global(s)%rank_ 
            tmp_send_sites(tmp_num_comm_sites(rank),rank)%label_ = global_linkorg_to_s(s)%sites_(k)
          endif
        endif
      endif
    enddo
  enddo
enddo
!! (3) 担当するサイトを始点とするlinkの終点
do rank=0,NPROCS-1
  do s=1,global_num_sites
  !write(*,*) s,global_linktip_from_s(s)%num_
    do k=1,global_linktip_from_s(s)%num_
      if( local_site_of_global(global_linktip_from_s(s)%sites_(k))%rank_  &
          /= local_site_of_global(s)%rank_ ) then
        if(local_site_of_global(global_linktip_from_s(s)%sites_(k))%rank_ == rank) then
        !! 重複チェック
          info=0
          do i=1,tmp_num_comm_sites(rank)
            if( tmp_send_sites(i,rank)%rank_ &
                == local_site_of_global(s)%rank_ &
              .and. &
                tmp_send_sites(i,rank)%label_  &
                == global_linktip_from_s(s)%sites_(k)) then
              info=1
              exit
            endif
          enddo
          if (info==0) then 
            tmp_num_comm_sites(rank)=tmp_num_comm_sites(rank)+1
            tmp_send_sites(tmp_num_comm_sites(rank),rank)%rank_ = local_site_of_global(s)%rank_ 
            tmp_send_sites(tmp_num_comm_sites(rank),rank)%label_ = global_linktip_from_s(s)%sites_(k)
          endif
        endif
      endif
    enddo
  enddo
enddo
num_comm_sites=0
do rank=0,NPROCS-1
  num_comm_sites=num_comm_sites+tmp_num_comm_sites(rank)
enddo


end subroutine set_local_sites


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! 計算に必要なlinkの情報を設定
subroutine set_local_links(num_comm_links)
use global
implicit none

integer, intent(out) :: num_comm_links 
integer :: tmp_num_comm_links(0:NPROCS-1)

integer s,l,f,i,j,k,ll,ll_label,rank
type(LOCAL_LABEL) :: tmp_send_links(1:NPROCS*global_num_links,0:NPROCS-1)
integer :: info


!! この段階で、num_necesarry_links=num_links
!! 計算に必要なのは、
!!   (1) 担当するsiteをtipとするlink
!!   (2) 担当するfaceを構成するlink
!!   (3) 担当するlinkを共有するfaceに含まれる全link
!! 重複があり得るので、逐一チェックしながら進める

tmp_num_comm_links=0
do rank=0,NPROCS-1
  do i=1,global_num_links
    tmp_send_links(i,rank)%rank_=0
    tmp_send_links(i,rank)%label_=0
  enddo
enddo
!!   (1) 担当するsiteをtipとするlink
do rank=0,NPROCS-1
  do l=1,global_num_links
    if( local_site_of_global(global_link_tip(l))%rank_ /= &
        local_site_of_global(global_link_org(l))%rank_ ) then
      if( local_site_of_global(global_link_org(l))%rank_ == rank ) then
        !! 重複チェック
        info=0
        do i=1,tmp_num_comm_links(rank)
          if( tmp_send_links(i,rank)%rank_ == local_site_of_global(global_link_tip(l))%rank_ &
            .and. &
              tmp_send_links(i,rank)%label_ == l ) then
            info=1
            exit
          endif
        enddo
        if( info == 0 ) then
          tmp_num_comm_links(rank)=tmp_num_comm_links(rank)+1
          tmp_send_links(tmp_num_comm_links(rank),rank)%rank_ = local_site_of_global(global_link_tip(l))%rank_
          tmp_send_links(tmp_num_comm_links(rank),rank)%label_ = l
        endif
      endif
    endif
  enddo
enddo

!!   (2) 担当するfaceを構成するlink
do rank=0,NPROCS-1
  do f=1,global_num_faces
    do j=1,global_links_in_f(f)%num_
      l=global_links_in_f(f)%link_labels_(j)
      if( local_face_of_global(f)%rank_ /= local_link_of_global(l)%rank_ ) then
        if( local_link_of_global(l)%rank_ == rank ) then
          !! 重複チェック
          info=0
          do i=1,tmp_num_comm_links(rank)
            if( tmp_send_links(i,rank)%rank_ == local_face_of_global(f)%rank_ &
              .and. &
                tmp_send_links(i,rank)%label_ == l ) then
              info=1
              exit
            endif
          enddo
          if( info == 0 ) then
            tmp_num_comm_links(rank)=tmp_num_comm_links(rank)+1
            tmp_send_links(tmp_num_comm_links(rank),rank)%rank_ = local_face_of_global(f)%rank_ 
            tmp_send_links(tmp_num_comm_links(rank),rank)%label_ = l
          endif
        endif
      endif
    enddo
  enddo
enddo

!!   (3) 担当するlinkを共有するfaceに含まれる全link
do rank=0,NPROCS-1
  do l=1,global_num_links
    do k=1,global_face_in_l(l)%num_
      f=global_face_in_l(l)%label_(k)
      do ll_label=1,global_links_in_f(f)%num_
        ll=global_links_in_f(f)%link_labels_(ll_label)
        if( local_link_of_global(l)%rank_ /= local_link_of_global(ll)%rank_ ) then 
          if( local_link_of_global(ll)%rank_ == rank ) then
            !! 重複チェック
            info=0
            do i=1,tmp_num_comm_links(rank)
              if( tmp_send_links(i,rank)%rank_ == local_link_of_global(l)%rank_ &
                  .and.  &
                  tmp_send_links(i,rank)%label_ == ll ) then
                info=1
                exit
              endif
            enddo
            if (info == 0 ) then
              tmp_num_comm_links(rank)=tmp_num_comm_links(rank)+1
              tmp_send_links(tmp_num_comm_links(rank),rank)%rank_ = local_link_of_global(l)%rank_ 
              tmp_send_links(tmp_num_comm_links(rank),rank)%label_ = ll
            endif
          endif
        endif
      enddo
    enddo
  enddo
enddo
num_comm_links=0
do rank=0,NPROCS-1
  num_comm_links=num_comm_links+tmp_num_comm_links(rank)
enddo



end subroutine set_local_links


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! 計算に必要なfaceの情報を設定
subroutine set_local_faces(num_comm_faces)
use global
implicit none

integer, intent(out) :: num_comm_faces 
integer :: tmp_num_comm_faces(0:NPROCS-1)

integer s,l,f,i,j,k,ll,ll_label,rank
type(LOCAL_LABEL) :: tmp_send_faces(1:NPROCS*global_num_faces,0:NPROCS-1)
integer :: info


tmp_num_comm_faces=0
!! 担当するlinkを共有するfaceが全て必要
do rank=0,NPROCS-1
  do l=1,global_num_links
    do k=1,global_face_in_l(l)%num_
      f=global_face_in_l(l)%label_(k)
      if( local_link_of_global(l)%rank_ &
          /= local_face_of_global(f)%rank_ ) then
        if( local_face_of_global(f)%rank_ == rank) then
          !重複チェック
          info=0
          do i=1,tmp_num_comm_faces(rank)
            if( tmp_send_faces(i,rank)%rank_ == local_link_of_global(l)%rank_ &
                .and. &
                tmp_send_faces(i,rank)%label_ == f ) then 
              info=1
              exit
            endif
          enddo
          if (info == 0 ) then
            tmp_num_comm_faces(rank)=tmp_num_comm_faces(rank)+1
            tmp_send_faces(tmp_num_comm_faces(rank),rank)%rank_ = local_link_of_global(l)%rank_ 
            tmp_send_faces(tmp_num_comm_faces(rank),rank)%label_ = f
          endif
        endif
      endif
    enddo
  enddo
enddo

num_comm_faces=0
do rank=0,NPROCS-1
  num_comm_faces=num_comm_faces+tmp_num_comm_faces(rank)
enddo

end subroutine set_local_faces



