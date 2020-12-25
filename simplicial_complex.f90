!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 単体的複体を表すmodule
!! 
!!  ver.01: linkの方向を表す数字を
!!    1: 順方向 -1: 逆方向 0: リンクなし
!!   に変更した
!!  ver.02: SmpComの構造を大幅に変更
!! 
!!  マニュアル
!!##############
!! 構造型 SmpCom （単体的複体）
!!  integer :: NumSites_  ! 複体を構成するサイトの数
!!  integer :: NumLinks_  ! 複体を構成するリンクの数
!!  integer :: NumFaces_  ! 複体を構成するfaceの数
!!  type(Link), allocatable :: LinkSet_(:)  ! 全linkの集合
!!  type(face), allocatable :: FaceSet_(:) ! 全faceの集合 
!!##############
!! init_smpcom(SC,Ns,Nl,Nf)
!!    simplical complex SC を初期化する。
!!    Ns: NumSites_
!!    Nl: NumLinks_
!!    Nf: NumFaces_
!!    allocate: LinkSet_(1:NumLinks_),  FaceSet_(1:NumFaces_)
!!    
!! put_link_sc(sc,l,origin,tip)
!!      SmpCom:scのl番目のlinkを<origin,tip>にセットする。
!!        1 <= l <= NumLinks_
!!
!! put_face_sc(sc,f,sites)
!!      SmpCom:scのf番目のfaceを設定する。
!!      ただし、事前にすべてのlinkをset_link_scで設定しておく必要がある。
!!        1 <= f <= NumFaces_
!!        sites: 大きさが FaceSize の整数配列
!!
!! get_numsite_sc(sc,numsite)
!!     scを構成する全サイト数を取得
!!
!! get_numlink_sc(sc,numlink)
!!     scを構成する全リンク数を取得
!!
!! get_numface_sc(sc,numface)
!!     scを構成する全face数を取得
!!
!! get_link_sc(sc,l,origin,tip)
!!     scのl番目のlinkのoriginとtipを取得
!!
!! get_linklabel_sc(sc,origin,tip,link_label,link_dir)
!!     scの中のlink:<origin,tip>のラベルと方向を取得
!!     もし該当するリンクがなければ 0 を返す
!!
!! get_links_from_s_sc(sc,s,link_labels,tips)
!!     sc のs番目のサイトから出ているリンクのラベルと行き先を取得
!!      link_labels : allocatable な integer
!!            tips  : allocatable な integer
!!     で、自動的にdeallocateされるので注意
!!
!! get_links_to_s_sc(sc,s,link_labels,origins)
!!     sc のs番目のサイトに入るリンクのラベルと出発点を取得
!!      link_labels   : allocatable な integer
!!            origins : allocatable な integer
!!     で、自動的にdeallocateされるので注意
!!
!! get_links_in_face_sc(sc,f,FaceSize,sites,link_labels,link_dirs)
!!     sc のf番目のfaceを構成する
!!      1) facesize: face fのサイズ（リンクの数）
!!      2) sites: サイトのラベル
!!      3) link_labels: link のラベル
!!      4) link_dirs: そのリンクの方向(+/-1)
!!     を取得する
!!     "sites","link_labels","link_dirs"はallocatableな整数配列でなければならず、
!!     仮にallocateされていたとしても、サイズも含めて上書きされる。
!!
!! get_face_components_sc(sc,f,sites)
!!     face f を構成するサイトのリストを取得する
!!     "sites"はallocatableな整数配列でなければならず、
!!     仮にallocateされていたとしても、サイズも含めて上書きされる。
!!
!! get_faces_in_site_sc(sc,s,faces_s)
!!     sc のサイト s を含むfaceのラベルを取得
!!      1) s:サイトのラベル
!!      2) faces_s: sを含むfaceのラベル
!!     "faces_s"はallocatableな整数配列でなければならず、
!!     仮にallocateされていたとしても、サイズも含めて上書きされる。
!! 
!! get_faces_in_link_sc(sc,l,faces_l)
!!     sc のリンク l を含むfaceのラベルを取得
!!      1) l:リンクのラベル
!!      2) faces_l: lを含むfaceのラベル
!!     "faces_l"はallocatableな整数配列でなければならず、
!!     仮にallocateされていたとしても、サイズも含めて上書きされる。
!!
!! ############
!! 構造型 Link
!! integer :: origin_
!! integer :: tip_
!! ############
!!  put_Link(link,s,t) 
!!    linkの始点と終点を s,t を代入する
!!  put_LinkOrigin(link,s) 
!!    linkの始点に s を代入する
!!  put_LinkTip(link,t) 
!!    linkの終点に t を代入する
!!  get_Link(link,s,t) 
!!    s,t にlinkの始点と終点を代入する
!!  get_LinkOrigin(link,s) 
!!    sにlinkの始点を代入する
!!  get_LinkTip(link,t) 
!!    tにlinkの終点を代入する
!!  
!!##############
!! 構造型 LinkSet
!!  リンク変数の集合 （「originがsであるようなリンクの集合」のように使う）
!!  integer :: NumLinks_
!!  type(Link), allocatable :: LinkSet_(:)
!!###############
!!  init_Linkset(listset,N)  
!!   linksetの大きさをNにセットし、その大きさのLinkの配列を宣言する
!!  put_LinkToLS(linkset,n,s,t)
!!   linksetのn番目のLinkの始点と終点をs,tに指定する
!!  put_OriginToLS(linkset,n,s)
!!   linksetのn番目のLinkの始点をsに指定する
!!  put_TipToLS(linkset,n,t)
!!   linksetのn番目のLinkの終点をtに指定する
!!  get_LinkFromLS(ls,n,s,t)
!!   linksetのn番目のLinkの始点と終点を(s,t)に代入する
!! 
!!###############
!! 構造型 Face
!!  integer :: FaceSize_ !faceを構成する=siteの数
!!  integer, allocatable :: SitesF_(:) faceを構成するsiteの集合
!!  integer, allocatable :: DirecL_(:) !リンクの方向（1:正, -1:逆)
!!  最初のリンクの原点が代表点。つまり、
!!     DirecL_S(1)=1 => SitesF_(1) 
!!     DirecL_S(1)=-1 => SitesF_(2) 
!!   が代表点
!!
!! init_face(FC, sites, dirs)
!!  type(face), intent(inout) :: FC
!!  integer, intent(in) :: sites(1:FaceSize)
!!  integer, intent(in) :: dirs(1:FaceSize)
!!  face:FCを構成する SitesF_, DirecL_ を大きさ numsites にallocateし、
!!   SitesF_ = sites
!!   DirecL_ = dirs
!!  としてデータも代入する。
!!  ただし、numsitesはsitesのサイズで、自動的に取得される
!!
!!###############
!!構造型 Face2
!!  integer :: RefPt_  !代表点
!!  integer :: NumLinks_ !faceを構成するlinkの数=siteの数
!!  type(Linkset) :: LinkS_ ! linkの集合
!!  integer, allocatable :: DirecS_(:) !リンクの方向（0:正, 1:逆)
!!###############
!!  init_face2(FC,Refpt, numl)
!!    face2:FC の 代表点(Refpt_)とLinksetのサイズ(NumLinks_)を設定し、
!!    linkset( LinkS_ )とリンクの方向( DirecS_(1:NumLinks_) )をメモリ上に定義する
!!


module simplicial_complex
    implicit none

!!!!!!!!! 構造型 !!!!!!!!!!!!!
!!!!! Link変数 !!!!
type Link
    integer :: origin_
    integer :: tip_
end type Link

!!!!! リンク変数の集合 （「originがsであるようなリンクの集合」のように使う）
type Linkset
    integer :: NumLinks_
    type(Link), allocatable :: LinkSet_(:)
end type Linkset

!!!!! Face変数 !!!!!
type Face
    integer :: FaceSize_ !faceを構成するsiteの数
    integer, allocatable :: SitesF_(:) ! faceを構成するsiteの集合
    integer, allocatable :: DirecL_(:) !リンクの方向（0:正, 1:逆)
end type Face

!!!!! Face変数(old) !!!!!
type Face2
    integer :: RefPt_ ! 代表点
    integer :: NumLinks_ !faceを構成するlinkの数=siteの数
    type(Linkset) :: LinkS_ ! faceを構成するlinkの集合
    integer, allocatable :: DirecS_(:) !リンクの方向（1:正, -1:逆)
end type Face2

!!!!! SimplicalComplex !!!!!
type SmpCom
    integer :: NumSites_  ! 複体を構成するサイトの数
    integer :: NumLinks_  ! 複体を構成するリンクの数
    integer :: NumFaces_  ! 複体を構成するフェイスの数
    type(link), allocatable :: LinkSet_(:)
    type(Face), allocatable :: FaceSet_(:)
end type SmpCom

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 制御用の関数/副プログラム
contains
!!!
subroutine put_Link(ln,s,t)
    implicit none
    integer, intent(in) :: s, t
    type(Link), intent(inout) :: ln

    ln%origin_=s
    ln%tip_=t
end subroutine put_Link

!!!
subroutine put_LinkOrigin(ln,s)
    implicit none
    integer, intent(in) :: s
    type(Link), intent(inout) :: ln

    ln%origin_=s
end subroutine put_LinkOrigin

!!!
subroutine put_LinkTip(ln,t)
    implicit none
    integer, intent(in) :: t
    type(Link), intent(inout) :: ln

    ln%tip_=t
end subroutine put_LinkTip
    
!!!
subroutine get_Link(ln,s,t)
    implicit none
    integer, intent(inout) :: s, t
    type(Link), intent(in) :: ln

    s=ln%origin_
    t=ln%tip_
end subroutine get_Link

!!!
subroutine get_LinkOrigin(ln,s)
    implicit none
    integer, intent(inout) :: s
    type(Link), intent(in) :: ln

    s=ln%origin_
end subroutine get_LinkOrigin

!!!
subroutine get_LinkTip(ln,t)
    implicit none
    integer, intent(inout) :: t
    type(Link), intent(in) :: ln

    t=ln%tip_
end subroutine get_LinkTip

!!!
subroutine init_Linkset(list,N)
    implicit none
    type(Linkset), intent(inout) :: list
    integer, intent(in) :: N

    list%NumLinks_=N
    allocate(list%LinkSet_(1:N))
end subroutine init_Linkset

!!! 
subroutine put_LinkToLS(list,n,s,t)
    implicit none
    type(Linkset), intent(inout) :: list
    integer, intent(in) :: n,s,t

    if ( n<1 .or. list%NumLinks_<n ) then 
        stop "put_LinkToLS: out of range"
    endif
    call put_Link(list%LinkSet_(n),s,t)
end subroutine put_LinkToLS 

!!! 
subroutine put_OriginToLS(list,n,s)
    implicit none
    type(Linkset), intent(inout) :: list
    integer, intent(in) :: n,s

    if ( n<1 .or. list%NumLinks_<n ) then 
        stop "put_OritinToLS: out of range"
    endif
    call put_LinkOrigin(list%LinkSet_(n),s)
end subroutine put_OriginToLS

!!! 
subroutine put_TipToLS(list,n,t)
    implicit none
    type(Linkset), intent(inout) :: list
    integer, intent(in) :: n,t

    if ( n<1 .or. list%NumLinks_<n ) then 
        stop "put_TipToLS: out of range"
    endif
    call put_LinkTip(list%LinkSet_(n),t)
end subroutine put_TipToLS

!!!
subroutine get_LinkFromLS(ls,n,s,t)
    implicit none
    type(Linkset) ls
    integer, intent(in) :: n
    integer, intent(inout) :: s,t

    call get_link(ls%LinkSet_(n),s,t)
end subroutine get_LinkFromLS    

!!!
subroutine init_face(FC, sites, dirs)
    implicit none
    type(face), intent(inout) :: FC
    integer, intent(in) :: sites(:)
    integer, intent(in) :: dirs(:)
    integer :: facesize

    facesize=size(sites)

    FC%FaceSize_=facesize
    allocate( FC%SitesF_(1:facesize) )
    allocate( FC%DirecL_(1:facesize) )
    FC%SitesF_=sites
    FC%DirecL_=dirs
end subroutine init_face


!!!
subroutine  init_face2(FC, numl)
    implicit none
    type(face2), intent(inout) :: FC
    integer, intent(in) :: numl

    call init_Linkset(FC%LinkS_, numl)
    allocate( FC%DirecS_(1:numl) )
    !write(*,*) "##",FC%LinkS_%Numlink_
end subroutine init_face2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! init_smpcom(SC,Ns,Nl,Nf)
!!    simplical complex SC を初期化する。
!!    Ns: NumSites_
!!    Nl: NumLinks_
!!    Nf: NumFaces_
!!    allocate: LinkSet_(1:NumLinks_),  FaceSet_(1:NumFaces_)
subroutine init_smpcom(sc,ns,nl,nf)
    implicit none
    type(SmpCom), intent(inout) :: sc
    integer, intent(in) :: ns,nl,nf

    sc%NumSites_=ns
    sc%NumLinks_=nl
    sc%NumFaces_=nf

    allocate(sc%LinkSet_(1:nl))
    allocate(sc%FaceSet_(1:nf))
end subroutine init_smpcom

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! put_link_sc(sc,l,origin,tip)
!!      SmpCom:scのl番目のlinkを<origin,tip>にセットする。
!!        1 <= l <= NumLinks_
subroutine put_link_sc(sc,l,origin,tip)
    implicit none
    type(SmpCom), intent(inout) :: sc
    integer, intent(in) :: origin,tip,l

    call put_link(sc%LinkSet_(l),origin,tip)
end subroutine put_link_sc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! put_face_sc(sc,f,sites)
!!      SmpCom:scのf番目のfaceを設定する。
!!      ただし、事前にすべてのlinkをset_link_scで設定しておく必要がある。
!!        1 <= f <= NumFaces_
!!        sites: 大きさが FaceSize の整数配列
subroutine put_face_sc(sc,f,sites)
    implicit none
    type(SmpCom), intent(inout) :: sc
    integer, intent(in) :: f 
    integer, intent(in) :: sites(:)
    integer, allocatable :: dirs(:)
    integer :: i,j,s,t,facesize

    facesize=size(sites)
    allocate( dirs(1:facesize))
! faceを構成するリンクの方向を決定
    dirs=0
    do i=1,facesize
        s=sites(i)
        if(i==facesize) then 
            t=sites(1)
        else
            t=sites(i+1)
        endif

        do j=1, sc%NumLinks_
! LS(s)の中に<st>があるか？
            if( sc%LinkSet_(j)%origin_ == s .and. sc%LinkSet_(j)%tip_ == t ) then
                dirs(i)=1
                exit
! LS(t)の中に<ts>があるか？
            elseif( sc%LinkSet_(j)%origin_ == t .and. sc%LinkSet_(j)%tip_ == s ) then
                dirs(i)=-1
                exit
            endif
        enddo
        if ( dirs(i) == 0 ) then
            write(*,*) "face",f,": no link <",s,t,"> or <",t,s,">"
            stop
        endif
    enddo

    call init_face(sc%FaceSet_(f),sites,dirs)

end subroutine put_face_sc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! get_numsite_sc(sc,numsite)
!!     scを構成する全サイト数を取得
subroutine get_numsite_sc(sc,numsite)
    implicit none
    type(SmpCom), intent(in) :: sc
    integer, intent(inout) :: numsite

    numsite=sc%NumSites_
end subroutine get_numsite_sc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! get_numlink_sc(sc,numlink)
!!     scを構成する全リンク数を取得
subroutine get_numlink_sc(sc,numlink)
    implicit none
    type(SmpCom), intent(in) :: sc
    integer, intent(inout) :: numlink

    numlink=sc%NumLinks_
end subroutine get_numlink_sc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! get_numface_sc(sc,numface)
!!     scを構成する全リンク数を取得
subroutine get_numface_sc(sc,numface)
    implicit none
    type(SmpCom), intent(in) :: sc
    integer, intent(inout) :: numface

    numface=sc%NumFaces_
end subroutine get_numface_sc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! get_link_sc(sc,l,origin,tip)
!!     scのl番目のlinkのoriginとtipを取得
!!
subroutine get_link_sc(sc,l,origin,tip)
    implicit none
    type(SmpCom), intent(in) :: sc
    integer, intent(in) :: l
    integer, intent(inout) :: origin,tip
    integer :: numlink

    call get_numlink_sc(sc,numlink)

    ! エラー処理
    if (l<0 .or. l>numlink) then 
        stop
    endif

    call get_link(sc%LinkSet_(l),origin,tip)

end subroutine get_link_sc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! get_linklabel_sc(sc,link_label,link_dir,origin,tip)
!!     scの中のlink:<origin,tip>のラベルと方向を取得
!!     もし該当するリンクがなければ 0 を返す
subroutine get_linklabel_sc(sc,origin,tip,link_label,link_dir)
    implicit none
    type(SmpCom), intent(in) :: sc
    integer, intent(in) :: origin,tip
    integer, intent(inout) :: link_label,link_dir
    integer :: l,s,t

    link_label=0
    link_dir=0
    do l=1,sc%NumLinks_
        call get_link_sc(sc,l,s,t)
        if( s==origin .and. t==tip) then
          link_label=l
          link_dir=1
          exit
        elseif( s==tip .and. t==origin ) then
          link_label=l
          link_dir=-1
          exit
        endif
    enddo
end subroutine get_linklabel_sc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! get_links_from_s_sc(sc,s,links,tips)
!!     sc のs番目のサイトから出ているリンクのラベルと行き先を取得
!!      link_labels : allocatable な integer
!!      tips  : allocatable な integer
!!     で、自動的にdeallocateされるので注意
subroutine get_links_from_s_sc(sc,s,link_labels,tips)
    implicit none
    type(SmpCom), intent(in) :: sc
    integer, intent(in) :: s
    integer, allocatable, intent(inout) :: tips(:)
    integer, allocatable, intent(inout) :: link_labels(:)
    integer :: l,label,dir,numsites,t
    integer :: num
    integer, allocatable :: stock_tips(:)
    integer, allocatable :: stock_l(:)

    allocate( stock_tips(1:sc%NumSites_) )
    allocate( stock_l(1:sc%NumSites_) )


    !! tipsとlink_labelsを確保（もし既にallocateされていたらそれを外す）
    if( allocated(tips) ) deallocate(tips)
    if( allocated(link_labels) ) deallocate(link_labels)

    num=0
    do t=1,sc%NumSites_
      call get_linklabel_sc(sc,s,t,label,dir)
      if (dir == 1) then
        num=num+1
        stock_tips(num)=t
        stock_l(num)=label
      endif
    enddo

    allocate( tips(1:num) )
    allocate( link_labels(1:num) )
    do l=1,num
        tips(l)=stock_tips(l)
        link_labels(l)=stock_l(l)
    enddo
    deallocate( stock_tips )
    deallocate( stock_l )
end subroutine get_links_from_s_sc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! get_links_to_s_sc(sc,s,links,origins)
!!     sc のs番目のサイトに入るリンクのラベルと出発点を取得
!!      link_labels   : allocatable な integer
!!      origins : allocatable な integer
!!     で、自動的にdeallocateされるので注意
subroutine get_links_to_s_sc(sc,t,link_labels,origins)
    implicit none
    type(SmpCom), intent(in) :: sc
    integer, intent(in) :: t
    integer, allocatable, intent(inout) :: origins(:)
    integer, allocatable, intent(inout) :: link_labels(:)
    integer :: l,s,label,dir
    integer :: num
    integer, allocatable :: stock_origins(:)
    integer, allocatable :: stock_l(:)

    allocate( stock_origins(1:sc%NumSites_) )
    allocate( stock_l(1:sc%NumSites_) )

    !! tipsとlink_labelsを確保（もし既にallocateされていたらそれを外す）
    if( allocated(origins) ) deallocate(origins)
    if( allocated(link_labels) ) deallocate(link_labels)

    num=0
    do s=1,sc%NumSites_
      call get_linklabel_sc(sc,s,t,label,dir)
      if (dir == 1) then
        num=num+1
        stock_origins(num)=s
        stock_l(num)=label
      endif
    enddo

    allocate( origins(1:num) )
    allocate( link_labels(1:num) )
    do l=1,num
        origins(l)=stock_origins(l)
        link_labels(l)=stock_l(l)
    enddo
    deallocate( stock_origins )
    deallocate( stock_l )
end subroutine get_links_to_s_sc




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! get_links_in_face_sc(sc,f,FaceSize,sites,link_labels,link_dirs)
!!     sc のf番目のfaceを構成する
!!      1) facesize: face fのサイズ（リンクの数）
!!      2) sites: サイトのラベル
!!      3) link_labels: link のラベル
!!      4) link_dirs: そのリンクの方向(+/-1)
!!     を取得する
!!     "sites","link_labels","link_dirs"はallocatableな整数配列でなければならず、
!!     仮にallocateされていたとしても、サイズも含めて上書きされる。
subroutine get_links_in_face_sc(sc,f,FaceSize,sites,link_labels,link_dirs)
    implicit none
    type(SmpCom), intent(in) :: sc
    integer, intent(in) :: f
    integer, intent(inout) :: FaceSize
    integer, allocatable, intent(inout) :: sites(:), link_labels(:), link_dirs(:)
    integer i,origin,tip

    !! サイトの数を取得
    FaceSize=sc%FaceSet_(f)%FaceSize_

    !! sitesとdirsを確保（もし既にallocateされていたらそれを外す）
    if( allocated(sites) ) deallocate(sites)
    if( allocated(link_labels) ) deallocate(link_labels)
    if( allocated(link_dirs) ) deallocate(link_dirs)
    allocate( sites(1:FaceSize), link_labels(1:FaceSize), link_dirs(1:FaceSize) )

    sites=sc%FaceSet_(f)%SitesF_
    link_dirs=sc%FaceSet_(f)%DirecL_
    do i=1,FaceSize
        origin=sites(i)
        if(i==FaceSize) then
            tip=sites(1)
        else
            tip=sites(i+1)
        endif
        call get_linklabel_sc(sc,origin,tip,link_labels(i),link_dirs(i))
    enddo
end subroutine get_links_in_face_sc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! get_face_components_sc(sc,f,sites)
!!     face f を構成するサイトのリストを取得する
!!     "sites"はallocatableな整数配列でなければならず、
!!     仮にallocateされていたとしても、サイズも含めて上書きされる。
subroutine get_face_components_sc(sc,f,sites)
    implicit none
    type(SmpCom), intent(in) :: sc
    integer, intent(in) :: f
    integer, allocatable, intent(inout) :: sites(:)
    integer fsize

    !! sitesとdirsを確保（もし既にallocateされていたらそれを外す）
    if( allocated(sites) ) deallocate(sites)

    fsize=sc%FaceSet_(f)%FaceSize_
    allocate( sites(fsize) )
    sites=sc%FaceSet_(f)%SitesF_


end subroutine get_face_components_sc



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! get_faces_in_site_sc(sc,s,faces_s)
!!     sc のサイト s を含むfaceのラベルを取得
!!      1) s:サイトのラベル
!!      2) faces_s: sを含むfaceのラベル
!!     "faces_s"はallocatableな整数配列でなければならず、
!!     仮にallocateされていたとしても、サイズも含めて上書きされる。
subroutine get_faces_in_site_sc(sc,s,faces_s)
    implicit none
    type(SmpCom), intent(in) :: sc
    integer, intent(in) :: s
    integer, allocatable :: faces_s(:)
    integer, allocatable :: stock_faces(:)
    integer f,i,num
    
    if( allocated(faces_s) ) deallocate(faces_s)

    allocate( stock_faces(1:sc%NumFaces_) )
    stock_faces=0
    num=0
    do f=1,sc%NumFaces_
        do i=1,sc%FaceSet_(f)%FaceSize_
            if( sc%FaceSet_(f)%SitesF_(i) == s ) then
                num=num+1
                stock_faces(num)=f
                exit
            endif
        enddo
    enddo

    allocate( faces_s(1:num) )
    do i=1,num
        faces_s(i)=stock_faces(i)
    enddo
    deallocate( stock_faces )
end subroutine get_faces_in_site_sc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! get_faces_in_link_sc(sc,l,faces_l)
!!     sc のリンク l を含むfaceのラベルを取得
!!      1) l:サイトのラベル
!!      2) faces_l: lを含むfaceのラベル
!!     "faces_l"はallocatableな整数配列でなければならず、
!!     仮にallocateされていたとしても、サイズも含めて上書きされる。
subroutine get_faces_in_link_sc(sc,l,faces_l)
    implicit none
    type(SmpCom), intent(in) :: sc
        integer, intent(in) :: l
    integer, allocatable :: faces_l(:)
    integer, allocatable :: stock_faces(:)
    integer tip,origin,s,t
    integer f,i,num
    
    if( allocated(faces_l) ) deallocate(faces_l)

    call get_link_sc(sc,l,origin,tip)

    allocate( stock_faces(1:sc%NumFaces_) )
    stock_faces=0
    num=0
    do f=1,sc%NumFaces_
        do i=1,sc%FaceSet_(f)%FaceSize_
            s=sc%FaceSet_(f)%SitesF_(i)
            if(i==sc%FaceSet_(f)%FaceSize_) then
                t=sc%FaceSet_(f)%SitesF_(1)
            else
                t=sc%FaceSet_(f)%SitesF_(i+1)
            endif

            if( (origin==s .and. tip==t) .or. (origin==t .and. tip==s) ) then
                num=num+1
                stock_faces(num)=f
                exit
            endif
        enddo
    enddo

    allocate( faces_l(1:num) )
    do i=1,num
        faces_l(i)=stock_faces(i)
    enddo
    deallocate( stock_faces )

end subroutine get_faces_in_link_sc

end module simplicial_complex

