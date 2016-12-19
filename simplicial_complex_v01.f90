!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 単体的複体を表すmodule
!! 
!!  ver.01: linkの方向を表す数字を
!!    1: 順方向 -1: 逆方向 0: リンクなし
!!   に変更した
!! 
!!  マニュアル
!!##############
!! 構造型 SmpCom （単体的複体）
!!  integer :: NumSite_  ! 複体を構成するサイトの数
!!  type(Linkset), allocatable :: LS_(:)
!!     ! 各サイトから出るlinkの集合(L_from_S(1:NumSite_))
!!  integer :: NumFace_ ! 複体を構成するfaceの数
!!  type(face), allocatable :: FS_(:)
!!##############
!!  init_smpcom(SC,Ns,Nf)
!!    simplical complex SC を初期化する。
!!    NumSite_: Ns, NumFace_: Nf
!!    allocate: LS_(1:Ns), FS_(1:Nf)
!!    
!! set_links_sc(sc,s,tips)
!!      SmpCom:scのs番目のlinkset( site s から出るリンクの集合 )を設定
!!        tips(1:numl) : sから出るリンクの行き先
!!      ただし、numlは自動的に読み込まれるので指定しなくて良い。 
!!
!! set_face_sc(sc,f,sites)
!!      SmpCom:scのf番目のfaceを設定
!!      （set_links_scで事前にlinkを設定しておく必要あり）
!!                       f: faceの識別番号
!!        sites(1:numsite): そのfaceを構成するサイトの番号
!!      ただし、numsite: f番目のfaceを構成するsiteの数 は
!!      sitesから自動的に読み込まれるので指定しなくて良い。
!!
!! get_numsite_sc(sc,numsite)
!!     scを構成する全サイト数を取得
!!
!! get_numlink_sc(sc,numlink)
!!     scを構成する全リンク数を取得
!!
!! get_link_sc(sc,l,origin,tip)
!!     scのl番目のlinkのoriginとtipを取得
!!
!! get_links_sc(sc,s,tips)
!!     sc のs番目のサイトから出ているリンクの行き先リストを取得
!!     "tips"はallocatableな整数配列で、予めdeallocateしておく必要がある。
!!
!! get_face_sc(sc,f,sizeF, sites,dirs)
!!     sc のf番目のfaceを構成する
!!       サイト数: sizeF
!!       faceを構成するサイト番号: sites(1:sizeF)
!!       各辺のリンクの向き: dires(1:sizeF)
!!     を取得する
!!     "f_sites"と"dir"はallocatableな整数配列でなければならず、
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
!!  integer :: NumLink_
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
!!  integer :: NumSites_ !faceを構成する=siteの数
!!  integer, allocatable :: SitesF_(:) faceを構成するsiteの集合
!!  integer, allocatable :: DirecL_(:) !リンクの方向（1:正, -1:逆)
!!  最初のリンクの原点が代表点。つまり、
!!     DirecL_S(1)=1 => SitesF_(1) 
!!     DirecL_S(1)=-1 => SitesF_(2) 
!!   が代表点
!!
!! init_face(FC, sites, dirs)
!!  type(face), intent(inout) :: FC
!!  integer, intent(in) :: sites(1:numsites)
!!  integer, intent(in) :: dirs(1:numsites)
!!  face:FCを構成する SitesF_, DirecL_ を大きさ numsites にallocateし、
!!   SitesF_ = sites
!!   DirecL_ = dirs
!!  としてデータも代入する。
!!  ただし、numsitesはsitesのサイズで、自動的に取得される
!!
!!###############
!!構造型 Face2
!!  integer :: RefPt_  !代表点
!!  integer :: NumLink_ !faceを構成するlinkの数=siteの数
!!  type(Linkset) :: LinkS_ ! linkの集合
!!  integer, allocatable :: DirecS_(:) !リンクの方向（0:正, 1:逆)
!!###############
!!  init_face2(FC,Refpt, numl)
!!    face2:FC の 代表点(Refpt_)とLinksetのサイズ(NumLink_)を設定し、
!!    linkset( LinkS_ )とリンクの方向( DirecS_(1:NumLink_) )をメモリ上に定義する
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
    integer :: NumLink_
    type(Link), allocatable :: LinkSet_(:)
end type Linkset

!!!!! Face変数 !!!!!
type Face
    integer :: NumSites_ !faceを構成する=siteの数
    integer, allocatable :: SitesF_(:) ! faceを構成するsiteの集合
    integer, allocatable :: DirecL_(:) !リンクの方向（0:正, 1:逆)
end type Face

!!!!! Face変数(old) !!!!!
type Face2
    integer :: RefPt_ ! 代表点
    integer :: NumLink_ !faceを構成するlinkの数=siteの数
    type(Linkset) :: LinkS_ ! faceを構成するlinkの集合
    integer, allocatable :: DirecS_(:) !リンクの方向（1:正, -1:逆)
end type Face2

!!!!! SimplicalComplex !!!!!
type SmpCom
    integer :: NumSite_  ! 複体を構成するサイトの数
    type(Linkset), allocatable :: LS_(:)
       ! 各サイトから出るlinkの集合(サイズは常にNumSite_)
    integer :: NumFace_ ! 複体を構成するfaceの数
    type(Face), allocatable :: FS_(:)
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

    list%NumLink_=N
    allocate(list%LinkSet_(1:N))
end subroutine init_Linkset

!!! 
subroutine put_LinkToLS(list,n,s,t)
    implicit none
    type(Linkset), intent(inout) :: list
    integer, intent(in) :: n,s,t

    if ( n<1 .or. list%NumLink_<n ) then 
        stop "put_LinkToLS: out of range"
    endif
    call put_Link(list%LinkSet_(n),s,t)
end subroutine put_LinkToLS 

!!! 
subroutine put_OriginToLS(list,n,s)
    implicit none
    type(Linkset), intent(inout) :: list
    integer, intent(in) :: n,s

    if ( n<1 .or. list%NumLink_<n ) then 
        stop "put_OritinToLS: out of range"
    endif
    call put_LinkOrigin(list%LinkSet_(n),s)
end subroutine put_OriginToLS

!!! 
subroutine put_TipToLS(list,n,t)
    implicit none
    type(Linkset), intent(inout) :: list
    integer, intent(in) :: n,t

    if ( n<1 .or. list%NumLink_<n ) then 
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
!subroutine init_face(FC, numsites, sites, dirs)
    implicit none
    type(face), intent(inout) :: FC
    !integer, intent(in) :: numsites
    !integer, intent(in) :: sites(1:numsites)
    !integer, intent(in) :: dirs(1:numsites)
    integer, intent(in) :: sites(:)
    integer, intent(in) :: dirs(:)
    integer :: numsites

    numsites=size(sites)

    FC%NumSites_=numsites
    allocate( FC%SitesF_(1:numsites) )
    allocate( FC%DirecL_(1:numsites) )
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


!!! 
subroutine init_smpcom(sc,ns,nf)
    implicit none
    type(SmpCom), intent(inout) :: sc
    integer, intent(in) :: ns,nf

    sc%NumSite_=ns
    sc%NumFace_=nf

    allocate(sc%LS_(1:ns))
    allocate(sc%FS_(1:nf))
end subroutine init_smpcom

!!!
subroutine set_links_sc(sc,s,tips)
!subroutine set_links_sc(sc,s,numl,tips)
    implicit none
    type(SmpCom), intent(inout) :: sc
    integer, intent(in) :: s, tips(:)
    !integer, intent(in) :: s, numl, tips(1:numl)
    integer :: i,numl

    if ( s<1 .or. s>sc%NumSite_ ) then
        stop "out of range: set_link_sc"
    endif

    numl=size(tips)
    call init_Linkset(sc%LS_(s),numl)

    do i=1,numl
        call put_LinkToLS(sc%LS_(s),i,s,tips(i))
    enddo
end subroutine set_links_sc

subroutine set_face_sc(sc,f,sites)
!subroutine set_face_sc(sc,f,numsite,sites)
    implicit none
    type(SmpCom), intent(inout) :: sc
    integer, intent(in) :: f 
    !integer, intent(in) :: numsite 
    integer, intent(in) :: sites(:)
    !integer, intent(in) :: sites(1:numsite)
    integer, allocatable :: dirs(:)
    !integer :: dirs(1:numsite)
    integer :: i,j,s,t,numsite

    numsite=size(sites)
    allocate( dirs(1:numsite))
! faceを構成するリンクの方向を決定
    !dirs=-1
    dirs=0
    do i=1,numsite
        s=sites(i)
        if(i==numsite) then 
            t=sites(1)
        else
            t=sites(i+1)
        endif

! LS(s)の中に<st>があるか？
        do j=1, sc%LS_(s)%NumLink_
            if( sc%LS_(s)%LinkSet_(j)%tip_ == t ) then
                !dirs(i)=0
                dirs(i)=1
                exit
            endif
        enddo
! LS(t)の中に<ts>があるか？
        !if (dirs(i) == -1) then 
        if (dirs(i) == 0) then 
            do j=1, sc%LS_(t)%NumLink_
             if( sc%LS_(t)%LinkSet_(j)%tip_ == s ) then
                 dirs(i)=-1
                 exit
             endif
            enddo
        endif
        !if ( dirs(i) == -1 ) then
        if ( dirs(i) == 0 ) then
            write(*,*) "no link <",s,t,"> or <",t,s,">"
            stop
        endif
    enddo

    call init_face(sc%FS_(f),sites,dirs)

end subroutine set_face_sc

subroutine get_numsite_sc(sc,numsite)
    implicit none
    type(SmpCom), intent(in) :: sc
    integer, intent(inout) :: numsite

    numsite=sc%NumSite_
end subroutine get_numsite_sc

subroutine get_numlink_sc(sc,numlink)
    implicit none
    type(SmpCom), intent(in) :: sc
    integer, intent(inout) :: numlink
    integer :: i,s,numsite

    call get_numsite_sc(sc,numsite) 

    numlink=0
    do s=1,numsite
        numlink=numlink+sc%LS_(s)%NumLink_
    enddo
end subroutine get_numlink_sc

subroutine get_link_sc(sc,l,origin,tip)
    implicit none
    type(SmpCom), intent(in) :: sc
    integer, intent(in) :: l
    integer, intent(inout) :: origin,tip
    integer :: i,s,numsite,numlink,tmp


    call get_numlink_sc(sc,numlink)
    call get_numsite_sc(sc,numsite) 

    ! エラー処理
    if (l<0 .or. l>numlink) then 
        stop
    endif

    tmp=0
    s=0
    do while ( tmp<l )
        s=s+1
        tmp=tmp+sc%LS_(s)%NumLink_
    enddo
    tmp=tmp-sc%LS_(s)%NumLink_
    do i=1,sc%LS_(s)%NumLink_
        tmp=tmp+1
        if( tmp==l ) then
            exit
        endif
    enddo

    origin=s
    call get_LinkTip(sc%LS_(s)%LinkSet_(i),tip)

end subroutine get_link_sc


! sc のs番目のサイトから出ているリンクの行き先リストを取得
! ※ "tips"はallocatableな整数配列である必要がある
subroutine get_links_sc(sc,s,tips)
    implicit none
    type(SmpCom), intent(in) :: sc
    integer, intent(in) :: s
    integer, allocatable, intent(inout) :: tips(:)
    integer :: i

    allocate( tips(1:sc%LS_(s)%NumLink_) )

    do i=1,sc%LS_(s)%NumLink_
        call get_LinkTip(sc%LS_(s)%LinkSet_(i),tips(i))
    enddo
end subroutine get_links_sc

subroutine get_face_sc(sc, f, sizeF, sites, dirs)
    implicit none
    type(SmpCom), intent(in) :: sc
    integer, intent(in) :: f
    integer, intent(inout) :: sizeF
    integer, allocatable, intent(inout) :: sites(:), dirs(:)
    integer i

    !! サイトの数を取得
    sizeF=sc%FS_(f)%NumSites_

    !! sitesとdirsを確保（もし既にallocateされていたらそれを外す）
    if( allocated(sites) ) deallocate(sites)
    if( allocated(dirs) ) deallocate(dirs)
    allocate( sites(1:sizeF), dirs(1:sizeF) )

    sites=sc%FS_(f)%SitesF_
    dirs=sc%FS_(f)%DirecL_
end subroutine get_face_sc

end module simplicial_complex

