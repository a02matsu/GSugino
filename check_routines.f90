module check_routines
use global_parameters
#ifdef PARALLEL
use parallel
#endif
implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! routines for check
subroutine test_module_simplicial_complex(sc)
use simplicial_complex
implicit none
type(SmpCom) sc
integer numsites,numlinks,numfaces
integer origin,tip
integer link_label, link_dir
integer,allocatable :: link_labels(:),tips(:),origins(:),faces(:),sites(:)
integer s,l,f
integer i,j,k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call get_numsite_sc(sc,numsites)
write(*,*) "test:get_numsite_sc: numsite=",numsites
write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call get_numlink_sc(sc,numlinks)
write(*,*) "test:get_numlink_sc: numlink=",numlinks
write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call get_numface_sc(sc,numfaces)
write(*,*) "test:get_numface_sc: numface=",numfaces
write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*) "test:get_link_sc"
do l=1,numlinks
call get_link_sc(sc,l,origin,tip)
  write(*,*) l,"'th link=(",origin,tip,")"
enddo
write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*) "test:get_linklabel_sc"
do i=1,numsites
do j=1,numsites
  if(i .ne. j) then
  call get_linklabel_sc(sc,i,j,link_label,link_dir)
  write(*,*) "(",i,j,"): label=",link_label,",dir=",link_dir
endif
enddo
enddo
write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*) "test:get_links_from_s_sc"
do s=1,numsites
  call get_links_from_s_sc(sc,s,link_labels,tips)
  write(*,*) "links from",s,"going to"
  do i=1,size(link_labels)
    write(*,*) "  ",tips(i),"is label=",link_labels(i)
  enddo
enddo
write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*) "test:get_links_to_s_sc"
do s=1,numsites
  call get_links_to_s_sc(sc,s,link_labels,origins)
  write(*,*) "links to",s,"coming from"
  do i=1,size(link_labels)
    write(*,*) "  ",origins(i),"is label=",link_labels(i)
  enddo
enddo
write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*) "test:get_face_components_sc"
do f=1,numfaces
  call get_face_components_sc(sc,f,sites)
  write(*,*) "face",f,"is made of", sites
enddo
write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*) "test:get_faces_in_site_sc"
do s=1,numsites
  call get_faces_in_site_sc(sc,s,faces)
  write(*,*) "site",s,"is included in the faces",faces 
enddo
write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*) "test:get_faces_in_link_sc"
do l=1,numlinks
  call get_faces_in_link_sc(sc,l,faces)
  write(*,*) "link",l,"is included in the faces",faces
enddo
write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

end subroutine test_module_simplicial_complex


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! check unitaryty
subroutine  check_unitary(UMAT)
use matrix_functions, only : matrix_norm
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_links)
integer l,i,j,k
double precision norm
complex(kind(0d0)) tmp(1:NMAT,1:NMAT)

write(*,*) "===== check unitarity of UMAT ========="
do l=1,num_links
  tmp=(0d0,0d0)
  do i=1,NMAT
    do j=1,NMAT
        do k=1,NMAT
            tmp(i,j)=tmp(i,j)+UMAT(i,k,l)*dconjg(UMAT(j,k,l))
        enddo
    enddo
  enddo
  do i=1,NMAT
    tmp(i,i)=tmp(i,i)-(1d0,0d0)
  enddo
  call matrix_norm(norm,tmp)
    write(*,*) "link:",l,":||U_l.U_l^\dagger - 1||=",norm
enddo

end subroutine check_unitary



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! test  Make_traceless_matrix_from_modes
subroutine  check_Ta_expansion
use mt95
use global_subroutines
use SUN_generators, only : make_traceless_matrix_from_modes, trace_MTa
implicit none

double precision MAT_re(1:NMAT,1:NMAT), MAT_im(1:NMAT,1:NMAT)
double precision MODES_re(1:dimG),MODES_im(1:dimG)
complex(kind(0d0)) MAT(1:NMAT,1:NMAT),tmp
complex(kind(0d0)) MAT2(1:NMAT,1:NMAT)
complex(kind(0d0)) MODES(1:dimG),MODES2(1:dimG)
integer i,j,k,a

call genrand_real3(MAT_re)
call genrand_real3(MAT_im)

MAT=MAT_re + (0d0,1d0)*MAT_im
tmp=(0d0,0d0)
do i=1,NMAT-1
    tmp=tmp+MAT(i,i)
enddo
MAT(NMAT,NMAT)=-tmp

do a=1,dimG
    call Trace_MTa(MODES(a),MAT,a,NMAT)
    write(*,*) MODES(a)
enddo

call Make_traceless_matrix_from_modes(MAT2,NMAT,MODES)

write(*,*) "### TEST1 ###"
do i=1,NMAT
do j=1,NMAT
write(*,*) i,j, MAT(i,j)-MAT2(i,j)
enddo
enddo


call genrand_real3(MODES_re)
call genrand_real3(MODES_im)
MODES=dcmplx(MODES_re) + (0d0,1d0)*dcmplx(MODES_im)
call Make_traceless_matrix_from_modes(MAT,NMAT,MODES)
do a=1,dimG
    call Trace_MTa(MODES2(a),MAT,a,NMAT)
enddo

write(*,*) "### TEST2 ###"
do a=1,dimG
write(*,*) MODES(a) - MODES2(a)
enddo



end subroutine check_Ta_expansion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! global SC の情報を出力する
#ifdef PARALLEL
subroutine check_global_sc
use parallel

integer :: rank, turn
integer :: l,f

turn=0
if ( turn .ne. MYRANK ) then
  call MPI_RECV(turn,1,MPI_INTEGER,MYRANK-1,MYRANK,MPI_COMM_WORLD,ISTATUS,IERR)
endif
write(*,*) "#####",MYRANK,"#####"

write(*,*) "global_num_sites,links,faces:", global_num_sites,global_num_links,global_num_faces

write(*,*) "### links ###"
do l=1,global_num_links
  write(*,'(a,I3,a,I3,a,I3,a)') "link",l,"=(",global_link_org(l),",",global_link_tip(l),")"
enddo


write(*,*)
write(*,*)
turn=MYRANK+1
if( MYRANK .ne. NPROCS-1 ) then
  call MPI_SEND(turn,1,MPI_INTEGER,MYRANK+1,MYRANK+1,MPI_COMM_WORLD,IERR)
endif




end subroutine check_global_sc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! local なラベルを出力する
subroutine check_local_sc
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
write(*,*) "### link ###"
do local=1,num_necessary_links
  write(*,*) "local link",local,"->","global",global_link_of_local(local)
enddo
write(*,*) "### face ###"
do local=1,num_necessary_faces
  write(*,*) "local face",local,"->","global",global_face_of_local(local)
enddo


write(*,*) "#!!!! SEND and RECV !!!!"
write(*,*) "### send site ###"
do local=1,size(send_sites,1)
  write(*,*) "send local site",send_sites(local)%label_,"to RANK",send_sites(local)%rank_
enddo
write(*,*) "### send link ###"
do local=1,size(send_links,1)
  write(*,*) "send local link",send_links(local)%label_,"to RANK",send_links(local)%rank_
enddo
write(*,*) "### send face ###"
do local=1,size(send_faces,1)
  write(*,*) "send local face",send_faces(local)%label_,"to RANK",send_faces(local)%rank_
enddo
write(*,*) "### recv site ###"
do local=1,size(recv_sites,1)
  write(*,*) "recv local site",recv_sites(local)%label_,"from RANK",recv_sites(local)%rank_
enddo
write(*,*) "### recv link ###"
do local=1,size(recv_links,1)
  write(*,*) "recv local link",recv_links(local)%label_,"from RANK",recv_links(local)%rank_
enddo
write(*,*) "### recv face ###"
do local=1,size(recv_faces,1)
  write(*,*) "recv local face",recv_faces(local)%label_,"from RANK",recv_faces(local)%rank_
enddo

write(*,*) "#!!!! SITE-LINK  !!!!"
do local = 1,num_links
  write(*,*) "local link",local,"is",link_org(local),link_tip(local)
enddo

write(*,*) "###"
do local=1,num_sites
  write(*,*) "local links from",local,";",linktip_from_s(local)%labels_
enddo
write(*,*) "###"
do local=1,num_sites
  write(*,*) "local links to",local,";",linkorg_to_s(local)%labels_
enddo


write(*,*)
write(*,*)
turn=MYRANK+1
if( MYRANK .ne. NPROCS-1 ) then
  call MPI_SEND(turn,1,MPI_INTEGER,MYRANK+1,MYRANK+1,MPI_COMM_WORLD,IERR)
endif

end subroutine check_local_sc
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! local なラベルを出力する
#ifdef PARALLEL
subroutine check_local_vals(PhiMat,UMat)
use parallel
implicit none

complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

integer :: rank, turn,i,j
integer :: local,global
complex(kind(0d0)) :: tmp

turn=0
if ( turn .ne. MYRANK ) then
  call MPI_RECV(turn,1,MPI_INTEGER,MYRANK-1,MYRANK,MPI_COMM_WORLD,ISTATUS,IERR)
endif
write(*,*) "#####",MYRANK,"#####"

write(*,*) "#!!!! local data to global data !!!!"
write(*,*) "### Tr|Phi(s)|^2 ###"
do local=1,num_necessary_sites
  tmp=(0d0,0d0)
  do i=1,NMAT
    do j=1,NMAT
      tmp=(PhiMat(i,j,local)*dconjg(PhiMat(i,j,local)))
    enddo
  enddo
  write(*,*) "s=",global_site_of_local(local),":",dble(tmp)
enddo
write(*,*) "### sum_{ij} U(l)_{ij} ###"
do local=1,num_necessary_links
  tmp=(0d0,0d0)
  do i=1,NMAT
    do j=1,NMAT
      tmp=tmp+UMat(i,j,local)
    enddo
  enddo
  write(*,*) "l=",global_link_of_local(local),":",tmp
enddo

write(*,*)
write(*,*)
turn=MYRANK+1
if( MYRANK .ne. NPROCS-1 ) then
  call MPI_SEND(turn,1,MPI_INTEGER,MYRANK+1,MYRANK+1,MPI_COMM_WORLD,IERR)
endif


end subroutine check_local_vals
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! force の出力
#ifdef PARALLEL
subroutine check_force(force,C)
use parallel
implicit none

integer, intent(in) :: C
complex(kind(0d0)), intent(in) :: force(:,:,:)
integer :: rank, turn,i,j,num
integer :: local,global
complex(kind(0d0)) :: tmp


num=size(force,3)
turn=0
if ( turn .ne. MYRANK ) then
  call MPI_RECV(turn,1,MPI_INTEGER,MYRANK-1,MYRANK,MPI_COMM_WORLD,ISTATUS,IERR)
endif

do local=1,num
  tmp=(0d0,0d0)
  do i=1,NMAT
    do j=1,NMAT
      tmp=(force(i,j,local)*dconjg(force(i,j,local)))
    enddo
  enddo
  if(C==1) write(*,*) "s=",global_site_of_local(local),":",dble(tmp)
  if(C==2) write(*,*) "l=",global_link_of_local(local),":",dble(tmp)
  if(C==3) write(*,*) "f=",global_face_of_local(local),":",dble(tmp)
enddo
turn=MYRANK+1
if( MYRANK .ne. NPROCS-1 ) then
  call MPI_SEND(turn,1,MPI_INTEGER,MYRANK+1,MYRANK+1,MPI_COMM_WORLD,IERR)
endif
end subroutine check_force
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! alpha beta の出力
subroutine check_alpha_beta
use parallel
implicit none

integer turn
integer :: ls,ll,lf
integer :: gs,gl,gf

!!!!!!! alpha_s !!!!!!!!!
call MPI_BARRIER(MPI_COMM_WORLD,IERR)
if(MYRANK==0) write (*,*) "### alpha_s ###"
turn=0
if ( turn .ne. MYRANK ) then
  call MPI_RECV(turn,1,MPI_INTEGER,MYRANK-1,MYRANK,MPI_COMM_WORLD,ISTATUS,IERR)
endif
do ls=1,num_necessary_sites
  gs=global_site_of_local(ls)
  write(*,*) MYRANK, ls, gs, alpha_s(ls)
enddo
turn=MYRANK+1
if( MYRANK .ne. NPROCS-1 ) then
  call MPI_SEND(turn,1,MPI_INTEGER,MYRANK+1,MYRANK+1,MPI_COMM_WORLD,IERR)
endif


!!!!!!! alpha_l !!!!!!!!!
call MPI_BARRIER(MPI_COMM_WORLD,IERR)
if(MYRANK==0) write (*,*) "### alpha_l ###"
turn=0
if ( turn .ne. MYRANK ) then
  call MPI_RECV(turn,1,MPI_INTEGER,MYRANK-1,MYRANK,MPI_COMM_WORLD,ISTATUS,IERR)
endif
do ll=1,num_necessary_links
  gl=global_link_of_local(ll)
  write(*,*) MYRANK, ll, gl, alpha_l(ll)
enddo
turn=MYRANK+1
if( MYRANK .ne. NPROCS-1 ) then
  call MPI_SEND(turn,1,MPI_INTEGER,MYRANK+1,MYRANK+1,MPI_COMM_WORLD,IERR)
endif

!!!!!!! alpha_f, beta_f !!!!!!!!!
call MPI_BARRIER(MPI_COMM_WORLD,IERR)
if(MYRANK==0) write (*,*) "### alpha_f and beta_f ###"
turn=0
if ( turn .ne. MYRANK ) then
  call MPI_RECV(turn,1,MPI_INTEGER,MYRANK-1,MYRANK,MPI_COMM_WORLD,ISTATUS,IERR)
endif
do lf=1,num_necessary_faces
  gf=global_face_of_local(lf)
  write(*,*) MYRANK, lf, gf, alpha_f(lf), beta_f(lf)
enddo
turn=MYRANK+1
if( MYRANK .ne. NPROCS-1 ) then
  call MPI_SEND(turn,1,MPI_INTEGER,MYRANK+1,MYRANK+1,MPI_COMM_WORLD,IERR)
endif

end subroutine check_alpha_beta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine to test Q-transformation of S
!subroutine check_QS(Umat,PhiMat)
!use matrix_functions, only : matrix_commutator, matrix_3_product
!use Dirac_operator, only : Prod_Dirac
!use simulation, only : make_bosonic_force_nomass
!#ifdef PARALLEL
!use parallel
!#endif
!implicit none
!
!complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
!
!complex(kind(0d0)) :: Qeta(1:NMAT,1:NMAT,1:num_necessary_sites)
!complex(kind(0d0)) :: Qlambda(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)) :: Qchi(1:NMAT,1:NMAT,1:num_necessary_faces)
!
!complex(kind(0d0)) :: tmp_eta(1:NMAT,1:NMAT,1:num_sites)
!complex(kind(0d0)) :: tmp_lambda(1:NMAT,1:NMAT,1:num_links)
!complex(kind(0d0)) :: tmp_chi(1:NMAT,1:NMAT,1:num_faces)
!complex(kind(0d0)) :: Omega(1:NMAT,1:NMAT,1:num_faces)
!
!complex(kind(0d0)) :: Bforce_s(1:NMAT,1:NMAT,1:num_sites)
!complex(kind(0d0)) :: Bforce_l(1:NMAT,1:NMAT,1:num_links)
!!complex(kind(0d0)) :: Fforce_s(1:NMAT,1:NMAT,1:num_sites)
!!complex(kind(0d0)) :: Fforce_l(1:NMAT,1:NMAT,1:num_links)
!
!complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
!double precision :: tmp,QS
!integer :: info,s,l,f,i,j,triger
!
!call make_Qfermion(Qeta,Qlambda,Qchi,Omega,Umat,PhiMat)
!call make_bosonic_force_nomass(Bforce_s,Bforce_l,Umat,PhiMat)
!
!call Prod_Dirac(tmp_eta,tmp_lambda,tmp_chi,Qeta,Qlambda,Qchi,UMAT,Phimat)
!! Q^2 \Omega を care する
!!do f=1,num_faces
!!  call matrix_commutator(tmpmat,PhiMat(:,:,sites_in_f(f)%label_(1)),Omega(:,:,f))
!!  tmp_chi(:,:,f)=tmp_chi(:,:,f)+(0d0,1d0)*dcmplx(beta_f(f))*tmpmat
!!enddo
!
!if(MYRANK==0) write(*,*) "# QS = 0 ?"
!!write(*,*) tmp_chi
!QS=0d0
!tmp=0d0
!do s=1,num_sites
!  !tmp=0d0
!  do i=1,NMAT
!    do j=1,NMAT
!      tmp_eta(i,j,s)=-tmp_eta(i,j,s)+dconjg(Bforce_s(j,i,s))
!      tmp=tmp+dble( tmp_eta(i,j,s)*dconjg(tmp_eta(i,j,s)) )
!    enddo
!  enddo
!  !write(*,*) "# site",global_site_of_local(s),tmp
!enddo
!call MPI_REDUCE(tmp,QS,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
!if(MYRANK==0) write(*,*) "#   site:",QS
!
!QS=0d0
!tmp=0d0
!do l=1,num_links
!  !tmp=0d0
!  do i=1,NMAT
!    do j=1,NMAT
!      tmp_lambda(i,j,l) = -tmp_lambda(i,j,l)+Bforce_l(i,j,l)
!      tmp=tmp+dble( tmp_lambda(i,j,l)*dconjg(tmp_lambda(i,j,l)) )
!    enddo
!  enddo
!  !write(*,*) "# link",global_link_of_local(l),tmp
!enddo
!call MPI_REDUCE(tmp,QS,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
!if(MYRANK==0) write(*,*) "#   link:",QS
!
!QS=0d0
!tmp=0d0
!do f=1,num_faces
!  !tmp=0d0
!  do i=1,NMAT
!    do j=1,NMAT
!      tmp_chi(i,j,f)=-tmp_chi(j,i,f)
!      tmp=tmp+dble( tmp_chi(i,j,f)*dconjg(tmp_chi(i,j,f)) )
!    enddo
!  enddo
!  !write(*,*) "# face",global_face_of_local(f),tmp
!enddo
!call MPI_REDUCE(tmp,QS,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
!if(MYRANK==0) write(*,*) "#   face:",QS
!
!end subroutine check_QS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to make Q\Psi
!subroutine make_Qfermion(Qeta,Qlambda,Qchi,Omega,Umat,PhiMat)
!use global_subroutines
!use matrix_functions, only : matrix_commutator, matrix_3_product
!use parallel
!implicit none
!
!complex(kind(0d0)), intent(out) :: Qeta(1:NMAT,1:NMAT,1:num_necessary_sites)
!complex(kind(0d0)), intent(out) :: Qlambda(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)), intent(out) :: Qchi(1:NMAT,1:NMAT,1:num_necessary_faces)
!complex(kind(0d0)), intent(out) :: Omega(1:NMAT,1:NMAT,1:num_faces)
!complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
!!complex(kind(0d0)), intent(in) :: PF_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
!!complex(kind(0d0)) :: comm(1:NMAT,1:NMAT)
!
!integer :: s,l,f,i,j
!complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT)
!
!do s=1,num_sites
!  call matrix_commutator(Qeta(:,:,s),PhiMat(:,:,s),PhiMat(:,:,s),'N','C')
!enddo
!
!do l=1,num_links
!  Qlambda(:,:,l)=(0d0,-1d0)*PhiMat(:,:,link_org(l))
!  call matrix_3_product(Qlambda(:,:,l),&
!    Umat(:,:,l),PhiMat(:,:,link_tip(l)),Umat(:,:,l),&
!    'N','N','C',(0d0,1d0),'ADD')
!enddo
!
!do f=1,num_faces
!  call Make_face_variable(Uf,f,UMAT)
!  if(m_omega == 0) then 
!    call Make_moment_map0(Omega(:,:,f),Uf)
!  elseif(m_omega == -1) then
!    call Make_moment_map_adm(Omega(:,:,f),Uf)
!  endif
!  Qchi(:,:,f)=(0d0,-0.5d0)*dcmplx(beta_f(f))*Omega(:,:,f)
!enddo
!
!#ifdef PARALLEL
!call syncronize_sites(Qeta)
!call syncronize_links(Qlambda)
!call syncronize_faces(Qchi)
!#endif
!
!end subroutine make_Qfermion




end module check_routines
