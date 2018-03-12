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
end module check_routines
