!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Subroutine to compute Dinv.PF directly by using Dinv_FILE
subroutine DinvPF_direct(DinvPF_eta, DinvPF_lambda, DinvPF_chi, &
    eta,lambda,chi,&
    Geta_eta, Glambda_eta, Gchi_eta, &
    Geta_lambda, Glambda_lambda, Gchi_lambda, &
    Geta_chi, Glambda_chi, Gchi_chi, &
    Umat,PhiMat)
use global_parameters
use parallel
implicit none

complex(kind(0d0)), intent(in) :: eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: chi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

complex(kind(0d0)), intent(in) ::  Geta_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_sites) 
complex(kind(0d0)) , intent(in) :: Glambda_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_sites) 
complex(kind(0d0)) , intent(in) :: Gchi_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_sites) 
!!!
complex(kind(0d0)) , intent(in) :: Geta_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_links) 
complex(kind(0d0)) , intent(in) :: Glambda_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_links) 
complex(kind(0d0)) , intent(in) :: Gchi_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_links) 
!!!
complex(kind(0d0)) , intent(in) :: Geta_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_faces)  
complex(kind(0d0)) , intent(in) :: Glambda_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_faces) 
complex(kind(0d0)) , intent(in) :: Gchi_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_faces) 

complex(kind(0d0)), intent(out) :: DinvPF_eta(1:NMAT,1:NMAT,1:num_sites) 
complex(kind(0d0)), intent(out) :: DinvPF_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(out) :: DinvPF_chi(1:NMAT,1:NMAT,1:num_faces)

complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: recvmat(1:NMAT,1:NMAT)
integer :: i,j,ls,ll,lf,gs,gl,gf
integer :: k,l
integer :: rank

DinvPF_eta=(0d0,0d0)
DinvPF_lambda=(0d0,0d0)
DinvPF_chi=(0d0,0d0)

do gs=1,global_num_sites
  tmpmat=(0d0,0d0)
  recvmat=(0d0,0d0)
  do ls=1,num_sites
    do j=1,NMAT
      do i=1,NMAT
        tmpmat=tmpmat+Geta_eta(:,:,i,j,gs,ls)*eta(j,i,ls)
      enddo
    enddo
  enddo
  !!
  do ll=1,num_links
    do j=1,NMAT
      do i=1,NMAT
        tmpmat=tmpmat+Geta_lambda(:,:,i,j,gs,ll)*lambda(j,i,ll)
      enddo
    enddo
  enddo
  !!
  do lf=1,num_faces
    do j=1,NMAT
      do i=1,NMAT
        tmpmat=tmpmat+Geta_chi(:,:,i,j,gs,lf)*chi(j,i,lf)
      enddo
    enddo
  enddo
  rank=local_site_of_global(gs)%rank_
  ls=local_site_of_global(gs)%label_
  call MPI_REDUCE(tmpmat,recvmat,NMAT*NMAT,MPI_DOUBLE_COMPLEX, &
    MPI_SUM,rank,MPI_COMM_WORLD,IERR)
  if( MYRANK==rank) then
    DinvPF_eta(:,:,ls)=recvmat
  endif
enddo
!!!!!!!
do gl=1,global_num_links
  tmpmat=(0d0,0d0)
  recvmat=(0d0,0d0)
  do ls=1,num_sites
    do j=1,NMAT
      do i=1,NMAT
        tmpmat=tmpmat+Glambda_eta(:,:,i,j,gl,ls)*eta(j,i,ls)
      enddo
    enddo
  enddo
  !!
  do ll=1,num_links
    do j=1,NMAT
      do i=1,NMAT
        tmpmat=tmpmat+Glambda_lambda(:,:,i,j,gl,ll)*lambda(j,i,ll)
      enddo
    enddo
  enddo
  !!
  do lf=1,num_faces
    do j=1,NMAT
      do i=1,NMAT
        tmpmat=tmpmat+Glambda_chi(:,:,i,j,gl,lf)*chi(j,i,lf)
      enddo
    enddo
  enddo
  rank=local_link_of_global(gl)%rank_
  ll=local_link_of_global(gl)%label_
  call MPI_REDUCE(tmpmat,recvmat,NMAT*NMAT,MPI_DOUBLE_COMPLEX, &
    MPI_SUM,rank,MPI_COMM_WORLD,IERR)
  if( MYRANK==rank ) then
    DinvPF_lambda(:,:,ll)=recvmat
  endif
enddo
!!!!!!!
do gf=1,global_num_faces
  tmpmat=(0d0,0d0)
  recvmat=(0d0,0d0)
  do ls=1,num_sites
    do j=1,NMAT
      do i=1,NMAT
        tmpmat=tmpmat+Gchi_eta(:,:,i,j,gf,ls)*eta(j,i,ls)
      enddo
    enddo
  enddo
  !!
  do ll=1,num_links
    do j=1,NMAT
      do i=1,NMAT
        tmpmat=tmpmat+Gchi_lambda(:,:,i,j,gf,ll)*lambda(j,i,ll)
      enddo
    enddo
  enddo
  !!
  do lf=1,num_faces
    do j=1,NMAT
      do i=1,NMAT
        tmpmat=tmpmat+Gchi_chi(:,:,i,j,gf,lf)*chi(j,i,lf)
      enddo
    enddo
  enddo
  rank=local_face_of_global(gf)%rank_
  lf=local_face_of_global(gf)%label_
  call MPI_REDUCE(tmpmat,recvmat,NMAT*NMAT,MPI_DOUBLE_COMPLEX, &
    MPI_SUM,rank,MPI_COMM_WORLD,IERR)
  if( MYRANK==rank ) then
    DinvPF_chi(:,:,lf)=recvmat
  endif

enddo
!!!!!!!

end subroutine DinvPF_direct

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_Dinv(ite, &
Geta_eta, Glambda_eta, Gchi_eta, &
Geta_lambda, Glambda_lambda, Gchi_lambda, &
Geta_chi, Glambda_chi, Gchi_chi, &
N_DinvFILE)
use global_parameters
use parallel
use SUN_generators, only : make_traceless_matrix_from_modes
implicit none

integer, intent(out) :: ite
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
integer, intent(in) :: N_DinvFILE

!!!
complex(kind(0d0)) :: Gany_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: Gany_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: Gany_chi(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)) :: modes_eta(1:dimG,1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: modes_lambda(1:dimG,1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)) :: modes_chi(1:dimG,1:NMAT,1:NMAT,1:num_faces)

integer :: gs,gl,gf
integer :: ls,ll,lf
integer :: a,b
integer :: i,j,k,l

complex(kind(0d0)) :: matrix(1:NMAT,1:NMAT)
integer :: tag, rank
integer :: err

if( MYRANK==0 ) then
  read(N_DinvFILE, '(I10)', advance='no', iostat=err) ite
endif
call MPI_BCAST(ite,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
call MPI_BCAST(err,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
if( err/=0 ) then 
  call stop_for_test
endif

Geta_eta=(0d0,0d0)
Geta_lambda=(0d0,0d0)
Geta_chi=(0d0,0d0)
Glambda_eta=(0d0,0d0)
Glambda_lambda=(0d0,0d0)
Glambda_chi=(0d0,0d0)
Gchi_eta=(0d0,0d0)
Gchi_lambda=(0d0,0d0)
Gchi_chi=(0d0,0d0)

do gs=1,global_num_sites
  do a=1,dimG
    call read_data(Gany_eta, Gany_lambda, Gany_chi, N_DinvFILE)
    do j=1,NMAT
      do i=1,NMAT
        do ls=1,num_sites
          modes_eta(a,i,j,ls)=Gany_eta(i,j,ls)
        enddo
        do ll=1,num_links
          modes_lambda(a,i,j,ll)=Gany_lambda(i,j,ll)
        enddo
        do lf=1,num_faces
          modes_chi(a,i,j,lf)=Gany_chi(i,j,lf)
        enddo
      enddo
    enddo
  enddo
  do j=1,NMAT
    do i=1,NMAT
      do ls=1,num_sites
        call Make_traceless_matrix_from_modes(&
          Geta_eta(:,:,i,j,gs,ls),NMAT,modes_eta(:,i,j,ls))
      enddo
      do ll=1,num_links
        call Make_traceless_matrix_from_modes(&
          Geta_lambda(:,:,i,j,gs,ll),NMAT,modes_lambda(:,i,j,ll))
      enddo
      do lf=1,num_faces
        call Make_traceless_matrix_from_modes(&
          Geta_chi(:,:,i,j,gs,lf),NMAT,modes_chi(:,i,j,lf))
      enddo
    enddo
  enddo
enddo
!!!!!!!!!!!!
do gl=1,global_num_links
  do a=1,dimG
    call read_data(Gany_eta, Gany_lambda, Gany_chi, N_DinvFILE)
    do j=1,NMAT
      do i=1,NMAT
        do ls=1,num_sites
          modes_eta(a,i,j,ls)=Gany_eta(i,j,ls)
        enddo
        do ll=1,num_links
          modes_lambda(a,i,j,ll)=Gany_lambda(i,j,ll)
        enddo
        do lf=1,num_faces
          modes_chi(a,i,j,lf)=Gany_chi(i,j,lf)
        enddo
      enddo
    enddo
  enddo
  do j=1,NMAT
    do i=1,NMAT
      do ls=1,num_sites
        call Make_traceless_matrix_from_modes(&
          Glambda_eta(:,:,i,j,gl,ls),NMAT,modes_eta(:,i,j,ls))
      enddo
      do ll=1,num_links
        call Make_traceless_matrix_from_modes(&
          Glambda_lambda(:,:,i,j,gl,ll),NMAT,modes_lambda(:,i,j,ll))
      enddo
      do lf=1,num_faces
        call Make_traceless_matrix_from_modes(&
          Glambda_chi(:,:,i,j,gl,lf),NMAT,modes_chi(:,i,j,lf))
      enddo
    enddo
  enddo
enddo
!!!!!!!!!!!!
do gf=1,global_num_faces
  do a=1,dimG
    call read_data(Gany_eta, Gany_lambda, Gany_chi, N_DinvFILE)
    do j=1,NMAT
      do i=1,NMAT
        do ls=1,num_sites
          modes_eta(a,i,j,ls)=Gany_eta(i,j,ls)
        enddo
        do ll=1,num_links
          modes_lambda(a,i,j,ll)=Gany_lambda(i,j,ll)
        enddo
        do lf=1,num_faces
          modes_chi(a,i,j,lf)=Gany_chi(i,j,lf)
        enddo
      enddo
    enddo
  enddo
  do j=1,NMAT
    do i=1,NMAT
      do ls=1,num_sites
        call Make_traceless_matrix_from_modes(&
          Gchi_eta(:,:,i,j,gf,ls),NMAT,modes_eta(:,i,j,ls))
      enddo
      do ll=1,num_links
        call Make_traceless_matrix_from_modes(&
          Gchi_lambda(:,:,i,j,gf,ll),NMAT,modes_lambda(:,i,j,ll))
      enddo
      do lf=1,num_faces
        call Make_traceless_matrix_from_modes(&
          Gchi_chi(:,:,i,j,gf,lf),NMAT,modes_chi(:,i,j,lf))
      enddo
    enddo
  enddo
enddo

if( MYRANK==0 ) then
  read(N_DinvFILE, *) 
endif

end subroutine read_Dinv

!!!!!!!!!!!!!!!!!!!!!!
subroutine read_data(&
    Gany_eta, Gany_lambda, Gany_chi,&
    N_DinvFILE)!,rank,tag)
!subroutine read_matrix(matrix,N_DinvFILE,rank,tag)
use global_parameters
use parallel
use SUN_generators, only : make_traceless_matrix_from_modes
implicit none
complex(kind(0d0)), intent(out) :: Gany_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(out) :: Gany_lambda(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(out) :: Gany_chi(1:NMAT,1:NMAT,1:num_faces)
!integer, intent(inout) :: tag
integer, intent(in) :: N_DinvFILE!, rank

double precision :: real, imag
character(50) :: C_real, C_imag
complex(kind(0d0)) :: ele
complex(kind(0d0)) :: modes(1:dimG)
complex(kind(0d0)) :: matrix(1:NMAT,1:NMAT)

integer :: rank, ls,ll,lf, tag
integer :: gs,gl,gf
integer :: a,b

tag=0
do gs=1,global_num_sites
  rank=local_site_of_global(gs)%rank_
  ls=local_site_of_global(gs)%label_
  tag=tag+1
  if( MYRANK==0 ) then
    do a=1,dimG
      read(N_DinvFILE,'(a17,a17)',advance='no') C_real, C_imag
      read(C_real,*) real
      read(C_imag,*) imag
      modes(a)=dcmplx(real)+(0d0,1d0)*dcmplx(imag)
    enddo
    if( rank /= 0 ) then 
      call MPI_SEND(modes,NMAT*NMAT-1,MPI_DOUBLE_COMPLEX,&
        rank,tag,MPI_COMM_WORLD,IERR)
    endif
  endif
  if( MYRANK==rank .and. rank/=0 ) then
    call MPI_RECV(modes,NMAT*NMAT-1,MPI_DOUBLE_COMPLEX,&
      0,tag,MPI_COMM_WORLD,ISTATUS,IERR)
  endif
  if( MYRANK == rank ) then
    call Make_traceless_matrix_from_modes(matrix,NMAT,MODES)
    Gany_eta(:,:,ls)=matrix
  endif
enddo
!!!
do gl=1,global_num_links
  rank=local_link_of_global(gl)%rank_
  ll=local_link_of_global(gl)%label_
  tag=tag+1
  if( MYRANK== 0 ) then
    do a=1,dimG
      read(N_DinvFILE,'(a17,a17)',advance='no') C_real, C_imag
      read(C_real,*) real
      read(C_imag,*) imag
      modes(a)=dcmplx(real)+(0d0,1d0)*dcmplx(imag)
    enddo
    if( rank /= 0 ) then 
      call MPI_SEND(modes,NMAT*NMAT-1,MPI_DOUBLE_COMPLEX,&
        rank,tag,MPI_COMM_WORLD,IERR)
    endif
  endif
  if( MYRANK==rank .and. rank/=0 ) then
    call MPI_RECV(modes,NMAT*NMAT-1,MPI_DOUBLE_COMPLEX,&
      0,tag,MPI_COMM_WORLD,ISTATUS,IERR)
  endif
  if( MYRANK == rank ) then
    call Make_traceless_matrix_from_modes(matrix,NMAT,MODES)
    Gany_lambda(:,:,ll)=matrix
  endif
enddo
!!!
do gf=1,global_num_faces
  rank=local_face_of_global(gf)%rank_
  lf=local_face_of_global(gf)%label_
  tag=tag+1
  if( MYRANK==0 ) then 
    do a=1,dimG
      read(N_DinvFILE,'(a17,a17)',advance='no') C_real, C_imag
      read(C_real,*) real
      read(C_imag,*) imag
      modes(a)=dcmplx(real)+(0d0,1d0)*dcmplx(imag)
    enddo
    if( rank /= 0 ) then 
      call MPI_SEND(modes,NMAT*NMAT-1,MPI_DOUBLE_COMPLEX,&
        rank,tag,MPI_COMM_WORLD,IERR)
    endif
  endif
  if( MYRANK==rank .and. rank/=0 ) then
    call MPI_RECV(modes,NMAT*NMAT-1,MPI_DOUBLE_COMPLEX,&
      0,tag,MPI_COMM_WORLD,ISTATUS,IERR)
  endif
  if( MYRANK == rank ) then
    call Make_traceless_matrix_from_modes(matrix,NMAT,MODES)
    Gany_chi(:,:,lf)=matrix
  endif
enddo

end subroutine read_data



