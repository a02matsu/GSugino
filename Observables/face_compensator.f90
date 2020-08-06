!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
subroutine calc_face_compensator(Acomp,Umat,PhiMat,Geta_chi)
use parallel
use global_parameters
use matrix_functions, only : hermitian_conjugate, matrix_power, trace_mm, make_unit_matrix, matrix_product, matrix_3_product, matrix_commutator
implicit none

complex(kind(0d0)), intent(out) :: Acomp
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: Geta_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_faces) 

complex(kind(0d0)) :: tmp_Acomp, tmp
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)), allocatable :: phibar_p(:,:,:)
complex(kind(0d0)) :: Omega(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT)
integer :: lf,ls,gs,gf
integer :: i,j,k,l,p
integer :: ratio

ratio = (NMAT*NMAT-1)*(global_num_sites-global_num_links+global_num_faces)/2
allocate( phibar_p(1:NMAT,1:NMAT,0:ratio) )

Acomp=(0d0,0d0)
tmp_Acomp=(0d0,0d0)
do lf=1, num_faces
  ls=sites_in_f(lf)%label_(1)
  gf=global_face_of_local(lf)
  gs=global_sites_in_f(gf)%label_(1)

  !! phibar_p = \bar(\PhiMat)^p
  call make_unit_matrix(phibar_p(:,:,0))
  call hermitian_conjugate(phibar_p(:,:,1), PhiMat(:,:,ls))
  do k=2,ratio
    call matrix_product(phibar_p(:,:,k),phibar_p(:,:,k-1),phibar_p(:,:,1))
  enddo
  !! Omega
  call Make_face_variable(Uf,lf,UMAT)
  call Make_moment_map_adm(Omega,Uf)

  tmp=(0d0,0d0)
  !! bosonic part
  do j=1,NMAT
    do i=1,NMAT
      tmp = tmp + (0d0,0.5d0)*beta_f(lf)*phibar_p(i,j,ratio)*Omega(j,i)
    enddo
  enddo
  !! fermionic part
  do p=0,ratio-1
    do l=1,NMAT
      do k=1,NMAT
        do j=1,NMAT
          do i=1,NMAT
            tmp = tmp + phibar_p(i,j,ratio-1-p)*phibar_p(k,l,p)&
              *Geta_chi(j,k,l,i,gs,lf)
          enddo
        enddo
      enddo
    enddo
  enddo

  tmp_Acomp=tmp_Acomp + alpha_f(lf)*tmp
enddo

call MPI_REDUCE(tmp_Acomp,Acomp,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)
  
Acomp=Acomp/dcmplx(dble(global_num_faces * NMAT))

end subroutine calc_face_compensator


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 4-fermio terms in C_face SF_site
subroutine calc_4fermi_in_CSFsite(CSF, Umat, Phimat, &
    Geta_eta, Gchi_eta )
use parallel
use global_parameters
use Dirac_operator , only : prod_Dirac_site
use matrix_functions, only : hermitian_conjugate, matrix_power, trace_mm, make_unit_matrix, matrix_product, matrix_3_product, matrix_commutator
implicit none

complex(kind(0d0)), intent(out) :: CSF
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: Geta_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_sites) 
complex(kind(0d0)), intent(in) :: Gchi_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_sites) 

complex(kind(0d0)), allocatable :: SMAT(:,:,:,:,:,:,:) 
complex(kind(0d0)), allocatable :: FMAT(:,:,:,:,:,:,:) 
complex(kind(0d0)), allocatable :: phibar_p(:,:,:)
complex(kind(0d0)), allocatable :: phibar_p2(:,:,:)
complex(kind(0d0)) :: DSmat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: DFmat(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)) :: Xi_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) :: trace1, trace2, trace3, trace4
complex(kind(0d0)) :: tmp1, tmp2, tmp3, tmp4

integer :: ls,lf,gs,gf
integer :: tag, rank
integer :: i,j,k,l,p,a,b
integer :: ratio
complex(kind(0d0)) :: tmp_CSF

ratio = (NMAT*NMAT-1)*(global_num_sites-global_num_links+global_num_faces)/2
allocate( SMAT(1:NMAT,1:NMAT,1:num_necessary_sites,1:NMAT,1:NMAT,0:ratio-1,1:global_num_faces) )
allocate( FMAT(1:NMAT,1:NMAT,1:num_necessary_sites,1:NMAT,1:NMAT,0:ratio-1,1:global_num_faces) )
allocate( phibar_p(1:NMAT,1:NMAT,0:ratio) )
allocate( phibar_p2(1:NMAT,1:NMAT,0:ratio) )

!! 
SMAT=(0d0,0d0)
FMAT=(0d0,0d0)
do gf=1,global_num_faces
  rank=local_face_of_global(gf)%rank_
  lf=local_face_of_global(gf)%label_
  ls=sites_in_f(lf)%label_(1)
  gs=global_sites_in_f(gf)%label_(1)

  !! phibar_p = phibar^{0...r}(:,:,gf)
  phibar_p=(0d0,0d0)
  phibar_p2=(0d0,0d0)
  call make_phibar_p(phibar_p2,PhiMat,ratio,gf)
  if( MYRANK==rank ) then
    !! phibar_p = \bar(\PhiMat)^p
    call make_unit_matrix(phibar_p(:,:,0))
    call hermitian_conjugate(phibar_p(:,:,1), PhiMat(:,:,ls))
    do k=2,ratio
      call matrix_product(phibar_p(:,:,k),phibar_p(:,:,k-1),phibar_p(:,:,1))
    enddo
  endif
  call MPI_BCAST(phibar_p,NMAT*NMAT*(ratio+1),MPI_DOUBLE_COMPLEX,rank,MPI_COMM_WORLD,IERR)
  write(*,*) phibar_p-phibar_p2

!! prepare SMAT and FMAT
  do p=0,ratio-1
    do j=1,NMAT
      do i=1,NMAT
        do ls=1,num_sites
          call matrix_product(SMAT(:,:,ls,i,j,p,gf),&
            phibar_p(:,:,p), Geta_eta(:,:,i,j,gs,ls))
          call matrix_product(FMAT(:,:,ls,i,j,p,gf),&
            phibar_p(:,:,p), Gchi_eta(:,:,i,j,gf,ls))
        enddo !ls
        call syncronize_sites(SMAT(:,:,:,i,j,p,gf))
        call syncronize_sites(FMAT(:,:,:,i,j,p,gf))
      enddo
    enddo
  enddo
enddo

tmp_CSF=(0d0,0d0)
CSF=(0d0,0d0)
!! (1) Dirac term
do gf=1,global_num_faces
  do p=0,ratio-1
    do i=1,NMAT
      do j=1,NMAT
        call prod_Dirac_site(DFmat,PhiMat,FMAT(:,:,:,j,i,ratio-p-1,gf))
        call prod_Dirac_site(DSmat,PhiMat,SMAT(:,:,:,j,i,p,gf))
        do ls=1,num_sites
          do b=1,NMAT
            do a=1,NMAT
              tmp_CSF=tmp_CSF &
                + (-0.5d0,0d0) * SMAT(a,b,ls,i,j,p,gf) * DFmat(b,a,ls) &
                + (0.5d0,0d0) * FMAT(a,b,ls,i,j,ratio-p-1,gf) * DSmat(b,a,ls) 
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
enddo
call MPI_REDUCE(tmp_CSF,CSF,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)

!! (2) mass term
call make_XiVec_site(Xi_eta,PhiMat)
do gf=1,global_num_faces
  do p=0,ratio-1
    do j=1,NMAT
      do i=1,NMAT
        trace1=(0d0,0d0)
        trace2=(0d0,0d0)
        trace3=(0d0,0d0)
        trace4=(0d0,0d0)
        do ls=1,num_sites
          tmp1=(0d0,0d0)
          tmp2=(0d0,0d0)
          tmp3=(0d0,0d0)
          tmp4=(0d0,0d0)
          do b=1,NMAT
            do a=1,NMAT
              tmp1=tmp1 + Xi_eta(a,b,ls)*Smat(b,a,ls,i,j,p,gf)
              tmp2=tmp2 + Phimat(a,b,ls)*Fmat(b,a,ls,j,i,ratio-p-1,gf)
              !!
              tmp3=tmp3 + Xi_eta(a,b,ls)*Fmat(b,a,ls,i,j,p,gf)
              tmp4=tmp4 + Phimat(a,b,ls)*Smat(b,a,ls,j,i,ratio-p-1,gf)
            enddo
          enddo
          call MPI_REDUCE(tmp1,trace1,1,MPI_DOUBLE_COMPLEX, &
            MPI_SUM,0,MPI_COMM_WORLD,IERR)
          call MPI_REDUCE(tmp2,trace2,1,MPI_DOUBLE_COMPLEX, &
            MPI_SUM,0,MPI_COMM_WORLD,IERR)
          call MPI_REDUCE(tmp3,trace3,1,MPI_DOUBLE_COMPLEX, &
            MPI_SUM,0,MPI_COMM_WORLD,IERR)
          call MPI_REDUCE(tmp4,trace4,1,MPI_DOUBLE_COMPLEX, &
            MPI_SUM,0,MPI_COMM_WORLD,IERR)
        enddo
        if( MYRANK==0 ) then
          CSF=CSF&
            -dcmplx( 0.5d0*mass_square_phi )*trace1*trace2 &
            +dcmplx( 0.5d0*mass_square_phi )*trace3*trace4
        endif
      enddo
    enddo
  enddo
enddo

if( MYRANK==0 ) then
  CSF=CSF / dcmplx( NMAT * global_num_faces )
endif

end subroutine calc_4fermi_in_CSFsite

!!!!
subroutine make_phibar_p(phibar_p,PhiMat,ratio,gf)
use global_parameters
use parallel
use matrix_functions, only : hermitian_conjugate, make_unit_matrix, matrix_product
implicit none

integer, intent(in) :: ratio,gf
complex(kind(0d0)), intent(out) :: phibar_p(1:NMAT,1:NMAT,0:ratio)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)

integer :: gs,lf,ls,rank,k

rank=local_face_of_global(gf)%rank_
lf=local_face_of_global(gf)%label_
ls=sites_in_f(lf)%label_(1)
gs=global_sites_in_f(gf)%label_(1)

!! phibar_p = phibar^{0...r}(:,:,gf)
if( MYRANK==rank ) then
  !! phibar_p = \bar(\PhiMat)^p
  call make_unit_matrix(phibar_p(:,:,0))
  call hermitian_conjugate(phibar_p(:,:,1), PhiMat(:,:,ls))
  do k=2,ratio
    call matrix_product(phibar_p(:,:,k),phibar_p(:,:,k-1),phibar_p(:,:,1))
  enddo
endif
call MPI_BCAST(phibar_p,NMAT*NMAT*(ratio+1),MPI_DOUBLE_COMPLEX,rank,MPI_COMM_WORLD,IERR)

end subroutine make_phibar_p
