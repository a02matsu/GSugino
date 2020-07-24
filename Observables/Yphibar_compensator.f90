!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! A = 1/NF sum_f ( 1/N Tr(\phibar^{dimG*\chi/2} )
subroutine calc_Yphibar_compensator(Acomp,PhiMat,UMAT)
use parallel
use global_parameters
use matrix_functions, only : matrix_product
implicit none

complex(kind(0d0)), intent(out) :: Acomp
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: UMat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) :: phibar(1:NMAT,1:NMAT)
complex(kind(0d0)) :: phibar_p(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Ymat(1:NMAT,1:NMAT)
integer :: ratio,eular
double precision :: radius, phase

complex(kind(0d0)) :: tmp_Acomp, tmp
integer :: ls,lf
integer :: i,j

eular=global_num_sites-global_num_links+global_num_faces 
ratio=(NMAT*NMAT-1)*eular/2
Acomp=(0d0,0d0)
tmp_Acomp=(0d0,0d0)
do lf=1,num_faces
  !! phi^{-dimG Eular/2}
  ls=sites_in_f(lf)%label_(1)
  call hermitian_conjugate(phibar,Phimat(:,:,ls))
  call matrix_power(phibar_p,phibar,ratio)
  !! Omega
  call Make_face_variable(Uf,lf,UMAT)
  call Make_moment_map_adm(Ymat,Uf)
  Ymat = Ymat * (0d0,0.5d0)*beta_f(lf)*Ymat
  !! tr( phi^{-dimG r} Y )
  call trace_MM(tmp,phibar_p, Ymat)
  !! 
  tmp_Acomp=tmp_Acomp + tmp/dcmplx(NMAT)
enddo

call MPI_REDUCE(tmp_Acomp,Acomp,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)
  
Acomp=Acomp/dcmplx(dble(global_num_faces))

end subroutine calc_Yphibar_compensator


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Q A_phibarY \Xi term 
subroutine calc_QCYphibar_Xi(&
    QC_Xi_site,QC_Xi_link,QC_Xi_face, &
    Geta_eta,Geta_lambda,Geta_chi,&
    Gchi_eta,Gchi_lambda,Gchi_chi,&
    Umat,PhiMat)
implicit none

complex(kind(0d0)), intent(out) :: QC_Xi_site
complex(kind(0d0)), intent(out) :: QC_Xi_link
complex(kind(0d0)), intent(out) :: QC_Xi_face
complex(kind(0d0)), intent(in) ::  Geta_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_sites) 
complex(kind(0d0)), intent(in) :: Geta_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_links) 
complex(kind(0d0)), intent(in) :: Geta_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_faces) 
complex(kind(0d0)), intent(in) :: Gchi_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_sites) 
complex(kind(0d0)), intent(in) :: Gchi_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_links) 
complex(kind(0d0)), intent(in) :: Gchi_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_faces) 
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)


complex(kind(0d0)) Xi_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)) Xi_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)) Xi_chi(1:NMAT,1:NMAT,1:num_necessary_faces)

complex(kind(0d0)), allocatable :: phibar_p(:,:,:)
complex(kind(0d0)) comm(1:NMAT,1:NMAT)
complex(kind(0d0)) Ymat(1:NMAT,1:NMAT)
complex(kind(0d0)) tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) DinvXi_site(1:NMAT,1:NMAT)
complex(kind(0d0)) DinvXi_link(1:NMAT,1:NMAT)
complex(kind(0d0)) DinvXi_face(1:NMAT,1:NMAT)
complex(kind(0d0)) trace
complex(kind(0d0)) tmp_QC_Xi_site
complex(kind(0d0)) tmp_QC_Xi_link
complex(kind(0d0)) tmp_QC_Xi_face
integer gf,ls,ll,lf,rank
integer i,j,k,l,kk

complex(kind(0d0)) tr_phibar2
integer :: ratio,eular
double precision :: radius, phase

eular=global_num_sites-global_num_links+global_num_faces 
ratio=(NMAT*NMAT-1)*eular/2

!! preparation 
allocate( phibar_p(1:NMAT,1:NMAT,0:ratio) )
call make_unit_matrix(phibar_p(:,:,0))


QC_Xi_site=(0d0,0d0)
QC_Xi_link=(0d0,0d0)
QC_Xi_face=(0d0,0d0)
tmp_QC_Xi_site=(0d0,0d0)
tmp_QC_Xi_link=(0d0,0d0)
tmp_QC_Xi_face=(0d0,0d0)
call make_XiVec(Xi_eta,Xi_lambda,Xi_chi,Umat,PhiMat)
!!!
!! gf=gs
do gf=1,global_num_faces
  gs=global_sites_in_f(gf)%label_(1)

  rank=local_site_of_global(gs)%rank_
  if( MYRANK==rank) then
    lf=local_face_of_global(gf)
    ls=local_site_of_global(gs)

    call hermitian_conjugate(phibar_p(:,:,1), PhiMat(:,:,ls))
    do k=2,ratio
      call matrix_product(phibar_p(:,:,k),phibar_p(:,:,k-1),phibar_p(:,:,1))
    enddo
    !! Y
    call Make_face_variable(tmpmat,lf,UMAT)
    call Make_moment_map_adm(Ymat,tmpmat)
    Ymat = Ymat * (0d0,0.5d0)*beta_f(lf)*Ymat
  endif
  call MPI_BCAST(phibar_p, NMAT*NMAT*(ratio+1), MPI_DOUBLE_COMPLEX,rank,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(Ymat, NMAT*NMAT, MPI_DOUBLE_COMPLEX,rank,MPI_COMM_WORLD,IERR)

  do kk=0, ratio-1
    call matrix_3_product(tmpmat, phibar_p(:,:,ratio-kk-1), Ymat, phibar_p(:,:,kk))
    !! [\phibar^r, \phi]
    do i=1,NMAT
      do j=1,NMAT
        do k=1,NMAT
          do l=1,NMAT
            !! site
            do ls=1,num_sites
              tmp_QC_Xi_site=tmp_QC_Xi_site&
                + tmpmat(i,j)*Xi_eta(k,l,ls)*Geta_eta(j,i,l,k,gs,ls)
            enddo
            !! link
            do ll=1,num_links
              tmp_QC_Xi_link=tmp_QC_Xi_link&
                + tmpmat(i,j)*Xi_lambda(k,l,ll)*Geta_lambda(j,i,l,k,gs,ll)
            enddo
            !! face
            do lf=1,num_faces
              tmp_QC_Xi_face=tmp_QC_Xi_face&
                + tmpmat(i,j)*Xi_chi(k,l,lf)*Geta_eta(j,i,l,k,gs,lf)
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
  !! 
  call matrix_commutator(comm,phibar_p(:,:,ratio),PhiMat(:,:,ls))
  do i=1,NMAT
    do j=1,NMAT
      do k=1,NMAT
        do l=1,NMAT
          !! site
          do ls=1,num_sites
            tmp_QC_Xi_site=tmp_QC_Xi_site&
              + comm(i,j)*Xi_eta(k,l,ls)*Gchi_eta(j,i,l,k,gf,ls)
          enddo
          !! link
          do ll=1,num_links
            tmp_QC_Xi_link=tmp_QC_Xi_link&
              + comm(i,j)*Xi_lambda(k,l,ll)*Gchi_lambda(j,i,l,k,gf,ll)
          enddo
          !! face
          do lf=1,num_faces
            tmp_QC_Xi_face=tmp_QC_Xi_face&
              + comm(i,j)*Xi_chi(k,l,lf)*Gchi_eta(j,i,l,k,gf,lf)
          enddo
        enddo
      enddo
    enddo
  enddo
  !!
enddo

call MPI_REDUCE(tmp_QC_Xi_site,QC_Xi_site,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)
call MPI_REDUCE(tmp_QC_Xi_link,QC_Xi_link,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)
call MPI_REDUCE(tmp_QC_Xi_face,QC_Xi_face,1,MPI_DOUBLE_COMPLEX, &
  MPI_SUM,0,MPI_COMM_WORLD,IERR)

if( MYRANK==0 ) then
  QC_Xi_site=QC_Xi_site/dcmplx( NMAT*global_num_faces ) 
  QC_Xi_link=QC_Xi_link/dcmplx( NMAT*global_num_faces ) 
  QC_Xi_face=QC_Xi_face/dcmplx( NMAT*global_num_faces ) 
endif

end subroutine calc_QCYphibar_Xi


