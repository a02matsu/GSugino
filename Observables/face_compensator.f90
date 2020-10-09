!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
subroutine calc_face_compensator(&
    Acomp,CSF_site,CSF_link,CSF_face,&
    Umat,PhiMat,&
    Geta_eta, Geta_lambda, Geta_chi, Gchi_eta, Gchi_lambda, Gchi_chi) 
use parallel
use global_parameters
implicit none

complex(kind(0d0)), intent(out) :: Acomp, CSF_site, CSF_link, CSF_face
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: Geta_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_sites) 
complex(kind(0d0)), intent(in) :: Geta_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_links) 
complex(kind(0d0)), intent(in) :: Geta_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_faces) 
complex(kind(0d0)), intent(in) :: Gchi_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_sites) 
complex(kind(0d0)), intent(in) :: Gchi_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_links) 
complex(kind(0d0)), intent(in) :: Gchi_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_faces) 
integer :: ratio

ratio = (NMAT*NMAT-1)*(global_num_sites-global_num_links+global_num_faces)/2
call calc_face_compensator_main(Acomp,Umat,PhiMat,Geta_chi)
call calc_4fermi_in_CSFsite(CSF_site, Umat, Phimat, Geta_eta, Gchi_eta, ratio)
call calc_4fermi_in_CSFlink(CSF_link, Umat, Phimat, Geta_eta, Gchi_eta, Geta_lambda, Gchi_lambda, ratio )
call calc_4fermi_in_CSFface(CSF_face, Umat, Phimat, Geta_eta, Gchi_eta, Geta_chi, Gchi_chi, Geta_lambda, Gchi_lambda ,ratio)


contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine calc_face_compensator_main(Acomp,Umat,PhiMat,Geta_chi)
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
  
  end subroutine calc_face_compensator_main
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! 4-fermio terms in C_face SF_site
  subroutine calc_4fermi_in_CSFsite(CSF, Umat, Phimat, &
      Geta_eta, Gchi_eta, ratio )
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
  integer, intent(in) :: ratio

  complex(kind(0d0)) :: SMAT(1:NMAT,1:NMAT,1:num_necessary_sites,1:NMAT,1:NMAT,0:ratio-1,1:global_num_faces)
  complex(kind(0d0)) :: FMAT(1:NMAT,1:NMAT,1:num_necessary_sites,1:NMAT,1:NMAT,0:ratio-1,1:global_num_faces)
  complex(kind(0d0)) :: phibar_p(1:NMAT,1:NMAT,0:ratio)

  complex(kind(0d0)) :: tmp1(1:NMAT,1:NMAT,0:ratio,1:global_num_faces)
  complex(kind(0d0)) :: tmp2(1:NMAT,1:NMAT,0:ratio,1:global_num_faces)
  complex(kind(0d0)) :: tmp3(1:NMAT,1:NMAT,0:ratio,1:global_num_faces)
  complex(kind(0d0)) :: tmp4(1:NMAT,1:NMAT,0:ratio,1:global_num_faces)
  complex(kind(0d0)) :: ttmp1(1:NMAT,1:NMAT,0:ratio,1:global_num_faces)
  complex(kind(0d0)) :: ttmp2(1:NMAT,1:NMAT,0:ratio,1:global_num_faces)
  complex(kind(0d0)) :: ttmp3(1:NMAT,1:NMAT,0:ratio,1:global_num_faces)
  complex(kind(0d0)) :: ttmp4(1:NMAT,1:NMAT,0:ratio,1:global_num_faces)
  
  !complex(kind(0d0)), allocatable :: SMAT(:,:,:,:,:,:,:) 
  !complex(kind(0d0)), allocatable :: FMAT(:,:,:,:,:,:,:) 
  !complex(kind(0d0)), allocatable :: phibar_p(:,:,:)
  complex(kind(0d0)) :: DSmat(1:NMAT,1:NMAT,1:num_sites)
  complex(kind(0d0)) :: DFmat(1:NMAT,1:NMAT,1:num_sites)
  complex(kind(0d0)) :: Xi_eta(1:NMAT,1:NMAT,1:num_necessary_sites)
  complex(kind(0d0)) :: trace,tmp
  !complex(kind(0d0)), allocatable :: tmp1(:,:,:,:)
  !complex(kind(0d0)), allocatable :: tmp2(:,:,:,:)
  !complex(kind(0d0)), allocatable :: tmp3(:,:,:,:)
  !complex(kind(0d0)), allocatable :: tmp4(:,:,:,:)
  !complex(kind(0d0)), allocatable :: ttmp1(:,:,:,:)
  !complex(kind(0d0)), allocatable :: ttmp2(:,:,:,:)
  !complex(kind(0d0)), allocatable :: ttmp3(:,:,:,:)
  !complex(kind(0d0)), allocatable :: ttmp4(:,:,:,:)
  
  integer :: ls,lf,gs,gf
  integer :: tag, rank
  integer :: i,j,k,l,p,a,b
  complex(kind(0d0)) :: tmp_CSF
  
  !ratio = (NMAT*NMAT-1)*(global_num_sites-global_num_links+global_num_faces)/2
  !allocate( SMAT(1:NMAT,1:NMAT,1:num_necessary_sites,1:NMAT,1:NMAT,0:ratio-1,1:global_num_faces) )
  !allocate( FMAT(1:NMAT,1:NMAT,1:num_necessary_sites,1:NMAT,1:NMAT,0:ratio-1,1:global_num_faces) )
  !allocate( phibar_p(1:NMAT,1:NMAT,0:ratio) )
!
  !allocate( tmp1(1:NMAT,1:NMAT,0:ratio,1:global_num_faces) )
  !allocate( tmp2(1:NMAT,1:NMAT,0:ratio,1:global_num_faces) )
  !allocate( tmp3(1:NMAT,1:NMAT,0:ratio,1:global_num_faces) )
  !allocate( tmp4(1:NMAT,1:NMAT,0:ratio,1:global_num_faces) )
  !allocate( ttmp1(1:NMAT,1:NMAT,0:ratio,1:global_num_faces) )
  !allocate( ttmp2(1:NMAT,1:NMAT,0:ratio,1:global_num_faces) )
  !allocate( ttmp3(1:NMAT,1:NMAT,0:ratio,1:global_num_faces) )
  !allocate( ttmp4(1:NMAT,1:NMAT,0:ratio,1:global_num_faces) )
  
  !! 
  SMAT=(0d0,0d0)
  FMAT=(0d0,0d0)
  do gf=1,global_num_faces
    gs=global_sites_in_f(gf)%label_(1)
    !rank=local_face_of_global(gf)%rank_
    !lf=local_face_of_global(gf)%label_
    !if( rank==MYRANK ) then
      !ls=sites_in_f(lf)%label_(1)
    !else 
      !ls=0
    !endif
  
    !! phibar_p = phibar^{0...r}(:,:,gf)
    call make_phibar_p(phibar_p,PhiMat,ratio,gf)
  
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
        enddo
      enddo
    enddo
  enddo
  call syncronize_sites_large(SMAT,ratio)
  call syncronize_sites_large(FMAT,ratio)
  
  tmp_CSF=(0d0,0d0)
  CSF=(0d0,0d0)
  !! (1) Dirac term
  do gf=1,global_num_faces
    do p=0,ratio-1
      do i=1,NMAT
        do j=1,NMAT
          DSmat=(0d0,0d0)
          DFmat=(0d0,0d0)
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
  tmp=(0d0,0d0)
  call MPI_REDUCE(tmp_CSF,tmp,1,MPI_DOUBLE_COMPLEX, &
    MPI_SUM,0,MPI_COMM_WORLD,IERR)
  if( MYRANK==0 ) then
    CSF=CSF+tmp
  endif
  
  !! (2) mass term
  call make_XiVec_site(Xi_eta,PhiMat)
  tmp1=(0d0,0d0)
  tmp2=(0d0,0d0)
  tmp3=(0d0,0d0)
  tmp4=(0d0,0d0)
  do gf=1,global_num_faces
    do p=0,ratio-1
      do j=1,NMAT
        do i=1,NMAT
          do ls=1,num_sites
            do b=1,NMAT
              do a=1,NMAT
                tmp1(i,j,p,gf)=tmp1(i,j,p,gf) + Xi_eta(a,b,ls)*Smat(b,a,ls,i,j,p,gf)
                tmp2(i,j,p,gf)=tmp2(i,j,p,gf) + Phimat(a,b,ls)*Fmat(b,a,ls,i,j,p,gf)
                !!
                tmp3(i,j,p,gf)=tmp3(i,j,p,gf) + Xi_eta(a,b,ls)*Fmat(b,a,ls,i,j,p,gf)
                tmp4(i,j,p,gf)=tmp4(i,j,p,gf) + Phimat(a,b,ls)*Smat(b,a,ls,i,j,p,gf)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
  call MPI_REDUCE(tmp1,ttmp1,global_num_faces*NMAT*NMAT*ratio,MPI_DOUBLE_COMPLEX, &
    MPI_SUM,0,MPI_COMM_WORLD,IERR)
  call MPI_REDUCE(tmp2,ttmp2,global_num_faces*NMAT*NMAT*ratio,MPI_DOUBLE_COMPLEX, &
    MPI_SUM,0,MPI_COMM_WORLD,IERR)
  call MPI_REDUCE(tmp3,ttmp3,global_num_faces*NMAT*NMAT*ratio,MPI_DOUBLE_COMPLEX, &
    MPI_SUM,0,MPI_COMM_WORLD,IERR)
  call MPI_REDUCE(tmp4,ttmp4,global_num_faces*NMAT*NMAT*ratio,MPI_DOUBLE_COMPLEX, &
    MPI_SUM,0,MPI_COMM_WORLD,IERR)

  if( MYRANK==0 ) then
    trace=(0d0,0d0)
    do gf=1,global_num_faces
      do p=0,ratio-1
        do j=1,NMAT
          do i=1,NMAT
            trace=trace&
              -ttmp1(i,j,p,gf)*ttmp2(j,i,ratio-p-1,gf) &
              +ttmp3(i,j,p,gf)*ttmp4(j,i,ratio-p-1,gf)
          enddo
        enddo
      enddo
    enddo
    CSF=CSF + dcmplx( 0.5d0*mass_square_phi )*trace

    CSF=CSF / dcmplx( NMAT * global_num_faces )
  endif
  end subroutine calc_4fermi_in_CSFsite
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! 4-fermio terms in C_face SF_link
  subroutine calc_4fermi_in_CSFlink(CSF, Umat, Phimat, &
      Geta_eta, Gchi_eta, Geta_lambda, Gchi_lambda, ratio)
  use parallel
  use global_parameters
  use Dirac_operator , only : prod_Dirac_link1, prod_dirac_link2
  use matrix_functions, only : hermitian_conjugate, matrix_power, trace_mm, make_unit_matrix, matrix_product, matrix_3_product, matrix_commutator
  implicit none
  
  complex(kind(0d0)), intent(out) :: CSF
  complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
  complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
  complex(kind(0d0)), intent(in) :: Geta_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_sites) 
  complex(kind(0d0)), intent(in) :: Gchi_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_sites) 
  complex(kind(0d0)), intent(in) :: Geta_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_links) 
  complex(kind(0d0)), intent(in) :: Gchi_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_links) 
  integer, intent(in) :: ratio
  
  complex(kind(0d0)) :: Seta(1:NMAT,1:NMAT,1:num_necessary_sites,1:NMAT,1:NMAT,0:ratio-1,1:global_num_faces)
  complex(kind(0d0)) :: Feta(1:NMAT,1:NMAT,1:num_necessary_sites,1:NMAT,1:NMAT,0:ratio-1,1:global_num_faces)
  complex(kind(0d0)) :: Slambda(1:NMAT,1:NMAT,1:num_necessary_links,1:NMAT,1:NMAT,0:ratio-1,1:global_num_faces)
  complex(kind(0d0)) :: Flambda(1:NMAT,1:NMAT,1:num_necessary_links,1:NMAT,1:NMAT,0:ratio-1,1:global_num_faces)
  complex(kind(0d0)) :: phibar_p(1:NMAT,1:NMAT,0:ratio)
 
  complex(kind(0d0)) :: tmp1(1:NMAT,1:NMAT,0:ratio,1:global_num_faces)
  complex(kind(0d0)) :: tmp2(1:NMAT,1:NMAT,0:ratio,1:global_num_faces)
  complex(kind(0d0)) :: tmp3(1:NMAT,1:NMAT,0:ratio,1:global_num_faces)
  complex(kind(0d0)) :: tmp4(1:NMAT,1:NMAT,0:ratio,1:global_num_faces)
  complex(kind(0d0)) :: ttmp1(1:NMAT,1:NMAT,0:ratio,1:global_num_faces)
  complex(kind(0d0)) :: ttmp2(1:NMAT,1:NMAT,0:ratio,1:global_num_faces)
  complex(kind(0d0)) :: ttmp3(1:NMAT,1:NMAT,0:ratio,1:global_num_faces)
  complex(kind(0d0)) :: ttmp4(1:NMAT,1:NMAT,0:ratio,1:global_num_faces)

  !complex(kind(0d0)) :: DSeta(1:NMAT,1:NMAT,1:num_sites)
  !complex(kind(0d0)) :: DFeta(1:NMAT,1:NMAT,1:num_sites)
  complex(kind(0d0)) :: D1Feta(1:NMAT,1:NMAT,1:num_links)
  complex(kind(0d0)) :: D1Slambda(1:NMAT,1:NMAT,1:num_sites)
  complex(kind(0d0)) :: D1Seta(1:NMAT,1:NMAT,1:num_links)
  complex(kind(0d0)) :: D1Flambda(1:NMAT,1:NMAT,1:num_sites)
  !!
  complex(kind(0d0)) :: D2Slambda(1:NMAT,1:NMAT,1:num_links)
  complex(kind(0d0)) :: D2Flambda(1:NMAT,1:NMAT,1:num_links)
  complex(kind(0d0)) :: Xi_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
  complex(kind(0d0)) :: trace,tmp

  !complex(kind(0d0)), allocatable :: Seta(:,:,:,:,:,:,:) 
  !complex(kind(0d0)), allocatable :: Feta(:,:,:,:,:,:,:) 
  !complex(kind(0d0)), allocatable :: Slambda(:,:,:,:,:,:,:) 
  !complex(kind(0d0)), allocatable :: Flambda(:,:,:,:,:,:,:) 
  !complex(kind(0d0)), allocatable :: phibar_p(:,:,:)
  !complex(kind(0d0)), allocatable :: tmp1(:,:,:,:)
  !complex(kind(0d0)), allocatable :: tmp2(:,:,:,:)
  !complex(kind(0d0)), allocatable :: tmp3(:,:,:,:)
  !complex(kind(0d0)), allocatable :: tmp4(:,:,:,:)
  !complex(kind(0d0)), allocatable :: ttmp1(:,:,:,:)
  !complex(kind(0d0)), allocatable :: ttmp2(:,:,:,:)
  !complex(kind(0d0)), allocatable :: ttmp3(:,:,:,:)
  !complex(kind(0d0)), allocatable :: ttmp4(:,:,:,:)
  
  integer :: ls,ll,lf,gs,gf
  integer :: tag, rank
  integer :: i,j,k,l,p,a,b
  complex(kind(0d0)) :: tmp_CSF
  
  !ratio = (NMAT*NMAT-1)*(global_num_sites-global_num_links+global_num_faces)/2
  !allocate( Seta(1:NMAT,1:NMAT,1:num_necessary_sites,1:NMAT,1:NMAT,0:ratio-1,1:global_num_faces) )
  !allocate( Feta(1:NMAT,1:NMAT,1:num_necessary_sites,1:NMAT,1:NMAT,0:ratio-1,1:global_num_faces) )
  !allocate( Slambda(1:NMAT,1:NMAT,1:num_necessary_links,1:NMAT,1:NMAT,0:ratio-1,1:global_num_faces) )
  !allocate( Flambda(1:NMAT,1:NMAT,1:num_necessary_links,1:NMAT,1:NMAT,0:ratio-1,1:global_num_faces) )
  !allocate( phibar_p(1:NMAT,1:NMAT,0:ratio) )
  !
  !allocate( tmp1(1:NMAT,1:NMAT,0:ratio,1:global_num_faces) )
  !allocate( tmp2(1:NMAT,1:NMAT,0:ratio,1:global_num_faces) )
  !allocate( tmp3(1:NMAT,1:NMAT,0:ratio,1:global_num_faces) )
  !allocate( tmp4(1:NMAT,1:NMAT,0:ratio,1:global_num_faces) )
  !allocate( ttmp1(1:NMAT,1:NMAT,0:ratio,1:global_num_faces) )
  !allocate( ttmp2(1:NMAT,1:NMAT,0:ratio,1:global_num_faces) )
  !allocate( ttmp3(1:NMAT,1:NMAT,0:ratio,1:global_num_faces) )
  !allocate( ttmp4(1:NMAT,1:NMAT,0:ratio,1:global_num_faces) )
  !! 
  Seta=(0d0,0d0)
  Feta=(0d0,0d0)
  Slambda=(0d0,0d0)
  Flambda=(0d0,0d0)
  do gf=1,global_num_faces
    !rank=local_face_of_global(gf)%rank_
    !lf=local_face_of_global(gf)%label_
    !ls=sites_in_f(lf)%label_(1)
    gs=global_sites_in_f(gf)%label_(1)
  
    !! phibar_p = phibar^{0...r}(:,:,gf)
    call make_phibar_p(phibar_p,PhiMat,ratio,gf)
  
    !! prepare SMAT and FMAT
    do p=0,ratio-1
      do j=1,NMAT
        do i=1,NMAT
          do ls=1,num_sites
            call matrix_product(Seta(:,:,ls,i,j,p,gf),&
              phibar_p(:,:,p), Geta_eta(:,:,i,j,gs,ls))
            call matrix_product(Feta(:,:,ls,i,j,p,gf),&
              phibar_p(:,:,p), Gchi_eta(:,:,i,j,gf,ls))
          enddo !ls
          do ll=1,num_links
            call matrix_product(Slambda(:,:,ll,i,j,p,gf),&
              phibar_p(:,:,p), Geta_lambda(:,:,i,j,gs,ll))
            call matrix_product(Flambda(:,:,ll,i,j,p,gf),&
              phibar_p(:,:,p), Gchi_lambda(:,:,i,j,gf,ll))
          enddo !ll
        enddo
      enddo
    enddo
  enddo
  call syncronize_sites_large(Seta,ratio)
  call syncronize_sites_large(Feta,ratio)
  call syncronize_links_large(Slambda,ratio)
  call syncronize_links_large(Flambda,ratio)
  
  tmp_CSF=(0d0,0d0)
  CSF=(0d0,0d0)
  !! (1) Dirac term
  do gf=1,global_num_faces
    do p=0,ratio-1
      do i=1,NMAT
        do j=1,NMAT
          !! link1
          D1Seta=(0d0,0d0)
          D1Feta=(0d0,0d0)
          D1Slambda=(0d0,0d0)
          D1Flambda=(0d0,0d0)
          call prod_Dirac_link1(D1Slambda,D1Feta,PhiMat,Umat,& 
            Feta(:,:,:,j,i,ratio-p-1,gf),Slambda(:,:,:,i,j,p,gf)) !! it does not affect but fixed here (j,i) -> (i,j) in Slambda
          call prod_Dirac_link1(D1Flambda,D1Seta,PhiMat,Umat,&
            Seta(:,:,:,j,i,p,gf),Flambda(:,:,:,i,j,ratio-p-1,gf)) !! it does not affect but fixed here (j,i) -> (i,j) in Flambda
          do ll=1,num_links
            do b=1,NMAT
              do a=1,NMAT
                tmp_CSF=tmp_CSF &
                  - Slambda(a,b,ll,i,j,p,gf) * D1Feta(b,a,ll) &
                  + Flambda(a,b,ll,i,j,ratio-p-1,gf) * D1Seta(b,a,ll) 
              enddo
            enddo
          enddo !ll
          !! link2
          D2Slambda=(0d0,0d0)
          D2Flambda=(0d0,0d0)
          call prod_Dirac_link2(D2Slambda,PhiMat,Umat,&
            Slambda(:,:,:,j,i,p,gf))
          call prod_Dirac_link2(D2Flambda,PhiMat,Umat,&
            Flambda(:,:,:,j,i,ratio-p-1,gf))
          do ll=1,num_links
            do b=1,NMAT
              do a=1,NMAT
                tmp_CSF=tmp_CSF &
                  - (0.5d0,0d0) * Slambda(a,b,ll,i,j,p,gf) * D2Flambda(b,a,ll) &
                  + (0.5d0,0d0) * Flambda(a,b,ll,i,j,ratio-p-1,gf) * D2Slambda(b,a,ll) 
              enddo
            enddo
          enddo !ll
        enddo
      enddo
    enddo
  enddo
  tmp=(0d0,0d0)
  call MPI_REDUCE(tmp_CSF,tmp,1,MPI_DOUBLE_COMPLEX, &
    MPI_SUM,0,MPI_COMM_WORLD,IERR)
  if( MYRANK==0 ) then
    CSF=CSF+tmp
  endif

  
  !! (2) mass term
  call make_XiVec_link(Xi_lambda,Umat,PhiMat)
  tmp1=(0d0,0d0)
  tmp2=(0d0,0d0)
  tmp3=(0d0,0d0)
  tmp4=(0d0,0d0)
  do gf=1,global_num_faces
    do p=0,ratio-1
      do j=1,NMAT
        do i=1,NMAT
          do ls=1,num_sites
            do b=1,NMAT
              do a=1,NMAT
                ! Phi.Feta
                tmp2(i,j,p,gf)=tmp2(i,j,p,gf) + Phimat(a,b,ls)*Feta(b,a,ls,i,j,p,gf)
                ! Phi.Seta
                tmp4(i,j,p,gf)=tmp4(i,j,p,gf) + Phimat(a,b,ls)*Seta(b,a,ls,i,j,p,gf)
              enddo
            enddo
          enddo
          do ll=1,num_links
            do b=1,NMAT
              do a=1,NMAT
                ! Xi.Slambda
                tmp1(i,j,p,gf)=tmp1(i,j,p,gf) + Xi_lambda(a,b,ll)*Slambda(b,a,ls,i,j,p,gf)
                ! Xi.Flambda
                tmp3(i,j,p,gf)=tmp3(i,j,p,gf) + Xi_lambda(a,b,ll)*Flambda(b,a,ls,i,j,p,gf)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  call MPI_REDUCE(tmp1,ttmp1,global_num_faces*NMAT*NMAT*ratio,MPI_DOUBLE_COMPLEX, &
    MPI_SUM,0,MPI_COMM_WORLD,IERR)
  call MPI_REDUCE(tmp2,ttmp2,global_num_faces*NMAT*NMAT*ratio,MPI_DOUBLE_COMPLEX, &
    MPI_SUM,0,MPI_COMM_WORLD,IERR)
  call MPI_REDUCE(tmp3,ttmp3,global_num_faces*NMAT*NMAT*ratio,MPI_DOUBLE_COMPLEX, &
    MPI_SUM,0,MPI_COMM_WORLD,IERR)
  call MPI_REDUCE(tmp4,ttmp4,global_num_faces*NMAT*NMAT*ratio,MPI_DOUBLE_COMPLEX, &
    MPI_SUM,0,MPI_COMM_WORLD,IERR)

  if( MYRANK==0 ) then
    trace=(0d0,0d0)
    do gf=1,global_num_faces
      do p=0,ratio-1
        do j=1,NMAT
          do i=1,NMAT
            trace=trace&
              -ttmp1(i,j,p,gf)*ttmp2(j,i,ratio-p-1,gf) &
              +ttmp3(i,j,p,gf)*ttmp4(j,i,ratio-p-1,gf)
          enddo
        enddo
      enddo
    enddo
    CSF=CSF + dcmplx( 0.5d0*mass_square_phi )*trace
  
    CSF=CSF / dcmplx( NMAT * global_num_faces )
  endif
  end subroutine calc_4fermi_in_CSFlink
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! 4-fermio terms in C_face SF_face
  subroutine calc_4fermi_in_CSFface(CSF, Umat, Phimat, &
      Geta_eta, Gchi_eta, Geta_chi, Gchi_chi, Geta_lambda, Gchi_lambda, ratio )
  use parallel
  use global_parameters
  use Dirac_operator , only : prod_Dirac_face1, Dirac_omega_adm
  use matrix_functions, only : hermitian_conjugate, matrix_power, trace_mm, make_unit_matrix, matrix_product, matrix_product2, matrix_3_product, matrix_commutator
  implicit none
  
  complex(kind(0d0)), intent(out) :: CSF
  complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
  complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
  complex(kind(0d0)), intent(in) :: Geta_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_sites) 
  complex(kind(0d0)), intent(in) :: Gchi_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_sites) 
  complex(kind(0d0)), intent(in) :: Geta_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_faces) 
  complex(kind(0d0)), intent(in) :: Gchi_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_faces) 
  complex(kind(0d0)), intent(in) :: Geta_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_links) 
  complex(kind(0d0)), intent(in) :: Gchi_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_links) 
  integer, intent(in) :: ratio
  
  complex(kind(0d0)) :: Seta(1:NMAT,1:NMAT,1:num_necessary_sites,1:NMAT,1:NMAT,0:ratio-1,1:global_num_faces)
  complex(kind(0d0)) :: Feta(1:NMAT,1:NMAT,1:num_necessary_sites,1:NMAT,1:NMAT,0:ratio-1,1:global_num_faces)
  complex(kind(0d0)) :: Schi(1:NMAT,1:NMAT,1:num_necessary_faces,1:NMAT,1:NMAT,0:ratio-1,1:global_num_faces)
  complex(kind(0d0)) :: Fchi(1:NMAT,1:NMAT,1:num_necessary_faces,1:NMAT,1:NMAT,0:ratio-1,1:global_num_faces)
  complex(kind(0d0)) :: Slambda(1:NMAT,1:NMAT,1:num_necessary_links,1:NMAT,1:NMAT,0:ratio-1,1:global_num_faces)
  complex(kind(0d0)) :: Flambda(1:NMAT,1:NMAT,1:num_necessary_links,1:NMAT,1:NMAT,0:ratio-1,1:global_num_faces)
  complex(kind(0d0)) :: phibar_p(1:NMAT,1:NMAT,0:ratio)
  
  complex(kind(0d0)) :: tmp1(1:NMAT,1:NMAT,0:ratio,1:global_num_faces)
  complex(kind(0d0)) :: tmp2(1:NMAT,1:NMAT,0:ratio,1:global_num_faces)
  complex(kind(0d0)) :: tmp3(1:NMAT,1:NMAT,0:ratio,1:global_num_faces)
  complex(kind(0d0)) :: tmp4(1:NMAT,1:NMAT,0:ratio,1:global_num_faces)
  complex(kind(0d0)) :: ttmp1(1:NMAT,1:NMAT,0:ratio,1:global_num_faces)
  complex(kind(0d0)) :: ttmp2(1:NMAT,1:NMAT,0:ratio,1:global_num_faces)
  complex(kind(0d0)) :: ttmp3(1:NMAT,1:NMAT,0:ratio,1:global_num_faces)
  complex(kind(0d0)) :: ttmp4(1:NMAT,1:NMAT,0:ratio,1:global_num_faces)
  !complex(kind(0d0)), allocatable :: Seta(:,:,:,:,:,:,:) 
  !complex(kind(0d0)), allocatable :: Feta(:,:,:,:,:,:,:) 
  !complex(kind(0d0)), allocatable :: Schi(:,:,:,:,:,:,:) 
  !complex(kind(0d0)), allocatable :: Fchi(:,:,:,:,:,:,:) 
  !complex(kind(0d0)), allocatable :: Slambda(:,:,:,:,:,:,:) 
  !complex(kind(0d0)), allocatable :: Flambda(:,:,:,:,:,:,:) 
  !complex(kind(0d0)), allocatable :: phibar_p(:,:,:)
  complex(kind(0d0)) :: DSchi(1:NMAT,1:NMAT,1:num_faces)
  complex(kind(0d0)) :: DFchi(1:NMAT,1:NMAT,1:num_faces)
  complex(kind(0d0)) :: DSlambda(1:NMAT,1:NMAT,1:num_links)
  complex(kind(0d0)) :: DFlambda(1:NMAT,1:NMAT,1:num_links)
  complex(kind(0d0)) :: Xi_chi(1:NMAT,1:NMAT,1:num_necessary_faces)
  complex(kind(0d0)) :: trace, tmp
  !complex(kind(0d0)), allocatable :: tmp1(:,:,:,:)
  !complex(kind(0d0)), allocatable :: tmp2(:,:,:,:)
  !complex(kind(0d0)), allocatable :: tmp3(:,:,:,:)
  !complex(kind(0d0)), allocatable :: tmp4(:,:,:,:)
  !complex(kind(0d0)), allocatable :: ttmp1(:,:,:,:)
  !complex(kind(0d0)), allocatable :: ttmp2(:,:,:,:)
  !complex(kind(0d0)), allocatable :: ttmp3(:,:,:,:)
  !complex(kind(0d0)), allocatable :: ttmp4(:,:,:,:)
  
  integer :: ls,ll,lf,gs,gf,lf2
  integer :: tag, rank
  integer :: i,j,k,l,p,a,b
  complex(kind(0d0)) :: tmp_CSF
  
  !ratio = (NMAT*NMAT-1)*(global_num_sites-global_num_links+global_num_faces)/2
  !allocate( Seta(1:NMAT,1:NMAT,1:num_necessary_sites,1:NMAT,1:NMAT,0:ratio-1,1:global_num_faces) )
  !allocate( Feta(1:NMAT,1:NMAT,1:num_necessary_sites,1:NMAT,1:NMAT,0:ratio-1,1:global_num_faces) )
  !allocate( Schi(1:NMAT,1:NMAT,1:num_necessary_faces,1:NMAT,1:NMAT,0:ratio-1,1:global_num_faces) )
  !allocate( Fchi(1:NMAT,1:NMAT,1:num_necessary_faces,1:NMAT,1:NMAT,0:ratio-1,1:global_num_faces) )
  !allocate( Slambda(1:NMAT,1:NMAT,1:num_necessary_links,1:NMAT,1:NMAT,0:ratio-1,1:global_num_faces) )
  !allocate( Flambda(1:NMAT,1:NMAT,1:num_necessary_links,1:NMAT,1:NMAT,0:ratio-1,1:global_num_faces) )
  !allocate( phibar_p(1:NMAT,1:NMAT,0:ratio) )
  
  !allocate( tmp1(1:NMAT,1:NMAT,0:ratio,1:global_num_faces) )
  !allocate( tmp2(1:NMAT,1:NMAT,0:ratio,1:global_num_faces) )
  !allocate( tmp3(1:NMAT,1:NMAT,0:ratio,1:global_num_faces) )
  !allocate( tmp4(1:NMAT,1:NMAT,0:ratio,1:global_num_faces) )
  !allocate( ttmp1(1:NMAT,1:NMAT,0:ratio,1:global_num_faces) )
  !allocate( ttmp2(1:NMAT,1:NMAT,0:ratio,1:global_num_faces) )
  !allocate( ttmp3(1:NMAT,1:NMAT,0:ratio,1:global_num_faces) )
  !allocate( ttmp4(1:NMAT,1:NMAT,0:ratio,1:global_num_faces) )
  !! 
  Seta=(0d0,0d0)
  Feta=(0d0,0d0)
  Slambda=(0d0,0d0)
  Flambda=(0d0,0d0)
  do gf=1,global_num_faces
    !rank=local_face_of_global(gf)%rank_
    !lf=local_face_of_global(gf)%label_
    !ls=sites_in_f(lf)%label_(1)
    gs=global_sites_in_f(gf)%label_(1)
  
    !! phibar_p = phibar^{0...r}(:,:,gf)
    call make_phibar_p(phibar_p,PhiMat,ratio,gf)
  
    !! prepare SMAT and FMAT
    do p=0,ratio-1
      do j=1,NMAT
        do i=1,NMAT
          do ls=1,num_sites
            call matrix_product2(Seta(:,:,ls,i,j,p,gf),&
              phibar_p(:,:,p), Geta_eta(:,:,i,j,gs,ls), NMAT)
            call matrix_product2(Feta(:,:,ls,i,j,p,gf),&
              phibar_p(:,:,p), Gchi_eta(:,:,i,j,gf,ls), NMAT)
          enddo !ls
          do lf2=1,num_faces
            call matrix_product2(Schi(:,:,lf2,i,j,p,gf),&
              phibar_p(:,:,p), Geta_chi(:,:,i,j,gs,lf2), NMAT)
            call matrix_product2(Fchi(:,:,lf2,i,j,p,gf),&
              phibar_p(:,:,p), Gchi_chi(:,:,i,j,gf,lf2), NMAT)
          enddo !lf2
          do ll=1,num_links
            call matrix_product2(Slambda(:,:,ll,i,j,p,gf),&
              phibar_p(:,:,p), Geta_lambda(:,:,i,j,gs,ll), NMAT)
            call matrix_product2(Flambda(:,:,ll,i,j,p,gf),&
              phibar_p(:,:,p), Gchi_lambda(:,:,i,j,gf,ll), NMAT)
          enddo !ll
        enddo
      enddo
    enddo
  enddo
  call syncronize_sites_large(Seta,ratio)
  call syncronize_sites_large(Feta,ratio)
  call syncronize_faces_large(Schi,ratio)
  call syncronize_faces_large(Fchi,ratio)
  call syncronize_links_large(Slambda,ratio)
  call syncronize_links_large(Flambda,ratio)
  
  tmp_CSF=(0d0,0d0)
  CSF=(0d0,0d0)
  !! (1) Dirac term
  do gf=1,global_num_faces
    do p=0,ratio-1
      do i=1,NMAT
        do j=1,NMAT
          !! face1
          DSchi=(0d0,0d0)
          DFchi=(0d0,0d0)
          call prod_Dirac_face1(DFchi,PhiMat,&
            Fchi(:,:,:,j,i,ratio-p-1,gf))
          call prod_Dirac_face1(DSchi,PhiMat,&
            Schi(:,:,:,j,i,p,gf))
          do lf=1,num_faces
            do b=1,NMAT
              do a=1,NMAT
                tmp_CSF=tmp_CSF &
                  - (0.5d0,0d0) * Schi(a,b,lf,i,j,p,gf) * DFchi(b,a,lf) &
                  + (0.5d0,0d0) * Fchi(a,b,lf,i,j,ratio-p-1,gf) * DSchi(b,a,lf) 
              enddo
            enddo
          enddo !lf
          !! face2
          DSchi=(0d0,0d0)
          DFchi=(0d0,0d0)
          DSlambda=(0d0,0d0)
          DFlambda=(0d0,0d0)
          call Dirac_Omega_adm(DFchi,DFlambda,Umat,&
            Schi(:,:,:,j,i,p,gf), Flambda(:,:,:,j,i,ratio-p-1,gf) )
          call Dirac_Omega_adm(DSchi,DSlambda,Umat,&
            Fchi(:,:,:,j,i,ratio-p-1,gf), Slambda(:,:,:,j,i,p,gf)) 
          do lf=1,num_faces
            do b=1,NMAT
              do a=1,NMAT
                tmp_CSF=tmp_CSF &
                  -  Schi(a,b,lf,i,j,p,gf) * DFchi(b,a,lf) &
                  +  Fchi(a,b,lf,i,j,ratio-p-1,gf) * DSchi(b,a,lf) 
              enddo
            enddo
          enddo !lf
        enddo
      enddo
    enddo
  enddo
  tmp=(0d0,0d0)
  call MPI_REDUCE(tmp_CSF,tmp,1,MPI_DOUBLE_COMPLEX, &
    MPI_SUM,0,MPI_COMM_WORLD,IERR)
  if( MYRANK==0 ) then
    CSF=CSF+tmp
  endif
  
  !! (2) mass term
  call make_XiVec_face(Xi_chi,Umat)
  tmp1=(0d0,0d0)
  tmp2=(0d0,0d0)
  tmp3=(0d0,0d0)
  tmp4=(0d0,0d0)
  do gf=1,global_num_faces
    do p=0,ratio-1
      do j=1,NMAT
        do i=1,NMAT
          do ls=1,num_sites
            do b=1,NMAT
              do a=1,NMAT
                ! Phi.Fchi
                tmp2(i,j,p,gf)=tmp2(i,j,p,gf) + Phimat(a,b,ls)*Feta(b,a,ls,i,j,p,gf)
                ! Phi.Schi
                tmp4(i,j,p,gf)=tmp4(i,j,p,gf) + Phimat(a,b,ls)*Seta(b,a,ls,i,j,p,gf)
              enddo
            enddo
          enddo
          !!!!!!!!!
          do lf=1,num_faces
            do b=1,NMAT
              do a=1,NMAT
                ! Xi.Schi
                tmp1(i,j,p,gf)=tmp1(i,j,p,gf) + Xi_chi(a,b,lf)*Schi(b,a,lf,i,j,p,gf)
                ! Xi.Fchi   
                tmp3(i,j,p,gf)=tmp3(i,j,p,gf) + Xi_chi(a,b,lf)*Fchi(b,a,lf,i,j,p,gf)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
  call MPI_REDUCE(tmp1,ttmp1,global_num_faces*NMAT*NMAT*ratio,MPI_DOUBLE_COMPLEX, &
    MPI_SUM,0,MPI_COMM_WORLD,IERR)
  call MPI_REDUCE(tmp2,ttmp2,global_num_faces*NMAT*NMAT*ratio,MPI_DOUBLE_COMPLEX, &
    MPI_SUM,0,MPI_COMM_WORLD,IERR)
  call MPI_REDUCE(tmp3,ttmp3,global_num_faces*NMAT*NMAT*ratio,MPI_DOUBLE_COMPLEX, &
    MPI_SUM,0,MPI_COMM_WORLD,IERR)
  call MPI_REDUCE(tmp4,ttmp4,global_num_faces*NMAT*NMAT*ratio,MPI_DOUBLE_COMPLEX, &
    MPI_SUM,0,MPI_COMM_WORLD,IERR)

  if( MYRANK==0 ) then
    trace=(0d0,0d0)
    do gf=1,global_num_faces
      do p=0,ratio-1
        do j=1,NMAT
          do i=1,NMAT
            trace=trace&
              -ttmp1(i,j,p,gf)*ttmp2(j,i,ratio-p-1,gf) &
              +ttmp3(i,j,p,gf)*ttmp4(j,i,ratio-p-1,gf)
          enddo
        enddo
      enddo
    enddo
    CSF=CSF + dcmplx( 0.5d0*mass_square_phi )*trace
  
    CSF=CSF / dcmplx( NMAT * global_num_faces )
  endif
  
  
  
  
  end subroutine calc_4fermi_in_CSFface
  
  
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine make_phibar_p(phibar_p,PhiMat,ratio,gf)
  use global_parameters
  use parallel
  use matrix_functions, only : hermitian_conjugate, make_unit_matrix, matrix_product
  implicit none
  
  integer, intent(in) :: ratio,gf
  complex(kind(0d0)), intent(out) :: phibar_p(1:NMAT,1:NMAT,0:ratio)
  complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
  
  integer :: gs,lf,ls,rank,k
  
  gs=global_sites_in_f(gf)%label_(1)
  rank=local_face_of_global(gf)%rank_
  if( rank == MYRANK) then
    lf=local_face_of_global(gf)%label_
    ls=sites_in_f(lf)%label_(1)
  endif
  
  phibar_p=(0d0,0d0)
  
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
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine syncronize_sites_large(eta,pekepeke)
  use global_parameters
  use parallel
  implicit none

  integer, intent(in) :: pekepeke
  complex(kind(0d0)), intent(inout) :: eta(1:NMAT,1:NMAT,1:num_necessary_sites,1:NMAT,1:NMAT,0:pekepeke-1,1:global_num_faces) 
  
  complex(kind(0d0)) :: tmp_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,0:pekepeke-1,1:global_num_faces,num_necessary_sites)
  integer :: s_send
  integer :: s_recv
  integer :: local, rank, tag
  !integer, allocatable :: ISEND(:), IRECV(:) ! for MPI_WAIT 
  integer :: ISEND(1:num_send_sites)
  integer :: IRECV(1:num_recv_sites)
  
  integer :: ls, gf, r, i,j, num
  
  num=NMAT**4*pekepeke*global_num_faces
  !! prepare tmp_lambda
  tmp_eta=(0d0,0d0)
  do ls=1,num_sites
    do gf=1,global_num_faces
      do r=0, pekepeke-1
        do j=1,NMAT
          do i=1,NMAT
            tmp_eta(:,:,i,j,r,gf,ls)=eta(:,:,ls,i,j,r,gf)
          enddo
        enddo
      enddo
    enddo
  enddo
  !!!!!!!!
  !allocate(ISEND(1:num_send_sites))
  !allocate(IRECV(1:num_recv_sites))
  do s_send=1,num_send_sites
    local=send_sites(s_send)%label_
    rank=send_sites(s_send)%rank_
    tag=10000*rank + global_site_of_local(local)
  
    call MPI_ISEND(tmp_eta(:,:,:,:,:,:,local),num,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,ISEND(s_send),IERR)
  enddo
  
  do s_recv=1,num_recv_sites
    local=recv_sites(s_recv)%label_
    rank=recv_sites(s_recv)%rank_
    tag=10000*MYRANK + global_site_of_local(local)
  
    call MPI_IRECV(tmp_eta(:,:,:,:,:,:,local),num,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,IRECV(s_recv),IERR)
  enddo
  
  do s_send=1,num_send_sites
    call MPI_WAIT(ISEND(s_send),ISTATUS,IERR)
  enddo
  do s_recv=1,num_recv_sites
    call MPI_WAIT(IRECV(s_recv),ISTATUS,IERR)
  enddo
  
  !deallocate(ISEND, IRECV)
  
  !! from tmp_eta to eta
  do ls=1,num_necessary_sites
    do gf=1,global_num_faces
      do r=0, pekepeke-1
        do j=1,NMAT
          do i=1,NMAT
            eta(:,:,ls,i,j,r,gf)=tmp_eta(:,:,i,j,r,gf,ls)
          enddo
        enddo
      enddo
    enddo
  enddo
  end subroutine syncronize_sites_large

  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine syncronize_links_large(lambda,pekepeke)
  use global_parameters
  use parallel
  implicit none

  integer, intent(in) :: pekepeke
  complex(kind(0d0)), intent(inout) :: lambda(1:NMAT,1:NMAT,1:num_necessary_links,1:NMAT,1:NMAT,0:pekepeke-1,1:global_num_faces) 
  
  complex(kind(0d0)) :: tmp_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,0:pekepeke-1,1:global_num_faces, 1:num_necessary_links)
  integer :: s_send
  integer :: s_recv
  integer :: local, rank, tag
  !integer, allocatable :: ISEND(:) ! for MPI_WAIT 
  !integer, allocatable :: IRECV(:) ! for MPI_WAIT 
  integer :: ISEND(1:num_send_links) ! for MPI_WAIT 
  integer :: IRECV(1:num_recv_links) ! for MPI_WAIT 
  complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
  
  integer :: ll, gf, r, i,j, num

  
  num=NMAT**4*pekepeke*global_num_faces
  !! prepare tmp_lambda
  tmp_lambda=(0d0,0d0)
  do ll=1,num_links
    do gf=1,global_num_faces
      do r=0, pekepeke-1
        do j=1,NMAT
          do i=1,NMAT
            tmp_lambda(:,:,i,j,r,gf,ll)=lambda(:,:,ll,i,j,r,gf)
          enddo
        enddo
      enddo
    enddo
  enddo
  !!!!!!!!
  !allocate(ISEND(1:num_send_links))
  !allocate(IRECV(1:num_recv_links))
  do s_send=1,num_send_links
    local=send_links(s_send)%label_
    rank=send_links(s_send)%rank_
    tag=10000*rank + global_link_of_local(local)
  
    call MPI_ISEND(tmp_lambda(:,:,:,:,:,:,local),num,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,ISEND(s_send),IERR)
  enddo
  
  do s_recv=1,num_recv_links
    local=recv_links(s_recv)%label_
    rank=recv_links(s_recv)%rank_
    tag=10000*MYRANK + global_link_of_local(local)
  
    !call MPI_RECV(lambda(:,:,local),NMAT*NMAT,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,ISTATUS,IERR)
    call MPI_IRECV(tmp_lambda(:,:,:,:,:,:,local),num,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,IRECV(s_recv),IERR)
  enddo
  
  do s_send=1,num_send_links
    call MPI_WAIT(ISEND(s_send),ISTATUS,IERR)
  enddo
  do s_recv=1,num_recv_links
    call MPI_WAIT(IRECV(s_recv),ISTATUS,IERR)
  enddo
  
  !deallocate(ISEND, IRECV)

  !! from tmp_lambda to tmp
  do ll=1,num_necessary_links
    do gf=1,global_num_faces
      do r=0, pekepeke-1
        do j=1,NMAT
          do i=1,NMAT
            lambda(:,:,ll,i,j,r,gf)=tmp_lambda(:,:,i,j,r,gf,ll)
          enddo
        enddo
      enddo
    enddo
  enddo
  end subroutine syncronize_links_large

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine syncronize_faces_large(chi,pekepeke)
  use global_parameters
  use parallel
  implicit none

  integer, intent(in) :: pekepeke
  complex(kind(0d0)), intent(inout) :: chi(1:NMAT,1:NMAT,1:num_necessary_faces,1:NMAT,1:NMAT,0:pekepeke-1,1:global_num_faces) 
  
  complex(kind(0d0)) :: tmp_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,0:pekepeke-1,1:global_num_faces, 1:num_necessary_faces)
  integer :: s_send
  integer :: s_recv
  integer :: local, rank, tag
  !integer, allocatable :: ISEND(:), IRECV(:) ! for MPI_WAIT 
  integer :: ISEND(1:num_send_faces)
  integer :: IRECV(1:num_recv_faces)
  complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
  
  integer :: lf, gf, r, i,j, num
  
  num=NMAT**4*pekepeke*global_num_faces
  !! prepare tmp_chi
  tmp_chi=(0d0,0d0)
  do lf=1,num_faces
    do gf=1,global_num_faces
      do r=0, pekepeke-1
        do j=1,NMAT
          do i=1,NMAT
            tmp_chi(:,:,i,j,r,gf,lf)=chi(:,:,lf,i,j,r,gf)
          enddo
        enddo
      enddo
    enddo
  enddo
  !!!!!!!!
  !allocate(ISEND(1:num_send_faces))
  !allocate(IRECV(1:num_recv_faces))
  do s_send=1,num_send_faces
    local=send_faces(s_send)%label_
    rank=send_faces(s_send)%rank_
    tag=10000*rank + global_face_of_local(local)
  
    call MPI_ISEND(tmp_chi(:,:,:,:,:,:,local),num,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,ISEND(s_send),IERR)
  enddo
  
  do s_recv=1,num_recv_faces
    local=recv_faces(s_recv)%label_
    rank=recv_faces(s_recv)%rank_
    tag=10000*MYRANK + global_face_of_local(local)
  
    !call MPI_RECV(chi(:,:,local),NMAT*NMAT,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,ISTATUS,IERR)
    call MPI_IRECV(tmp_chi(:,:,:,:,:,:,local),num,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,IRECV(s_recv),IERR)
  enddo
  
  do s_send=1,num_send_faces
    call MPI_WAIT(ISEND(s_send),ISTATUS,IERR)
  enddo
  do s_recv=1,num_recv_faces
    call MPI_WAIT(IRECV(s_recv),ISTATUS,IERR)
  enddo
  
  !deallocate(ISEND, IRECV)

  !! from tmp_chi to tmp
  do lf=1,num_necessary_faces
    do gf=1,global_num_faces
      do r=0, pekepeke-1
        do j=1,NMAT
          do i=1,NMAT
            chi(:,:,lf,i,j,r,gf)=tmp_chi(:,:,i,j,r,gf,lf)
          enddo
        enddo
      enddo
    enddo
  enddo
  end subroutine syncronize_faces_large
end subroutine calc_face_compensator
