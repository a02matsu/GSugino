!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
subroutine calc_Qexact_operators(&
    opS1, opL1, opL2, opF1, opF2, &
    Umat,PhiMat,&
    Geta_eta, Geta_lambda, Geta_chi, Gchi_eta, Gchi_lambda, Gchi_chi) 
use parallel
use global_parameters
use matrix_functions, only : matrix_product, make_unit_matrix
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
integer :: ratio, k
complex(kind(0d0)), allocatable :: phibar_p(:,:,:,:)
integer :: ls

ratio = (NMAT*NMAT-1)*(global_num_sites-global_num_links+global_num_faces)/2

!! preparation
allocate( phibar_p(1:NMAT,1:NMAT,0:ratio+1,1:necessary_num_sites) )

do ls=1,num_sites
  call make_unit_matrix(phibar_p(:,:,0,ls))
  do k=1,ratio+1
    call matrix_product(phibar_p(:,:,k,ls), phibar_p(:,:,k-1,ls), phimat(:,;,ls), 'N', 'C')
  enddo
enddo
call syncronize_phibar(phibar_p,ratio)

call opS1(op, PhiMat, Geta_eta, phibar_p, ratio)

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine opS1(op, PhiMat, Geta_eta, phibar_p, ratio)
  use parallel
  use global_parameters
  use matrix_functions, only : matrix_commutator, matrix_power, trace_mm, matrix_product
  implicit none

  complex(kind(0d0)), intent(out) :: op
  complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
  complex(kind(0d0)), intent(in) :: Geta_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_sites) 
  integer, intent(in) :: ratio
  complex(kind(0d0)), intent(in) :: phibar_p(1:NMAT, 1:NMAT,0:ratio+1,1:num_sites)

  complex(kind(0d0)) :: comm1(1:NMAT,1:NMAT), comm2(1:NMAT,1:NMAT), tmpmat(1:NMAT,1:NMAT)
  complex(Kind(0d0)) :: tmp, trace
  integer :: ls,gs,k,i,j,l,r

  tmp=(0d0,0d0)
  do ls=1,num_sites
    gs=global_site_of_local(ls)

    !! bosonic part
    call matrix_commutator(comm1, phimat(:,:,ls), phimat(:,:,ls), 'N','C')
    call matrix_commutator(comm2, phimat(:,:,ls), phibar_p(:,:,ratio+1,ls))
    trace=(0d0,0d0)
    call trace_mm(trace, comm1, comm2)
    tmp=tmp+trace

    !! fermionic part
    do r=0, ratio
      call matrix_product(tmpmat, Phimat(:,:,ls), phibar_p(:,:,r,ls))

      do i=1,NMAT
        do j=1,NMAT
          do k=1,NMAT
            do l=1,NMAT
              tmp = tmp - (2d0,0d0) * tmpmat(i,j) * phibar_p(k,l,ratio-p,ls) &
                * Geta_eta(l,i,j,k,gs,ls) 
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo !! ls
  
  call MPI_REDUCE(tmp,op,1,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,IERR)

  end subroutine opS1



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine opL1(op, PhiMat, Umat, Geta_lambda, Glambda_lambda, phibar_p, ratio)
  use parallel
  use global_parameters
  use matrix_functions, only : matrix_commutator, matrix_power, trace_mm, matrix_product
  implicit none

  complex(kind(0d0)), intent(out) :: op
  complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
  complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
  complex(kind(0d0)), intent(in) :: Geta_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_links) 
  complex(kind(0d0)), intent(in) :: Glambda_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_links) 
  integer, intent(in) :: ratio
  complex(kind(0d0)), intent(in) :: phibar_p(1:NMAT, 1:NMAT,0:ratio+1,1:num_necessary_sites)

  complex(Kind(0d0)) :: tmp, trace
  complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
  complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
  integer :: ls,lt,gs,ll,gl
  integer :: k,i,j,l,r

  tmp=(0d0,0d0)
  do ll=1,num_links
    gl=global_link_of_local(ll)
    ls=link_org(ll)
    lt=link_tip(ll)
    gs=global_site_of_local(ls)

    !! bosonic part
    trace=(0d0,0d0)
    call matrix_product(tmpmat1, Umat(:,:,ll), phimat(:,:,lt) )
    call matrix_product(tmpmat2, Umat(:,:,ll), phibar_p(:,:,ratio+1,ls), 'C', 'N')
    call trace_mm(trace, tmpmat1, tmpmat2)
    tmp = tmp + trace

    trace=(0d0,0d0)
    call trace_mm(trace, phimat(:,:,ls), phibar_p(:,:,ratio+1,ls))
    tmp = tmp - trace

    !! fermionic part
    trace=(0d0,0d0)
    do k=1,NMAT
      do j=1,NMAT
        do i=1,NMAT
          trace = trace + phibar_p(k,i,ratio+1,ls) &
            * Glambda_lambda(i,j,j,k,gl,ll)
        enddo
      enddo
    enddo
    tmp = tmp + trace

    trace=(0d0,0d0)
    do r=0, ratio
      do i=1,NMAT
        do j=1,NMAT
          do k=1,NMAT
            do l=1,NMAT
              trace = trace + phibar_p(i,j,r,ls) * phibar_p(k,l,ratio-r,ls) &
                * Geta_lambda(j,k,l,i,gs,ll)
            enddo
          enddo
        enddo
      enddo
    enddo
    tmp = tmp + (0d0,-1d0)*trace
  enddo

  call MPI_REDUCE(tmp,op,1,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,IERR)

  end subroutine opL1


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine opL2(op, PhiMat, Umat, Geta_lambda, Glambda_lambda, phibar_p, ratio)
  use parallel
  use global_parameters
  use matrix_functions, only : matrix_commutator, matrix_power, trace_mm, matrix_product, matrix_3_product
  implicit none

  complex(kind(0d0)), intent(out) :: op
  complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
  complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
  complex(kind(0d0)), intent(in) :: Geta_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_links) 
  complex(kind(0d0)), intent(in) :: Glambda_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_links) 
  integer, intent(in) :: ratio
  complex(kind(0d0)), intent(in) :: phibar_p(1:NMAT, 1:NMAT,0:ratio+1,1:num_necessary_sites)

  complex(Kind(0d0)) :: tmp, trace
  complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
  complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
  integer :: ls,lt,gs,ll,gl
  integer :: k,i,j,l,r


  tmp=(0d0,0d0)
  do ll=1,num_links
    gl=global_link_of_local(ll)
    ls=link_org(ll)
    lt=link_tip(ll)
    gs=global_site_of_local(ls)


    !! bosonic part
    trace=(0d0,0d0)
    call trace_mm(trace, phimat(:,:,lt), phibar_p(:,:,ratio+1,lt))
    tmp = tmp + trace

    trace=(0d0,0d0)
    call matrix_product(tmpmat1, Umat(:,:,ll), phibar_p(:,:,ratio+1,lt))
    call matrix_product(tmpmat2, Umat(:,:,ll), phimat(:,:,ls), 'C','N' )
    call trace_mm(trace, tmpmat1, tmpmat2)
    tmp = tmp - trace

    !! fermionic part
    call matrix_3_product(tmpmat1,Umat(:,:,ll),phibar_p(:,:,ratio+1,lt),Umat(:,:,ll),&
      'N','N','C')
    trace=(0d0,0d0)
    do i=1,NMAT
      do j=1,NMAT
        do k=1,NMAT
          trace = trace - Glambda_lambda(i,j,j,k,gl,ll) * tmpmat1(k,i)
        enddo
      enddo
    enddo
    tmp = tmp + trace

    trace=(0d0,0d0)
    do r=0,ratio
      call matrix_product(tmpmat1, Umat(:,:,ll), phibar_p(:,:,r,lt))
      call matrix_product(tmpmat2, phimat_p(:,:,ratio-k,lt),Umat(:,:,ll), 'N','C' )
      do i=1,NMAT
        do j=1,NMAT
          do k=1,NMAT
            do l=1,NMAT
              trace = trace &
                -(0d0,1d0)*tmpmat1(i,j)*tmpmat2(k,l)&
                 * Geta_lambda(j,k,l,i,gs,ll) 
               ! take care for the order of eta and lambda
            enddo
          enddo
        enddo
      enddo
    enddo
    tmp = tmp + trace
  enddo

  call MPI_REDUCE(tmp,op,1,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,IERR)

  end subroutine opL2

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine opF1(op, PhiMat, Umat, Geta_chi, Gchi_chi, phibar_p, ratio)
  use parallel
  use global_parameters
  use matrix_functions, only : matrix_commutator, matrix_power, trace_mm, matrix_product, matrix_3_product
  implicit none

  complex(kind(0d0)), intent(out) :: op
  complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
  complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
  complex(kind(0d0)), intent(in) :: Geta_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_faces) 
  complex(kind(0d0)), intent(in) :: Gchi_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_faces) 
  integer, intent(in) :: ratio
  complex(kind(0d0)), intent(in) :: phibar_p(1:NMAT, 1:NMAT,0:ratio+1,1:num_necessary_sites)

  complex(Kind(0d0)) :: tmp, trace
  complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
  complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
  complex(kind(0d0)) :: Ymat(1:NMAT,1:NMAT)
  complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT)

  integer :: lf,ls,lt,gs,ll,gl,gf
  integer :: k,i,j,l,r


  tmp=(0d0,0d0)
  do lf=1,num_links
    gf=global_face_of_local(lf)
    ls=sites_in_f(lf)%label_(1)
    gs=global_site_of_local(ls)

    !! Ymat
    call Make_face_variable(Uf,lf,UMAT)
    call Make_moment_map_adm(Ymat,Uf)
    Ymat = Ymat * (0d0,0.5d0)*beta_f(lf)

    !! bosonic part
    trace=(0d0,0d0)
    call matrix_product(tmpmat1,Ymat,Ymat)
    call trace_mm(trace,tmpmat1,phibar_p(:,:,ratio,ls))
    tmp = tmp + trace

    !! fermionic part
    trace=(0d0,0d0)
    do k=1,NMAT
      do j=1,NMAT
        do i=1,NMAT
          do l=1,NMAT
            trace = trace - Phimat(i,j,ls) * phibar_p(k,l,ratio,ls) &
              * Gchi_chi(l,i,j,k,gf,lf)
          enddo
        enddo
      enddo
    enddo
    tmp = tmp + trace

    trace=(0d0,0d0)
    call matrix_product(tmpmat1, phimat(:,:,ls), phibar_p(:,:,ratio,ls))
    do i=1,NMAT
      do k=1,NMAT
        do j=1,NMAT
          trace=trace+tmpmat(i,j)*Gchi_chi(j,k,k,i,gf,lf)
        enddo
      enddo
    enddo
    tmp = tmp + trace

    trace=(0d0,0d0)
    do r=0, ratio-1
      call matrix_product(tmpmat1, Ymat, phibar_p(:,:,r))
      do i=1,NMAT
        do l=1,NMAT
          do k=1,NMAT
            do j=1,NMAT
              trace=trace+tmpmat(i,j)*phibar_p(k,l,ratio-r-1,ls)&
                * Geta_chi(j,k,l,i,gs,lf)
            enddo
          enddo
        enddo
      enddo
    enddo
    tmp = tmp + trace
  enddo

  call MPI_REDUCE(tmp,op,1,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,IERR)
  end subroutine opF1


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine opF2(op, Umat, Geta_chi, Glambda_chi, phibar_p, ratio)
  use global_parameters
  use matrix_functions, only : matrix_product
  implicit none
  
  complex(kind(0d0)), intent(out) :: op
  complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
  complex(kind(0d0)), intent(in) :: Geta_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_faces) 
  complex(kind(0d0)), intent(in) :: Gchi_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_faces) 
  integer, intent(in) :: ratio
  complex(kind(0d0)), intent(in) :: phibar_p(1:NMAT, 1:NMAT,0:ratio+1,1:num_necessary_sites)
  
  complex(kind(0d0)) :: Xmat(1:NMAT,1:NMAT)
  complex(kind(0d0)) :: Ymat(1:NMAT,1:NMAT)
  complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT)
  complex(kind(0d0)) :: Usin(1:NMAT,1:NMAT)
  complex(kind(0d0)) :: Tmat(1:NMAT,1:NMAT)
  complex(kind(0d0)) :: Omega(1:NMAT,1:NMAT)
  complex(kind(0d0)) :: tmpmat1(1:NMAT,1:NMAT)
  complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
  complex(kind(0d0)) :: Bval
  complex(kind(0d0)) :: tmp
  complex(kind(0d0)) :: trace

  integer :: l_place ! place of the link in lf
  integer :: dir ! direction of ll in lf
  integer :: gs,gl,ls,lf,ll
  integer :: i,j,k,l,p
  

  tmp=(0d0,0d0)
  do lf=1,num_faces
    ls=sites_in_f(lf)%label_(1)
    gs=global_site_of_local(ls)

    !! Omega and B
    call Make_face_variable(Uf,lf,UMAT)
    call Make_moment_map_adm(Omega,Uf)
    call calc_Bval(Bval,Uf)

    !! bosonic part
    trace=(0d0,0d0)
    call matrix_product(tmpmat1, Omega,Omega)
    call trace_mm(trace, tmpmat1, phibar_p(:,:,ratio,ls))
    tmp=tmp+(0d0,0.5d0)*beta_f(lf)*trace

  
    !! fermionic part
    trace=(0d0,0d0)
    do p=0,ratio-1
      call matrix_product(tmpmat1,Omega,phibar_p(:,:,p,ls))
      do l=1,NMAT
        do k=1,NMAT
          do j=1,NMAT
            do i=1,NMAT
              trace=trace + &
                Geta_chi(i,j,k,l,gs,lf)*tmpmat1(l,i)*phibar_p(j,k,ratio-p-1,ls)
            enddo
          enddo
        enddo
      enddo
    enddo
    tmp = tmp + trace

    do l_place=1,links_in_f(lf)%num_
      ll=links_in_f(lf)%link_labels_(l_place)
      dir=links_in_f(lf)%link_dirs_(l_place)
      gl = global_link_of_local( links_in_f(lf)%link_labels_(l_place) )
      call calc_XYmat(Xmat,Ymat,lf,l_place,UMAT)
      
      call matrix_product(tmpmat1,Ymat,phibar_p(:,:,ratio,ls))
      call matrix_product(tmpmat2,Xmat,phibar_p(:,:,ratio,ls),'C','N')
      trace=(0d0,0d0)
      do l=1,NMAT
        do k=1,NMAT
          do j=1,NMAT
            do i=1,NMAT
              trace=trace+&
                Glambda_chi(i,j,i,l,gl,lf) &
                * (Xmat(l,i)*tmpmat1(j,k) + dconjg(Ymat(i,l))*tmpmat2(j,k)) &
            enddo
          enddo
        enddo
      enddo
      tmp = tmp + dcmplx(dir)/Bval*trace
      
      !! Usin
      do i=1,NMAT
        do j=1,NMAT
          tmpmat1(i,j) = Uf(i,j) - dconjg(Uf(j,i))
        enddo
      enddo
      call matrix_product(Usin,tmpmat1,phibar_p(:,:,ratio,ls))
      !! Tmat
      Tmat=(0d0,0d0)
      call matrix_product(Tmat,Ymat,Xmat)
      call matrix_product(Tmat,Xmat,Ymat,'C','C',(-1d0,0d0),'ADD')
      !! Tmat.Dinv
      tmpmat1=(0d0,0d0)
      do i=1,NMAT
        do j=1,NMAT
          do k=1,NMAT
            do l=1,NMAT
              tmpmat1(i,j)=tmpmat1(i,j) &
                +Glambda_chi(i,j,k,l,gl,lf)*Usin(l,k)
            enddo
          enddo
        enddo
      enddo
      trace=(0d0,0d0)
      call trace_mm(trace,Tmat,tmpmat1)
      tmp = tmp + dcmplx(dir)/(Bval*Bval*dcmplx(e_max*e_max))
    enddo

  call MPI_REDUCE(tmp,op,1,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,IERR)
  end subroutine opF2


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine syncronize_phibar(eta,ratio)
  use global_parameters
  use parallel
  implicit none

  integer, intent(in) :: ratio
  complex(kind(0d0)), intent(inout) :: eta(1:NMAT,1:NMAT,0:ratio+1,1:num_necessary_sites)
  
  integer :: s_send
  integer :: s_recv
  integer :: local, rank, tag
  integer :: ISEND(1:num_send_sites)
  integer :: IRECV(1:num_recv_sites)
  
  integer :: ls, gf, r, i,j, num
  
  num = NMAT*NMAT*(ratio+2)
  !!!!!!!!
  do s_send=1,num_send_sites
    local=send_sites(s_send)%label_
    rank=send_sites(s_send)%rank_
    tag=10000*rank + global_site_of_local(local)
  
    call MPI_ISEND(eta(:,:,:,local),num,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,ISEND(s_send),IERR)
  enddo
  
  do s_recv=1,num_recv_sites
    local=recv_sites(s_recv)%label_
    rank=recv_sites(s_recv)%rank_
    tag=10000*MYRANK + global_site_of_local(local)
  
    call MPI_IRECV(eta(:,:,:,local),num,MPI_DOUBLE_COMPLEX,rank,tag,MPI_COMM_WORLD,IRECV(s_recv),IERR)
  enddo
  
  do s_send=1,num_send_sites
    call MPI_WAIT(ISEND(s_send),ISTATUS,IERR)
  enddo
  do s_recv=1,num_recv_sites
    call MPI_WAIT(IRECV(s_recv),ISTATUS,IERR)
  enddo
  
  end subroutine syncronize_phibar


