!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! calculate D^\mu J_\mu for U(1)_R current on the lattice 
!module global_exact_U1R
!implicit none
!
!type A_in_B
!  integer :: num_ ! number of elements
!  integer, allocatable :: label_(:) ! global link label
!  complex(kind(0d0)), allocatable :: val_(:) ! 1~num_
!end type A_in_B
!
!end module global_exact_U1R


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate D^\mu J_\mu for U(1)_V current
subroutine calc_exact_U1R(divJ,Geta_lambda,Glambda_eta,Glambda_lambda,Gchi_lambda,Glambda_chi,PhiMat,UMAT)
use global_parameters
!use global_exact_U1R
!use initialization_calcobs
use parallel
!use matrix_functions, only : matrix_3_product, matrix_product
implicit none


complex(kind(0d0)), intent(out) :: divJ(1:num_sites)
complex(kind(0d0)), intent(in) :: Geta_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_links) 
complex(kind(0d0)), intent(in) :: Glambda_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_sites) 
complex(kind(0d0)), intent(in) :: Glambda_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_links) 
complex(kind(0d0)), intent(in) :: Gchi_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_faces,1:num_links) 
complex(kind(0d0)), intent(in) :: Glambda_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_faces) 
complex(kind(0d0)), intent(in) :: PhiMAT(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)


complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Xmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Ymat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Bval
complex(kind(0d0)) :: trace
complex(kind(0d0)) :: PUPU(1:num_necessary_links)
complex(kind(0d0)) :: LUEU(1:num_necessary_links)
complex(kind(0d0)) :: LamLam(1:NMAT,1:NMAT,1:num_links)

!type A_in_B
!  integer :: num_ ! number of elements
!  integer, allocatable :: label_(:) ! global link label
!  complex(kind(0d0)), allocatable :: val_(:) ! 1~num_
!end type A_in_B
type(A_in_B) :: CXLY_Fbase(1:num_faces)
type(A_in_B) :: CXLY_Lbase(1:num_links)


integer :: ls, ll, lf, src, tip, l_place
integer :: gl, gf
integer :: i,j,k,l,kk,jj

!!! preparation
divJ=(0d0,0d0)
! for bosonic part
call calc_Jboson(PUPU,PhiMat,Umat)
! for lambda-lambda part
call lambdalambda(LamLam,Glambda_lambda)
! for fermionic site-link and link-link part
call calc_Jfermion1(LUEU,LamLam,PhiMat,Umat,Geta_lambda)
! for fermionic face-link part
call calc_Jfermion2(CXLY_Fbase,CXLY_Lbase,Umat,Glambda_chi)


do ls=1,num_sites
  !!! bosonic part
  do k=1,linkorg_to_s(ls)%num_
    ll=linkorg_to_s(ls)%labels_(k)
    divJ(ls)=divJ(ls) + PUPU(ll) 
  enddo
  do k=1,linktip_from_s(ls)%num_
    ll=linktip_from_s(ls)%labels_(k)
    divJ(ls)=divJ(ls) - PUPU(ll)
  enddo

  !!! eta-lambda and lambda-lambda part
  do k=1,linkorg_to_s(ls)%num_
    ll=linkorg_to_s(ls)%labels_(k)
    divJ(ls) = divJ(ls) + LUEU(ll)
  enddo
  do k=1,linktip_from_s(ls)%num_
    ll=linktip_from_s(ls)%labels_(k)
    divJ(ls) = divJ(ls) - LUEU(ll)
  enddo

  !!! chi-lambda part 
  do kk=1,linktip_from_s(ls)%num_ 
    ll=linktip_from_s(ls)%labels_(kk)
    do jj=1,face_in_l(ll)%num_
      lf=face_in_l(ll)%label_(jj)

      do i=1,CXLY_Lbase(ll)%num_
        if( CXLY_Lbase(ll)%label_(i) == lf ) then
          divJ(ls) = divJ(ls) + CXLY_Lbase(ll)%val_(i)
          exit
        endif
      enddo
    enddo
  enddo
  !!
  do lf=1,num_faces
    if( sites_in_f(lf)%label_(1) == ls ) then
      do i=1, CXLY_Fbase(lf)%num_ 
        ll=CXLY_Fbase(lf)%label_(i)
        divJ(ls) = divJ(ls) - CXLY_Fbase(ll)%val_(i)
      enddo
    endif
  enddo
enddo

end subroutine calc_exact_U1R

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! (i*alpha_l)*Tr( Phi U_ll \bar{Phi} U_ll\dag - \bar{Phi} U_ll Phi U_ll\dag )
subroutine calc_Jboson(PUPU,PhiMat,Umat)
use global_parameters
use matrix_functions, only : matrix_3_product, matrix_product
use global_subroutines, only : syncronize_linkval
implicit none
complex(kind(0d0)), intent(out) :: PUPU(1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMAT(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)

complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
integer :: src,tip
integer :: i,j,ll

PUPU=(0d0,0d0)
do ll=1,num_links
  src=link_org(ll)
  tip=link_tip(ll)
  call matrix_3_product(tmpmat,Umat(:,:,ll),PhiMat(:,:,tip),Umat(:,:,ll),&
    'N','C','C')
  do j=1,NMAT
    do i=1,NMAT
      PUPU(ll)=PUPU(ll)+PhiMat(i,j,src)*tmpmat(j,i)
    enddo
  enddo
  call matrix_3_product(tmpmat,Umat(:,:,ll),PhiMat(:,:,tip),Umat(:,:,ll),&
    'N','N','C')
  do j=1,NMAT
    do i=1,NMAT
      PUPU(ll)=PUPU(ll)+dconjg(PhiMat(j,i,src))*tmpmat(j,i)
    enddo
  enddo
  PUPU(ll) = PUPU(ll) * (0d0,1d0)*dcmplx(alpha_l(ll))
enddo

call syncronize_linkval(PUPU)

end subroutine calc_Jboson


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compute lambda_{l,i,j} lambda_{l,j,k}
!! 
subroutine lambdalambda(LamLam,Glambda_lambda)
use global_parameters
!use initialization_calcobs
use global_subroutines, only : syncronize_links
implicit none

complex(kind(0d0)), intent(out) :: LamLam(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: Glambda_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_links) 
integer :: i,j,k,ll,gl

LamLam=(0d0,0d0)
do ll=1,num_links
  gl=global_link_of_local(ll)
  do i=1,NMAT
    do j=1,NMAT
      do k=1,NMAT
        LamLam(i,j,ll)=LamLam(i,j,ll)+Glambda_lambda(i,k,k,j,gl,ll)
      enddo
    enddo
  enddo
enddo

end subroutine lambdalambda

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! (i*alpha_l) * Tr( i/2 \lambda U \eta Uinv +\lambda\lambda U \bar{\Phi} Uinv
subroutine calc_Jfermion1(LUEU,LamLam,PhiMat,Umat,Geta_lambda)
use global_parameters
use matrix_functions, only : matrix_3_product
implicit none

complex(kind(0d0)), intent(out) :: LUEU(1:num_necessary_links)
complex(kind(0d0)), intent(in) :: LamLam(1:NMAT,1:NMAT,1:num_links)
complex(kind(0d0)), intent(in) :: PhiMAT(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: Geta_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_links) 

complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
integer :: src,tip
integer :: ll
integer :: i,j,k,l

LUEU=(0d0,0d0)
do ll=1,num_links
  src=global_link_of_local(link_org(ll))
  do l=1,NMAT
    do k=1,NMAT
      do j=1,NMAT
        do i=1,NMAT
          LUEU(ll)=LUEU(ll) - dconjg(UMAT(k,j,ll))*(Umat(l,i,ll)) &
                 * Geta_lambda(i,j,k,l,src,ll)
        enddo
      enddo
    enddo
  enddo
  LUEU(ll) = LUEU(ll) * (0d0,0.5d0)

  tip=link_tip(ll)
  call matrix_3_product(tmpmat,Umat(:,:,ll),PhiMat(:,:,tip),Umat(:,:,ll),&
    'N','C','C')
  do j=1,NMAT
    do i=1,NMAT
      LUEU(ll) = LUEU(ll) + LamLam(i,j,ll)*tmpmat(j,i)
    enddo
  enddo
  LUEU(ll) = LUEU(ll) * (0d0,1d0)*dcmplx( alpha_l(ll) )
enddo

call syncronize_linkval(LUEU)


end subroutine calc_Jfermion1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Xmat = U_1 ... U_l
!! Ymat = U_{l+1} ... U_n
!! Bval = 1 - (2-Uf-Uf^\dagger)/e_max^2
subroutine calc_XYB(Xmat,Ymat,Bval,lf, l_place,Umat)
use global_parameters
use global_subroutines, only : calc_XYmat, make_face_variable
implicit none

complex(kind(0d0)), intent(out) :: Xmat(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(out) :: Ymat(1:NMAT,1:NMAT)
complex(kind(0d0)), intent(out) :: Bval
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
integer, intent(in) :: lf,l_place
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT)
integer :: i

call calc_XYmat(Xmat,Ymat,lf,l_place,UMAT)

!!!!!!!!!!!!!!!
!! Bval = 1 - (2-Uf-Uf^\dagger)/e_max^2
Bval=(1d0,0d0)
call make_face_variable(Uf(:,:),lf,Umat)
do i=1,NMAT
  Bval=Bval - ((2d0,0d0)-Uf(i,i)-dconjg(Uf(i,i)))/(e_max*e_max) 
enddo

end subroutine calc_XYB


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! alpha_f*beta_f*epsilon_{f,l} 
!!  * ( -1/2B Tr(\chi X \lambda Y + \chi Ydag \lambda Xdag)
!!      + 1/2B^2\e_max^2 Tr( chi(Uf-Ufinv) ) * Tr( X \lambda Y - Ydag \lambda Xdag )
subroutine calc_Jfermion2(CXLY_Fbase,CXLY_Lbase,Umat,Glambda_chi)
use global_parameters
!use global_exact_U1R
!use initialization_calcobs
use global_subroutines, only : make_face_variable
use matrix_functions, only : matrix_product
implicit none


!type A_in_B
!  integer :: num_ ! number of elements
!  integer, allocatable :: label_(:) ! global link label
!  complex(kind(0d0)), allocatable :: val_(:) ! 1~num_
!end type A_in_B
type(A_in_B), intent(out) :: CXLY_Fbase(1:num_faces)
type(A_in_B), intent(out) :: CXLY_Lbase(1:num_links)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: Glambda_chi(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_faces) 

complex(kind(0d0)) :: Xmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Ymat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Usin(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Tmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Bval
complex(kind(0d0)) :: trace,tmp
integer :: lf,ll,gl,gf
integer :: l_place
integer :: i,j,k,l
integer :: rank_recv, rank_send,tag,info

do lf=1,num_faces
  CXLY_Fbase(lf)%num_= links_in_f(lf)%num_
  allocate( CXLY_Fbase(lf)%label_(1:links_in_f(lf)%num_) )
  allocate( CXLY_Fbase(lf)%val_(1:links_in_f(lf)%num_) )
  do i=1,CXLY_Fbase(lf)%num_
    CXLY_Fbase(lf)%label_(i) = links_in_f(lf)%link_labels_(i)
  enddo
  CXLY_Fbase(lf)%val_ = (0d0,0d0)
enddo

do lf=1,num_faces
  do l_place=1,links_in_f(lf)%num_
    gl = global_link_of_local( links_in_f(lf)%link_labels_(l_place) )
    call calc_XYB(Xmat,Ymat,Bval,lf,l_place,Umat)

    trace=(0d0,0d0)
    do l=1,NMAT
      do k=1,NMAT
        do j=1,NMAT
          do i=1,NMAT
            trace=trace + (0.5d0,0d0)/Bval &
              * Glambda_chi(i,j,k,l,gl,lf) &
              * (Ymat(j,k)*Xmat(l,i) + dconjg(Xmat(k,j)*Ymat(i,l)))
          enddo
        enddo
      enddo
    enddo
    !! Usin
    call make_face_variable(Uf,lf,Umat)
    do i=1,NMAT
      do j=1,NMAT
        Usin(i,j) = Uf(i,j) - dconjg(Uf(j,i))
      enddo
    enddo
    !! Tmat
    Tmat=(0d0,0d0)
    call matrix_product(Tmat,Ymat,Xmat)
    call matrix_product(Tmat,Xmat,Ymat,'C','C',(-1d0,0d0),'ADD')
    !do i=1,NMAT
    !  do j=1,NMAT
    !    do k=1,NMAT
    !      Tmat(i,j) = Tmat(i,j) &
    !        + Ymat(i,k)*Xmat(k,j) - dconjg(Xmat(k,i)*Ymat(j,k))
    !    enddo
    !  enddo
    !enddo
    !!
    do i=1,NMAT
      do j=1,NMAT
        tmp=(0d0,0d0)
        do k=1,NMAT
          do l=1,NMAT
            tmp=tmp+Glambda_chi(i,j,k,l,gl,lf)*Usin(l,k)
          enddo
        enddo
        trace = trace - (0.5d0,0d0) / (e_max*e_max*Bval*Bval) &
          * tmp * Tmat(j,i)
      enddo
    enddo
    CXLY_Fbase(lf)%val_ &
      = dcmplx( alpha_f(lf)*beta_f(lf)*dble(links_in_f(lf)%link_dirs_(l_place)) ) &
          * trace
  enddo
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! fransfer the information to CXLY_L
do ll=1,num_links
  CXLY_Lbase(ll)%num_= face_in_l(ll)%num_
  allocate( CXLY_Lbase(ll)%label_(1:face_in_l(ll)%num_) )
  allocate( CXLY_Lbase(ll)%val_(1:face_in_l(ll)%num_) )
  do i=1,CXLY_Lbase(ll)%num_
    CXLY_Lbase(ll)%label_(i) = face_in_l(ll)%label_(i)
  enddo
  CXLY_Lbase(lf)%val_ = (0d0,0d0)
enddo

do gf=1,global_num_faces
  rank_send=local_face_of_global(gf)%rank_
  lf=local_face_of_global(gf)%label_

  do k=1,global_links_in_f(gf)%num_
    gl=global_links_in_f(gf)%link_labels_(k)
    rank_recv=local_link_of_global(gl)%rank_
    ll=local_link_of_global(gl)%label_

    tag=gf*(global_num_faces-1)+gl-1

    if( MYRANK == rank_send ) then
      info=1
      do i=1,CXLY_Fbase(lf)%num_
        if( links_in_f(lf)%link_labels_(i) == ll ) then 
          info=0
          exit
        endif
      enddo
      if( info==1 ) then 
        write(*,*) "there is no ll in lf" 
        stop
      endif
      if( MYRANK /= rank_recv ) then 
        call MPI_SEND(CXLY_Fbase(lf)%val_(i),1,MPI_DOUBLE_COMPLEX,rank_recv,tag,MPI_COMM_WORLD,IERR)
      endif
    elseif( MYRANK == rank_recv ) then
      info=1
      do j=1,CXLY_Lbase(ll)%num_
        if( face_in_l(ll)%label_(j) == lf ) then
          info=0
          exit
        endif
      enddo
      if( info==1 ) then 
        write(*,*) "there is no lf in ll" 
        stop
      endif
      if( MYRANK /= rank_send ) then
        call MPI_RECV(CXLY_Lbase(ll)%val_(j),1,MPI_DOUBLE_COMPLEX,rank_send,tag,MPI_COMM_WORLD,ISTATUS,IERR)
      else
        CXLY_Lbase(ll)%val_(j) = CXLY_Fbase(lf)%val_(i)
      endif
    endif
  enddo
enddo

end subroutine calc_Jfermion2


