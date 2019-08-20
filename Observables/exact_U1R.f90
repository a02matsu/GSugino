!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! calculate D^\mu J_\mu for U(1)_R current on the lattice 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate D^\mu J_\mu for U(1)_V current
subroutine calc_exact_U1R(divJ,Geta_lambda,Glambda_eta,Glambda_lambda,Gchi_lambda,Glambda_chi,PhiMat,UMAT)
use global_parameters
!use initialization_calcobs
use parallel
use matrix_functions, only : matrix_3_product, matrix_product
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
complex(kind(0d0)) :: Xmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Ymat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Bval
complex(kind(0d0)) :: trace
complex(kind(0d0)) :: LamLam(1:NMAT,1:NMAT,1:num_necessary_links)
integer :: ls, ll, lf, src, tip, l_place
integer :: gl, gf
integer :: i,j,k,l,kk,jj

divJ=(0d0,0d0)
do ls=1,num_sites
  !!! bosonic part
  do k=1,linkorg_to_s(ls)%num_
    ll=linkorg_to_s(ls)%labels_(k)
    call PUPU(trace,ll,PhiMat,Umat,1)
    divJ(ls)=divJ(ls) + (0d0,1d0)*dcmplx(alpha_l(ll))*trace
  enddo
  do k=1,linktip_from_s(ls)%num_
    ll=linktip_from_s(ls)%labels_(k)
    call PUPU(trace,ll,PhiMat,Umat,2)
    divJ(ls)=divJ(ls) - (0d0,1d0)*dcmplx(alpha_l(ll))*trace
  enddo
  !write(*,*) MYRANK, "test1"
  !!! eta-lambda part 1
  do kk=1,linkorg_to_s(ls)%num_
    ll=linkorg_to_s(ls)%labels_(kk)
    trace=(0d0,0d0)
    do l=1,NMAT
      do k=1,NMAT
        do j=1,NMAT
          do i=1,NMAT
            trace=trace + UMAT(j,k,ll)*dconjg(Umat(i,l,ll)) &
                   * Glambda_eta(i,j,k,l,global_link_of_local(ll),ls)
          enddo
        enddo
      enddo
    enddo
   divJ(ls) = divJ(ls) + (-0.5d0,0d0)*dcmplx( alpha_l(ll) ) * trace
  enddo
  !write(*,*) MYRANK, "test2"
  !!! eta-lambda part 2
  do kk=1,linktip_from_s(ls)%num_
    ll=linktip_from_s(ls)%labels_(kk)
    tip=global_site_of_local(link_tip(ll))
    trace=(0d0,0d0)
     do l=1,NMAT
       do k=1,NMAT
         do j=1,NMAT
           do i=1,NMAT
             trace = trace + dconjg(UMAT(k,j,ll))*Umat(l,i,ll) &
                    * Geta_lambda(i,j,k,l,tip,ll)
           enddo
         enddo
       enddo
     enddo
     divJ(ls) = divJ(ls) + (0.5d0,0d0)*dcmplx( alpha_l(ll) ) * trace
   enddo
  !write(*,*) MYRANK, "test3"
  !!! lambda-lambda part 1
  call lambdalambda(LamLam,Glambda_lambda)
  do k=1,linkorg_to_s(ls)%num_
    ll=linkorg_to_s(ls)%labels_(k)
    call matrix_3_product(tmpmat,Umat(:,:,ll),PhiMat(:,:,ls),Umat(:,:,ll),&
      'N','C','C')
    trace=(0d0,0d0)
    do j=1,NMAT
      do i=1,NMAT
        trace = trace + LamLam(i,j,ll)*tmpmat(j,i)
      enddo
    enddo
  divJ(ls) = divJ(ls) + (0d0,1d0)*dcmplx( alpha_l(ll) ) * trace
  enddo
  !write(*,*) MYRANK, "test4"
  !!! lambda-lambda part 2
  do k=1,linktip_from_s(ls)%num_
    ll=linktip_from_s(ls)%labels_(k)
    tip=link_tip(ll)
    !write(*,*) tip
    call matrix_3_product(tmpmat,Umat(:,:,ll),PhiMat(:,:,tip),Umat(:,:,ll),&
      'N','C','C')
    trace=(0d0,0d0)
    do j=1,NMAT
      do i=1,NMAT
        trace = trace + LamLam(i,j,ll)*tmpmat(j,i)
      enddo
    enddo
  divJ(ls) = divJ(ls) - (0d0,1d0)*dcmplx( alpha_l(ll) ) * trace
  enddo
  !write(*,*) MYRANK, "test5"
  !!! chi-lambda part 1
  do jj=1,linktip_from_s(ls)%num_
    ll=linktip_from_s(ls)%labels_(jj)
    do kk=1,face_in_l(ll)%num_
      lf=face_in_l(ll)%label_(kk)
      gf=global_face_of_local(lf)
      ! find the place of ll in lf
      do l_place=1,links_in_f(lf)%num_
        if( links_in_f(lf)%link_labels_(l_place) == ll ) then 
          exit
        endif
      enddo
      call calc_XYB(Xmat,Ymat,Bval,lf,l_place,Umat)
      trace=(0d0,0d0)
      do l=1,NMAT
        do k=1,NMAT
          do j=1,NMAT
            do i=1,NMAT
              trace=trace+Gchi_lambda(i,j,k,l,gf,ll)&
                *Xmat(j,k)*Ymat(l,i)
            enddo
          enddo
        enddo
      enddo
      divJ(ls) = divJ(ls) &
        + (0d0,0.5d0)*dcmplx( alpha_f(lf)*beta_f(lf) ) &
          *(0d0,1d0)*dcmplx(dble(links_in_f(lf)%link_dirs_(l_place))) / Bval 
    enddo
  enddo
  !write(*,*) MYRANK, "test6"
  !!! chi-lambda part 2
  do lf=1,num_faces
    if( sites_in_f(lf)%label_(1) == ls ) then
      do l_place=1,links_in_f(lf)%num_
        ll=links_in_f(lf)%link_labels_(l_place)
        gl=global_link_of_local(ll)

        call calc_XYB(Xmat,Ymat,Bval,lf,l_place,Umat)
        trace=(0d0,0d0)
        do l=1,NMAT
          do k=1,NMAT
            do j=1,NMAT
              do i=1,NMAT
                trace=trace - Glambda_chi(i,j,k,l,gl,lf)&
                  *Ymat(j,k)*Xmat(l,i)
              enddo
            enddo
          enddo
        enddo
        divJ(ls) = divJ(ls) &
          - (0d0,0.5d0)*dcmplx( alpha_f(lf)*beta_f(lf) ) &
            *(0d0,1d0)*dcmplx(dble(links_in_f(lf)%link_dirs_(l_place))) / Bval 
      enddo
    endif
  enddo
  !write(*,*) MYRANK, "test7"
enddo

end subroutine calc_exact_U1R

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! C=1: trace = Tr( Phi U_ll \bar{Phi} U_ll\dag )
!! C=2: trace = Tr( \bar{Phi} U_ll Phi U_ll\dag )
subroutine PUPU(trace,ll,PhiMat,Umat,C)
use global_parameters
use matrix_functions, only : matrix_3_product, matrix_product
implicit none
complex(kind(0d0)), intent(out) :: trace
complex(kind(0d0)), intent(in) :: PhiMAT(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: UMAT(1:NMAT,1:NMAT,1:num_necessary_links)
integer, intent(in) :: ll
integer, intent(in) :: C

complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
integer :: src,tip
integer :: i,j

trace=(0d0,0d0)
src=link_org(ll)
tip=link_tip(ll)
if( C==1 ) then 
  call matrix_3_product(tmpmat,Umat(:,:,ll),PhiMat(:,:,tip),Umat(:,:,ll),&
    'N','C','C')
  do j=1,NMAT
    do i=1,NMAT
      trace=trace+PhiMat(i,j,src)*tmpmat(j,i)
    enddo
  enddo
elseif( C==2 ) then
  call matrix_3_product(tmpmat,Umat(:,:,ll),PhiMat(:,:,tip),Umat(:,:,ll),&
    'N','N','C')
  do j=1,NMAT
    do i=1,NMAT
      trace=trace+dconjg(PhiMat(j,i,src))*tmpmat(j,i)
    enddo
  enddo
endif
end subroutine PUPU


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compute lambda_{l,i,j} lambda_{l,j,k}
!! 
subroutine lambdalambda(LamLam,Glambda_lambda)
use global_parameters
use global_subroutines, only : syncronize_links
implicit none

complex(kind(0d0)), intent(out) :: LamLam(1:NMAT,1:NMAT,1:num_necessary_links)
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

call syncronize_links(LamLam)

end subroutine lambdalambda

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
