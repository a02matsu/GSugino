subroutine calc_QS_3fermion_link(QS_eta,QS_lambda,&
    Glambda_eta, Geta_lambda, Glambda_lambda, &
    Phimat, Umat)
use matrix_functions, only : make_unit_matrix, hermitian_conjugate,matrix_3_product
implicit none

complex(kind(0d0)), intent(out) :: QS_eta(1:NMAT,1:NMAT,1:num_sites)
complex(kind(0d0)), intent(out) :: QS_lambda(1:NMAT,1:NMAT,1:num_links)

complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: Glambda_eta(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_sites)
complex(kind(0d0)), intent(in) :: Geta_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_sites,1:num_links)
complex(kind(0d0)), intent(in) :: Glambda_lambda(1:NMAT,1:NMAT,1:NMAT,1:NMAT,1:global_num_links,1:num_links)


complex(kind(0d0)) :: factor, qs, qt, ql
complex(kind(0d0)) :: tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) :: tmpmat2(1:NMAT,1:NMAT)
integer :: ll,ls,gl,gs
integer :: i,j,k,l,a,c,d

QS_eta=(0d0,0d0)
do ls=1,num_sites
  do a=1,linkorg_to_s(ls)%num_
    ll=linkorg_to_s(ls)%labels_(a)
    !! factor
    qs=U1Rfactor_site(link_org(ll))
    qt=U1Rfactor_site(link_tip(ll))
    ql=U1Rfactor_link(ll)
    factor= qs/ql*( (1d0,0d0) - qt/qs/ql )*overall_factor
    !!
    if( cdabs(factor) > 1d-8 ) then
      gl=global_link_of_local(ll)
      tmpmat=(0d0,0d0)
      do j=1,NMAT
        do i=1,NMAT
          do k=1,NMAT
            tmpmat(i,j)=tmpmat(i,j)+Glambda_lambda(i,k,k,j,gl,ll)
          enddo
        enddo
      enddo
      call matrix_3_product( tmpmat2, Umat(:,:,ll),tmpmat,Umat(:,:,ll),'C','N','N')
      QS_eta(:,:,ls)=QS_eta(:,:,ls) + factor*tmpmat2
    endif
  enddo
enddo

QS_lambda=(0d0,0d0)
do ll=1,num_links
  qs=U1Rfactor_site(link_org(ll))
  qt=U1Rfactor_site(link_tip(ll))
  ql=U1Rfactor_link(ll)
  factor= qs/ql*( (1d0,0d0) - qt/qs/ql )*overall_factor
  if( cdabs(factor) > 1d-8 ) then
    gl=global_link_of_local(ll)
    ls=link_tip(ll)
    gs=global_site_of_local(ls)
    !!
    tmpmat=(0d0,0d0)
    do k=1,NMAT
    do l=1,NMAT
      do c=1,NMAT
        do d=1,NMAT
          do a=1,NMAT
            tmpmat(l,k)=tmpmat(l,k)&
              +Glambda_eta(l,a,c,d,gl,ls)*Umat(a,c,ll)*dconjg(Umat(k,d,ll)) &
              +Geta_lambda(c,d,a,k,gs,ll)*Umat(l,c,ll)*dconjg(Umat(a,d,ll))
          enddo
        enddo
      enddo
    enddo
    enddo
    QS_lambda(:,:,ll)=QS_lambda(:,:,ll) + factor*tmpmat
  endif
enddo


end subroutine calc_QS_3fermion_link

