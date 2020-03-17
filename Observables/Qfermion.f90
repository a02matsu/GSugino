!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutine to make Q\Psi
subroutine make_Qfermion(Qeta,Qlambda,Qchi,Omega,Umat,PhiMat)
use global_subroutines
use matrix_functions, only : matrix_commutator, matrix_3_product,matrix_power
use parallel
implicit none

complex(kind(0d0)), intent(out) :: Qeta(1:NMAT,1:NMAT,1:num_necessary_sites)
complex(kind(0d0)), intent(out) :: Qlambda(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(out) :: Qchi(1:NMAT,1:NMAT,1:num_necessary_faces)
complex(kind(0d0)), intent(out) :: Omega(1:NMAT,1:NMAT,1:num_faces)
complex(kind(0d0)), intent(in) :: Umat(1:NMAT,1:NMAT,1:num_necessary_links)
complex(kind(0d0)), intent(in) :: PhiMat(1:NMAT,1:NMAT,1:num_necessary_sites)
!complex(kind(0d0)), intent(in) :: PF_lambda(1:NMAT,1:NMAT,1:num_necessary_links)
!complex(kind(0d0)) :: comm(1:NMAT,1:NMAT)

integer :: s,l,f,i,j
complex(kind(0d0)) :: Uf(1:NMAT,1:NMAT)
complex(kind(0d0)) :: Ufm(1:NMAT,1:NMAT)

Qeta=(0d0,0d0)
do s=1,num_sites
  call matrix_commutator(Qeta(:,:,s),PhiMat(:,:,s),PhiMat(:,:,s),'N','C')
  Qeta(:,:,s)=Qeta(:,:,s)*U1Rfactor_site(s)
enddo

Qlambda=(0d0,0d0)
do l=1,num_links
  Qlambda(:,:,l)=(0d0,-1d0)*PhiMat(:,:,link_org(l))
  call matrix_3_product(Qlambda(:,:,l),&
    Umat(:,:,l),PhiMat(:,:,link_tip(l)),Umat(:,:,l),&
    'N','N','C',(0d0,1d0)*U1Rfactor_link(l)**2d0*U1R_ratio(l)**2d0,'ADD')
  Qlambda(:,:,l)=Qlambda(:,:,l)*U1Rfactor_site( link_org(l) )
enddo

Qchi=(0d0,0d0)
do f=1,num_faces
  call Make_face_variable(Uf,f,UMAT)
  if(m_omega == 0) then 
    call Make_moment_map0(Omega(:,:,f),Uf)
  elseif(m_omega == -1) then
    call Make_moment_map_adm(Omega(:,:,f),Uf)
  else
    call matrix_power(Ufm,Uf,m_omega)
    call Make_moment_map(Omega(:,:,f),Ufm)
  endif
  Qchi(:,:,f)=(0d0,-0.5d0)*dcmplx(beta_f(f))*Omega(:,:,f)
  Qchi(:,:,f)=Qchi(:,:,f)*U1Rfactor_site(sites_in_f(f)%label_(1))
enddo

#ifdef PARALLEL
call syncronize_sites(Qeta)
call syncronize_links(Qlambda)
call syncronize_faces(Qchi)
#endif

end subroutine make_Qfermion


