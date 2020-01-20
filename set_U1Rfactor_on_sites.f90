subroutine set_U1Rfactor_on_sites
use global_parameters
implicit none

!integer :: U1Rfactor_site(1:global_num_sites)
integer :: groupS(1:global_num_sites)
integer :: groupT(1:global_num_sites)
integer :: groupU(1:global_num_sites)

integer :: s,t,l
integer :: i,j,k
integer :: numS,numT,numU
integer :: info

!!!!!!!!!!!
groupS=0
numS=0
!!!!!!!!!!!
groupT=0
numT=1
T(1)=1
!!!!!!!!!!!
groupU=0
numU=0

do while ( numS + numT < global_num_sites ) 
  do i=1,numT
    s=groupT(i)
    do j=1,global_linktip_from_s(s)%num_
      t=global_linktip_from_s(s)%sites_(j)
      l=global_linktip_from_s(s)%labels_(j)

      info=1
      do k=1,numS
        if( t==groupS(k) ) then
          info=0
          exit
        endif
      enddo
      do k=1,numT
        if( t==groupT(k) ) then
          info=0
          exit
        endif
      enddo
      if( info==1 ) then

        U1Rfactor_site(t) = U1Rfactor_site(s)*U1Rfactor(


      endif





end subroutine set_U1Rfactor_on_sites
