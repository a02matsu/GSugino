!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! write Dirac matrix
subroutine writeout_Dirac(Dirac)
use global_parameters
implicit none

complex(kind(0d0)), intent(in) :: Dirac(:,:)

integer :: i,s,l,f,a

i=0
do s=1,global_num_sites
  do a=1,dimG
    i=i+1
    call write_line(Dirac,i,'(S)',s,a)
  enddo
enddo
!!
do l=1,global_num_links
  do a=1,dimG
    i=i+1
    call write_line(Dirac,i,'(L)',l,a)
  enddo
enddo
!!
do f=1,global_num_faces
  do a=1,dimG
    i=i+1
    call write_line(Dirac,i,'(F)',f,a)
  enddo
enddo


contains
  subroutine write_line(Dirac,i,C,label,a)
  implicit none
  
  complex(kind(0d0)), intent(in) :: Dirac(:,:)
  integer, intent(in) :: i,label,a
  character(3), intent(in) :: C
  
  integer s,l,f,j,b
  
  j=0
  do s=1,global_num_sites
    do b=1,dimG
      j=j+1
      if( cdabs(Dirac(i,j)) > 1d-8) then
        write(*,'(A,2X,i5,2X,i3,2X,A,2X,i5,2X,i3,2X,E15.8,2X,E15.8,2X)') &
          C,label,a,"; (S)",s,b, dble(Dirac(i,j)),dble((0d0,-1d0)*Dirac(i,j))
      endif
    enddo
  enddo
  !!
  do l=1,global_num_links
    do b=1,dimG
      j=j+1
      if( cdabs(Dirac(i,j)) > 1d-8) then
        write(*,'(A,2X,i5,2X,i3,2X,A,2X,i5,2X,i3,2X,E15.8,2X,E15.8,2X) ') &
          C,label,a,"; (L)",l,b, dble(Dirac(i,j)),dble((0d0,-1d0)*Dirac(i,j))
      endif
    enddo
  enddo
  !!
  do f=1,global_num_faces
    do b=1,dimG
      j=j+1
      if( cdabs(Dirac(i,j)) > 1d-8) then
        write(*,'(A,2X,i5,2X,i3,2X,A,2X,i5,2X,i3,2X,E15.8,2X,E15.8,2X) ') &
          C,label,a,"; (F)",f,b, dble(Dirac(i,j)),dble((0d0,-1d0)*Dirac(i,j))
      endif
    enddo
  enddo
  end subroutine write_line

end subroutine writeout_Dirac
