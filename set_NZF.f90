!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! set nonzero values of the structure constant
subroutine set_NZF
use SUN_generators
implicit none
integer :: a,b,c,i
integer :: num(1:dimG)
 
call get_num_nonzerof(NZF,NMAT)
allocate( NZF_index(1:3,1:NZF) )
allocate( NZF_value(1:NZF) )
call get_nonzerof(NZF_index,NZF_value,NZF,NMAT)

allocate( NZF_a(1:dimG) )

!! for f(a,b,e)*f(c,d,e)
call get_num_nonzero_ff(NZFF,NMAT)
allocate( NZFF_index(1:4,1:NZFF) )
allocate( NZFF_value(1:NZFF) )
call get_nonzero_ff( NZFF_index, NZFF_value, NZFF, NMAT)

!! count the number of nonzero components of f(a,b,c) with a fixed a
num=0
do i=1,NZF
  a=NZF_index(1,i)
  num(a)=num(a)+1
enddo

do a=1,dimG
!write(*,*) a,num(a)
NZF_a(a)%num_=num(a)
allocate( NZF_a(a)%b_(1:num(a)) )
allocate( NZF_a(a)%c_(1:num(a)) )
allocate( NZF_a(a)%value_(1:num(a)) )
enddo

num=0
do i=1,NZF
  a=NZF_index(1,i)
  b=NZF_index(2,i)
  c=NZF_index(3,i)
  num(a)=num(a)+1
  NZF_a(a)%b_(num(a))=b
  NZF_a(a)%c_(num(a))=c
  NZF_a(a)%value_(num(a))=NZF_value(i)
enddo
end subroutine set_NZF

