program main
implicit none

integer :: sizeM
integer :: sizeN
double precision :: length1=1d0
double precision :: length2=1d0
!character(20) :: filename="S2MisumiM32N32R1.dat"
integer, parameter :: OUT=6

double precision, parameter :: PI=dacos(-1d0)
integer :: num_sites, num_links, num_faces
double precision :: a 
double precision, allocatable :: alpha_s(:)
double precision, allocatable :: alpha_l(:)
double precision, allocatable :: alpha_f(:)
double precision, allocatable :: beta_f(:)

double precision :: tmp
integer :: s,l,f
integer :: i,j,k,n
integer :: f_end
integer :: s1, s2, s3, s4

integer :: iarg
character(50) :: CM,CN
double precision :: U1R_ratio ! exp(2\pi i U1R_ratio) is the U(1) monodoromy of the cycle (common to the two cycles)
character(50) :: C_U1R_ratio
double precision :: U1R_link1 ! U1R_ratio/sizeM 
double precision :: U1R_link2 ! U1R_ratio/sizeN 

iarg=iargc()
call getarg(1,CM)
call getarg(2,CN)
read(CM,*) sizeM
read(CN,*) sizeN
if( iarg >= 3 ) then 
  call getarg(3,C_U1R_ratio)
  read(C_U1R_ratio,*) U1R_ratio
else
  U1R_ratio=0d0
endif

num_sites=sizeM*sizeN
num_links=2*sizeM*sizeN
num_faces=sizeM*sizeN

a=dsqrt( length1 * length2 / dble(sizeM*sizeN) )

U1R_link1 = U1R_ratio / dble(sizeM) 
U1R_link2 = U1R_ratio / dble(sizeN) 


!write(*,*) radius,num_sites,num_links,num_faces,PI

allocate( alpha_s(1:num_sites) )
allocate( alpha_l(1:num_links) )
allocate( alpha_f(1:num_faces) )
allocate( beta_f(1:num_faces) )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! initialization 
alpha_s=1d0
alpha_l=1d0
alpha_f=1d0
beta_f=1d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! alpha_s 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! alpha_f, beta_f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! alpha_l

!open(unit=10, status='replace',file=filename,action='write')
write(OUT,'(a,I2,a,I2,X,a,F8.4,a,F8.4)') "# Torus, (M,N)=(",sizeM,",",sizeN,"), L1=",length1, ", L2=",length2
!!!!!!!!!!!
write(OUT,'(I6,A)') num_sites, "  # number of sites"
write(OUT,'(I6,A)') num_links, "  # number of links"
write(OUT,'(I6,A)') num_faces, "  # number of faces"
write(OUT,'(E15.8,A)') a, "   #lattice spacing"
!!!!!!!!!!!!
write(OUT,'(A)') "# alpha_s: s, alpha_s(s)"
do s=1,num_sites
  write(OUT,'(I6,E15.8)') s, alpha_s(s)
enddo
!!!!!!!!!!!!
write(OUT,'(A)') "# alpha_l: l, org(l), tip(l), alpha_l(l)"
l=0
do n=0, sizeN-1
  do i=1,sizeM
    if( i /= sizeM ) then
      l=l+1
      write(OUT,'(I6,X,I6,X,I6,X,E15.8,X,E15.8)') l, n*sizeM+i, n*sizeM+i+1, alpha_l(l), U1R_link1
    else
      l=l+1
      write(OUT,'(I6,X,I6,X,I6,X,E15.8,X,E15.8)') l, n*sizeM+i, n*sizeM+1, alpha_l(l), U1R_link1
    endif
    !!
    if( n /= sizeN-1 ) then
      l=l+1
      write(OUT,'(I6,X,I6,X,I6,X,E15.8,X,E15.8)') l, n*sizeM+i, (n+1)*sizeM+i, alpha_l(l), U1R_link2
    else
      l=l+1
      write(OUT,'(I6,X,I6,X,I6,X,E15.8,X,E15.8)') l, n*sizeM+i, i, alpha_l(l), U1R_link2
    endif
  enddo
enddo
!!!!!!!!!!!!
write(OUT,'(A)') "# alpha_f: f, size(f), s_1, ..., s_Sf, alpha(f), beta(f)"
f=0
do n=0, sizeN-1
  do i=1, sizeM
    f=f+1
    !!
    s1=n*sizeN + i
    !!
    if( i/=sizeM ) then
      s2=n*sizeN + i + 1
    else
      s2=n*sizeN + 1
    endif
    !!
    if( n /= sizeN-1 ) then
      if( i==sizeM ) then 
        s3=(n+1)*sizeN + 1
      else
        s3=(n+1)*sizeN + i + 1
      endif
    else
      if( i==sizeN ) then
        s3=1
      else
        s3=i+1
      endif
    endif
    !!
    if( n /= sizeN-1 ) then
      s4=(n+1)*sizeN + i
    else
      s4 = i
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(OUT,'(I6,X,I6,X,I6,X,I6,X,I6,X,I6,X,E15.8,X,E15.8)') f, 4, &
      s1, s2, s3, s4, alpha_f(f),beta_f(f)
  enddo
enddo

!close(OUT)

end program main

