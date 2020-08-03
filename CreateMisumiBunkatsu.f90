program main
implicit none

integer :: sizeM
integer :: sizeN
double precision :: radius=1d0
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

integer :: iarg
character(50) :: CM,CN

iarg=iargc()
call getarg(1,CM)
call getarg(2,CN)
read(CM,*) sizeM
read(CN,*) sizeN




if( mod(sizeN,2) == 1 ) then 
  write(*,*) "Set even N"
  stop
endif

num_sites=sizeM*sizeN/2
num_links=sizeM*(sizeN-1)
num_faces=sizeM*(sizeN-2)/2+2

a=dsqrt( 8d0*PI*radius**2 / dble( sizeM*(sizeN-2)+4 ) )

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
tmp=PI*radius**2/(dble(sizeM)*a**2)&
  *( 2d0 - dcos(Pi/dble(sizeN)) - dcos(3d0*Pi/dble(sizeN)) )
alpha_s(1:sizeM)=tmp
alpha_s( sizeM*(sizeN/2-1)+1:sizeM*sizeN/2 )=tmp
!!
do n=1,sizeN/2-2
  alpha_s(n*sizeM+1:(n+1)*sizeM)=PI*radius**2/( dble(sizeM)*a**2 )&
    *( dcos( dble(2*n-1)*Pi/dble(sizeN) ) - dcos( dble(2*n+3)*Pi/dble(sizeN) ) )
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! alpha_f, beta_f
do n=1,sizeN/2-1 
  tmp=2d0*Pi*radius**2/(dble(sizeM)*a**2) &
    *( dcos( dble(2*n-1)*Pi/dble(sizeN) ) - dcos( dble(2*n+1)*Pi/dble(sizeN) ) )
  alpha_f( (n-1)*sizeM+2 : n*sizeM +1) = tmp
  beta_f( (n-1)*sizeM+2 : n*sizeM +1) = 1d0/tmp
enddo
tmp=2d0*Pi*radius**2/a**2*(1d0 - dcos(Pi/dble(sizeN))) 
alpha_f( 1 ) = tmp
alpha_f( sizeM*(sizeN/2-1)+2 ) = tmp
beta_f( 1 ) = 1d0/tmp
beta_f( sizeM*(sizeN/2-1)+2 ) = 1d0/tmp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! alpha_l
if( mod(sizeN,4) == 0 ) then
  do n=1,sizeN/2
    tmp=0d0
    do k=1,n-1
      tmp=tmp + dble((-1)**(n+k+1))/dsin( dble(2*k)*Pi/dble(sizeN) )
    enddo
    do k=n,sizeN/2-1
      tmp=tmp - dble((-1)**(n+k+1))/dsin( dble(2*k)*Pi/dble(sizeN) )
    enddo
    tmp=tmp * dble(sizeM)*dsin(Pi/dble(sizeN))/Pi  

    alpha_l((n-1)*sizeM+1:n*sizeM)=tmp
  enddo
elseif( sizeN == 6 ) then
  do n=1,3
    alpha_l((n-1)*sizeM+1:n*sizeM)=dble(sizeM)/dsqrt(3d0)/Pi
  enddo
else
  tmp=0d0
  do k=1, (sizeN-2)/4-2
    tmp=tmp + 4d0*dble( (-1)**( (sizeN-2)/4+k ) )/dsin( dble(2*k)*Pi/dble(sizeN) ) 
  enddo
  tmp=tmp - 1d0/dsin( 2d0*Pi/dble(sizeN)*dble( (sizeN-2)/4 -1 ) ) &
          - 1d0/dsin( 2d0*Pi/dble(sizeN)*dble( (sizeN-2)/4 ) )
  tmp = dabs(tmp) * dble(sizeM)*dsin(Pi/dble(sizeN))/(2d0*Pi) 

  alpha_l(1:sizeM) = tmp
  do n=2, sizeN/2
    tmp = 2d0*sizeM*dsin(Pi/sizeN)/Pi * 1d0/dsin( dble(2*(n-1))*Pi/dble(sizeN) ) - tmp
    alpha_l( (n-1)*sizeM+1 : n*sizeM ) = tmp
  enddo
endif

do n=1, sizeN/2-1
  alpha_l( sizeM*sizeN/2 + (n-1)*sizeM + 1 : sizeM*sizeN/2 + n*sizeM ) &
    = dble( sizeN**2 )/(dble(2*sizeM)*Pi) &
      *( dcos( dble(2*n-1)*Pi/dble(sizeN) ) - dcos( dble(2*n+1)*Pi/dble(sizeN) ) )
enddo



!open(unit=10, status='replace',file=filename,action='write')
write(OUT,'(a,I2,a,I2,X,a,F8.4)') "# Misumi-bunkatsu, (M,N)=(",sizeM,",",sizeN,"), R=",radius
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
do n=1,sizeN/2
  do i=1,sizeM-1
    l=(n-1)*sizeM+i
    write(OUT,'(I6,X,I6,X,I6,X,E15.8,X,A)') l, (n-1)*sizeM+i, (n-1)*sizeM+i+1, alpha_l(l), "0d0"
  enddo
  l=n*sizeM
  write(OUT,'(I6,X,I6,X,I6,X,E15.8,X,A)' ) l, n*sizeM, (n-1)*sizeM+1, alpha_l(l), "0d0"
enddo
do n=1,sizeN/2-1
  do i=1,sizeM
    l= sizeM*sizeN/2 + (n-1)*sizeM + i
    write(OUT,'(I6,X,I6,X,I6,X,E15.8,X,A)') l , (n-1)*sizeM+i, n*sizeM+i, alpha_l(l), "0d0"
  enddo
enddo
!!!!!!!!!!!!
write(OUT,'(A)') "# alpha_f: f, size(f), s_1, ..., s_Sf, alpha(f), beta(f)"
f=1
write(OUT,'(I6,X,I6,X)',advance='no') f, sizeM
do i=1,sizeM
  write(OUT,'(I6,X)',advance='no') i
enddo
write(OUT,'(E15.8,X,E15.8)') alpha_f(f), beta_f(f)
!!
do n=1,sizeN/2-1
  do i=1,sizeM-1
    f=(n-1)*sizeM+i+1
    write(OUT,'(I6,X,I6,X,I6,X,I6,X,I6,X,I6,X,E15.8,X,E15.8)') f, 4, &
      (n-1)*sizeM+i, (n-1)*sizeM+i+1, n*sizeM+i+1, n*sizeM+i, &
      alpha_f(f),beta_f(f)
  enddo
    f=n*sizeM+1
    write(OUT,'(I6,X,I6,X,I6,X,I6,X,I6,X,I6,X,E15.8,X,E15.8)') f, 4, &
      n*sizeM, (n-1)*sizeM+1, n*sizeM+1, (n+1)*sizeM, &
      alpha_f(f),beta_f(f)
enddo
!!
f=(sizeN/2-1)*sizeM+2
write(OUT,'(I6,X,I6,X)',advance='no') f, sizeM
do i=1,sizeM
  write(OUT,'(I6,X)',advance='no') sizeM*(sizeN/2-1)+i
enddo
write(OUT,'(E15.8,X,E15.8)') alpha_f(f),beta_f(f)

!close(OUT)

end program main
