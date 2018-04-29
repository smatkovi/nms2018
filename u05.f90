program wannier
  implicit none
  integer :: i
  real :: Eo=0.01d0
  real, dimension(1001) :: v1,v2,V,d
  real, dimension(1000) :: e
  real, dimension(2002) :: w
  real, dimension(1001,1001) :: H
  do i=0, 1000
    V(i)=0d0
    if(mod(i, 100) .le. 50) V(i)=Eo
  enddo
  H(1,1)=-1.d0*V(1)
  H(2,1)=0.5

  do i=2, 1000 
    H(i,i-1) = 0.5
    H(i,i) = -1.d0*V(i)
    H(i,i+1) = 0.5
  enddo

  do i=1, 1000
    v2(i)=1.d0
    v2(i)=cos( ((i-500))*1.d0/1000.d0*4.d0*atan(1.d0) )
  enddo
  v2=v2/norm2(v2)

  

  do i=1, 1000
    v1(i)=0.d0
  enddo
  
  write(*,*) sum(abs(v1-v2))
  i=0
  do while(sum(abs(v1-v2)) .gt. 0.027d0)
    v1=v2
    v2=matmul(H,v1)
    v2=v2/norm2(v2)
    i=i+1
    if(mod(i,1000) .eq. 0) write(*,*) i, sum(abs(v1-v2))
  enddo
  write(*,*) i

  do i=1, 1000
    d(i)=H(i,i)
  enddo

  
  do i=1, 1000
    e(i)=H(i+1,i)
    w(i)=H(i+1,i)
  enddo

  call ssteqr('T', 1001, d, e, H, 1001, w, i)
  write(*,*) i
  
  open(13, file='u05.txt')

  do i=1, 1000
    write (13,*) i, v1(i)
  enddo

  close(13)
end program
