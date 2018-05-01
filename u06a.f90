program laplace
  implicit none
  integer :: i,j,k=1000
  real, dimension(151,151) :: u1, u2  
  
  do  i=1, 150
    do j=1, 150
      u2(i,j)=0
      if(i .ge. 50 .and. i .le. 100 .and. j .ge. 50 .and. j .le. 100) u2(i,j)=100
    enddo
  enddo
  write(*,13) abs(sum(u1-u2))
  13 format (e10.3)
  do while(abs(sum(u1-u2)) .ge. 1e-6)!(k .ge. 1)!abs(sum(u1-u2)) .ge. 1) 
    u1=u2
      
    do  i=1, 150
      do j=1, 150
        if(i .ge. 50 .and. i .le. 100 .and. j .ge. 50 .and. j .le. 100) then
          u2(i,j) = 1000.d0
        else if(i .eq. 1 .or. i .eq. 150 .or. j .eq. 1 .or. j .eq. 150) then
          u2(i,j) = 0.d0
        else
          u2(i, j) = (u2(i-1, j) + u2(i, j-1) + u1(i+1, j) + u1(i, j+1))/4.d0
        endif
      enddo
    enddo
    k=k-1
  enddo
  write(*,13) abs(sum(u1-u2))

  open(unit=13, file='u06a.txt')
  do i=1, 150
      write(13,*) (u2(i,j), j=1,150)
  enddo
  close(13)
end program
