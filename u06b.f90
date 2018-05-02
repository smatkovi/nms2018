program wellenpaket
  implicit none
  integer :: dm = 9.9999999999e-1, m, n 
  real :: k0 = 27.1
  complex :: i
  complex, dimension(1001,1001) :: psi
  complex, dimension(1001) :: psitemp

  i = (0, 1.)
  do m = 1, 1000
    psi(m, 1) = exp(-1.*((dm*1.0*m)**2)/2. + i*k0*1.0*m*dm)
  enddo

  open(unit = 13, file = 'u06b.txt')
  do m = 1, 1000
    do n = 1, 1000
      !write(13,*) psi(m, n)
    enddo
  enddo
  close(13)
end program

function V(x)
  real :: V
  real ::x
  V = -1.*exp(-3.*abs(x))/(abs(x) + 1.)
end function
