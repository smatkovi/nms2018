program helium
  real :: y(1001)
  integer :: n
  
  y(1)=1.0
  y(2)=0.0
  do n=3, 1000
    y(n) = 2.0*y(n-1)*(1.0 - 5.0/12.0*g(n*1.0-1.0) ) - y(n-2)*(1.0 + (1.0/12.0)*g(n*1.0-2.0) )
    y(n) = y(n)/(1.0 + g(n*1.0)/12.0)
  end do
  write(*,*) y

  contains
  real function g(r)
    real :: r
    g = 1.0/r 
  end function g
end program helium
