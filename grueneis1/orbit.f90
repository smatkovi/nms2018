PROGRAM ORBIT
   !USE constants
   IMPLICIT NONE
   
   !real :: force_sun
!  Define position and velocity Vector/Arrays
   REAL, DIMENSION(2) :: r,v,rn,vn

!  Time and time step
   REAL :: t, dt=1
   INTEGER :: N, MAXN=3650

!  Intermediate quantities for Runge-Kutta (RK4)
   REAL, DIMENSION(2) :: r1, r2, r3, r4, v1, v2, v3, v4
   !real, dimension(2) :: f1, f2
!  Define simulation length and time step size first

   WRITE(*,*)'!Welcome to Orbit calculator!'
!   WRITE(*,*)'!Time step is set to:       ',dt,' days.'
!   WRITE(*,*)'!Number of steps is set to: ',MAXN,'.'

!  Define initial paramters
   t=0
   CALL init_earth(r,v)
   
   OPEN(unit=11,file="EarthOrbit_Euler.dat")
   OPEN(unit=12,file="EarthOrbit_Energy_Euler.dat")
   WRITE(*,*)' '
   WRITE(*,*)'!Calculating EarthOrbit using first order Euler method.'
   WRITE(*,*)'!Results will be written to EarthOrbit_Euler.dat and EarthOrbit_Energy_Euler.dat files.'
   DO N=1,MAXN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  IMPLEMENT EULER INTEGRATION HERE  !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      rn = r + dt*N*v
      vn = v - dt*N*13.3159162*333054.253*r/(norm2(r)**3)
      WRITE(11,*)N*dt, r
!      WRITE(12,*)t, ENERGY(Me,r,v)
      r=rn
      v=vn

   ENDDO

   CLOSE(11)
   CLOSE(12)


!  Define initial paramters for earth again
   t=0
   CALL init_earth(r,v)

   OPEN(unit=11,file="EarthOrbit_RK4.dat")
!   OPEN(unit=12,file="EarthOrbit_Energy_RK4.dat")
   WRITE(*,*)' '
   WRITE(*,*)'!Calculating EarthOrbit using RK4 method.'
   WRITE(*,*)'!Results will be written to EarthOrbit_RK4.dat and EarthOrbit_Energy_RK4.dat files.'

   DO N=1,MAXN
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  IMPLEMENT RK4   INTEGRATION HERE  !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!     Integrate r and v using RK4
      r1 =  f1(v)
      r2 =  f1(v + dt/2.0*v1)
      r3 =  f1(v + dt/2.0*v2)
      r4 =  f1(v + dt*v3) 
      rn = r + dt*(v1 + 2.0*v2 + 2.0*v3 + v4)/6.0

      v1 = f2(r)
      v2 = f2(r + dt/2.0*r1)
      v3 = f2(r + dt/2.0*r2)
      v4 = f2(r + dt*r3) 
      vn = v + dt*(v1 + 2.0*v2 + 2.0*v3 + v4)/6.0
      
      
      WRITE(11,*)t, r
      r=rn
      v=vn
!      WRITE(12,*)t, ENERGY(Me,r,v)
   ENDDO
   
   CLOSE(11)
   CLOSE(12)




CONTAINS

!                                                                                           ->
! This function returns the gravitational force of the sun acting on a mass m at coordinate  r  
!

!real function FORCE_SUN(m,r)
!   USE constants
!   IMPLICIT NONE
!   REAL, DIMENSION(2) :: FORCE_SUN
!   REAL, DIMENSION(2) :: r      ! position of mass
!   REAL :: m                    ! and mass m

 !  FORCE_SUN
!END


!                                                                                  ->
! This function returns the kinetic and potential energy of a mass m at coordinate  r  in the grav. field of the sun
!

!FUNCTION ENERGY(m,r,v)
 !  USE constants
  ! IMPLICIT NONE
   !REAL :: ENERGY
   !REAL :: m
   !REAL, DIMENSION(2) :: r,v

!
!   GESAMTENERGIE MUSS HIER RICHTIG BERECHNET WERDEN.
!

  ! RETURN
!END

!                                 ->    ->
! This routine defines the initial r and v values for the earth.
real, dimension(2) function f1(v)
  real, dimension(2) :: v
  f1 = v
end function f1
!contains
real, dimension(2) function f2(r)
  real, dimension(2) :: r
  f2 = -13.3159162*1.989e30*r/(norm2(r)**3)
end function f2
subroutine init_earth(r,v)
   IMPLICIT NONE
   REAL, DIMENSION(2) :: r,v

!
!   ANFANGSWERTE MUESSEN HIER EINGESETZT WERDEN
!
   r = (/ 1.0 , 0.0 /)
   v = (/ 0.0 , -0.017326 /)

end subroutine init_earth
end program
