PROGRAM heat_equation_solver
IMPLICIT NONE

! Parameters
INTEGER, PARAMETER :: nx = 30 ! Number of spatial points
INTEGER, PARAMETER :: nt = 100 ! Number of time steps
REAL, PARAMETER :: L = 1.0 ! Length of domain [0,L]
REAL, PARAMETER :: T_final = 0.5 ! Final time
REAL, PARAMETER :: alpha = 0.02 ! Thermal diffusivity
INTEGER, PARAMETER :: nout = 10 ! Number of output snapshots

! Variables
REAL :: dx, dt
REAL :: r ! r = alpha*dt/dx^2
REAL :: u(0:nx), u_new(0:nx)
REAL :: x(0:nx), t
INTEGER :: i, n, iout

! Variables for tridiagonal solver
REAL :: a(nx - 1), b(nx - 1), c(nx - 1), d(nx - 1)

! Data output
CHARACTER(LEN=100) :: filename

! Initialize parameters
dx = L / REAL(nx)
dt = T_final / REAL(nt)
r = alpha * dt / (dx * dx)

! Initialize spatial grid
DO i = 0, nx
  x(i) = i * dx
END DO

! Initial condition: u(x,0) = sin(pi*x)
DO i = 0, nx
  u(i) = SIN(3.14159265359 * x(i))
END DO

! Output initial condition to CSV file
WRITE (filename, '(A,I0.3,A)') 'heat_solution_t', 0, '.csv'
OPEN (UNIT=10, FILE=TRIM(filename), STATUS='REPLACE')

! Write header
WRITE (10, '(A)') 'x,u'

! Write data
DO i = 0, nx
  WRITE (10, '(F10.6,A,F10.6)') x(i), ',', u(i)
END DO

CLOSE (10)

! Set up tridiagonal system coefficients for backward Euler
! The system is:  -r*u_{i-1} + (1+2r)*u_i - r*u_{i+1} = u^n_i
! Where u^n_i is the value from the previous time step

! Coefficient of u_{i-1}: a[i] = -r
! Coefficient of u_i: b[i] = (1+2r)
! Coefficient of u_{i+1}: c[i] = -r
DO i = 1, nx - 1
  a(i) = -r
  b(i) = 1.0 + 2.0 * r
  c(i) = -r
END DO

! Time stepping loop
t = 0.0
iout = 0
DO n = 1, nt
  ! Update time
  t = t + dt

  ! Set up the right-hand side for the tridiagonal system: d[i] = u^n_i
  DO i = 1, nx - 1
    d(i) = u(i)
  END DO

  ! Adjust for boundary conditions
  d(1) = d(1) + r * u(0) ! Left boundary contribution
  d(nx - 1) = d(nx - 1) + r * u(nx) ! Right boundary contribution

  ! Solve the tridiagonal system using Thomas algorithm
  CALL thomas_algorithm(a, b, c, d, nx - 1, u_new(1:nx - 1))

  ! Apply boundary conditions (fixed at 0)
  u_new(0) = 0.0
  u_new(nx) = 0.0

  ! Update solution for next time step
  u = u_new

  ! Output solution at specific time steps
  IF (n == NINT(REAL(iout + 1) * REAL(nt) / REAL(nout))) THEN
    iout = iout + 1

    ! Create new CSV file for this time step
    WRITE (filename, '(A,I0.3,A)') 'heat_solution_t', iout, '.csv'
    OPEN (UNIT=10, FILE=TRIM(filename), STATUS='REPLACE')

    ! Write header
    WRITE (10, '(A)') 'x,u'

    ! Write data with comma delimiter
    DO i = 0, nx
      WRITE (10, '(F10.6,A,F10.6)') x(i), ',', u(i)
    END DO

    CLOSE (10)

    ! Print progress
            PRINT '(A,F8.6,A,A)', 'Output snapshot at t = ', t, ' to file: ', TRIM(filename)
  END IF
END DO

! Print information to console
PRINT *, ""
PRINT *, "Simulation completed!"
PRINT *, "Parameters:"
PRINT *, "  dx =", dx
PRINT *, "  dt =", dt
PRINT *, "  r =", r
PRINT *, "  Method: Backward Euler (implicit) - Unconditionally stable"
PRINT *, ""
PRINT *, "Results saved to CSV files: heat_solution_t*.csv"
PRINT *, "  - Total files:", nout + 1, " (including initial condition)"
PRINT *, "  - Format: comma-delimited CSV with header 'x,u'"

CONTAINS
! Thomas algorithm for solving tridiagonal systems
! a: lower diagonal, b: main diagonal, c: upper diagonal, d: right-hand side
! n: system size, x: solution vector
SUBROUTINE thomas_algorithm(a, b, c, d, n, x)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  REAL, INTENT(IN) :: a(n), b(n), c(n), d(n)
  REAL, INTENT(OUT) :: x(n)
  REAL :: c_prime(n), d_prime(n)
  INTEGER :: i

  ! Forward elimination
  c_prime(1) = c(1) / b(1)
  d_prime(1) = d(1) / b(1)

  DO i = 2, n
    c_prime(i) = c(i) / (b(i) - a(i) * c_prime(i - 1))
  d_prime(i) = (d(i) - a(i) * d_prime(i - 1)) / (b(i) - a(i) * c_prime(i - 1))
  END DO

  ! Back substitution
  x(n) = d_prime(n)

  DO i = n - 1, 1, -1
    x(i) = d_prime(i) - c_prime(i) * x(i + 1)
  END DO
END SUBROUTINE thomas_algorithm

END PROGRAM heat_equation_solver
