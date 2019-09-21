! This program is a driver program to test the subroutines developed
! to compute the 35 objective functions and their first-, second- and
! third-order derivatives. It uses the classical routines proposed by
! Moré, Garbow, and Hillstrom (OBJFCN and GRDFCN, 1981), the HESFCN
! proposed by Averbukh, Figueroa, and Schlick (1994), and a new THRFCN
! to compute third-order derivatives.
!
! Just like the driver program that is contained in Algorithm 566 and
! its remark, this driver program uses Taylor expansions in order to
! verify the resulting error from third-order expansion, which
! indicates whether the first-, second- and third-order derivative of
! the objective function are correct. The explanation given below is
! very similar to the one given by the authors in the Algorithm 566.
!
! It is employed a Taylor expansion of the objective function around
! the given initial point X0 for each problem. The Taylor series
! is expanded at x0 + eps*y, where eps is a scalar and y is a random
! permutation vector.
!
! Denoting by g, H, and T the first-, second-, and third-order
! derivative of the objective function, respectively, and by Mv^a the
! result of applying the a-order tensor M to the vector v (e.g. gy is
! the dot product between the gradient and the vector y and Hy^2 is
! the dot product between the vector y and the matrix-vector product
! Hy - this is always a scalar), the Taylor expression of order 3
! around x0 is given by
!
! f(x0+eps*y) = f(x0) + eps*gy + (1/2) eps^2 Hy^2 + (1/6) eps^3 Ty^3
!                     + O(eps^4).                                     (***)
!
! Notice that, if first-order Taylor expansion is employed, the error
! is O(eps^2) and we are able to check the correctness of the
! first-order derivative, and if second-order Taylor is employed, the
! error is O(eps^3) and we are able to check the correctness of the
! first- and second-order derivatives.
!
! The test is iterative and consist in dividing eps by 2 at each
! iteration. Then, the objective function is evaluated at x0+eps*y and
! the Taylor expansion (***) is computed for the given eps. The
! difference between them gives the error. Then, if the error for EPS
! is E1, the error for EPS/2 should be E1/4 if the first-order
! derivative is correct, E1/8 if the second-order derivative is
! correct and E1/16 if the third-order derivative is correct. The
! variable RATIO contains the OLD/NEW error, at each iteration.
!
! Therefore, this driver gives an output of RATIO until EPS and/or the
! difference between the old and new objective function value is very
! small.
!
! If RATIO tends to 4, 8, or 16, the first-order derivative is
! correct, the first- and second-order derivatives are correct, or the
! first-, second-, and third-order derivatives are correct.

program driver2

  implicit none

  ! -----------------------------------------------------------------------
  ! DESCRIPTION OF LOCAL VARIABLES
  ! -----------------------------------------------------------------------
  ! macheps      machine precision
  ! rk           real number precision (change it according to set_parameter)
  ! -----------------------------------------------------------------------
  ! allocstat    the status variable for allocate
  ! i,j          counters
  ! n            number of variables for the given problem
  ! ntries       number of tries to test the derivatives
  ! problem      the number of the problem to be considered
  ! sizehu       the size of Hessian array
  ! sizetu       the size of third derivative tensor array
  ! testder      the order of the Taylor expansion to be used
  ! -----------------------------------------------------------------------
  ! diff         the current difference between the real objective
  !              function value and the Taylor expansion
  ! diffold      the previous difference between the real objective
  !              function value and the Taylor expansion
  ! eps          the size of the step to give from the current point
  !              along the direction y
  ! epslim       a small scalar to measure the "small variations"
  ! factor       the factor to which the initial point should be scaled
  ! f0           the objective function value at the initial point
  ! fc           the objective function value at x0+eps*y
  ! fold         the objective function value at the previous point
  ! gy           the first-order term of the Taylor expansion
  ! rand         a random number obtained using drand
  ! ratio        the ratio between old/new errors
  ! seed         seed for random number generator
  ! t3y          the third-order term of the Taylor expansion
  ! taylor       the Taylor expansion of f at x0+eps*y
  ! yhy          the second-order term of the Taylor expansion
  ! -----------------------------------------------------------------------
  ! pdim         an integer array that contains n for each problem
  ! g0           the gradient of the objective function at x0
  ! hd0,hu0      the Hessian of the objective function at x0
  ! td0,tu0      the third-order derivative of the objective function at x0
  ! x0           the initial point of the problem
  ! xc           x0+eps*y
  ! y            a permutation vector
  ! -----------------------------------------------------------------------
  
  ! LOCAL PARAMETERS
  integer,       parameter :: rk      = kind( 0.0d0 )
  real(kind=rk), parameter :: macheps = epsilon( 0.0d0 )

  ! LOCAL SCALARS
  integer       :: allocstat, i, j, n, ntries, problem, sizehu,      &
                   sizetu, testder
  real(kind=rk) :: diff, diffold, eps, epslim, factor, f0, fc, fold, &
                   gy, ratio, rand, seed, t3y, taylor, yhy

  ! LOCAL ARRAYS
  integer,                    dimension(18) :: pdim
  real(kind=rk), allocatable, dimension(:)  :: g0, hd0, hu0, td0, tu0, &
                                               x0, xc, y

  ! Change ntries and testder as your need.
  ntries  =  2
  testder =  3

  ! Dimension of the 18 unconstrained minimization problems
  pdim = [ 3, 6, 3,  2,  3, 10, 6, 4, 4, 2, 4, 3, 10, 10, 12, 2, 4, 8 ]

  ! Seed for drand
  seed = 123456.0_rk

  open ( 50, file = "driver2.out" )

  do problem = 1, 18
     write ( 50, 110 ) problem

     n = pdim(problem)

     sizehu = n*(n-1)/2
     sizetu = ( (n-1)*( (n-2)*(n-3) + 9*(n-2) + 12 ) )/6

     ! Allocate arrays
     allocate( g0(n), hd0(n), hu0(sizehu), td0(n), tu0(sizetu), x0(n), &
          xc(n), y(n), stat=allocstat )
     if ( allocstat .ne. 0 ) then
        write ( *, * ) "Allocation error in DRIVER2."
        stop
     end if

     epslim = macheps * real( n*n ) * 1.d+2

     ! The factor by which the initial point is scaled
     factor = 1.0d0

     do i = 1, ntries
        ! Get the initial point for the given problem
        call initpt( n, x0(1:n), problem, factor )

        ! Compute the objective function and its first-, second- and
        ! third-order derivatives for the initial point
        call objfcn( n, x0(1:n), f0, problem )
        call grdfcn( n, x0(1:n), g0(1:n), problem )
        call hesfcn( n, x0(1:n), hd0(1:n), hu0(1:sizehu), problem )
        call trdfcn( n, x0(1:n), td0(1:n), tu0(1:sizetu), problem )

        ! Obtain a random perturbation vector
        do j = 1, n
           ! Generates a random number between -1 and 1
           rand = 2.0_rk*drand( seed ) - 1.0_rk

           ! Defines the permutation vector
           if ( x0(j) .ne. 0.0d0 ) then
              y(j) = rand*x0(j)
           else
              y(j) = rand
           end if
        end do

        ! Print x0 in the screen
        write ( 50, 130, advance='no' )
        write ( 50, 100 ) ( x0(j), j = 1, n )

        ! Print the perturbation vector y in the screen
        write ( 50, 140, advance='no' )
        write ( 50, 100 ) ( y(j), j = 1, n )

        ! In the following code, we compute the first-, second-, and
        ! third-order terms of Taylor expansion. Depending on the value
        ! of ** testder **, which defines the order of the Taylor
        ! expansion chosen by the user, the respective term of the Taylor
        ! expansion is set to zero.

        ! Computes the first-order term of the Taylor expansion
        gy  = dot_product( g0(1:n), y(1:n) )
        if ( testder .ge. 2 ) then
           ! Computes the second-order term of the Taylor expansion
           call vmvprod( n, hd0(1:n), hu0(1:sizehu), y(1:n), yhy )

           if ( testder .ge. 3 ) then
              ! Computes the third-order term of the Taylor expansion
              call tvprod( n, td0(1:n), tu0(1:sizetu), y(1:n), t3y )
           else
              t3y = 0
           end if
        else
           yhy = 0
           t3y = 0
        end if

        eps  = 0.5d0
        diff = 0.0d0
        fc   = f0

        write ( 50, 940 )

        do
           xc(1:n) = x0(1:n) + eps * y(1:n)

           fold = fc
           call objfcn( n, xc(1:n), fc, problem )

           taylor = f0                 &
                + eps * gy             &
                + 0.5d0 * eps**2 * yhy &
                + eps**3/6.0d0 * t3y

           diffold = diff
           diff    = fc - taylor

           if ( abs( diff ) .lt. abs( epslim * fc ) ) then
              write ( 50, 920 ) abs( diff ), abs( epslim * fc )
              exit
           end if

           if ( abs( fc - fold ) .lt. abs( epslim * fold ) ) then
              write ( 50, 930 ) abs( fc - fold ), abs( epslim * fold )
              exit
           end if

           if ( diffold .eq. 0.0d0 .or. diff .eq. 0.0d0 ) then
              write ( 50, 950 ) eps,fc,taylor,diff
           else
              ratio = diffold / diff
              write ( 50, 950 ) eps,fc,taylor,diff,ratio
           end if

           eps = 0.5d0 * eps
           if (eps .le. macheps) exit
        end do

        factor = 5.0d0 * factor
     end do

     ! Deallocate arrays
     deallocate( g0, hd0, hu0, td0, tu0, x0, xc, y, stat=allocstat )
     if ( allocstat .ne. 0 ) then
        write ( *, * ) "Deallocation error in DRIVER2."
        stop
     end if
  end do

  close ( 50 )
  write ( *, * ) "Check file 'driver2.out' for driver2 output."

  ! NONEXECUTABLE STATEMENTS
100 format(    6( F10.4, 1X ) )
110 format( /, /, 90('='), /, 40X, 'PROBLEM ', I0, / )
130 format( /, 1X, 'X0 = ' )
140 format( /, 1X, 'Y  = ' )

920 format(/T5,'diff is small (', 1PE16.8, ', less than ', 1PE16.8, &
               ' in absolute value)'/)
930 format(/T5,'change in function value is very small (', 1PE16.8, &
            ', less than ', 1PE16.8,' in absolute value)'/)
940 format(/,1X,74('-'),&
           /,1X,7X,'EPS',15X,'F',10X,'TAYLOR',12X,'DIFF',11X,'RATIO', &
           /,1X,74('-'))
950 format(1X,1PE10.4,1PE16.8,1PE16.8,1PE16.8,1PE16.8)

contains

  ! --------------------------------------------------------------------

  subroutine vmvprod(n, diaga, uppera, x, y)

    ! vmvprod performs the vector-matrix-vector product x'*a*x, and
    ! stores the result in the scalar y. A is a symmetric nxn matrix,
    ! with diagonal elements stored in diaga and the strict upper
    ! triangular part stored by colums in uppera. x is a vector of
    ! length n.

    ! SCALAR ARGUMENTS
    integer,       intent(in) :: n
    real(kind=rk), intent(out) :: y

    ! ARRAY ARGUMENTS
    real(kind=rk), dimension(n*(n-1)/2), intent(in) :: uppera
    real(kind=rk), dimension(n),         intent(in) :: diaga, x

    ! LOCAL SCALARS
    integer :: i, j, l

    y = 0.0d0

    do i = 1, n
       y = y + diaga(i) * x(i) ** 2
    end do

    l = 1
    do j = 1, n
       do i = 1, j-1
          y = y + 2.0d0 * uppera(l)*x(i)*x(j)
          l = l + 1
       end do
    end do

  end subroutine vmvprod

  ! --------------------------------------------------------------------

  subroutine tvprod(n, diagt, uppert, x, y)

    ! tvprod performs the application of a three-dimensional array
    ! three times over x, and stores the result in the scalar y. The
    ! three-dimensional array is a symmetric n x n x n tensor, with
    ! diagonal elements stored in diagt and the strict upper part
    ! stored by columns in uppert. x is a vector of length n.

    ! SCALAR ARGUMENT
    integer,       intent(in)  :: n
    real(kind=rk), intent(out) :: y

    ! ARRAY ARGUMENTS
    real(kind=rk), dimension(*), intent(in)  :: uppert
    real(kind=rk), dimension(n), intent(in)  :: diagt, x

    ! LOCAL SCALARS
    integer :: i, j, k, l

    y = 0.0d0
    l = 1
    do k = 1, n
       do j = 1, k
          do i = 1, j
             if ( i .eq. j .and. j .eq. k ) then
                y = y + diagt(i) * x(i) ** 3
             elseif ( i .ne. j .and. j .ne. k ) then
                y = y + 6.0d0 * uppert(l)*x(i)*x(j)*x(k)
                l = l + 1
             else
                y = y + 3.0d0 * uppert(l)*x(i)*x(j)*x(k)
                l = l + 1
             end if
          end do
       end do
    end do

  end subroutine tvprod

  ! --------------------------------------------------------------------

  real(kind=rk) function drand(ix)
    
    ! This is the random number generator of Schrage:
    !
    ! L. Schrage, A more portable Fortran random number generator, ACM
    ! Transactions on Mathematical Software 5 (1979), 132-138.

    ! SCALAR ARGUMENT
    real(kind=rk), intent(inout) :: ix

    ! LOCAL SCALARS
    real(kind=rk) :: a,p,b15,b16,xhi,xalo,leftlo,fhi,k

    data a/16807.0_rk/,b15/32768.0_rk/,b16/65536.0_rk/,p/2147483647.0_rk/

    xhi= ix/b16
    xhi= xhi - dmod(xhi,1.0_rk)
    xalo= (ix-xhi*b16)*a
    leftlo= xalo/b16
    leftlo= leftlo - dmod(leftlo,1.0_rk)
    fhi= xhi*a + leftlo
    k= fhi/b15
    k= k - dmod(k,1.0_rk)
    ix= (((xalo-leftlo*b16)-p)+(fhi-k*b15)*b16)+k
    if (ix.lt.0) ix= ix + p
    drand= ix*4.656612875e-10_rk

    return

  end function drand

end program driver2
