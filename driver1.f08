! This program is a driver program to test the subroutines developed
! to compute the 35 objective functions and their first-, second- and
! third-order derivatives. It uses the new Fortran 2008 routines
! developed.
!
! Just like the driver program that is contained in Algorithm 566 and
! its remark, this driver program uses Taylor expansions in order to
! verify the resulting error from third-order expansion, which
! indicates whether the first-, second- and third-order derivative of
! the objective function are correct. The explanation given below is
! very similar to the one given by the authors in the Algorithm 566.
!
! It is employed a Taylor expansion of the objective function around
! the given initial point XC for each problem. The Taylor series
! is expanded at xc + eps*y, where eps is a scalar and y is a random
! permutation vector.
!
! Denoting by g, H, and T the first-, second-, and third-order
! derivative of the objective function, respectively, and by Mv^a the
! result of applying the a-order tensor M to the vector v (e.g. gy is
! the dot product between the gradient and the vector y and Hy^2 is
! the dot product between the vector y and the matrix-vector product
! Hy - this is always a scalar), the Taylor expression of order 3
! around xc is given by
!
! f(xc+eps*y) = f(xc) + eps*gy + (1/2) eps^2 Hy^2 + (1/6) eps^3 Ty^3
!                     + O(eps^4).                                     (***)
!
! Notice that, if first-order Taylor expansion is employed, the error
! is O(eps^2) and we are able to check the correctness of the
! first-order derivative, and if second-order Taylor is employed, the
! error is O(eps^3) and we are able to check the correctness of the
! first- and second-order derivatives.
!
! The test is iterative and consist in dividing eps by 2 at each
! iteration. Then, the objective function is evaluated at xc+eps*y and
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

program driver1

  use mgh
  
  implicit none

  ! -----------------------------------------------------------------------
  ! DESCRIPTION OF LOCAL VARIABLES
  ! -----------------------------------------------------------------------
  ! macheps      machine precision
  ! -----------------------------------------------------------------------
  ! allocstat    the status variable for allocate
  ! flag         the status variable for the routines
  ! i,j,k,l      counters
  ! m            number of equations for the given problem
  ! n            number of variables for the given problem
  ! ntries       number of tries to test the derivatives
  ! problem      the number of the problem to be considered
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
  ! fc           the objective function value at the initial point
  ! fnew         the objective function value at xc+eps*y
  ! fold         the objective function value at the previous point
  ! gy           the first-order term of the Taylor expansion
  ! rand         a random number obtained using drand
  ! ratio        the ratio between old/new errors
  ! seed         the seed used in random number generator
  ! t3y          the third-order term of the Taylor expansion
  ! taylor       the Taylor expansion of f at xc+eps*y
  ! yhy          the second-order term of the Taylor expansion
  ! -----------------------------------------------------------------------
  ! gc           the gradient of the objective function at xc
  ! hc           the Hessian of the objective function at xc
  ! hy           the Hessian applied over the array y
  ! tc           the third-order derivative of the objective function at xc
  ! xc           the initial point of the problem
  ! xnew         xc+eps*y
  ! y            a permutation vector
  ! -----------------------------------------------------------------------
  
  ! LOCAL PARAMETERS
  real(kind=rk), parameter :: macheps = epsilon( 0.0_rk )

  ! LOCAL SCALARS
  character(len=60) :: name
  integer           :: allocstat, flag, i, j, k, l, m, n, ntries, problem, &
                       testder
  real(kind=rk)     :: diff, diffold, eps, epslim, factor, fc, fnew, fold, &
                       gy, ratio, rand, seed, t3y, taylor, yhy

  ! LOCAL ARRAYS
  real(kind=rk), allocatable, dimension(:)     :: gc, hy, xc, xnew, y
  real(kind=rk), allocatable, dimension(:,:)   :: hc
  real(kind=rk), allocatable, dimension(:,:,:) :: tc

  ! Change ntries and testder as your need.
  ntries     =  4
  testder    =  3

  ! Seed for drand
  seed = 123456.0_rk

  open ( 50, file = "driver1.out" )

  do problem = 1, 35
     ! Set the number of the problem given by the user
     call mgh_set_problem( problem, flag )

     ! Set the number of variables and equations for the chosen
     ! problem. The subroutine mgh_get_dims validates the input m and n,
     ! and returns the correct values for m and n
     call mgh_get_dims( n, m )

     ! Allocate arrays
     allocate( gc(n), hc(n,n), hy(n), tc(n,n,n), xc(n), xnew(n), y(n), &
          stat=allocstat )
     if ( allocstat .ne. 0 ) then
        write ( *, * ) "Allocation memory error on DRIVER1."
        stop
     end if

     ! Show preliminar label
     call mgh_get_name( name )
     write ( 50, 110 ) trim( name ), n, m

     epslim = macheps*real( n*n )*1.e+2_rk

     ! The factor by which the initial point is scaled
     factor = 1.0_rk

     do i = 1, ntries
        ! Get the initial point for the given problem
        call mgh_get_x0( xc(1:n), factor )

        ! Compute the objective function and its first-, second- and
        ! third-order derivatives for the initial point
        call mgh_evalf( xc(1:n), fc, flag )
        if ( flag .ne. 0 ) then
           write ( *, * ) "Error during evalf, flag ", flag
           stop
        end if

        call mgh_evalg( xc(1:n), gc(1:n), flag )
        if ( flag .ne. 0 ) then
           write ( *, * ) "Error during evalg, flag ", flag
           stop
        end if

        call mgh_evalh( xc(1:n), hc(1:n,1:n), flag )
        if ( flag .ne. 0 ) then
           write ( *, * ) "Error during evalh, flag ", flag
           stop
        end if

        call mgh_evalt( xc(1:n), tc(1:n,1:n,1:n), flag )
        if ( flag .ne. 0 ) then
           write ( *, * ) "Error during evalt, flag ", flag
           stop
        end if

        ! Obtain a random perturbation vector
        do j = 1, n
           ! Generates a random number between -1 and 1
           rand = 2.0_rk*drand( seed ) - 1.0_rk

           ! Defines the permutation vector
           if ( xc(j) .ne. 0.0_rk ) then
              y(j) = rand*xc(j)
           else
              y(j) = rand
           end if
        end do

        ! Print xc in the screen
        write ( 50, 130, advance='no' )
        write ( 50, 100 ) ( xc(j), j = 1, n )

        ! Print the perturbation vector y in the screen
        write ( 50, 140, advance='no' )
        write ( 50, 100 ) ( y(j), j = 1, n )

        ! In the following code, we compute the first-, second-, and
        ! third-order terms of Taylor expansion. Depending on the value
        ! of ** testder **, which defines the order of the Taylor
        ! expansion chosen by the user, the respective term of the Taylor
        ! expansion is set to zero.

        ! Computes the first-order term of the Taylor expansion
        gy  = dot_product( gc(1:n), y(1:n) )
        if ( testder .ge. 2 ) then
           ! Computes the second-order term of the Taylor expansion
           hy = 0.0_rk
           do l = 1, n
              hy(l) = hy(l) + hc(l,l)*y(l)
              do j = 1, l-1
                 hy(l) = hy(l) + hc(j,l)*y(j)
                 hy(j) = hy(j) + hc(j,l)*y(l)
              end do
           end do
           yhy = dot_product( y(1:n), hy(1:n) )

           t3y = 0.0_rk
           if ( testder .ge. 3 ) then
              ! Computes the third-order term of the Taylor expansion
              do l = 1, n
                 do k = 1, l
                    do j = 1, k
                       if ( j .eq. k .and. k .eq. l ) then
                          t3y = t3y + tc(j,k,l)*y(j)*y(k)*y(l)
                       elseif ( j .eq. k .or. k .eq. l ) then
                          t3y = t3y + 3.0_rk*tc(j,k,l)*y(j)*y(k)*y(l)
                       else
                          t3y = t3y + 6.0_rk*tc(j,k,l)*y(j)*y(k)*y(l)
                       end if
                    end do
                 end do
              end do
           end if
        else
           yhy = 0
           t3y = 0
        end if

        eps  = 0.5_rk
        diff = 0.0_rk
        fnew = fc

        write ( 50, 940 )

        do
           xnew(1:n) = xc(1:n) + eps * y(1:n)

           fold = fnew
           call mgh_evalf( xnew, fnew, flag )
           if ( flag .ne. 0 ) then
              write ( *, * ) "Error during evalf, flag = ", flag
              stop
           end if

           taylor = fc                  &
                + eps * gy              &
                + 0.5_rk * eps**2 * yhy &
                + eps**3/6.0_rk * t3y

           diffold = diff
           diff    = fnew - taylor

           if ( abs( diff ) .lt. abs( epslim * fnew ) ) then
              write ( 50, 920 ) abs( diff ), abs( epslim * fnew )
              exit
           end if

           if ( abs( fnew - fold ) .lt. abs( epslim * fold ) ) then
              write ( 50, 930 ) abs( fnew - fold ), abs( epslim * fold )
              exit
           end if

           if ( diffold .eq. 0.0_rk .or. diff .eq. 0.0_rk ) then
              write ( 50, 950 ) eps,fnew,taylor,diff
           else
              ratio = diffold / diff
              write ( 50, 950 ) eps,fnew,taylor,diff,ratio
           end if

           eps = 0.5_rk * eps
           if ( eps .le. macheps ) exit
        end do

        factor = 5.0_rk * factor
     end do

     ! Deallocate arrays
     deallocate( gc, hc, hy, tc, xc, xnew, y, stat=allocstat )
     if ( allocstat .ne. 0 ) then
        write ( *, * ) "Deallocation error on DRIVER1."
        stop
     end if
  end do

  close ( 50 )
  write ( *, * ) "Check file 'driver1.out' for driver1 output."

  ! NONEXECUTABLE STATEMENTS
100 format(    6( F10.4, 1X ) )
110 format(  /, /, 90('='),              &
             /, 1X, 'PROBLEM NAME: ', A, &
             /, /,  '        n = ', I0,  &
                /,  '        m = ', I0, / )
130 format( /, 1X, 'XC = ' )
140 format( /, 1X, 'Y  = ' )

920 format(/T5,'diff is small (', 1PE16.8, ', less than ', 1PE16.8, &
               ' in absolute value)'/)
930 format(/T5,'change in function value is very small (', 1PE16.8, ', less than ', &
            1PE16.8,' in absolute value)'/)
940 format(/,1X,74('-'),&
           /,1X,7X,'EPS',15X,'F',10X,'TAYLOR',12X,'DIFF',11X,'RATIO', &
           /,1X,74('-'))
950 format(1X,1PE10.4,1PE16.8,1PE16.8,1PE16.8,1PE16.8)

contains

  function drand(ix)
    
    ! This is the random number generator of Schrage:
    !
    ! L. Schrage, A more portable Fortran random number generator, ACM
    ! Transactions on Mathematical Software 5 (1979), 132-138.
    
    ! FUNCTION TYPE
    real(kind=rk) :: drand

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

end program driver1
