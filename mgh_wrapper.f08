! --------------------------------------------------------------------

subroutine initpt(n, x, nprob, factor)

  use mgh

  implicit none

  ! ------------------------------------------------------------------
  ! This subroutine gives, for each problem, the corresponding initial
  ! point.
  ! ------------------------------------------------------------------
  ! DESCRIPTION OF THE ARGUMENTS
  ! ------------------------------------------------------------------
  ! n        number of variables (problem-specific)
  ! nprob    number of the problem (1..18)
  ! factor   a scalar by which the initial point will be scaled
  ! ------------------------------------------------------------------
  ! x        the output initial point
  ! ------------------------------------------------------------------
  ! DESCRIPTION OF LOCAL VARIABLES
  ! ------------------------------------------------------------------
  ! problem  maps the original problem number (1..18) to the generic
  !          problem number (1..35)
  ! ------------------------------------------------------------------
  ! flag     output flag of mgh_* routines
  ! n_mgh    the default n given by mgh_get_dims (if input n is not
  !          valid)
  ! ------------------------------------------------------------------
  
  ! SCALAR ARGUMENTS
  integer,       intent(in) :: n, nprob
  real(kind=rk), intent(in) :: factor

  ! ARRAY ARGUMENT
  real(kind=rk), dimension(n), intent(out) :: x

  ! LOCAL PARAMETER
  integer, dimension(18) :: problem = [   &
       7, 18,  9,  3, 12, 25, 20, 23, 24, &
       4, 16, 11, 26, 21, 22,  5, 14, 35 ]

  ! LOCAL SCALAR
  integer :: flag, n_mgh

  call mgh_set_problem( problem(nprob), flag )
  if ( flag .ne. 0 ) then
     write ( *, * ) "ERROR (INITPT): Problem number is invalid."
     stop
  end if

  call mgh_set_dims( n=n, flag=flag )
  if ( flag .eq. -1 ) then
     call mgh_get_dims( n_mgh )
     write ( *, 100 ) n, n_mgh
     
     if ( n .lt. n_mgh ) return
  end if

  call mgh_get_x0( x, factor )

  ! NONEXECUTABLE STATEMENTS
100 format(    1X, 'WARNING (INITPT): n given (', I0, ') is not valid.', &
            /, 1X, '                  A possible value is ', I0, '.' ) 

end subroutine initpt

! --------------------------------------------------------------------

subroutine objfcn(n, x, f, nprob)

  use mgh

  implicit none

  ! ------------------------------------------------------------------
  ! This subroutine gives, for each problem, the corresponding
  ! objective function value at x.
  ! ------------------------------------------------------------------
  ! DESCRIPTION OF THE ARGUMENTS
  ! ------------------------------------------------------------------
  ! n        number of variables (problem-specific)
  ! nprob    number of the problem (1..18)
  ! f        the output objective function value at x
  ! ------------------------------------------------------------------
  ! x        the point in which the objective function should be
  !          evaluated
  ! ------------------------------------------------------------------
  ! DESCRIPTION OF LOCAL VARIABLES
  ! ------------------------------------------------------------------
  ! problem  maps the original problem number (1..18) to the generic
  !          problem number (1..35)
  ! ------------------------------------------------------------------
  ! flag     output flag of mgh_* routines
  ! n_mgh    the default n given by mgh_get_dims (if input n is not
  !          valid)
  ! ------------------------------------------------------------------

  ! SCALAR ARGUMENTS
  integer,       intent(in)  :: n, nprob
  real(kind=rk), intent(out) :: f

  ! ARRAY ARGUMENT
  real(kind=rk), dimension(n), intent(in) :: x

  ! LOCAL PARAMETER
  integer, dimension(18) :: problem = [   &
       7, 18,  9,  3, 12, 25, 20, 23, 24, &
       4, 16, 11, 26, 21, 22,  5, 14, 35 ]

  ! LOCAL SCALAR
  integer :: flag, n_mgh

  call mgh_set_problem( problem(nprob), flag )
  if ( flag .ne. 0 ) then
     write ( *, * ) "ERROR (OBJFCN): Invalid problem number."
     stop
  end if

  call mgh_set_dims( n=n, flag=flag )
  if ( flag .eq. -1 ) then
     call mgh_get_dims( n_mgh )
     write ( *, 100 ) n, n_mgh
     
     if ( n .lt. n_mgh ) return
  end if

  call mgh_evalf( x, f, flag )
  if ( flag .ne. 0 ) then
     write ( *, * ) "WARNING (OBJFCN): mgh_evalf returned flag ", flag
  end if

    ! NONEXECUTABLE STATEMENTS
100 format(    1X, 'WARNING (OBJFCN): n given (', I0, ') is not valid.', &
            /, 1X, '                  The correct is ', I0, '.' ) 

end subroutine objfcn

! --------------------------------------------------------------------

subroutine grdfcn(n, x, g, nprob)

  use mgh

  implicit none

  ! ------------------------------------------------------------------
  ! This subroutine gives, for each problem, the corresponding
  ! objective function gradient value at x.
  ! ------------------------------------------------------------------
  ! DESCRIPTION OF THE ARGUMENTS
  ! ------------------------------------------------------------------
  ! n        number of variables (problem-specific)
  ! nprob    number of the problem (1..18)
  ! ------------------------------------------------------------------
  ! g        the output objective function gradient at x
  ! x        the point in which the gradient should be evaluated
  ! ------------------------------------------------------------------
  ! DESCRIPTION OF LOCAL VARIABLES
  ! ------------------------------------------------------------------
  ! problem  maps the original problem number (1..18) to the generic
  !          problem number (1..35)
  ! ------------------------------------------------------------------
  ! flag     output flag of mgh_* routines
  ! n_mgh    the default n given by mgh_get_dims (if input n is not
  !          valid)
  ! ------------------------------------------------------------------

  ! SCALAR ARGUMENTS
  integer, intent(in) :: n, nprob

  ! ARRAY ARGUMENTS
  real(kind=rk), dimension(n), intent(in)  :: x
  real(kind=rk), dimension(n), intent(out) :: g

  ! LOCAL PARAMETER
  integer, dimension(18) :: problem = [   &
       7, 18,  9,  3, 12, 25, 20, 23, 24, &
       4, 16, 11, 26, 21, 22,  5, 14, 35 ]

  ! LOCAL SCALAR
  integer :: flag, n_mgh

  call mgh_set_problem( problem(nprob), flag )
  if ( flag .ne. 0 ) then
     write ( *, * ) "ERROR (GRDFCN): Invalid problem number."
     stop
  end if

  call mgh_set_dims( n=n, flag=flag )
  if ( flag .eq. -1 ) then
     call mgh_get_dims( n_mgh )
     write ( *, 100 ) n, n_mgh

     if ( n .lt. n_mgh ) return
  end if

  call mgh_evalg( x, g, flag )
  if ( flag .ne. 0 ) then
     write ( *, * ) "WARNING (GRDFCN): mgh_evalg returned flag ", flag
  end if

    ! NONEXECUTABLE STATEMENTS
100 format(    1X, 'WARNING (GRDFCN): n given (', I0, ') is not valid.', &
            /, 1X, '                  The correct is ', I0, '.' ) 

end subroutine grdfcn

! --------------------------------------------------------------------

subroutine hesfcn(n, x, hesd, hesu, nprob)

  use mgh

  implicit none

  ! ------------------------------------------------------------------
  ! This subroutine gives, for each problem, the corresponding
  ! objective function Hessian evaluated at x.
  ! ------------------------------------------------------------------
  ! DESCRIPTION OF THE ARGUMENTS
  ! ------------------------------------------------------------------
  ! n        number of variables (problem-specific)
  ! nprob    number of the problem (1..18)
  ! ------------------------------------------------------------------
  ! hesd     the diagonal of the Hessian
  ! hesu     the upper triangle of the Hessian stored by columns
  ! x        the point in which the Hessian should be evaluated
  ! ------------------------------------------------------------------
  ! DESCRIPTION OF LOCAL VARIABLES
  ! ------------------------------------------------------------------
  ! problem  maps the original problem number (1..18) to the generic
  !          problem number (1..35)
  ! ------------------------------------------------------------------
  ! i,j,l    counters
  ! flag     output flag of mgh_* routines
  ! h        the Hessian computed by mgh_evalh
  ! n_mgh    the default n given by mgh_get_dims (if input n is not
  !          valid)
  ! ------------------------------------------------------------------

  ! SCALAR ARGUMENTS
  integer :: n, nprob

  ! ARRAY ARGUMENTS
  real(kind=rk), dimension(n),         intent(in)  :: x
  real(kind=rk), dimension(n),         intent(out) :: hesd
  real(kind=rk), dimension(n*(n-1)/2), intent(out) :: hesu

  ! LOCAL PARAMETER
  integer, dimension(18) :: problem = [   &
       7, 18,  9,  3, 12, 25, 20, 23, 24, &
       4, 16, 11, 26, 21, 22,  5, 14, 35 ]

  ! LOCAL SCALARS
  integer :: flag, i, j, l, n_mgh

  ! LOCAL ARRAY
  real(kind=rk), dimension(n,n) :: h

  call mgh_set_problem( problem(nprob), flag )
  if ( flag .ne. 0 ) then
     write ( *, * ) "ERROR (HESFCN): Invalid problem number."
     stop
  end if

  call mgh_set_dims( n=n, flag=flag )
  if ( flag .eq. -1 ) then
     call mgh_get_dims( n_mgh )
     write ( *, 100 ) n, n_mgh
     
     if ( n .lt. n_mgh ) return
  end if

  call mgh_evalh( x, h, flag )
  if ( flag .ne. 0 ) then
     write ( *, * ) "WARNING (HESFCN): mgh_evalh returned flag ", flag
  end if

  l = 1
  do j = 1, n
     do i = 1, j-1
        hesu(l) = h(i,j)
        l = l + 1
     end do
     hesd(j) = h(j,j)
  end do

    ! NONEXECUTABLE STATEMENTS
100 format(    1X, 'WARNING (HESFCN): n given (', I0, ') is not valid.', &
            /, 1X, '                  The correct is ', I0, '.' ) 

end subroutine hesfcn

! --------------------------------------------------------------------

subroutine trdfcn(n, x, td, tu, nprob)

  use mgh

  implicit none

  ! ------------------------------------------------------------------
  ! This subroutine gives, for each problem, the corresponding
  ! objective function third-derivative tensor evaluated at x.
  ! ------------------------------------------------------------------
  ! DESCRIPTION OF THE ARGUMENTS
  ! ------------------------------------------------------------------
  ! n        number of variables (problem-specific)
  ! nprob    number of the problem (1..18)
  ! ------------------------------------------------------------------
  ! td       the diagonal of the tensor
  ! tu       the upper part of the tensor stored by columns
  ! x        the point in which the tensor should be evaluated
  ! ------------------------------------------------------------------
  ! DESCRIPTION OF LOCAL VARIABLES
  ! ------------------------------------------------------------------
  ! problem  maps the original problem number (1..18) to the generic
  !          problem number (1..35)
  ! ------------------------------------------------------------------
  ! i,j,k,l  counters
  ! flag     output flag of mgh_* routines
  ! n_mgh    the default n given by mgh_get_dims (if input n is not
  !          valid)
  ! t        the tensor computed by mgh_evalt
  ! ------------------------------------------------------------------

  ! SCALAR ARGUMENTS
  integer, intent(in) :: n, nprob

  ! ARRAY ARGUMENTS
  real(kind=rk), dimension(n), intent(in)  :: x
  real(kind=rk), dimension(*), intent(out) :: tu
  real(kind=rk), dimension(n), intent(out) :: td

  ! LOCAL PARAMETER
  integer, dimension(18) :: problem = [   &
       7, 18,  9,  3, 12, 25, 20, 23, 24, &
       4, 16, 11, 26, 21, 22,  5, 14, 35 ]

  ! LOCAL SCALARS
  integer :: flag, i, j, k, l, n_mgh

  ! LOCAL ARRAY
  real(kind=rk), dimension(n,n,n) :: t

  call mgh_set_problem( problem(nprob), flag )
  if ( flag .ne. 0 ) then
     write ( *, * ) "ERROR (TRDFCN): Invalid problem number."
     stop
  end if

  call mgh_set_dims( n=n, flag=flag )
  if ( flag .eq. -1 ) then
     call mgh_get_dims( n_mgh )
     write ( *, 100 ) n, n_mgh
     
     if ( n .lt. n_mgh ) return
  end if

  call mgh_evalt( x, t, flag )
  if ( flag .ne. 0 ) then
     write ( *, * ) "WARNING (TRDFCN): mgh_evalt returned flag ", flag
  end if

  l = 1
  do k = 1, n
     do j = 1, k
        do i = 1, j
           if ( i .eq. j .and. j .eq. k ) then
              td(i) = t(i,j,k)
           else
              tu(l) = t(i,j,k)
              l = l + 1
           end if
        end do
     end do
  end do

  ! NONEXECUTABLE STATEMENTS
100 format(    1X, 'WARNING (TRDFCN): n given (', I0, ') is not valid.', &
            /, 1X, '                  The correct is ', I0, '.' ) 

end subroutine trdfcn
