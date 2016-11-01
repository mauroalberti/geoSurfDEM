! Invert intersection points
! to derive a best-fit-plane

! By Mauro Alberti - alberti.m65@gmail.com
! October-November 2016

! Contains modified code originally from:
!  - NSIDC web site for GLAS files reading
!  - John Burkardt and released with LGPL license, for Lapack testing

!   Licensing:
!   This code is distributed under the GNU GPL v. 3 license.


module kinds

    integer, parameter :: i1b = selected_int_kind(1)
    integer, parameter :: i2b = selected_int_kind(4)
    integer, parameter :: i4b = selected_int_kind(9)
    integer, parameter :: r4b = selected_real_kind( 6, 37)
    integer, parameter :: r8b = selected_real_kind(15, 307)

end module kinds


subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows in A.
!
!    Input, integer N, the number of columns in A.
!
!    Input, real ( kind = r8b ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end


subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = r8b ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

    end do

  end do

  write ( *, '(a)' ) ' '

  return
end


subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, real ( kind = r8b ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,2x,g16.8)' ) i, a(i)
  end do

  return
end


subroutine svd(m, n, a, s, u, v)

    use kinds

    ! modified from: TEST07 tests DGESVD
    !                 by John Burkardt

    integer ( kind = i4b ), intent(in) :: m, n
    real ( kind = r8b ), intent(in) :: a(m, n)

    real ( kind = r8b ), intent(out) :: s(min(m, n))
    real ( kind = r8b ), intent(out) :: u(m, m)
    real ( kind = r8b ), intent(out) :: v(n, n)

    ! internal variables

    character :: jobu
    character :: jobvt
    integer ( kind = i4b ) :: lwork
    integer ( kind = i4b ) :: i
    integer ( kind = i4b ) :: info
    integer ( kind = i4b ) :: lda
    integer ( kind = i4b ) :: ldu
    integer ( kind = i4b ) :: ldvt
    real ( kind = r8b ) :: vt(n,n)
    real ( kind = r8b ) :: sigma(m, n)
    real ( kind = r8b ) :: b(m, n) ! debug variable
    real ( kind = r8b ), allocatable :: work(:)

    ! debug printouts

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  For a double precision real matrix (D)'
    write ( *, '(a)' ) '  in general storage mode (GE):'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  DGESVD computes the singular value decomposition:'
    write ( *, '(a)' ) '    A = U * S * V'''

    call r8mat_print(m, n, a, '  The matrix A:')

    !  Compute the singular values and singular vectors.

    lwork = 3*min(m,n) + max ( max(m,n), 2*min(m,n) )
    allocate (work(lwork))

    jobu = 'A'
    jobvt = 'A'
    lda = m
    ldu = m
    ldvt = n

    call dgesvd ( jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, &
                  lwork, info )

    if (info /= 0) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  DGESVD returned nonzero INFO = ', info
    end if

    call r8vec_print ( min ( m, n ), s, '  Singular values' )
    call r8mat_print ( m, m, u, '  Left singular vectors U:' )
    call r8mat_print ( n, n, vt, '  Right singular vectors V'':' )

    sigma(1:m, 1:n) = 0.0D+00
    do i = 1, min ( m, n )
        sigma(i,i) = s(i)
    end do

    b(1:m, 1:n) = matmul ( u(1:m, 1:m), matmul ( sigma(1:m, 1:n), vt(1:n, 1:n) ) )

    call r8mat_print ( m, n, b, '  The product U * S * V'':' )

    v = transpose(vt)

 end subroutine svd


! write program header on screen

subroutine WelcomeScreen()

    write (*,"(A)") '*****************************************'
    write (*,*)
    write (*,"(16X,A)") 'geoSurfDEM - Invert'
    write (*,*)
    write (*,"(16X,A)")	'vers. 0.0.1'
    write (*,"(17X,A)") '2016-11-31'
    write (*,*)
    write (*,"(5X,A)") 'Fortran program for determining Best-Fit-Plane to surface traces'
    write (*,*)
    write (*,"(A)") '*****************************************'

    write (*,*)
    write (*,*)
    write (*,"(4X,A)") 'Input parameters definition'
    write (*,*)

end subroutine WelcomeScreen


! choice and opening of input file, that stores the list of data sources

subroutine InputFileDef()

    character (len=50) :: InputFileName
    integer :: ios
    logical :: exists

    do
        write(*,"(A)", ADVANCE="no") 'Enter the name of input file: '
        read(*,*) InputFileName
        inquire(file=trim(InputFileName),exist=exists)
        if (.NOT.exists) then
            write(*,*) 'INPUT FILE: not found'
            cycle
        end if
        ! open the input file with sequential access
        open(unit=15,file=trim(InputFileName),status='old',access='sequential'  &
            ,form='formatted',iostat=ios)
        if (ios /= 0) then
            write(*,*) 'INPUT FILE: not opened'
            cycle
        end if
        exit
    end do

	write (*,*)

end subroutine InputFileDef


! define output files storing results

subroutine OutputFileDef()

	character (len= 37) :: OutputFileName

	do
		write (*,"(A)", ADVANCE="no") 'Enter output filename: '
		read (*,*) OutputFileName
		! create output file
		open (unit=16, file=trim(OutputFileName), status='NEW' &
  		, access='sequential', form='formatted', iostat=ios)
		if (ios /= 0) then
  			write (*,"(A,/,A)") 'Error with output file creation.','Change name'
  			cycle
		end if
        exit
	end do

    write (*,*)

end subroutine OutputFileDef


program invert_attitudes

    use kinds

    implicit none

    !type(point) :: pt1, pt2, pt3
    !pt1 = point(0.0, 0.0, 0.0)
    !pt2 = point(1.0, 0.0, 0.0)
    !pt3 = point(0.0, 1.0, 0.0)

    integer ( kind = i4b ) :: m
    integer ( kind = i4b ), parameter :: n = 3
    real ( kind = r8b ) :: x, y, z
    real ( kind = r8b ), dimension(:, :), allocatable :: a
    real ( kind = r8b ), dimension(:), allocatable :: s
    real ( kind = r8b ), dimension(:, :), allocatable :: u
    real ( kind = r8b ), dimension(n, n) :: v
    integer ( kind = i2b ) :: num_points, i
    integer ( kind = i2b ) :: ios

    call WelcomeScreen()

    call InputFileDef()

    call OutputFileDef()

    num_points = 0
    do
        read (15, *, iostat=ios) x, y, z
        ! if finished reading, program stops
        if (ios/=0) then
            exit
        else
            num_points = num_points + 1
        end if
    end do
    m = num_points
    allocate(a(m, n), s(min(m, n)), u(m, m))

    rewind 15

    do i = 1, m
        read (15, *, iostat=ios) a(i, 1), a(i, 2), a(i, 3)
    end do

    !a = reshape([0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0], shape(a)) ! column-major

    call svd(m, n, a, s, u, v)

    call r8mat_print ( n, n, v, '  The matrix V:' )

    write ( *, '(a)' ) 'Program completed'


end program
