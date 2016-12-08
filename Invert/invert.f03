! Invert intersection points
! to derive a best-fit-geol_plane

! By Mauro Alberti - alberti.m65@gmail.com
! October-November 2016

! Contains modified code originally from:
!  - NSIDC web site for GLAS files reading
!  - John Burkardt and released with LGPL license, for Lapack testing

!   Licensing:
!   This code is distributed under the GNU GPL v. 3 license.


module var_types

    implicit none

    integer, parameter :: i1b = selected_int_kind(1)
    integer, parameter :: i2b = selected_int_kind(4)
    integer, parameter :: i4b = selected_int_kind(9)
    integer, parameter :: r4b = selected_real_kind( 6, 37)
    integer, parameter :: r8b = selected_real_kind(15, 307)

    type :: real_opt
        real (kind=r8b) :: real_val
        logical :: valid
    end type real_opt

end module var_types


module math_constants

    use var_types

    implicit none

    real (kind=r8b), parameter :: pi = 3.1415926535897931_r8b  ! pi radians
    real (kind=r8b), parameter :: rad2degr = 180.0_r8b/pi, degr2rad = pi/180.0_r8b  ! for conversion from radians to degrees

end module math_constants


module vector_processing

    use var_types
    use math_constants

    implicit none

    type :: cart_vect

        sequence
        real (kind=r8b) :: x, y, z

    end type cart_vect

    type cart_vect_opt
        type(cart_vect) :: vect
        logical :: valid
    end type cart_vect_opt

    type :: orthonorm_triad
        type(cart_vect) :: X, Y, Z
    end type orthonorm_triad

    real (kind=r8b), parameter :: vect_magn_min_thresh = 1.0e-7
    real (kind=r8b), parameter :: vect_div_min_thresh = 1.0e-7

    type(orthonorm_triad), parameter :: frame0 = orthonorm_triad(cart_vect(1.0_r8b,0.0_r8b,0.0_r8b), &
                                                                 cart_vect(0.0_r8b,1.0_r8b,0.0_r8b), &
                                                                 cart_vect(0.0_r8b,0.0_r8b,1.0_r8b))

    contains

    ! cart_vect magnitude

    real (kind=r8b) function vect_magn(vect1)

        type(cart_vect), intent(in) :: vect1

        vect_magn = sqrt((vect1%x)**2 + (vect1%y)**2 + (vect1%z)**2)

    end function vect_magn

    ! test if is zero vector

    logical function almost_zero_vect(vect1)

        type(cart_vect), intent(in) :: vect1

        almost_zero_vect = vect_magn(vect1) < vect_magn_min_thresh

    end function almost_zero_vect

    ! test if a cart_vect has magnitude = 1

    logical function almost_unit_vect(vect1)

        type(cart_vect), intent(in) :: vect1

        almost_unit_vect = dabs(1.0 - vect_magn(vect1)) < vect_magn_min_thresh

    end function almost_unit_vect

    ! calculates the product of a cartesian vector by a scalar

    type (cart_vect) function vect_by_scal(vect1, scalar1)

        type(cart_vect), intent(in) :: vect1
        real (kind=r8b), intent(in) :: scalar1

        vect_by_scal%x = scalar1 * vect1%x
        vect_by_scal%y = scalar1 * vect1%y
        vect_by_scal%z = scalar1 * vect1%z

    end function vect_by_scal

    ! divides a cartesian vector by a scalar

    type (cart_vect_opt) function vect_scal_div_opt(vect1, scalar1)

        type(cart_vect), intent(in) :: vect1
        real (kind=r8b), intent(in) :: scalar1

        vect_scal_div_opt%valid = scalar1 > vect_div_min_thresh
        if (vect_scal_div_opt%valid) then
            vect_scal_div_opt%vect = vect_by_scal(vect1, 1.0/scalar1)
        endif

    end function

    ! cart_vect normalization

    type(cart_vect_opt) function unit_vect_opt(vect1)

        type(cart_vect), intent(in) :: vect1

        unit_vect_opt = vect_scal_div_opt(vect1, vect_magn(vect1))

    end function unit_vect_opt

    ! calculates the sum of two vectors

    type(cart_vect) function vect_sum(vect1, vect2)

        type(cart_vect), intent(in) :: vect1, vect2

        vect_sum%x = vect1%x + vect2%x
        vect_sum%y = vect1%y + vect2%y
        vect_sum%z = vect1%z + vect2%z

    end function vect_sum

    ! calculates the difference between two vectors

    type(cart_vect) function vect_diff(vect1,vect2)

        type(cart_vect), intent(in) :: vect1, vect2

        vect_diff%x = vect1%x - vect2%x
        vect_diff%y = vect1%y - vect2%y
        vect_diff%z = vect1%z - vect2%z

    end function vect_diff

    ! angle (in radians) between two vectors (radiants, 0-pi)

    type(real_opt) function vect_angl_rad_opt(vect1, vect2) result(angle_val)

        type(cart_vect), intent(in) :: vect1, vect2
        real (kind=r8b) :: vect1_magn, vect2_magn, scaled_scal_prod

        vect1_magn = vect_magn(vect1)
        vect2_magn = vect_magn(vect2)
        angle_val%valid = vect1_magn * vect2_magn >= vect_magn_min_thresh
        if (angle_val%valid) then
            scaled_scal_prod = scal_prod(vect1, vect2)/(vect1_magn*vect2_magn)
            if (scaled_scal_prod < -1.) then
                angle_val%real_val = pi
            else if (scaled_scal_prod > 1.) then
                angle_val%real_val = 0.0
            else
                angle_val%real_val = acos(scaled_scal_prod)
            end if
        end if

    end  function

    ! angle (in radians) between two axes (radiants, 0-pi/2)

    type(real_opt) function axes_angle_rad_opt(vect1, vect2)

        type(cart_vect), intent(in) :: vect1, vect2
        type(real_opt) :: vect_angle_rad

        ! angle between vectors (in radians)
        vect_angle_rad = vect_angl_rad_opt(vect1, vect2)
        axes_angle_rad_opt%valid = vect_angle_rad%valid
        if (axes_angle_rad_opt%valid) then
            axes_angle_rad_opt%real_val = min(vect_angle_rad%real_val, pi - vect_angle_rad%real_val)
        end if

    end function

    ! scalar product of two vectors (given as their cartesian coordinates)

    real (kind=r8b) function scal_prod(vect1, vect2)

        type(cart_vect), intent(in) :: vect1, vect2

        scal_prod = vect1%x * vect2%x + vect1%y * vect2%y + vect1%z * vect2%z

    end function scal_prod

    ! vectorial product of two vectors (given as their cartesian coordinates)

    type(cart_vect) function vect_prod(vect1,vect2)

        type(cart_vect), intent(in) :: vect1, vect2

        vect_prod%x = (vect1%y * vect2%z)-(vect1%z * vect2%y)
        vect_prod%y = (vect1%z * vect2%x)-(vect1%x * vect2%z)
        vect_prod%z = (vect1%x * vect2%y)-(vect1%y * vect2%x)

    end function vect_prod

    ! vect1 projection on vect2

    type (cart_vect_opt) function vect_project_opt(vect1, vect2)

        type(cart_vect), intent(in) :: vect1, vect2
        real (kind=r8b) :: scalprod
        type(cart_vect_opt) :: vect_2_opt

        vect_2_opt = unit_vect_opt(vect2)
        vect_project_opt%valid = vect_2_opt%valid
        if (vect_project_opt%valid) then
            scalprod = scal_prod(vect1, vect_2_opt%vect)
            vect_project_opt%vect = vect_by_scal(vect2, scalprod)
        end if

    end function

    ! converts a 3D cart_vect to a 3x1 array

    function vect2arr(vect1) result(array1)

        type(cart_vect), intent(in) :: vect1
        real(kind=r8b) :: array1(3)

        array1 = [vect1%x, vect1%y, vect1%z]

    end function vect2arr

    ! converts a 3x1 array to a 3D cart_vect

    function arr2vect(array1) result(vect1)

        real(kind=r8b), intent(in) :: array1(3)
        type(cart_vect):: vect1

        vect1%x = array1(1)
        vect1%y = array1(2)
        vect1%z = array1(3)

    end function arr2vect


end module vector_processing


module geologic_processing

    use vector_processing

    implicit none

    type :: geol_axis
        real (kind=r8b) :: trend, plunge  ! as degrees
    end type geol_axis

    type :: geol_axis_opt
        type(geol_axis) :: geoax
        logical :: valid
    end type

    type :: geol_plane
        real (kind=r8b) :: strike_rhr, dip_dir, dip_angle ! as degrees
    end type geol_plane

    type :: geol_plane_opt
        type(geol_plane) :: geolplane
        logical :: valid
    end type

    contains

    ! calculates the polar components from the cartesian ones

    type(cart_vect) function pole2cartes(axis1) result(vect1)

        type(geol_axis), intent(in) :: axis1

        vect1%x = cos(degr2rad*axis1%plunge) * cos(degr2rad*axis1%trend)
        vect1%y = cos(degr2rad*axis1%plunge) * sin(degr2rad*axis1%trend)
        vect1%z = sin(degr2rad*axis1%plunge)

    end function pole2cartes

    ! calculates polar components from cartesian ones

    type(geol_axis_opt) function cartes2pole_opt(vect1) result(axis1_opt)

        type(cart_vect), intent(in) :: vect1
        type(cart_vect_opt) :: vect_unit_opt
        type(cart_vect) :: vect_unit
        type(geol_axis) :: geoaxis

        vect_unit_opt = unit_vect_opt(vect1)
        axis1_opt%valid = vect_unit_opt%valid

        if (axis1_opt%valid) then ! polar coordinates calculation
            vect_unit = vect_unit_opt%vect
            if (vect_unit%z > 1.0_r8b) then
                vect_unit%z = 1.0_r8b
            elseif (vect_unit%z < -1.0_r8b) then
                vect_unit%z = -1.0_r8b
            endif
            geoaxis%plunge = rad2degr*dasin(vect_unit%z)
            if (dabs(geoaxis%plunge)> 89.9_r8b) then
                geoaxis%trend = 0.0_r8b
            else
                geoaxis%trend = modulo(rad2degr*atan2(vect_unit%x,vect_unit%y), 360.0_r8b)
            endif
            axis1_opt%geoax = geoaxis
        endif

    end function cartes2pole_opt

    ! calculates the down geol_axis from the geol_axis

    type(geol_axis) function axis2downaxis(axis1) result(axis2)

        type(geol_axis), intent(in) :: axis1

        if (axis1%plunge < 0.0) then
          axis2%plunge = -axis1%plunge
          axis2%trend = modulo(axis1%trend + 180.0_r8b, 360.0_r8b)
        else
          axis2 = axis1
        endif

    end function axis2downaxis

    ! calculates dip direction of a fault

    subroutine dipdir_calc(geolplane1)

        type(geol_plane), intent(inout) :: geolplane1

        geolplane1%dip_dir = modulo(geolplane1%strike_rhr + 90.0_r8b, 360.0_r8b)

    end subroutine dipdir_calc

    ! calculates the geological plane normal to a dow-axis (downward-oriented axis)

    type(geol_plane) function downaxis2geolplane(downaxis1) result(geolplane1)

        type(geol_axis), intent(in) :: downaxis1

        geolplane1%strike_rhr = modulo(downaxis1%trend + 90.0_r8b, 360.0_r8b)
        call dipdir_calc(geolplane1)
        geolplane1%dip_angle = 90.0_r8b - downaxis1%plunge

    end function downaxis2geolplane

     ! calculates geological plane normal

    type(cart_vect) function geolplane2normalvect(geoplane1) result(geolplanenorm1)

        type(geol_plane), intent(in) :: geoplane1

        ! formulas from Aki and Richards, 1980
        geolplanenorm1%x = -sin(degr2rad * geoplane1%dip_angle) * sin(degr2rad * geoplane1%strike_rhr)
        geolplanenorm1%y =  sin(degr2rad * geoplane1%dip_angle) * cos(degr2rad * geoplane1%strike_rhr)
        geolplanenorm1%z = -cos(degr2rad * geoplane1%dip_angle)

    end function

    ! calculates the geological plane given its normal

    type(geol_plane_opt) function normalvect2geolplane(normvect1) result(geolplaneopt1)

        !implicit none

        type(cart_vect), intent(in) :: normvect1
        type(geol_axis_opt) :: geolaxisopt1

        geolaxisopt1 = cartes2pole_opt(normvect1)
        geolplaneopt1%valid = geolaxisopt1%valid
        if (geolplaneopt1%valid) then
            geolplaneopt1%geolplane = downaxis2geolplane(axis2downaxis(geolaxisopt1%geoax))
        end if

    end function


end module geologic_processing


module bestfitplane

    use vector_processing

    implicit none

    contains

    ! calculates SVD solution (s, sigma, u, v) for a given array a with m rows and n columns

    type (cart_vect_opt) function svd(m, pts) result(res_vect)
     ! modified from: TEST07 tests DGESVD
     !                 by John Burkardt

        implicit none

        integer (kind = i4b), intent(in) :: m
        integer (kind= i4b), parameter :: n = 3
        real (kind = r8b), intent(in) :: pts(m, n)

        real (kind = r8b) :: s(min(m, n))
        real (kind = r8b) :: u(m, m)

        ! internal variables

        character :: jobu
        character :: jobvt
        integer ( kind = i4b ) :: lwork
        integer ( kind = i4b ) :: info
        integer ( kind = i4b ) :: lda
        integer ( kind = i4b ) :: ldu
        integer ( kind = i4b ) :: ldvt
        real ( kind = r8b ) :: vt(n,n)
        real ( kind = r8b ), allocatable :: work(:)

        !  Compute the singular values and singular vectors

        lwork = 3*min(m,n) + max ( max(m,n), 2*min(m,n) )
        allocate (work(lwork))

        jobu = 'A'
        jobvt = 'A'
        lda = m
        ldu = m
        ldvt = n

        call dgesvd ( jobu, jobvt, m, n, pts, lda, s, u, ldu, vt, ldvt, work, &
                      lwork, info )

        if (info /= 0) then
            res_vect%valid = .FALSE.
        else
             res_vect%valid = .TRUE.
             res_vect%vect = cart_vect(vt(3,1),vt(3,2),vt(3,3))
        end if

     end function

end module bestfitplane


module io

    implicit none

    integer :: ios ! status of input/output connection

    contains

    ! write program header on screen

    subroutine welcome_screen()


        write (*,"(A)") '*****************************************'
        write (*,*)
        write (*,"(16X,A)") 'geoSurfDEM - Invert'
        write (*,*)
        write (*,"(16X,A)") 'vers. 0.0.1'
        write (*,"(17X,A)") '2016-11-31'
        write (*,*)
        write (*,"(5X,A)") 'Fortran program for determining Best-Fit Geological Planes to surface traces'
        write (*,*)
        write (*,"(A)") '*****************************************'

        write (*,*)
        write (*,*)
        write (*,"(4X,A)") 'Input parameters definition'
        write (*,*)

    end subroutine welcome_screen

    ! choice and opening of input file, that stores the list of data sources

    subroutine input_file_def()


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

    end subroutine input_file_def

    ! define output files storing results

    subroutine output_file_def()


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

    end subroutine output_file_def


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
    !    Input, integer N, the number of components of the cart_vect.
    !
    !    Input, real ( kind = r8b ) A(N), the cart_vect to be printed.
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


end module io


program invert_attitudes

    use var_types
    use geologic_processing
    use bestfitplane
    use io

    implicit none

    integer ( kind = i4b ), parameter :: n = 3
    real ( kind = r8b ) :: x, y, z
    integer ( kind = i4b ) :: num_points, i
    real ( kind = r8b ), dimension(:, :), allocatable :: points

    type(cart_vect_opt) :: svd_res
    type(geol_plane_opt) :: geolplaneopt_res
    type(geol_plane) :: geoplane_res

    ! printout of welcome screen

    call welcome_screen()

    ! user-definition of source file

    call input_file_def()

    ! user-definition of output file

    call output_file_def()

    ! defines the number of input points

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

    ! allocate arrays

    allocate(points(num_points, n))

    ! read source point data for single solution calculation

    rewind 15
    do i = 1, num_points
        read (15, *, iostat=ios) points(i, 1), points(i, 2), points(i, 3)
    end do

    ! write source point data in output file

    write (16, *) 'point source data'
    do i = 1, num_points
        write (16, *) points(i, 1), points(i, 2), points(i, 3)
    end do
    write (16, *) ''

    ! calculates the SVD solution and the geological plane

    geolplaneopt_res%valid = .false.
    svd_res = svd(num_points, points)
    if (svd_res%valid) then
        write (16, *) svd_res%vect%x, svd_res%vect%y, svd_res%vect%z
        geolplaneopt_res = normalvect2geolplane(svd_res%vect)
    end if

    ! output results

    if (.not. geolplaneopt_res%valid) then
        write (16, *) 'invalid result'
    else
        geoplane_res = geolplaneopt_res%geolplane

        write (16, *) geoplane_res%strike_rhr, geoplane_res%dip_dir, geoplane_res%dip_angle
    end if

    ! write end message

    write ( *, '(a)' ) 'Program completed'


end program
