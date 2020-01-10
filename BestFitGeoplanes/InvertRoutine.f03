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
            geoaxis%plunge = -rad2degr*dasin(vect_unit%z)
            if (dabs(geoaxis%plunge)> 89.9_r8b) then
                geoaxis%trend = 0.0_r8b
            elseif (dabs(geoaxis%plunge) < -89.9_r8b) then
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

    ! calculates the geological plane normal to a down-axis (downward-oriented axis)

    type(geol_plane) function downaxis2geolplane(downaxis1) result(geolplane1)

        type(geol_axis), intent(in) :: downaxis1

        geolplane1%strike_rhr = modulo(downaxis1%trend + 90.0_r8b, 360.0_r8b)
        call dipdir_calc(geolplane1)
        geolplane1%dip_angle = 90.0_r8b - downaxis1%plunge

    end function downaxis2geolplane

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

        call dgesvd (jobu, jobvt, m, n, pts, lda, s, u, ldu, vt, ldvt, work, &
                    lwork, info)

        if (info /= 0) then
            res_vect%valid = .false.
        else
             res_vect%valid = .true.
             res_vect%vect = cart_vect(vt(3,1), vt(3,2), vt(3,3))
        end if

     end function


    function array_minus_mean(num_points, pts_array) result(array_mm)

        use var_types
        implicit none

        integer num_points
        real (kind = r8b), dimension(3, num_points) :: pts_array, array_mm
        real (kind = r8b), dimension(3) :: sums

        sums = SUM(pts_array, DIM=2)
        array_mm(1, 1:num_points) = pts_array(1, 1:num_points) - (sums(1) / num_points)
        array_mm(2, 1:num_points) = pts_array(2, 1:num_points) - (sums(2) / num_points)
        array_mm(3, 1:num_points) = pts_array(3, 1:num_points) - (sums(3) / num_points)

    end function


end module bestfitplane


subroutine invert_attitudes(num_points, pts_array_cform, success, dipdir, dipang) bind(c)

    use iso_c_binding

    use var_types
    use geologic_processing
    use bestfitplane

    implicit none

    ! C-connection variables

    integer(c_int) :: num_points
    real(c_double) :: pts_array_cform(3, num_points)
    logical(c_bool) :: success
    real(c_double) :: dipdir, dipang

    ! Fortran internal variables

    real (kind = r8b), dimension(3, num_points) :: points_less_mean
    real (kind = r8b), dimension(num_points, 3) :: points
    type(cart_vect_opt) :: svd_res, unvect_opt
    type(geol_plane_opt) :: geolplaneopt_res
    integer :: i


    ! define point array (Fortran format)

    points_less_mean = array_minus_mean(num_points, pts_array_cform)

    points = transpose(points_less_mean)

    ! calculates the SVD solution and the geological plane

    geolplaneopt_res%valid = .false.
    svd_res = svd(num_points, points)
    if (svd_res%valid) then
        unvect_opt = unit_vect_opt(svd_res%vect)
        if (unvect_opt%valid) then
            geolplaneopt_res = normalvect2geolplane(unvect_opt%vect)
        end if
    end if

    ! define results

    if (.not. geolplaneopt_res%valid) then
        success = .false.
        dipdir = -99.9
        dipang = -99.9
    else
        success = .true.
        dipdir = geolplaneopt_res%geolplane%dip_dir
        dipang = geolplaneopt_res%geolplane%dip_angle
    end if

end subroutine

