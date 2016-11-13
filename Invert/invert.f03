! Invert intersection points
! to derive a best-fit-geol_plane

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


module data_types

    use kinds

    ! variable declarations

    implicit none

    type :: cartes_vect
        sequence
        real (kind=r8b)	:: x, y, z
    end type cartes_vect

    type :: geol_axis
        real (kind=r8b)	:: trend, plunge
    end type geol_axis

    type :: orthonorm_triad
        type(cartes_vect) :: X, Y, Z
    end type orthonorm_triad

    type :: geol_plane
        real (kind=r8b)	:: strike_rhr, dip_dir, dip_angle
    end type geol_plane

    type :: fault_datum
        type (geol_plane) :: flt_pln
        type (geol_axis) :: slickln_ax
        real (kind=r8b)	:: slickln_rake_aki
        type(cartes_vect) :: flt_norm_vctr, slickenl_vctr
        real(kind=r8b) :: x, y, z
    end type fault_datum

    type :: stress_prcomp
        type (geol_axis) :: S1_ax, S2_ax, S3_ax
        type (cartes_vect) :: S1_vctr, S2_vctr, S3_vctr
        real (kind=r8b) :: sigma_1, sigma_2, sigma_3, phi
        real (kind=r8b) :: fluid_press, cohes_streng, int_frict_res
    end type stress_prcomp

    type :: stress_solution
        type(cartes_vect) :: tract_vctr, norm_stress_vctr, shear_stress_vctr
        real (kind=r8b)	:: theor_rake, tract_stress_magn, norm_stress_magn, shear_stress_magn
        real (kind=r8b)	:: mod_slip_tend, deform_index
        real (kind=r8b) :: dilat_tend, fract_stab, leakage_factor
        type (geol_axis) :: theor_slicknln
    end type stress_solution

    integer :: ios !status of input/output connection
    real (kind=r8b) :: pi  	! pi radians
    real (kind=r8b) :: rad2degr, degr2rad	! for conversion from radians to degrees

    type(orthonorm_triad), parameter :: frame0 = orthonorm_triad(cartes_vect(1.0,0.0,0.0), &
                                                                 cartes_vect(0.0,1.0,0.0), &
                                                                 cartes_vect(0.0,0.0,1.0))
    real (kind=r8b),parameter :: unit_vect_tol = 1.0e-6 !maximum accepted value for unitary cartes_vect magnitude difference with unit
    real (kind=r8b),parameter :: vect_norm_tol = 1.0e-5 !minimum accepted value for cartes_vect magnitude to apply normalization
    real (kind=r8b), parameter :: shear_magn_thresh = 1.0e-5 !minimum accepted value for shear stress to be considered meaningful

end module data_types


module vector_processing

    use data_types

    contains

    ! calculates the sum of two vectors

    type(cartes_vect) function vect_sum(vect1, vect2)

        type(cartes_vect), intent(in) :: vect1, vect2

        vect_sum%x = vect1%x + vect2%x
        vect_sum%y = vect1%y + vect2%y
        vect_sum%z = vect1%z + vect2%z

    end function vect_sum

    ! calculates the difference between two vectors

    type(cartes_vect) function vect_diff(vect1,vect2)

        type(cartes_vect), intent(in) :: vect1, vect2

        vect_diff%x = vect1%x - vect2%x
        vect_diff%y = vect1%y - vect2%y
        vect_diff%z = vect1%z - vect2%z

    end function vect_diff

    ! angle (in radians) between two vectors (radiants, 0-pi)

    real (kind=r8b) function vect_angl_rad(vect1,vect2)

        type(cartes_vect), intent(in) :: vect1, vect2

        vector1_magn = vect_magn(vect1)
        vector2_magn = vect_magn(vect2)
        if ((vector1_magn < vect_norm_tol).or.(vector2_magn < vect_norm_tol)) then
            write(*,*) 'Error in cartes_vect magnitude (function vect_angl_rad). Hit any key to stop'
            read(*,*)
            stop	! STOP PROGRAM FOR ERROR IN DATA INPUT
        end if

        ! scalar product between two vectors
        scaled_scal_prod = scal_prod(vect1, vect2)/(vector1_magn*vector2_magn)

        ! angle between vectors (in radians)
        if (scaled_scal_prod < -1.) then
            vect_angl_rad = pi
        else if (scaled_scal_prod > 1.) then
            vect_angl_rad = 0.0
        else
            vect_angl_rad = acos(scaled_scal_prod)
        end if

    end  function vect_angl_rad

    ! angle (in radians) between two axes (radiants, 0-pi/2)

    real (kind=r8b) function axes_angle_rad(vect1, vect2)

        type(cartes_vect), intent(in) :: vect1, vect2

        ! angle between vectors (in radians)
        axes_angle_rad = vect_angl_rad(vect1,vect2)
        axes_angle_rad = min(axes_angle_rad, pi-axes_angle_rad)

    end function axes_angle_rad

    ! cartes_vect normalization

    type(cartes_vect) function vect_normaliz(vect1) result(vect2)

        type(cartes_vect), intent(in) :: vect1

        vector1_magn = vect_magn(vect1)

        if (vector1_magn < vect_norm_tol) then
          write(*,*) 'Error in cartes_vect magnitude processing. Hit any key to stop'
          read (*,*)
          stop ! STOP PROGRAM FOR ERROR IN DATA INPUT
        end if

        vect2%x = vect1%x / vector1_magn
        vect2%y = vect1%y / vector1_magn
        vect2%z = vect1%z / vector1_magn


    end function vect_normaliz

    ! cartes_vect magnitude

    real (kind=r8b) function vect_magn(vect1)

        type(cartes_vect), intent(in) :: vect1

        vect_magn = sqrt((vect1%x)**2 + (vect1%y)**2 + (vect1%z)**2)

    end function vect_magn

    ! calculates the product of a cartes_vect by a scalar

    type (cartes_vect) function vect_by_scal(vect1, scalar1)

        type(cartes_vect), intent(in) :: vect1
        real (kind=r8b), intent(in) :: scalar1

        vect_by_scal%x = scalar1 * vect1%x
        vect_by_scal%y = scalar1 * vect1%y
        vect_by_scal%z = scalar1 * vect1%z

    end function vect_by_scal

    ! scalar product of two vectors (given as their cartesian coordinates)

    real (kind=r8b) function scal_prod(vect1, vect2)

        type(cartes_vect), intent(in) :: vect1, vect2

        scal_prod = vect1%x * vect2%x + vect1%y * vect2%y + vect1%z * vect2%z

    end function scal_prod

    ! vectorial product of two vectors (given as their cartesian coordinates)

    type(cartes_vect) function vect_prod(vect1,vect2)

        type(cartes_vect), intent(in) :: vect1, vect2

        vect_prod%x = (vect1%y * vect2%z)-(vect1%z * vect2%y)
        vect_prod%y = (vect1%z * vect2%x)-(vect1%x * vect2%z)
        vect_prod%z = (vect1%x * vect2%y)-(vect1%y * vect2%x)

    end function vect_prod

    ! vect1 projection on vect2

    type (cartes_vect) function vect_project(vect1,vect2)

        type(cartes_vect), intent(in) :: vect1, vect2
        real (kind=r8b) :: scalprod

        scalprod = scal_prod(vect1, vect_normaliz(vect2))
        vect_project = vect_by_scal(vect2, scalprod)

    end function vect_project

    ! test if a cartes_vect has magnitude = 1

    logical function is_unit_vect(vect1) result(vect_is_normalized)

        type(cartes_vect), intent(in) :: vect1
        real (kind=r8b) :: vector1_magn

        vector1_magn = vect_magn(vect1)

        abs_diff = dabs(1-vector1_magn)

        if (abs_diff > unit_vect_tol) then
            vect_is_normalized = .false.
        else
            vect_is_normalized = .true.
        endif

    end function is_unit_vect

    ! converts a 3D cartes_vect to a 3x1 array

    function vect2arr(vect1) result(array1)

        type(cartes_vect), intent(in) :: vect1
        real(kind=r8b) :: array1(3)

        array1 = (/vect1%x, vect1%y, vect1%z/)

    end function vect2arr

    ! converts a 3x1 array to a 3D cartes_vect

    function arr2vect(array1) result(vect1)

        real(kind=r8b), intent(in) :: array1(3)
        type(cartes_vect):: vect1

        vect1%x = array1(1)
        vect1%y = array1(2)
        vect1%z = array1(3)

    end function arr2vect


end module vector_processing


module geometric_processing

    use vector_processing

    contains

    ! calculates the polar components from the cartesian ones

    type(cartes_vect) function pole2cartes(axis1) result(vect1)

        type(geol_axis), intent(in) :: axis1


        vect1%x = cos(degr2rad*axis1%plunge) * cos(degr2rad*axis1%trend)
        vect1%y  = cos(degr2rad*axis1%plunge) * sin(degr2rad*axis1%trend)
        vect1%z  =  sin(degr2rad*axis1%plunge)


    end function pole2cartes

    ! calculates polar components from cartesian ones

    type(geol_axis) function cartes2pole(vect1) result(axis1)

        type(cartes_vect):: vect1
        logical :: vect_is_normalized

        vect_is_normalized = is_unit_vect(vect1)

        if (.not.vect_is_normalized) then
          vect1 = vect_normaliz(vect1)
        endif

        ! polar coordinates calculation

        if (vect1%z > 1.0_r8b) then
            vect1%z = 1.0_r8b
        elseif (vect1%z < -1.0_r8b) then
            vect1%z = -1.0_r8b
        endif

        axis1%plunge = rad2degr*dasin(vect1%z)

        if (dabs(axis1%plunge)> 89.5) then
                axis1%trend = 0.0
        else
            axis1%trend = rad2degr*atan2(vect1%y,vect1%x)
            if (axis1%trend < 0.0) then
                axis1%trend = 360.0 + axis1%trend
            endif
        endif

    end function cartes2pole

    ! calculates the down geol_axis from the geol_axis

    type(geol_axis) function axis2downaxis(axis1) result(axis2)

        type(geol_axis), intent(in) :: axis1

        if (axis1%plunge < 0.0) then
          axis2%plunge = -axis1%plunge
          axis2%trend = axis1%trend + 180.0
          if (axis2%trend >= 360.0) then
            axis2%trend = axis2%trend - 360.0
          endif
        else
          axis2%trend = axis1%trend
          axis2%plunge = axis1%plunge
        endif

    end function axis2downaxis


end module geometric_processing


module fault_processing

    use geometric_processing

    implicit none

    contains

    ! calculates dip direction of a fault

    subroutine dipdir_calc(faultplane1)

        type(geol_plane), intent(inout) :: faultplane1

        faultplane1%dip_dir = faultplane1%strike_rhr + 90.0
        if (faultplane1%dip_dir >= 360.0) then
          faultplane1%dip_dir = faultplane1%dip_dir - 360.0
        endif

    end subroutine dipdir_calc

    ! calculates fault normal

    type(cartes_vect) function faultplanenorm_calc(faultplane1) result(faultnorm)

        type(geol_plane), intent(in) :: faultplane1

        ! Fault Normal cartesian coordinates
        ! formulas from Aki and Richards, 1980
        faultnorm%x = -sin(degr2rad*faultplane1%dip_angle) * sin(degr2rad*faultplane1%strike_rhr)
        faultnorm%y = sin(degr2rad*faultplane1%dip_angle) * cos(degr2rad*faultplane1%strike_rhr)
        faultnorm%z = -cos(degr2rad*faultplane1%dip_angle)

    end function faultplanenorm_calc

    ! calculates cartes_vect (cartesian) components of fault record

    subroutine fault_cartcmp(faultplane1, rake1, flt_norm_vctr, slickenl_vctr)

        ! QUATf_faultpole2faultvector

        type(geol_plane), intent(in) :: faultplane1
        real (kind=r8b)	:: rake1
        type (cartes_vect) , intent(out) :: flt_norm_vctr, slickenl_vctr

        ! Fault Normal and slickln_ax cartesian coordinates
        ! formulas from Aki and Richards, 1980
        flt_norm_vctr%x = -sin(degr2rad*faultplane1%dip_angle)*sin(degr2rad*faultplane1%strike_rhr)
        flt_norm_vctr%y = sin(degr2rad*faultplane1%dip_angle)*cos(degr2rad*faultplane1%strike_rhr)
        flt_norm_vctr%z = -cos(degr2rad*faultplane1%dip_angle)

        slickenl_vctr%x = cos(degr2rad*rake1) * cos(degr2rad*faultplane1%strike_rhr)   &
         + sin(degr2rad*rake1) * cos(degr2rad*faultplane1%dip_angle) * sin(degr2rad*faultplane1%strike_rhr)
        slickenl_vctr%y = cos(degr2rad*rake1) * sin(degr2rad*faultplane1%strike_rhr)  &
         - sin(degr2rad*rake1) * cos(degr2rad*faultplane1%dip_angle) * cos(degr2rad*faultplane1%strike_rhr)
        slickenl_vctr%z = -sin(degr2rad*rake1) * sin(degr2rad*faultplane1%dip_angle)

    end subroutine fault_cartcmp

    ! calculate the slickln_ax (trend and plunge) from the rake angle

    type(geol_axis) function rake2slickenline(strike, dip, rake) result(slicken_pole)

        implicit none

        real (kind=r8b), intent(in) :: strike, dip, rake
        real (kind=r8b) :: strike_rd, dip_rd, rake_rd
        type (cartes_vect) :: slick_vect

        strike_rd = degr2rad*strike
        dip_rd = degr2rad*dip
        rake_rd = degr2rad*rake

        slick_vect%x = cos(rake_rd) * cos(strike_rd) + sin(rake_rd) * sin(strike_rd) * cos(dip_rd)
        slick_vect%y = cos(rake_rd) * sin(strike_rd) - sin(rake_rd) * cos(strike_rd) * cos(dip_rd)
        slick_vect%z = -sin(rake_rd) * sin(dip_rd)

        ! determination of trend and plunge of slickln_ax
        slicken_pole = cartes2pole(slick_vect)

    end function rake2slickenline


end module fault_processing


module stress_processing

    use fault_processing

    implicit none

    contains

    ! calculation of S2_ax as cartes_vect product of S3_ax and S1_ax

    subroutine S2_calc(stress_pc_1)

        type(stress_prcomp), intent(inout) :: stress_pc_1

        stress_pc_1%S2_vctr = vect_prod(stress_pc_1%S3_vctr,stress_pc_1%S1_vctr)
        stress_pc_1%S2_ax = cartes2pole(stress_pc_1%S2_vctr)
        stress_pc_1%S2_ax = axis2downaxis(stress_pc_1%S2_ax)

    end subroutine S2_calc

    ! calculation of sigma_2 from S1_ax, S3_ax and PHI values

    subroutine sigma2_calc(stress_pc_1)

        type(stress_prcomp), intent(inout) :: stress_pc_1

        stress_pc_1%sigma_2 = stress_pc_1%phi*stress_pc_1%sigma_1 				&
                                    + (1-stress_pc_1%phi)*stress_pc_1%sigma_3


    end subroutine sigma2_calc

    ! define rotation matrix based on stress principal axes

    function rotmatr(stress_pc_1) result(rot_matrix)

        ! based on Kuipers, 2002, p.161, eqs. 7.8

        type(stress_prcomp), intent(in) :: stress_pc_1
        real (kind=r8b) :: rot_matrix(3,3)

        rot_matrix(1,1) = scal_prod(stress_pc_1%S1_vctr, frame0%X)
        rot_matrix(1,2) = scal_prod(stress_pc_1%S2_vctr, frame0%X)
        rot_matrix(1,3) = scal_prod(stress_pc_1%S3_vctr, frame0%X)

        rot_matrix(2,1) = scal_prod(stress_pc_1%S1_vctr, frame0%Y)
        rot_matrix(2,2) = scal_prod(stress_pc_1%S2_vctr, frame0%Y)
        rot_matrix(2,3) = scal_prod(stress_pc_1%S3_vctr, frame0%Y)

        rot_matrix(3,1) = scal_prod(stress_pc_1%S1_vctr, frame0%Z)
        rot_matrix(3,2) = scal_prod(stress_pc_1%S2_vctr, frame0%Z)
        rot_matrix(3,3) = scal_prod(stress_pc_1%S3_vctr, frame0%Z)


    end function rotmatr

    ! calculate stress tensor expressed in frame components

    function stresstensorcalc(stress_pc_1,rot_matrix) result(tens)

        ! from Kagan and Knopoff, 1985a, p. 433

        type(stress_prcomp), intent(in) :: stress_pc_1
        real (kind=r8b), intent(in)  :: rot_matrix(3,3)

        real (kind=r8b) :: stresstens0(3,3), tens(3,3)

        ! stress eigentensor
        stresstens0 = reshape((/0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/),shape=(/3,3/),order=(/2,1/))
        stresstens0(1,1) = stress_pc_1%sigma_1
        stresstens0(2,2) = stress_pc_1%sigma_2
        stresstens0(3,3) = stress_pc_1%sigma_3

        tens = matmul(rot_matrix,matmul(stresstens0,transpose(rot_matrix)))

    end function stresstensorcalc

    ! calculation of stress solution for the given geol_plane and stress tensor

    type(stress_solution) function stresssol_calc(tens1,stress_pc_1,faultrec1) result(stress_sol_1)

        ! based on Xu, 2004 (Geoph. J. Int., 157,1316-1330) and references within.

        ! QUATstress_stress&fault2rake
        !last modified: 2008/01/05

        ! INPUT
        !    SELF with 2 objects:
        !      (0): fault geol_plane orientation: strike rhr, dip direction, dip angle
        !      (1): stress field tensor components: x,y,z (sigma 1,2,3)

        !  OUTPUT
        !    normal and tangential stresses magnitudes
        !    slip tendency
        !    theoretical rake angle of tangential stress

        !  STEPS
        !    input of data:
        !      fault geol_plane orientation: strike rhr, dip direction, dip angle     stress field tensor components: x,y,z (sigma 1,2,3)
        !    fault normal calculation: n (1,2,3)
        !    traction cartes_vect c.: sigma (1,2,3)
        !    normal stress c.: sigma n (1,2,3)
        !    shear stress c: tau s (1,2,3)
        !    slip tendency c.: Ts
        !    theoretical rake angle c: lambda
        !    output of results:
        !      magnitudes of normal and shear stresses,
        !      slip tendency,
        !      theoretical rake angle

        !    input of data:
        !      fault geol_plane orientation: strike rhr, dip direction, dip angle
        !      stress field tensor components: x,y,z (sigma 1,2,3)

        implicit none

        real (kind=r8b), intent(in)  :: tens1(3,3)
        type(stress_prcomp), intent(in) :: stress_pc_1
        type(fault_datum), intent(in)  :: faultrec1

        type(cartes_vect) :: shearstress_unitvect
        type(geol_axis) :: strike_pole, dipdir_pole
        type(cartes_vect) :: strike_vect, dipdir_vect, flt_norm_vctr
        real (kind=r8b) :: faultnormal_array(3), scalprod_shearstress_strike, scalprod_shearstress_dipdir
        real (kind=r8b) :: tau

        ! calculation of fault normal calculation: n (1,2,3)
        flt_norm_vctr = faultplanenorm_calc(faultrec1%flt_pln)
        faultnormal_array = vect2arr(flt_norm_vctr)

        !    traction cartes_vect c.: sigma (1,2,3)
        stress_sol_1%tract_vctr = arr2vect(- matmul(tens1,faultnormal_array))
        stress_sol_1%tract_stress_magn = vect_magn(stress_sol_1%tract_vctr)
        if (scal_prod(stress_sol_1%tract_vctr, flt_norm_vctr)<0.0_r8b) then
            stress_sol_1%tract_stress_magn = -stress_sol_1%tract_stress_magn
        endif

        !    normal stress c.: SigmaN
        stress_sol_1%norm_stress_vctr = vect_project(stress_sol_1%tract_vctr, flt_norm_vctr)
        stress_sol_1%norm_stress_magn = vect_magn(stress_sol_1%norm_stress_vctr)
        if (scal_prod(stress_sol_1%norm_stress_vctr, flt_norm_vctr)<0.0_r8b) then
            stress_sol_1%norm_stress_magn = -stress_sol_1%norm_stress_magn
        endif

        !    shear stress c: tau s (1,2,3)
        stress_sol_1%shear_stress_vctr = vect_diff(stress_sol_1%tract_vctr,stress_sol_1%norm_stress_vctr)
        stress_sol_1%shear_stress_magn = vect_magn(stress_sol_1%shear_stress_vctr)


        if (stress_sol_1%shear_stress_magn > shear_magn_thresh) then

            shearstress_unitvect = vect_normaliz(stress_sol_1%shear_stress_vctr)

            ! theoretical rake angle c: lambda
            strike_pole = geol_axis(faultrec1%flt_pln%strike_rhr,0)
            dipdir_pole = geol_axis(faultrec1%flt_pln%dip_dir,faultrec1%flt_pln%dip_angle)

            strike_vect = vect_normaliz(pole2cartes(strike_pole))
            dipdir_vect = vect_normaliz(pole2cartes(dipdir_pole))

            scalprod_shearstress_strike = scal_prod(shearstress_unitvect, strike_vect)

            scalprod_shearstress_dipdir = scal_prod(shearstress_unitvect, dipdir_vect)

            if (scalprod_shearstress_strike>1) then
                scalprod_shearstress_strike=1
            elseif (scalprod_shearstress_strike<-1) then
                scalprod_shearstress_strike=-1
            endif

            stress_sol_1%theor_rake = rad2degr*acos(scalprod_shearstress_strike)

            if (scalprod_shearstress_dipdir > 0.0_r8b) then
                stress_sol_1%theor_rake = - stress_sol_1%theor_rake
            endif

            stress_sol_1%theor_slicknln = rake2slickenline(   &
                faultrec1%flt_pln%strike_rhr,faultrec1%flt_pln%dip_angle,stress_sol_1%theor_rake)

            stress_sol_1%mod_slip_tend = stress_sol_1%shear_stress_magn/dabs(stress_sol_1%tract_stress_magn)

            stress_sol_1%deform_index = &
            (dabs(stress_sol_1%tract_stress_magn) - stress_sol_1%shear_stress_magn)/stress_sol_1%tract_stress_magn

            stress_sol_1%dilat_tend = &
            (stress_pc_1%sigma_1 - dabs(stress_sol_1%norm_stress_magn))/(stress_pc_1%sigma_1 - stress_pc_1%sigma_3)

            tau = stress_pc_1%cohes_streng+stress_pc_1%int_frict_res*dabs(stress_sol_1%norm_stress_magn)
            stress_sol_1%fract_stab = &
            (dabs(stress_sol_1%norm_stress_magn) - tau)/ stress_pc_1%int_frict_res

            stress_sol_1%leakage_factor = stress_pc_1%fluid_press/(dabs(stress_sol_1%norm_stress_magn) - tau)

        endif

    end function stresssol_calc


end module stress_processing


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
!    Input, integer N, the number of components of the cartes_vect.
!
!    Input, real ( kind = r8b ) A(N), the cartes_vect to be printed.
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

! calculates SVD solution (s, sigma, u, v) for a given array a with m rows and n columns

subroutine svd(m, n, a, s, sigma, u, v)

    use kinds

    ! modified from: TEST07 tests DGESVD
    !                 by John Burkardt

    integer ( kind = i4b ), intent(in) :: m, n
    real ( kind = r8b ), intent(in) :: a(m, n)

    real ( kind = r8b ), intent(out) :: s(min(m, n))
    real ( kind = r8b ), intent(out) :: sigma(m, n)
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
    real ( kind = r8b ) :: b(m, n)
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

subroutine welcome_screen()

    write (*,"(A)") '*****************************************'
    write (*,*)
    write (*,"(16X,A)") 'geoSurfDEM - Invert'
    write (*,*)
    write (*,"(16X,A)")	'vers. 0.0.1'
    write (*,"(17X,A)") '2016-11-31'
    write (*,*)
    write (*,"(5X,A)") 'Fortran program for determining Best-Fit-geol_plane to surface traces'
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


program invert_attitudes

    use kinds

    implicit none

    integer ( kind = i4b ) :: m
    integer ( kind = i4b ), parameter :: n = 3
    real ( kind = r8b ) :: x, y, z
    real ( kind = r8b ), dimension(:, :), allocatable :: a
    real ( kind = r8b ), dimension(:), allocatable :: s
    real ( kind = r8b ), dimension(:, :), allocatable :: sigma
    real ( kind = r8b ), dimension(:, :), allocatable :: u
    real ( kind = r8b ), dimension(n, n) :: v
    integer ( kind = i4b ) :: num_points, i
    integer ( kind = i2b ) :: ios

    real ( kind = r8b ), dimension(:, :), allocatable :: b ! debug variable

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
    m = num_points

    ! allocate arrays

    allocate(a(m, n), s(min(m, n)), sigma(m, n), u(m, m), b(m, n))

    ! read source point data for single solution calculation

    rewind 15
    do i = 1, m
        read (15, *, iostat=ios) a(i, 1), a(i, 2), a(i, 3)
    end do

    ! write source point data in output file

    write (16, *) 'point source data'
    do i = 1, m
        write (16, *) a(i, 1), a(i, 2), a(i, 3)
    end do
    write (16, *) ''

    ! calculates the SVD solution

    call svd(m, n, a, s, sigma, u, v)

    ! debug

    call r8mat_print ( n, n, v, '  The matrix V:' )

    ! write singular values in output file

    write (16, *) 'Singular values'
    do i = 1, m
        write (16, *) s(i)
    end do
    write (16, *) ''

    ! write left singular vectors U in output file

    write (16, *) 'Left singular vectors U'
    do i = 1, m
        write (16, *) u(i,1), u(i,2), u(i,3)
    end do
    write (16, *) ''

    ! calculates resulting product cartes_vect

    b = matmul ( u, matmul ( sigma, transpose(v) ) )

    ! writes resulting product cartes_vect in output file

    write (16, *) 'The product U * S * Vt'
    do i = 1, m
        write (16, *) b(i,1), b(i,2), b(i,3)
    end do
    write (16, *) ''

    ! writes V array solution in output file

    write (16, *) 'V array solution'
    do i = 1, 3
        write (16, *) v(i, 1), v(i, 2), v(i, 3)
    end do
    write (16, *) ''

    ! write end message

    write ( *, '(a)' ) 'Program completed'


end program
