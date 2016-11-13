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


module data_types

    use kinds

    ! variable declarations

    implicit none

    type :: range_1D
        real (kind=r8b) :: min, max
    end type range_1D

    type :: range_3D
        type (range_1D) :: r1D(3)
    end type range_3D

    type :: vector
        sequence
        real (kind=r8b)	:: x, y, z
    end type vector

    type :: axis
        real (kind=r8b)	:: trend, plunge
    end type axis

    type :: orthonormal_triad
        type(vector) :: X, Y, Z
    end type orthonormal_triad

    type :: plane
        real (kind=r8b)	:: strike_rhr, dipdirection, dipangle
    end type plane

    type :: fault_datum
        integer (kind=i4b) :: id
        type (plane) :: fltplane
        type (axis) :: slickenline
        real (kind=r8b)	:: rake_aki
        type(vector) :: faultnorm_vect, slickenl_vect
        real(kind=r8b) :: x, y, z, t
    end type fault_datum

    type :: stress_princcompon
        type (axis) :: S1_pole, S2_pole, S3_pole
        type (vector) :: S1_vector, S2_vector, S3_vector
        real (kind=r8b) :: sigma1, sigma2, sigma3, phi
        real (kind=r8b) :: fluid_pressure, cohesive_strength, intern_friction_resist
    end type stress_princcompon

    type :: stress_solution
        logical :: valid_solution
        type(vector) :: traction_vect, normalstress_vect, shearstress_vect
        real (kind=r8b)	:: theor_rake, tractionvect_magn, normalstress_magn, shearstress_magn
        real (kind=r8b)	:: modif_sliptendency, deformation_index
        real (kind=r8b) :: dilation_tendency, fracture_stability, leakage_factor
        type (axis) :: theor_slickenline
    end type stress_solution

    type :: stressanalysis_param
        type(stress_princcompon) :: stressprinccomp1
        real (kind=r8b) :: tens(3,3)
        real (kind=r8b) :: rot_matrix(3,3)
    end type stressanalysis_param

    type :: record_variables
        logical :: spatial_var,time_var
        integer (kind=i1b) :: spat_coord_system, readcase
    end type record_variables

    type :: faultsimulation_param
        integer (kind=i4b) :: faults_totnumb
        type (range_3D) :: spat_bound
        type (range_1D) :: temp_bound
    end type faultsimulation_param

    type :: analysis_parameters
        integer (kind=i4b) :: time_initial(8)
        type (record_variables) :: fault_loctimevar
        integer (kind=i1b) :: input_case
        character (len=50) :: inputfile_name
        integer (kind=i2b) :: inputfile_headerrownumb
        type (faultsimulation_param) :: faultsimulpar
        integer (kind=i4b) :: faults_totnumb
    end type analysis_parameters


    integer :: ios !status of input/output connection
    real (kind=r8b) :: pi  	! pi radians
    real (kind=r8b) :: r2d, d2r	! for conversion from radians to degrees

    type(orthonormal_triad), parameter :: frame0 = orthonormal_triad(vector(1.0,0.0,0.0),vector(0.0,1.0,0.0),vector(0.0,0.0,1.0))
    real (kind=r8b),parameter :: unitary_vect_tolerance = 1.0e-6 !maximum accepted value for unitary vector magnitude difference with unit
    real (kind=r8b),parameter :: vect_normaliz_tolerance = 1.0e-5 !minimum accepted value for vector magnitude to apply normalization
    real (kind=r8b), parameter :: shearmagnitude_minthresh = 1.0e-5 !minimum accepted value for shear stress to be considered meaningful

end module data_types


module vector_processing

    use data_types

    contains

    ! calculates the sum of two vectors

    type(vector) function vector_sum(vector1, vector2)

        type(vector), intent(in) :: vector1, vector2

        vector_sum%x = vector1%x + vector2%x
        vector_sum%y = vector1%y + vector2%y
        vector_sum%z = vector1%z + vector2%z

    end function vector_sum

    ! calculates the difference between two vectors

    type(vector) function vector_diff(vector1,vector2)

        type(vector), intent(in) :: vector1, vector2

        vector_diff%x = vector1%x - vector2%x
        vector_diff%y = vector1%y - vector2%y
        vector_diff%z = vector1%z - vector2%z

    end function vector_diff

    ! angle (in radians) between two vectors (radiants, 0-pi)

    real (kind=r8b) function vector_angle_rad(vector1,vector2)

        type(vector), intent(in) :: vector1, vector2

        vector1_magn = vector_magn(vector1)
        vector2_magn = vector_magn(vector2)
        if ((vector1_magn < vect_normaliz_tolerance).or.(vector2_magn < vect_normaliz_tolerance)) then
            write(*,*) 'Error in vector magnitude (function vector_angle_rad). Hit any key to stop'
            read(*,*)
            stop	! STOP PROGRAM FOR ERROR IN DATA INPUT
        end if

        ! scalar product between two vectors
        scaledscalarproduct = vector_scalprod(vector1, vector2)/(vector1_magn*vector2_magn)

        ! angle between vectors (in radians)
        if (scaledscalarproduct < -1.) then
            vector_angle_rad = pi
        else if (scaledscalarproduct > 1.) then
            vector_angle_rad = 0.0
        else
            vector_angle_rad = acos(scaledscalarproduct)
        end if

    end  function vector_angle_rad

    ! angle (in radians) between two axes (radiants, 0-pi/2)

    real (kind=r8b) function axes_angle_rad(vector1,vector2)

        type(vector), intent(in) :: vector1, vector2

        ! angle between vectors (in radians)
        axes_angle_rad = vector_angle_rad(vector1,vector2)
        axes_angle_rad = min(axes_angle_rad, pi-axes_angle_rad)

    end function axes_angle_rad

    ! vector normalization

    type(vector) function vector_normalization(vector1) result(vector2)

        type(vector), intent(in) :: vector1

        vector1_magn = vector_magn(vector1)

        if (vector1_magn < vect_normaliz_tolerance) then
          write(*,*) 'Error in vector magnitude processing. Hit any key to stop'
          read (*,*)
          stop ! STOP PROGRAM FOR ERROR IN DATA INPUT
        end if

        vector2%x = vector1%x/vector1_magn
        vector2%y = vector1%y/vector1_magn
        vector2%z = vector1%z/vector1_magn


    end function vector_normalization

    ! vector magnitude

    real (kind=r8b) function vector_magn(vector1)

        type(vector), intent(in) :: vector1

        vector_magn = sqrt((vector1%x)**2 + (vector1%y)**2 + (vector1%z)**2)

    end function vector_magn

    ! calculates the product of a vector by a scalar

    type (vector) function vectorbyscalar(vector1, scalar1)

        type(vector), intent(in) :: vector1
        real (kind=r8b), intent(in) :: scalar1

        vectorbyscalar%x = scalar1 * vector1%x
        vectorbyscalar%y = scalar1 * vector1%y
        vectorbyscalar%z = scalar1 * vector1%z

    end function vectorbyscalar

    ! scalar product of two vectors (given as their cartesian coordinates)

    real (kind=r8b) function vector_scalprod(vector1, vector2)

        type(vector), intent(in) :: vector1, vector2

        vector_scalprod = vector1%x * vector2%x + vector1%y * vector2%y + vector1%z * vector2%z

    end function vector_scalprod

    ! vectorial product of two vectors (given as their cartesian coordinates)

    type(vector) function vector_vectprod(vector1,vector2)

        type(vector), intent(in) :: vector1, vector2

        vector_vectprod%x=(vector1%y*vector2%z)-(vector1%z*vector2%y)
        vector_vectprod%y=(vector1%z*vector2%x)-(vector1%x*vector2%z)
        vector_vectprod%z=(vector1%x*vector2%y)-(vector1%y*vector2%x)

    end function vector_vectprod

    ! vector1 projection on vector2

    type (vector) function vector_projection(vector1,vector2)

        type(vector), intent(in) :: vector1, vector2
        real (kind=r8b) :: scalprod

        scalprod = vector_scalprod(vector1, vector_normalization(vector2))
        vector_projection = vectorbyscalar(vector2, scalprod)

    end function vector_projection

    ! test if a vector has magnitude = 1

    logical function vect_normaliztest(vector1) result(vect_is_normalized)

        type(vector), intent(in) :: vector1
        real (kind=r8b) :: vector1_magn

        vector1_magn = vector_magn(vector1)

        abs_diff = dabs(1-vector1_magn)

        if (abs_diff > unitary_vect_tolerance) then
            vect_is_normalized = .false.
        else
            vect_is_normalized = .true.
        endif

    end function vect_normaliztest

    ! converts a 3D vector to a 3x1 array

    function vector2array(vector1) result(array1)

        type(vector), intent(in) :: vector1
        real(kind=r8b) :: array1(3)

        array1 = (/vector1%x, vector1%y, vector1%z/)

    end function vector2array

    ! converts a 3x1 array to a 3D vector

    function array2vector(array1) result(vector1)

        real(kind=r8b), intent(in) :: array1(3)
        type(vector):: vector1

        vector1%x = array1(1)
        vector1%y = array1(2)
        vector1%z = array1(3)

    end function array2vector


end module vector_processing


module geometric_processing

    use vector_processing

    contains

    ! calculates the polar components from the cartesian ones

    type(vector) function pole2cartesian(axis1) result(vector1)

        type(axis), intent(in) :: axis1


        vector1%x = cos(d2r*axis1%plunge) * cos(d2r*axis1%trend)
        vector1%y  = cos(d2r*axis1%plunge) * sin(d2r*axis1%trend)
        vector1%z  =  sin(d2r*axis1%plunge)


    end function pole2cartesian

    ! calculates polar components from cartesian ones

    type(axis) function cartesian2pole(vector1) result(axis1)

        type(vector):: vector1
        logical :: vect_is_normalized

        vect_is_normalized = vect_normaliztest(vector1)

        if (.not.vect_is_normalized) then
          vector1 = vector_normalization(vector1)
        endif

        ! polar coordinates calculation

        if (vector1%z > 1.0_r8b) then
            vector1%z = 1.0_r8b
        elseif (vector1%z < -1.0_r8b) then
            vector1%z = -1.0_r8b
        endif

        axis1%plunge = r2d*dasin(vector1%z)

        if (dabs(axis1%plunge)> 89.5) then
                axis1%trend = 0.0
        else
            axis1%trend = r2d*atan2(vector1%y,vector1%x)
            if (axis1%trend < 0.0) then
                axis1%trend = 360.0 + axis1%trend
            endif
        endif

    end function cartesian2pole

    ! calculates the down axis from the axis

    type(axis) function axis2downaxis(axis1) result(axis2)

        type(axis), intent(in) :: axis1

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

        type(plane), intent(inout) :: faultplane1

        faultplane1%dipdirection = faultplane1%strike_rhr + 90.0
        if (faultplane1%dipdirection >= 360.0) then
          faultplane1%dipdirection = faultplane1%dipdirection - 360.0
        endif

    end subroutine dipdir_calc

    ! calculates fault normal

    type(vector) function faultplanenormal_calc(faultplane1) result(faultnorm)

        type(plane), intent(in) :: faultplane1

        ! Fault Normal cartesian coordinates
        ! formulas from Aki and Richards, 1980
        faultnorm%x = -sin(d2r*faultplane1%dipangle) * sin(d2r*faultplane1%strike_rhr)
        faultnorm%y = sin(d2r*faultplane1%dipangle) * cos(d2r*faultplane1%strike_rhr)
        faultnorm%z = -cos(d2r*faultplane1%dipangle)

    end function faultplanenormal_calc

    ! calculates vector (cartesian) components of fault record

    subroutine fault_cartcmp(faultplane1, rake1, faultnorm_vect, slickenl_vect)

        ! QUATf_faultpole2faultvector

        type(plane), intent(in) :: faultplane1
        real (kind=r8b)	:: rake1
        type (vector) , intent(out) :: faultnorm_vect, slickenl_vect

        ! Fault Normal and Slickenline cartesian coordinates
        ! formulas from Aki and Richards, 1980
        faultnorm_vect%x = -sin(d2r*faultplane1%dipangle)*sin(d2r*faultplane1%strike_rhr)
        faultnorm_vect%y = sin(d2r*faultplane1%dipangle)*cos(d2r*faultplane1%strike_rhr)
        faultnorm_vect%z = -cos(d2r*faultplane1%dipangle)

        slickenl_vect%x = cos(d2r*rake1) * cos(d2r*faultplane1%strike_rhr)   &
         + sin(d2r*rake1) * cos(d2r*faultplane1%dipangle) * sin(d2r*faultplane1%strike_rhr)
        slickenl_vect%y = cos(d2r*rake1) * sin(d2r*faultplane1%strike_rhr)  &
         - sin(d2r*rake1) * cos(d2r*faultplane1%dipangle) * cos(d2r*faultplane1%strike_rhr)
        slickenl_vect%z = -sin(d2r*rake1) * sin(d2r*faultplane1%dipangle)

    end subroutine fault_cartcmp

    ! calculate the slickenline (trend and plunge) from the rake angle

    type(axis) function rake2slickenline(strike, dip, rake) result(slicken_pole)

        implicit none

        real (kind=r8b), intent(in) :: strike, dip, rake
        real (kind=r8b) :: strike_rd, dip_rd, rake_rd
        type (vector) :: slick_vect

        strike_rd = d2r*strike
        dip_rd = d2r*dip
        rake_rd = d2r*rake

        slick_vect%x = cos(rake_rd) * cos(strike_rd) + sin(rake_rd) * sin(strike_rd) * cos(dip_rd)
        slick_vect%y = cos(rake_rd) * sin(strike_rd) - sin(rake_rd) * cos(strike_rd) * cos(dip_rd)
        slick_vect%z = -sin(rake_rd) * sin(dip_rd)

        ! determination of trend and plunge of Slickenline
        slicken_pole = cartesian2pole(slick_vect)

    end function rake2slickenline


end module fault_processing


module stress_processing

    use fault_processing

    implicit none

    contains

    ! calculation of S2_pole as vector product of S3_pole and S1_pole

    subroutine S2_calc(stressprinccomp1)

        type(stress_princcompon), intent(inout) :: stressprinccomp1

        stressprinccomp1%S2_vector = vector_vectprod(stressprinccomp1%S3_vector,stressprinccomp1%S1_vector)
        stressprinccomp1%S2_pole = cartesian2pole(stressprinccomp1%S2_vector)
        stressprinccomp1%S2_pole = axis2downaxis(stressprinccomp1%S2_pole)

    end subroutine S2_calc

    ! calculation of sigma2 from S1_pole, S3_pole and PHI values

    subroutine sigma2_calc(stressprinccomp1)

        type(stress_princcompon), intent(inout) :: stressprinccomp1

        stressprinccomp1%sigma2 = stressprinccomp1%phi*stressprinccomp1%sigma1 				&
                                    + (1-stressprinccomp1%phi)*stressprinccomp1%sigma3


    end subroutine sigma2_calc

    ! define rotation matrix based on stress principal axes

    function rotmatr(stressprinccomp1) result(rot_matrix)

        ! based on Kuipers, 2002, p.161, eqs. 7.8

        type(stress_princcompon), intent(in) :: stressprinccomp1
        real (kind=r8b) :: rot_matrix(3,3)

        rot_matrix(1,1) = vector_scalprod(stressprinccomp1%S1_vector, frame0%X)
        rot_matrix(1,2) = vector_scalprod(stressprinccomp1%S2_vector, frame0%X)
        rot_matrix(1,3) = vector_scalprod(stressprinccomp1%S3_vector, frame0%X)

        rot_matrix(2,1) = vector_scalprod(stressprinccomp1%S1_vector, frame0%Y)
        rot_matrix(2,2) = vector_scalprod(stressprinccomp1%S2_vector, frame0%Y)
        rot_matrix(2,3) = vector_scalprod(stressprinccomp1%S3_vector, frame0%Y)

        rot_matrix(3,1) = vector_scalprod(stressprinccomp1%S1_vector, frame0%Z)
        rot_matrix(3,2) = vector_scalprod(stressprinccomp1%S2_vector, frame0%Z)
        rot_matrix(3,3) = vector_scalprod(stressprinccomp1%S3_vector, frame0%Z)


    end function rotmatr

    ! calculate stress tensor expressed in frame components

    function stresstensorcalc(stressprinccomp1,rot_matrix) result(tens)

        ! from Kagan and Knopoff, 1985a, p. 433

        type(stress_princcompon), intent(in) :: stressprinccomp1
        real (kind=r8b), intent(in)  :: rot_matrix(3,3)

        real (kind=r8b) :: stresstens0(3,3), tens(3,3)

        ! stress eigentensor
        stresstens0 = reshape((/0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/),shape=(/3,3/),order=(/2,1/))
        stresstens0(1,1) = stressprinccomp1%sigma1
        stresstens0(2,2) = stressprinccomp1%sigma2
        stresstens0(3,3) = stressprinccomp1%sigma3

        tens = matmul(rot_matrix,matmul(stresstens0,transpose(rot_matrix)))

    end function stresstensorcalc

    ! calculation of stress solution for the given plane and stress tensor

    type(stress_solution) function stresssolution_calc(tens1,stressprinccomp1,faultrec1) result(stresssolution1)

        ! based on Xu, 2004 (Geoph. J. Int., 157,1316-1330) and references within.

        ! QUATstress_stress&fault2rake
        !last modified: 2008/01/05

        ! INPUT
        !    SELF with 2 objects:
        !      (0): fault plane orientation: strike rhr, dip direction, dip angle
        !      (1): stress field tensor components: x,y,z (sigma 1,2,3)

        !  OUTPUT
        !    normal and tangential stresses magnitudes
        !    slip tendency
        !    theoretical rake angle of tangential stress

        !  STEPS
        !    input of data:
        !      fault plane orientation: strike rhr, dip direction, dip angle     stress field tensor components: x,y,z (sigma 1,2,3)
        !    fault normal calculation: n (1,2,3)
        !    traction vector c.: sigma (1,2,3)
        !    normal stress c.: sigma n (1,2,3)
        !    shear stress c: tau s (1,2,3)
        !    slip tendency c.: Ts
        !    theoretical rake angle c: lambda
        !    output of results:
        !      magnitudes of normal and shear stresses,
        !      slip tendency,
        !      theoretical rake angle

        !    input of data:
        !      fault plane orientation: strike rhr, dip direction, dip angle
        !      stress field tensor components: x,y,z (sigma 1,2,3)

        implicit none

        real (kind=r8b), intent(in)  :: tens1(3,3)
        type(stress_princcompon), intent(in) :: stressprinccomp1
        type(fault_datum), intent(in)  :: faultrec1

        type(vector) :: shearstress_unitvect
        type(axis) :: strike_pole, dipdir_pole
        type(vector) :: strike_vect, dipdir_vect, faultnorm_vect
        real (kind=r8b) :: faultnormal_array(3), scalprod_shearstress_strike, scalprod_shearstress_dipdir
        real (kind=r8b) :: tau

        ! calculation of fault normal calculation: n (1,2,3)
        faultnorm_vect = faultplanenormal_calc(faultrec1%fltplane)
        faultnormal_array = vector2array(faultnorm_vect)

        !    traction vector c.: sigma (1,2,3)
        stresssolution1%traction_vect = array2vector(- matmul(tens1,faultnormal_array))
        stresssolution1%tractionvect_magn = vector_magn(stresssolution1%traction_vect)
        if (vector_scalprod(stresssolution1%traction_vect, faultnorm_vect)<0.0_r8b) then
            stresssolution1%tractionvect_magn = -stresssolution1%tractionvect_magn
        endif

        !    normal stress c.: SigmaN
        stresssolution1%normalstress_vect = vector_projection(stresssolution1%traction_vect, faultnorm_vect)
        stresssolution1%normalstress_magn = vector_magn(stresssolution1%normalstress_vect)
        if (vector_scalprod(stresssolution1%normalstress_vect, faultnorm_vect)<0.0_r8b) then
            stresssolution1%normalstress_magn = -stresssolution1%normalstress_magn
        endif

        !    shear stress c: tau s (1,2,3)
        stresssolution1%shearstress_vect = vector_diff(stresssolution1%traction_vect,stresssolution1%normalstress_vect)
        stresssolution1%shearstress_magn = vector_magn(stresssolution1%shearstress_vect)


        if (stresssolution1%shearstress_magn <= shearmagnitude_minthresh) then
          stresssolution1%valid_solution = .false.
        else
          stresssolution1%valid_solution = .true.
        end if

        if (stresssolution1%valid_solution) then

            shearstress_unitvect = vector_normalization(stresssolution1%shearstress_vect)

            ! theoretical rake angle c: lambda
            strike_pole = axis(faultrec1%fltplane%strike_rhr,0)
            dipdir_pole = axis(faultrec1%fltplane%dipdirection,faultrec1%fltplane%dipangle)

            strike_vect = vector_normalization(pole2cartesian(strike_pole))
            dipdir_vect = vector_normalization(pole2cartesian(dipdir_pole))

            scalprod_shearstress_strike = vector_scalprod(shearstress_unitvect, strike_vect)

            scalprod_shearstress_dipdir = vector_scalprod(shearstress_unitvect, dipdir_vect)

            if (scalprod_shearstress_strike>1) then
                scalprod_shearstress_strike=1
            elseif (scalprod_shearstress_strike<-1) then
                scalprod_shearstress_strike=-1
            endif

            stresssolution1%theor_rake = r2d*acos(scalprod_shearstress_strike)

            if (scalprod_shearstress_dipdir > 0.0_r8b) then
                stresssolution1%theor_rake = - stresssolution1%theor_rake
            endif

            stresssolution1%theor_slickenline = rake2slickenline(   &
                faultrec1%fltplane%strike_rhr,faultrec1%fltplane%dipangle,stresssolution1%theor_rake)

            stresssolution1%modif_sliptendency = stresssolution1%shearstress_magn/dabs(stresssolution1%tractionvect_magn)

            stresssolution1%deformation_index = &
            (dabs(stresssolution1%tractionvect_magn) - stresssolution1%shearstress_magn)/stresssolution1%tractionvect_magn

            stresssolution1%dilation_tendency = &
            (stressprinccomp1%sigma1 - dabs(stresssolution1%normalstress_magn))/(stressprinccomp1%sigma1 - stressprinccomp1%sigma3)

            tau = stressprinccomp1%cohesive_strength+stressprinccomp1%intern_friction_resist*dabs(stresssolution1%normalstress_magn)
            stresssolution1%fracture_stability = &
            (dabs(stresssolution1%normalstress_magn) - tau)/ stressprinccomp1%intern_friction_resist

            stresssolution1%leakage_factor = stressprinccomp1%fluid_pressure/(dabs(stresssolution1%normalstress_magn) - tau)

        endif

    end function stresssolution_calc


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

    call WelcomeScreen()

    ! user-definition of source file

    call InputFileDef()

    ! user-definition of output file

    call OutputFileDef()

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

    ! calculates resulting product vector

    b = matmul ( u, matmul ( sigma, transpose(v) ) )

    ! writes resulting product vector in output file

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
