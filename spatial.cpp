#include "cmath"
#include "spatial.hpp"
#include <cassert>


double plane_parallel_threshold_degr = 0.1; // in degrees
double plane_equidistance_threshold = 1.0e-10; // data distance units

Point2D::Point2D() {

}


Point2D::Point2D(double x,  double y) {

    pt[0] = x;
    pt[1] = y;

};


Point2D::~Point2D() {

}


std::array<double, 2> Point2D::arr_pt()  // C++11
{
    return pt;
}



double Point2D::x() {

    return pt[0];

};


double Point2D::y() {

    return pt[1];

};


Point2D Point2D::operator=(Point2D p2) {

    return Point2D(p2.x(), p2.y());

};


double Point2D::distance(Point2D another) {

    double dx = x() - another.x();
    double dy = y() - another.y();

    return sqrt( dx*dx + dy*dy);

};


Point2D Point2D::moveby(double move_x, double move_y) {

    double new_x = x() + move_x;
    double new_y = y() + move_y;

    return Point2D(new_x, new_y);

};


Point2D Point2D::moveby(Point2D shift_pt) {

    double new_x = x() + shift_pt.x();
    double new_y = y() + shift_pt.y();

    return Point2D(new_x, new_y);

};


Point2D Point2D::rotateby(Matrix2 rm) {

    double new_x = rm.m11()*x() + rm.m12()*y();
    double new_y = rm.m21()*x() + rm.m22()*y();

    return Point2D(new_x, new_y);

};


Point3D::Point3D() {

    valid = false;

    _x = 0.0;
    _y = 0.0;
    _z = 0.0;

};


Point3D::Point3D(double x, double y, double z, bool valid_) {

    valid =  valid_;

    _x = x;
    _y = y;
    _z = z;

};


Point3D::Point3D(double x, double y, double z) {

    valid =  true;

    _x = x;
    _y = y;
    _z = z;

};



Point3D::Point3D(Point2D pt2d, double z, bool valid_) {

    valid = valid_;

    _x = pt2d.x();
    _y = pt2d.y();
    _z = z;

};


Point3D::Point3D(Point2D pt2d, double z) {

    valid = true;

    _x = pt2d.x();
    _y = pt2d.y();
    _z = z;

};


Point3D::~Point3D() {

}


bool Point3D::is_valid() {

    return valid;

};


double Point3D::x() {

    return _x;
};


double Point3D::y() {

    return _y;

};


double Point3D::z() {

    return _z;

};


Point3D Point3D::moveby(double dx, double dy, double dz) {

    return Point3D(x() + dx,
                   y() + dy,
                   z() + dz);
};


double Point3D::distance(Point3D another) {

    double dx = x() - another.x();
    double dy = y() - another.y();
    double dz = z() - another.z();

    return sqrt( dx*dx + dy*dy + dz*dz);

};


bool Point3D::is_coincident(Point3D another) {

    if (distance(another) > 1e-12) {
        return false; }
    else {
        return true;
    }

};


Segment3D::Segment3D(Point3D start_pt, Point3D end_pt) :
    _start_pt(start_pt), _end_pt(end_pt) {

};


Segment3D::~Segment3D() {

};


Point3D Segment3D::start_pt(){

    return _start_pt;

};


Point3D Segment3D::end_pt(){

    return _end_pt;

};


double Segment3D::dx() {

    return _end_pt.x() - _start_pt.x();

};


double Segment3D::dy() {

    return _end_pt.y() - _start_pt.y();

};


double Segment3D::dz() {

    return _end_pt.z() - _start_pt.z();

};


Vector3D Segment3D::as_vector() {

    return Vector3D(dx(), dy(), dz() );

};


Vector3D Segment3D::as_versor() {

    return as_vector().versor();
};


double Segment3D::length() {

    double d_x = dx();
    double d_y = dy();
    double d_z = dz();

    return sqrt( d_x*d_x + d_y*d_y + d_z*d_z );

};


bool Segment3D::is_point_projection_in_segment(Point3D pt_3d) {

    Vector3D pt_vector = Segment3D( start_pt(), pt_3d ).as_vector();
    double scal_prod = as_vector().scalar_prod( pt_vector );
    double segm_length = length();
    if (0 <= scal_prod and scal_prod <= segm_length*segm_length) {
        return true; }
    else {
        return false; };

};


Line3D Segment3D::as_line() {

    Point3D start_pt3d = start_pt();
    Vector3D line_vect = as_versor();
    return Line3D(start_pt3d, line_vect);

};


Vector3D::Vector3D(double x, double y, double z) :
    _x(x), _y(y), _z(z) {
    };


Vector3D::Vector3D( Point3D pt3d_a, Point3D pt3d_b) {

    _x = pt3d_b.x() - pt3d_a.x();
    _y = pt3d_b.y() - pt3d_a.y();
    _z = pt3d_b.z() - pt3d_a.z();

};


Vector3D::~Vector3D() {

};


double Vector3D::x() {

    return _x;

};


double Vector3D::y() {

    return _y;

};


double Vector3D::z() {

    return _z;

};


double Vector3D::length() {

    return sqrt(_x*_x + _y*_y + _z*_z);

};


Vector3D Vector3D::scale(double scale_factor) {

    return Vector3D(_x*scale_factor, _y*scale_factor, _z*scale_factor);

};


Vector3D Vector3D::versor() {

    return scale( 1.0 / length() );

};


double Vector3D::scalar_prod(Vector3D another) {

    return _x*another._x + _y*another._y + _z*another._z;

};


Vector3D Vector3D::vector_prod(Vector3D another) {

    double x = _y * another.z() - _z * another.y();
    double y = _z * another.x() - _x * another.z();
    double z = _x * another.y() - _y * another.x();

    return Vector3D(x, y, z);

};


double Vector3D::angle(Vector3D another) {

    double cos_angle = scalar_prod(another) / (length()*another.length());
    if (cos_angle > 1.0) {
        cos_angle = 1.0;}
    else if (cos_angle < -1.0) {
        cos_angle = -1.0;
    };

    double angle_degr = degrees(acos(cos_angle));

    return angle_degr;

};


bool Vector3D::isodirection(Vector3D another) {

    Vector3D vect_prod = vector_prod(another);
    double norm = vect_prod.length();
    if (norm < 1e-12) {
        return true; }
    else {
        return false; }

};


Point3D Vector3D::move_pt(Point3D pt3d) {

    double x0 = pt3d.x();
    double y0 = pt3d.y();
    double z0 = pt3d.z();

    double dx = _x;
    double dy = _y;
    double dz = _z;

    return Point3D(x0+dx, y0+dy, z0+dz);
};


Line3D::Line3D(Point3D orig_pt3d, Vector3D vers) :
    _orig_pt3d(orig_pt3d), _vers(vers) {

    };


Line3D::~Line3D() {

};


Point3D Line3D::orig_pt() {

    return _orig_pt3d;

};


Vector3D Line3D::versor() {

    return _vers;
};


bool Line3D::isparallel(Line3D another) {

    Vector3D vers_a = versor();
    Vector3D vers_b = another.versor();
    if (vers_a.isodirection(vers_b)) {
        return true; }
    else {
        return false;};
    };


bool Line3D::iscoincident(Line3D another) {

    if (not isparallel(another)) {
        return false; };

    Point3D orig_pt_a = orig_pt();
    Point3D orig_pt_b = another.orig_pt();

    if (orig_pt_a.is_coincident(orig_pt_b)) {
        return true; }
    else {
        Vector3D test_line = Vector3D(orig_pt_a, orig_pt_b);
        if (test_line.isodirection(versor())) {
            return true;}
        else {
            return false;}
        }
};


Point3D Line3D::intersect_coplanar(Line3D another) {

    Point3D line_a_start_pt = orig_pt();
    Vector3D vers_a = versor();

    std::cout << "\nIntersect_coplanar\n";
    std::cout << "line a\n";
    std::cout << " - pt: " << line_a_start_pt.x() << " " << line_a_start_pt.y() << " " << line_a_start_pt.z() << "\n";
    std::cout << " - versor: " << vers_a.x() << " " << vers_a.y() << " " << vers_a.z() << "\n";

    Point3D line_b_start_pt = another.orig_pt();
    Vector3D vers_b = another.versor();

    std::cout << "line b\n";
    std::cout << " - pt: " << line_b_start_pt.x() << " " << line_b_start_pt.y() << " " << line_b_start_pt.z() << "\n";
    std::cout << " - versor: " << vers_b.x() << " " << vers_b.y() << " " << vers_b.z() << "\n";

    Point3D point_a = Point3D(line_a_start_pt.x(), line_a_start_pt.y(), line_a_start_pt.z());
    Point3D point_b = Point3D(line_b_start_pt.x(), line_b_start_pt.y(), line_b_start_pt.z());

    std::cout << "\point_a\n";
    std::cout << "x: " << point_a.x() << " y: " << point_a.y() << " z: " << point_a.z() << "\n";
    std::cout << "\point_b\n";
    std::cout << "x: " << point_b.x() << " y: " << point_b.y() << " z: " << point_b.z() << "\n";

    double delta_distance = 100.0;
    Vector3D displ_vector_a = vers_a.scale(delta_distance);
    Vector3D displ_vector_b = vers_b.scale(delta_distance);

    if (point_b.distance(point_a) < delta_distance) {
        point_b = displ_vector_b.move_pt(point_b);
    };

    Point3D point_c = displ_vector_b.move_pt(point_b);
    point_c = displ_vector_a.move_pt(point_c);

    assert( not point_a.is_coincident(point_b));
    assert( not point_c.is_coincident(point_a));
    assert( not point_c.is_coincident(point_b));

    Vector3D vect1 = Vector3D(point_a, point_b);
    Vector3D vect2 = Vector3D(point_a, point_c);

    assert( not vect1.isodirection(vect2));

    std::cout << "\point_a\n";
    std::cout << "x: " << point_a.x() << " y: " << point_a.y() << " z: " << point_a.z() << "\n";
    std::cout << "\point_b\n";
    std::cout << "x: " << point_b.x() << " y: " << point_b.y() << " z: " << point_b.z() << "\n";
    std::cout << "\point_c\n";
    std::cout << "x: " << point_c.x() << " y: " << point_c.y() << " z: " << point_c.z() << "\n";

    CartesianPlane colinear_plane = CartesianPlane(point_a, point_b, point_c);

    std::cout << "\ncolinear plane\n";
    std::cout << " a: " << colinear_plane.a() << " b: " << colinear_plane.b() << " c: " << colinear_plane.c() << " d: " << colinear_plane.d() << "\n";

    //code inspired to: http://geomalgorithms.com/a05-_intersect-1.html#intersect2D_2Segments()
    Vector3D w_vect = Vector3D( line_a_start_pt, line_b_start_pt );
    Vector3D vers_a_perp = colinear_plane.perp_versor_in_plane(vers_a);

    double factor_numerator = - vers_a_perp.scalar_prod(w_vect);
    double factor_denominator = vers_a_perp.scalar_prod(vers_b);

    std::cout << "\nfactor_numerator: " << factor_numerator << "\n";
    std::cout << "\nfactor_denominator: " << factor_denominator << "\n";

    double factor_scaling = factor_numerator / factor_denominator;

    std::cout << "\nfactor_scaling: " << factor_scaling << "\n";

    Point3D intersection_pt3d = vers_b.scale(factor_scaling).move_pt(line_b_start_pt);

    return intersection_pt3d;

};


/*
CartesianPlane::CartesianPlane(double a, double b, double c, double d) {

    std::cout << "** " << a << " " << b << " " << c << " " << d << "\n";

    double division_factor = sqrt(a*a + b*b + c*c);
    std::cout << "*** " << division_factor << "\n";
    if (d > 0.0) {
        division_factor = - division_factor;
    }
    std::cout << "**** " << c / division_factor << "\n";
    _a = a / division_factor;
    _b = b / division_factor;
    _c = c / division_factor;
    _d = d / division_factor;
};
*/


std::vector<double> normalize(std::vector<double> coeff) {

    double a = coeff[0];
    double b = coeff[1];
    double c = coeff[2];
    double d = coeff[3];

    std::vector<double> norm_coeff(4);

    double division_factor = sqrt(a*a + b*b + c*c);
    if (d > 0.0) {
        division_factor = - division_factor;
    }

    norm_coeff[0] = a / division_factor;
    norm_coeff[1] = b / division_factor;
    norm_coeff[2] = c / division_factor;
    norm_coeff[3] = d / division_factor;

    return norm_coeff;
};


CartesianPlane::CartesianPlane(Point3D pt_a, Point3D pt_b, Point3D pt_c) {

    Matrix3 matr_a = Matrix3(pt_a.y(), pt_a.z(), 1,
                             pt_b.y(), pt_b.z(), 1,
                             pt_c.y(), pt_c.z(), 1);

    Matrix3 matr_b = Matrix3(pt_a.x(), pt_a.z(), 1,
                             pt_b.x(), pt_b.z(), 1,
                             pt_c.x(), pt_c.z(), 1);

    Matrix3 matr_c = Matrix3(pt_a.x(), pt_a.y(), 1,
                             pt_b.x(), pt_b.y(), 1,
                             pt_c.x(), pt_c.y(), 1);

    Matrix3 matr_d = Matrix3(pt_a.x(), pt_a.y(), pt_a.z(),
                             pt_b.x(), pt_b.y(), pt_b.z(),
                             pt_c.x(), pt_c.y(), pt_c.z());

    std::vector<double> coeff(4);
    coeff[0] = matr_a.determinant();
    coeff[1] = - matr_b.determinant();
    coeff[2] = matr_c.determinant();
    coeff[3] = - matr_d.determinant();

    //std::cout << "* " << coeff[0] << " " << coeff[1] << " " << coeff[2] << " " << coeff[3] << "\n";

    std::vector<double> normalized_coeff = normalize(coeff);

    _a = normalized_coeff[0];
    _b = normalized_coeff[1];
    _c = normalized_coeff[2];
    _d = normalized_coeff[3];

};


CartesianPlane::~CartesianPlane() {
};


double CartesianPlane::a() {
    return _a;
};

double CartesianPlane::b() {
    return _b;
};

double CartesianPlane::c() {
    return _c;
};

double CartesianPlane::d() {
    return _d;
};


Vector3D CartesianPlane::normal_versor() {

    return Vector3D( _a, _b, _c ).versor();

};


Vector3D CartesianPlane::perp_versor_in_plane(Vector3D inplane_vect3d) {

    return inplane_vect3d.vector_prod( normal_versor() ).versor();
};


Vector3D CartesianPlane::intersect_versor(CartesianPlane another) {

    // TODO: consider the special case of two parallel planes
    Vector3D vers_a = normal_versor();
    Vector3D vers_b = another.normal_versor();

    return vers_a.vector_prod(vers_b).versor();

};


Point3D CartesianPlane::intersect_point3d(CartesianPlane another) {

    Matrix2 matrix_ab = Matrix2(_a, _b, another._a, another._b);
    Matrix2 matrix_ac = Matrix2(_a, _c, another._a, another._c);
    Matrix2 matrix_bc = Matrix2(_b, _c, another._b, another._c);

    double determinant_ab = matrix_ab.determinant();
    double determinant_ac = matrix_ac.determinant();
    double determinant_bc = matrix_bc.determinant();

    Matrix2 matrix_ad = Matrix2(_a, - _d, another._a, - another._d);
    Matrix2 matrix_bd = Matrix2(_b, - _d, another._b, - another._d);
    Matrix2 matrix_dc = Matrix2(- _d, _c, - another._d, another._c);
    Matrix2 matrix_db = Matrix2(- _d, _b, - another._d, another._b);

    double x, y, z;
    if ( fabs( determinant_ab ) > 1e-12 ) {
        z = 0.0;
        x = matrix_db.determinant() / determinant_ab;
        y = matrix_ad.determinant() / determinant_ab;
    }
    else if ( fabs( determinant_ac ) > 1e-12 ) {
        y = 0.0;
        x = matrix_dc.determinant() / determinant_ac;
        z = matrix_ad.determinant() / determinant_ac;
    }
    else if ( fabs( determinant_bc ) > 1e-12 ) {
        x = 0.0;
        y = matrix_dc.determinant() / determinant_bc;
        z = matrix_bd.determinant() / determinant_bc;
    }

    return Point3D(x, y, z);
};


Line3D CartesianPlane::intersect(CartesianPlane another) {

    Vector3D vers = intersect_versor(another);
    Point3D pt3d = intersect_point3d(another);

    return Line3D( pt3d, vers);

};


bool CartesianPlane::point_in_plane(Point3D pt) {

    double val = std::abs(_a*pt.x() + _b*pt.y() + _c*pt.z() + _d);
    //std::cout << val << "\n";
    if (val < 1e-12) {
        return true; }
    else {
        return false;
    };

};


double CartesianPlane::angle(CartesianPlane another) {

Vector3D normal1 = normal_versor();
Vector3D normal2 = another.normal_versor();

double angle_degr = normal1.angle(normal2);

return angle_degr;

};


bool CartesianPlane::isparallel(CartesianPlane another) {

    double angle_degr = angle(another);

    if (angle_degr < plane_parallel_threshold_degr) {
        return true;}
    else if (angle_degr > 180.0 - plane_parallel_threshold_degr) {
        return true;}
    else {
        return false;};

};

bool CartesianPlane::isequidistant(CartesianPlane another) {

    if (fabs(d() - another.d()) < plane_equidistance_threshold) {
        return true; }
    else {
        return false; };

};






