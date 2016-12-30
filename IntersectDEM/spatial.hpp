#include <string>
#include <vector>

#include <array>

#include "algebr.hpp"

#ifndef SPATIAL_HPP
#define SPATIAL_HPP


std::vector<double> normalize(std::vector<double>);


class Point2D {

    std::array<double, 2> pt;

public:
    //Point2D();
    Point2D();
    Point2D(double, double);
    ~Point2D();
    //double arr_pt [2];
    std::array<double, 2> arr_pt();
    double x();
    double y();
    Point2D operator=(Point2D);
    double distance(Point2D);
    Point2D moveby(double, double);
    Point2D moveby(Point2D);
    Point2D rotateby(Matrix2);

};


class Point3D {
    bool _valid;
    double _x, _y, _z;
public:
    Point3D();
    Point3D(double, double, double, bool);
    Point3D(Point2D, double, bool);
    ~Point3D();
    bool is_valid();
    double x();
    double y();
    double z();
    Point3D moveby(double, double, double);
    double distance(Point3D);
    bool is_coincident(Point3D);
};


class Vector3D {

    double _x, _y, _z;

public:

    Vector3D(double, double, double);
    Vector3D(Point3D, Point3D);
    ~Vector3D();

    double x();
    double y();
    double z();

    double length();
    Vector3D scale(double);
    Vector3D versor();
    double scalar_prod(Vector3D);
    Vector3D vector_prod(Vector3D);
    double angle(Vector3D);
    bool isodirection(Vector3D);
    //Point3D intersect_coplanar(Vector3D);
    Point3D move_pt(Point3D);

};


class Line3D {

    Point3D _orig_pt3d;
    Vector3D _vers;

public:

    Line3D(Point3D, Vector3D);
    ~Line3D();
    bool isparallel(Line3D);
    bool iscoincident(Line3D);
    Point3D intersect_coplanar(Line3D);
    Point3D orig_pt();
    Vector3D versor();

};


class Segment3D {

    Point3D _start_pt, _end_pt;

public:

    Segment3D(Point3D, Point3D);
    ~Segment3D();
    Point3D start_pt();
    Point3D end_pt();
    double dx();
    double dy();
    double dz();
    Vector3D as_vector();
    Vector3D as_versor();
    double length();
    bool is_point_projection_in_segment(Point3D);
    Line3D as_line();

};


class CartesianPlane {
    /*
    Cartesian plane, expressed by normal equation:
    ax + by + cz + d = 0,
    where a, b and c are the normal vector directional cosines
    and abs(d) is the distance from reference frame center
    */

    double _a, _b, _c, _d;

public:

    //CartesianPlane(double, double, double, double);
    CartesianPlane();

    std::vector<double> define_params(Point3D, Point3D, Point3D);

    CartesianPlane(Point3D, Point3D, Point3D);
    ~CartesianPlane();
    void set_params(Point3D, Point3D, Point3D);

    double a();
    double b();
    double c();
    double d();

    Vector3D normal_versor();
    Vector3D perp_versor_in_plane(Vector3D);
    Vector3D intersect_versor(CartesianPlane);
    Point3D intersect_point3d(CartesianPlane);
    Line3D intersect(CartesianPlane);

    double point_distance(Point3D);
    bool point_in_plane(Point3D);
    double angle(CartesianPlane);
    bool isparallel(CartesianPlane);
    bool isequidistant(CartesianPlane);

};


class GeologicalPlane {

    double _dipdir, _dipangle;

public:

    GeologicalPlane();
    GeologicalPlane(double, double);
    ~GeologicalPlane();

    double dipdir();
    double dipangle();

};


class GeolAxis {

    /*
    Structural axis,
    defined by trend and plunge (both in degrees)
    Trend range: [0.0, 360.0[ clockwise, from 0 (North)
    Plunge: [-90.0, 90.0], negative value: upward axis, positive values: downward axis
    */

    double _trend, _plunge;

public:

    GeolAxis(double, double);
    ~GeolAxis();

    double trend();
    double plunge();

    GeolAxis to_down_axis();
    GeologicalPlane normal_geolplane();

};


GeologicalPlane to_geolplane(CartesianPlane);


GeolAxis to_geol_axis(Vector3D);


#endif
