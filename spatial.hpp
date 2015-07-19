#include <vector>
#include <cmath>


#include "algebr.hpp"


#ifndef SPATIAL_HPP
#define SPATIAL_HPP


class Range1D {
    double r_min, r_max;

public:
    bool within(double val);
    bool intersects(Range1D);

};


class Point2D {
    double _x, _y;
public:
    //Point2D();
    Point2D();
    Point2D(double, double);
    ~Point2D();
    double x();
    double y();
    Point2D operator=(Point2D);
    double distance(Point2D);
};


class Point3D {
    double _x, _y, _z;
public:
    Point3D(double, double, double);
    ~Point3D();
    double x();
    double y();
    double z();
    double distance(Point3D);
    bool is_coincident(Point3D);
};


class Vector3D {

    double _x, _y, _z;

public:

    Vector3D( double, double, double);
    Vector3D( Point3D, Point3D);
    ~Vector3D();

    double x();
    double y();
    double z();

    double length();
    Vector3D scale(double);
    Vector3D versor();
    double scalar_prod(Vector3D);
    Vector3D vector_prod(Vector3D);
    bool is_colinear(Vector3D);
    Point3D intersect_coplanar(Vector3D);
    Point3D move_pt(Point3D);


};


class Line3D {

    Point3D _orig_pt3d;
    Vector3D _vers;

public:

    Line3D(Point3D, Vector3D);
    ~Line3D();
    Point3D intersect_coplanar(Line3D);

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
    double length();
    bool is_point_projection_in_segment(Point3D);
    Line3D as_line();

};


class CartesianPlane {
    /*
    Cartesian plane, expressed by equation:
    ax + by + cz + d = 0
    */

    double _a, _b, _c, _d;

public:

    CartesianPlane(Point3D, Point3D, Point3D);
    CartesianPlane(double, double, double, double);
    ~CartesianPlane();

    Vector3D normal_versor();
    Vector3D perp_versor_in_plane(Vector3D);
    Vector3D intersect_versor(CartesianPlane);
    Point3D intersect_point3d(CartesianPlane);
    Line3D intersect(CartesianPlane);

};


class Triangle3D {

    Point3D _pt_a, _pt_b, _pt_c;

public:

    Triangle3D(Point3D, Point3D, Point3D);
    ~Triangle3D();
    CartesianPlane to_cartes_plane();

};


#endif
