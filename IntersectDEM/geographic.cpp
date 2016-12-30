#include <iostream>

#include "geographic.hpp"



Space3DPartition::Space3DPartition() {

};


Space3DPartition::Space3DPartition(const Range1D& r_x, const Range1D& r_y, const Range1D& r_z) :
    range_x(r_x), range_y(r_y), range_z(r_z)  {
    };


Space3DPartition::~Space3DPartition() {

};


Range1D Space3DPartition::x_range() {

    return range_x;

};


Range1D Space3DPartition::y_range() {

    return range_y;

};


Range1D Space3DPartition::z_range() {

    return range_z;

};


bool Space3DPartition::intersects(Space3DPartition another) {

    if (x_range().intersects(another.x_range()) and
        y_range().intersects(another.y_range()) and
        z_range().intersects(another.z_range()) ) {
            return true; }
    else {
        return false; };

};


Triangle3D::Triangle3D() {

};


Triangle3D::Triangle3D(Point3D pt_a, Point3D pt_b, Point3D pt_c) {

    _pts[0] = pt_a;
    _pts[1] = pt_b;
    _pts[2] = pt_c;

};


Triangle3D::~Triangle3D() {

};


Point3D Triangle3D::pt( unsigned int i) {

    return _pts[i];

};


Triangle3D Triangle3D::shift(double dx, double dy, double dz) {

    Point3D orig_pt1 = pt(0);
    Point3D orig_pt2 = pt(1);
    Point3D orig_pt3 = pt(2);

    Point3D shifted_pt1 = orig_pt1.moveby(dx, dy,dz);
    Point3D shifted_pt2 = orig_pt2.moveby(dx, dy,dz);
    Point3D shifted_pt3 = orig_pt3.moveby(dx, dy,dz);

    return Triangle3D(shifted_pt1,
                      shifted_pt2,
                      shifted_pt3);

};


double Triangle3D::area() {

    Vector3D vect1 = Vector3D(pt(0), pt(1));
    Vector3D vect2 = Vector3D(pt(0), pt(2));

    Vector3D vect_prod = vect1.vector_prod(vect2);

    return vect_prod.length()/2.0;

};


CartesianPlane Triangle3D::to_cartes_plane() {

    CartesianPlane cart_plane = CartesianPlane(pt(0), pt(1), pt(2));
    //std::cout << "to cartesian plane: " << cart_plane.a() << "," << cart_plane.b() << "," << cart_plane.c() << "," << cart_plane.d() << "\n";
    return cart_plane;

};


Space3DPartition Triangle3D::space_volume() {

    double xs[] = {pt(0).x(), pt(1).x(), pt(2).x()};
    double ys[] = {pt(0).y(), pt(1).y(), pt(2).y()};
    double zs[] = {pt(0).z(), pt(1).z(), pt(2).z()};

    double min_x = *std::min_element(xs,xs+3);
    double max_x = *std::max_element(xs,xs+3);

    double min_y = *std::min_element(ys,ys+3);
    double max_y = *std::max_element(ys,ys+3);

    double min_z = *std::min_element(zs,zs+3);
    double max_z = *std::max_element(zs,zs+3);

    Range1D range_x = Range1D(min_x, max_x);
    Range1D range_y = Range1D(min_y, max_y);
    Range1D range_z = Range1D(min_z, max_z);

    return Space3DPartition(range_x, range_y, range_z);

};


double Triangle3D::point_distance(Point3D pt) {

    CartesianPlane cartes_plane = to_cartes_plane();
    return cartes_plane.point_distance(pt);
};


RectangularDomain::RectangularDomain() {
};


RectangularDomain::RectangularDomain(const Point2D& pt2d_, const double& rot_alpha_degr_, const double& l_size_, const double& m_size_, const unsigned int& ncols_, const unsigned int& nrows_) :
pt2d(pt2d_), rot_alpha_degr(rot_alpha_degr_), l_size(l_size_), m_size(m_size_), n_cols(ncols_), n_rows(nrows_) {
};


RectangularDomain::~RectangularDomain() {

};


Point2D RectangularDomain::pt() {

    return pt2d;

};


double RectangularDomain::rot_angle() {

    return rot_alpha_degr;

};


double RectangularDomain::l() {

    return l_size;

};


double RectangularDomain::m() {

    return m_size;

};


unsigned int RectangularDomain::ncols() {

    return n_cols;

};


unsigned int RectangularDomain::nrows() {

    return n_rows;

};


Range1D RectangularDomain::range_x() {

    double x_0 = pt().x();
    double alpha_rad = radians(rot_alpha_degr);
    double x_min = x_0 - nrows()*m()*sin(alpha_rad);
    double x_max = x_0 + ncols()*l()*cos(alpha_rad);

    return Range1D(x_min, x_max);

};


Range1D RectangularDomain::range_y() {

    double y_0 = pt().y();
    double alpha_rad = radians(rot_alpha_degr);
    double y_max = y_0 + ncols()*l()*sin(alpha_rad) + nrows()*m()*cos(alpha_rad);


    return Range1D(y_0, y_max);

};



DataRRGrid::DataRRGrid() {

};


DataRRGrid::DataRRGrid(RectangularDomain _rdomain, std::vector<double> _data, double _nodata_value):
        rdomain(_rdomain), data_vals(NumericData(_data, _nodata_value))
    {

};


DataRRGrid::~DataRRGrid() {

};


RectangularDomain DataRRGrid::rect_domain() {

    return rdomain;

};


NumericData DataRRGrid::data() {

    return data_vals;

};


Space3DPartition DataRRGrid::space_partition() {

    Range1D x_range = rect_domain().range_x();
    Range1D y_range = rect_domain().range_y();
    Range1D z_range = data().range();

    return Space3DPartition(x_range,
                            y_range,
                            z_range);

};


MeshTriangleStrip::MeshTriangleStrip() {

    };


MeshTriangleStrip::MeshTriangleStrip(const std::vector<Point3D>& pts, const std::vector<std::vector<unsigned int> >& trianglestrips) :
    _pts(pts), _trianglestrips(trianglestrips) {
    };


MeshTriangleStrip::~MeshTriangleStrip() {

};


std::vector<Point3D> MeshTriangleStrip::pts() {

    return _pts;

};


std::vector<std::vector<unsigned int> > MeshTriangleStrip::trianglestrips(){

    return _trianglestrips;

};


PointTriplet::PointTriplet() {

    triplet[0] = Point3D();
    triplet[1] = Point3D();
    triplet[2] = Point3D();

};


PointTriplet::PointTriplet(Point3D pt1, Point3D pt2, Point3D pt3) {

    triplet[0] = pt1;
    triplet[1] = pt2;
    triplet[2] = pt3;

};


Point3D PointTriplet::get(unsigned int i) {

    return triplet[i];

};


unsigned int PointTriplet::valid_pts() {

    unsigned int num_valid_pts = 0;

    for (int i = 0; i < 3; i++ ) {
        Point3D curr_pt = get(i);
        if (curr_pt.is_valid()) {
            num_valid_pts++;
        }
    };

    return num_valid_pts;

};


bool PointTriplet::is_valid() {

    if (valid_pts() == 3) {
        return true;
    } else {
        return false;
    };
};


PointTriplet PointTriplet::update(Point3D pt4) {

    return PointTriplet(get(1), get(2), pt4);

};


Triangle3D PointTriplet::to_triangle() {

    return Triangle3D(triplet[0], triplet[1], triplet[2]);

};



