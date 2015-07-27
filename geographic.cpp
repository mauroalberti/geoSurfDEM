#include <iostream>

#include "geographic.hpp"



bool Space3DPartition::intersects(const Space3DPartition& another) {

    if (range_x.intersects(another.range_x) and
        range_y.intersects(another.range_y) and
        range_z.intersects(another.range_z) ) {
            return true;
        }
    else {
        return false;
        };

};


UnrotatedRectangularDomain::UnrotatedRectangularDomain() {
};


UnrotatedRectangularDomain::UnrotatedRectangularDomain(const Point2D& pt2d_, const double& l_size_, const double& m_size_) :
pt2d(pt2d_), l_size(l_size_), m_size(m_size_) {
};


UnrotatedRectangularDomain::~UnrotatedRectangularDomain() {
};


Range1D UnrotatedRectangularDomain::get_x_range() {

    return Range1D(pt2d.x(), pt2d.x() + l_size);
};


Range1D UnrotatedRectangularDomain::get_y_range() {

    return Range1D(pt2d.y(), pt2d.y() + m_size);
};


RectangularDomain::RectangularDomain() {
};


RectangularDomain::RectangularDomain(const Point2D& pt2d_, const double& rot_alpha_degr_, const double& l_size_, const double& m_size_) :
pt2d(pt2d_), rot_alpha_degr(rot_alpha_degr_), l_size(l_size_), m_size(m_size_) {
};


RectangularDomain::~RectangularDomain() {
};


RectRegularGrid::RectRegularGrid() {
};


RectRegularGrid::RectRegularGrid(const RectangularDomain& domain_, const unsigned int& ncols_, const unsigned int& nrows_) :
domain(domain_), ncols(ncols_), nrows(nrows_) {
};


RectRegularGrid::~RectRegularGrid() {
};


DataRRGrid::DataRRGrid() {
};


DataRRGrid::DataRRGrid(const RectRegularGrid& rrgrid, const std::vector<double> data, const double nodata_value) :
    _rrgrid(rrgrid), _data(data), _nodata_value(nodata_value) {
};


DataRRGrid::~DataRRGrid() {
};


RectRegularGrid DataRRGrid::rrgrid() {
    return _rrgrid;
};



std::vector<double> DataRRGrid::data() {
    return _data;
};


double DataRRGrid::nodata_value() {
    return _nodata_value;
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



