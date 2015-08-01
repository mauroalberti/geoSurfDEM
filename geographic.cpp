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


RectangularDomain::RectangularDomain() {
};


RectangularDomain::RectangularDomain(const Point2D& pt2d_, const double& rot_alpha_degr_, const double& l_size_, const double& m_size_) :
pt2d(pt2d_), rot_alpha_degr(rot_alpha_degr_), l_size(l_size_), m_size(m_size_) {
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


Range1D RectangularDomain::range_x() {

    double x_0 = pt().x();
    //std::cout << "range_x x_0 " << x_0 << "\n";
    double alpha_rad = radians(rot_alpha_degr);
    //std::cout << "range_x alpha_rad " << alpha_rad << "\n";
    double x_min = x_0 - m()*sin(alpha_rad);
    //std::cout << "range_x x_min " << x_min << "\n";
    double x_max = x_0 + l()*cos(alpha_rad);
    //std::cout << "range_x x_max " << x_max << "\n";

    return Range1D(x_min, x_max);

};


Range1D RectangularDomain::range_y() {

    double y_0 = pt().y();
    double alpha_rad = radians(rot_alpha_degr);
    double y_max = y_0 + l()*sin(alpha_rad) + m()*cos(alpha_rad);

    return Range1D(y_0, y_max);

};



RectRegularGrid::RectRegularGrid() {
};


RectRegularGrid::RectRegularGrid(const RectangularDomain& domain_, const unsigned int& ncols_, const unsigned int& nrows_) :
domain(domain_), ncols(ncols_), nrows(nrows_) {

};


RectRegularGrid::~RectRegularGrid() {

};


RectangularDomain RectRegularGrid::rr_domain() {

    return domain;

};


unsigned int RectRegularGrid::cols() {

    return ncols;

};


unsigned int RectRegularGrid::rows() {

    return nrows;

};



DataRRGrid::DataRRGrid() {

};


DataRRGrid::DataRRGrid(RectRegularGrid _rrgrid, std::vector<double> _data, double _nodata_value):
        rrgrid(_rrgrid), data_vals(NumericData(_data, _nodata_value))
    {

    //rrgrid = RectRegularGrid(_rrgrid.rr_domain(), _rrgrid.cols(), _rrgrid.rows() );
    //data_vals = NumericData(_data, _nodata_value);

    //std::cout << "DataRRGrid _rrgrid domain pt2d " << _rrgrid.rr_domain().pt().x() << " " << _rrgrid.rr_domain().pt().y() << "\n";
    //std::cout << "DataRRGrid rrgrid domain pt2d " << rrgrid.rr_domain().pt().x() << " " << rrgrid.rr_domain().pt().y() << "\n";
    // OK

};


DataRRGrid::~DataRRGrid() {

};


RectRegularGrid DataRRGrid::rr_grid() {

    //std::cout << "DataRRGrid::rr_grid() domain pt2d " << rrgrid.rr_domain().pt().x() << " " << rrgrid.rr_domain().pt().y() << "\n";

    return rrgrid;

};



NumericData DataRRGrid::data() {

    return data_vals;

};


Space3DPartition DataRRGrid::space_partition() {

    Range1D x_range = rr_grid().rr_domain().range_x();
    Range1D y_range = rr_grid().rr_domain().range_y();
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



