
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




