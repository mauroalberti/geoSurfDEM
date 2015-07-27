
#include "spatial.hpp"


class Space3DPartition {
    Range1D range_x;
    Range1D range_y;
    Range1D range_z;

public:
    bool intersects(const Space3DPartition&);

};


class UnrotatedRectangularDomain {
    Point2D pt2d;  // generation point defined by x and y value
    double l_size, m_size;  // l, m domain sizes
public:
    UnrotatedRectangularDomain();
    UnrotatedRectangularDomain(const Point2D&, const double&, const double&);
    ~UnrotatedRectangularDomain();
    Range1D get_x_range();
    Range1D get_y_range();
};


class RectangularDomain {
    Point2D pt2d;  // generation point defined by x and y value
    double rot_alpha_degr;  // rotation angle Alpha, in degrees - with respect to x-axis
    double l_size, m_size;  // l, m domain sizes
public:
    RectangularDomain();
    RectangularDomain(const Point2D&, const double&, const double&, const double&);
    ~RectangularDomain();

};


class RectRegularGrid {
    // ? planar CRS - starting with undefined CRS
    RectangularDomain domain; // rectangular domain
    unsigned int ncols, nrows;  // cel size/number, along l and m directions
    //array // 2D array
public:
    RectRegularGrid();
    RectRegularGrid(const RectangularDomain&, const unsigned int&, const unsigned int&);
    ~RectRegularGrid();
};


class DataRRGrid {

    RectRegularGrid _rrgrid;
    std::vector<double> _data;
    double _nodata_value;

public:
    DataRRGrid();
    DataRRGrid(const RectRegularGrid&, const std::vector<double>, const double);
    ~DataRRGrid();
    RectRegularGrid rrgrid();
    std::vector<double> data();
    double nodata_value();

};


class MeshTriangleStrip {

    std::vector<Point3D> _pts;
    std::vector<std::vector<unsigned int> > _trianglestrips;

public:
    MeshTriangleStrip();
    MeshTriangleStrip(const std::vector<Point3D>&, const std::vector<std::vector<unsigned int> >& );
    ~MeshTriangleStrip();

    std::vector<Point3D> pts();
    std::vector<std::vector<unsigned int> > trianglestrips();

};


class PointTriplet {

    Point3D triplet[3];

public:
    PointTriplet();
    PointTriplet(Point3D, Point3D, Point3D);
    PointTriplet update(Point3D);
    Point3D get(unsigned int i);
    unsigned int valid_pts();
    bool is_valid();
    Triangle3D to_triangle();

};
