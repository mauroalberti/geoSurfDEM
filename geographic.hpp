
#include "spatial.hpp"


class Space3DPartition {

    Range1D range_x;
    Range1D range_y;
    Range1D range_z;

public:
    Space3DPartition();
    Space3DPartition(const Range1D&, const Range1D&, const Range1D&);
    ~Space3DPartition();
    Range1D x_range();
    Range1D y_range();
    Range1D z_range();
    bool intersects(const Space3DPartition&);

};


class RectangularDomain {

    Point2D pt2d;  // generation point defined by x and y value
    double rot_alpha_degr;  // rotation angle Alpha, in degrees - with respect to x-axis
    double l_size, m_size;  // l, m domain sizes

public:

    RectangularDomain();
    RectangularDomain(const Point2D&, const double&, const double&, const double&);
    ~RectangularDomain();
    Point2D pt();
    double rot_angle();
    double l();
    double m();
    Range1D range_x();
    Range1D range_y();

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
    RectangularDomain rr_domain();
    unsigned int cols();
    unsigned int rows();

};


class DataRRGrid {

    RectRegularGrid rrgrid;
    NumericData data_vals;

public:
    DataRRGrid();
    DataRRGrid(RectRegularGrid, std::vector<double>, double);
    ~DataRRGrid();
    RectRegularGrid rr_grid();
    NumericData data();
    Space3DPartition space_partition();

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
