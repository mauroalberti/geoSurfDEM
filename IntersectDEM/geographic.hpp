
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
    bool intersects(Space3DPartition);

};


class Triangle3D {

    Point3D _pts [3];

public:

    Triangle3D();
    Triangle3D(Point3D, Point3D, Point3D);
    ~Triangle3D();
    Point3D pt( unsigned int);
    Triangle3D shift(double, double, double);
    double area();
    CartesianPlane to_cartes_plane();
    Space3DPartition space_volume();
    double point_distance(Point3D);

};


class RectangularDomain {

    Point2D pt2d;  // generation point defined by x and y value
    double rot_alpha_degr;  // rotation angle Alpha, in degrees - with respect to x-axis
    double l_size, m_size; // l, m domain sizes
    int n_cols, n_rows;

public:

    RectangularDomain();
    RectangularDomain(const Point2D&, const double&, const double&, const double&, const unsigned int&, const unsigned int&);
    ~RectangularDomain();
    Point2D pt();
    double rot_angle();
    double l();
    double m();
    unsigned int ncols();
    unsigned int nrows();
    Range1D range_x();
    Range1D range_y();

};


class DataRRGrid {

    RectangularDomain rdomain;
    NumericData data_vals;

public:
    DataRRGrid();
    DataRRGrid(RectangularDomain, std::vector<double>, double);
    ~DataRRGrid();
    RectangularDomain rect_domain();
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
