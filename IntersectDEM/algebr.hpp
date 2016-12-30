#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>


double radians(const double);
double degrees(const double);


class Range1D {

    double r_min, r_max;

public:
    Range1D();
    Range1D(double, double);
    ~Range1D();
    double start();
    double end();
    bool within(double val);
    bool intersects(Range1D);

};


class Matrix2 {

    double _matr[2][2];

public:

    Matrix2(double, double, double, double);
    ~Matrix2();
    double m11();
    double m12();
    double m21();
    double m22();
    double determinant();
};


Matrix2 angle_to_rot_matrix( double);


class Matrix3 {

    double _matr[3][3];

public:

    Matrix3(double, double, double, double, double, double, double, double, double);
    ~Matrix3();
    double determinant();
};


class NumericData {

    std::vector<double> _values;
    double _nodata_value;

public:
    NumericData();
    NumericData(std::vector<double>, double);
    ~NumericData();
    std::vector<double> values();
    double nodata_value();
    double min();
    double max();
    Range1D range();

};
