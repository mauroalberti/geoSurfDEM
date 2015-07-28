#include <vector>
#include <algorithm>

#include <iostream>


class Matrix2 {

    double _matr[2][2];

public:

    Matrix2(double, double, double, double);
    ~Matrix2();
    double determinant();
};


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

};
