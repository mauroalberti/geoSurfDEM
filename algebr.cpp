
#include "algebr.hpp"


Matrix2::Matrix2(double m11, double m12,double m21, double m22) {

    _matr[0][0] = m11;
    _matr[0][1] = m12;

    _matr[1][0] = m21;
    _matr[1][1] = m22;

};


Matrix2::~Matrix2() {
};


double Matrix2::determinant() {

    double m11 = _matr[0][0];
    double m12 = _matr[0][1];

    double m21 = _matr[1][0];
    double m22 = _matr[1][1];

    return m11*m22 - m12*m21;

};


Matrix3::Matrix3(double m11, double m12, double m13, double m21, double m22, double m23, double m31, double m32, double m33) {

    _matr[0][0] = m11;
    _matr[0][1] = m12;
    _matr[0][2] = m13;

    _matr[1][0] = m21;
    _matr[1][1] = m22;
    _matr[1][2] = m23;

    _matr[2][0] = m31;
    _matr[2][1] = m32;
    _matr[2][2] = m33;

};


Matrix3::~Matrix3() {
};


double Matrix3::determinant() {

    double m11 = _matr[0][0];
    double m12 = _matr[0][1];
    double m13 = _matr[0][2];

    double m21 = _matr[1][0];
    double m22 = _matr[1][1];
    double m23 = _matr[1][2];

    double m31 = _matr[2][0];
    double m32 = _matr[2][1];
    double m33 = _matr[2][2];

    double minor_11 = m22*m33 - m23*m32;
    double minor_12 = m21*m33 - m23*m31;
    double minor_13 = m21*m32 - m22*m31;

    double cofactor_11 = m11;
    double cofactor_12 = -m12;
    double cofactor_13 = m13;

    return cofactor_11*minor_11 + cofactor_12*minor_12 + cofactor_13*minor_13;

};



