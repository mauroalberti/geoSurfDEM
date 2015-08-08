
#include "algebr.hpp"


const double PI  = 3.141592653589793238463;


double radians(const double degrees) {

    return degrees*PI/180.0;

};


double degrees(const double radians) {

    return radians*180.0/PI;

};


Range1D::Range1D() {

};


Range1D::Range1D(double val_start, double val_end) : r_min(val_start), r_max(val_end) {

};


Range1D::~Range1D() {

};


double Range1D::start() {

    return r_min;

};


double Range1D::end() {

    return r_max;

};


bool Range1D::within(double val) {

    //std::cout << "within " << start() << " - " <<  val << " - " << end();

    if (start() <= val and val <= end()) {
        //std::cout << "within " << start() << " - " <<  val << " - " << end() << "  -> YS\n";
        return true; }
    else {
        return false; }

};


bool Range1D::intersects(Range1D another) {

    if (within(another.start()) or within(another.end()) or
        another.within( start()) or another.within( end() ) ) {
        return true; }
    else {
        return false; };

};


Matrix2::Matrix2(double m11, double m12,double m21, double m22) {

    _matr[0][0] = m11;
    _matr[0][1] = m12;

    _matr[1][0] = m21;
    _matr[1][1] = m22;

};


Matrix2::~Matrix2() {

};


double Matrix2::m11() {

    return _matr[0][0];

};


double Matrix2::m12() {

    return _matr[0][1];

};


double Matrix2::m21() {

    return _matr[1][0];

};


double Matrix2::m22() {

    return _matr[1][1];

};


double Matrix2::determinant() {

    double m11 = _matr[0][0];
    double m12 = _matr[0][1];

    double m21 = _matr[1][0];
    double m22 = _matr[1][1];

    return m11*m22 - m12*m21;

};


Matrix2 angle_to_rot_matrix(double rot_angle_degr) {

    double omega = radians(rot_angle_degr);

    double sin_omega = std::sin(omega);
    double cos_omega = std::cos(omega);

    double rm11 = cos_omega;
    double rm12 = -sin_omega;
    double rm21 = sin_omega;
    double rm22 = cos_omega;

    return Matrix2(rm11, rm12, rm21, rm22);

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


NumericData::NumericData() {

};


NumericData::NumericData(std::vector<double> vals, double nodata_val) :
    _values(vals), _nodata_value(nodata_val) {

};


NumericData::~NumericData() {

};


std::vector<double> NumericData::values() {

    return _values;

};


double NumericData::nodata_value() {

    return _nodata_value;

};


std::vector<double> filter_data(std::vector<double> indata, double filter_out_value) {

    std::vector<double> bar (indata.size());

    auto it = std::copy_if(indata.begin(), indata.end(), bar.begin(), [filter_out_value](double d){return (fabs(d-filter_out_value) > 1.0e-6);} );
    bar.resize(std::distance(bar.begin(),it));  // shrink container to new size

    return bar;

};


double NumericData::min() {

    std::vector<double> foo = filter_data(values(), nodata_value());
    return *std::min_element(foo.begin(), foo.end());

};


double NumericData::max() {

    std::vector<double> foo = filter_data(values(), nodata_value());
    //std::cout << "foo size " << foo.size() << "\n";
    return *std::max_element(foo.begin(), foo.end());

};


Range1D NumericData::range() {

    return Range1D(min(), max());

};



