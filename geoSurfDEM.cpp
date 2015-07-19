#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <vector>

#include <stdlib.h>     /* atof */



#include "geographic.hpp"



int get_dem_param_int(const std::string& rec_line, const std::string& par_name) {
    // read single key-uint value from row - ESRI ASCII DEM format
    unsigned int par_val;
    std::string str_raw, str;
    std::istringstream instr(rec_line);
    instr >> str_raw >> par_val;
    transform(str_raw.begin(), str_raw.end(), back_inserter(str), ::toupper);
    if (str != par_name)
        throw -1;
    return par_val; }

double get_dem_param_double(const std::string& rec_line, const std::string& par_name) {
    // read single key-double value from row - ESRI ASCII DEM format
    double par_val;
    std::string str_raw, str;
    std::istringstream instr(rec_line);
    instr >> str_raw >> par_val;
    transform(str_raw.begin(), str_raw.end(), back_inserter(str), ::toupper);
    if (str != par_name)
        throw -1;
    return par_val; }

DataRRGrid read_esri_ascii_dem( std::string dem_filepath ) {
    // read DEM in ESRI Arc/Info ASCII format

    unsigned int ncols, nrows;
    float xllcorner, yllcorner, cellsize, nodata_value;

    std::ifstream infile;
    infile.open(dem_filepath.c_str(), std::ios::binary);

    if (infile.fail() ) {
        std::cout << "\nInput file " << dem_filepath << " not available\n";
        throw -1; }

    std::string rec_line;

    // get/check NCOLS
    getline(infile, rec_line);
    try {
        ncols = get_dem_param_int(rec_line, "NCOLS"); }
    catch (int e) {
        std::cout << "\nError with number of columns in input file first line\n";
        throw -1; }

    // get/check NROWS
    getline(infile, rec_line);
    try {
        nrows = get_dem_param_int(rec_line, "NROWS"); }
    catch (int e) {
        std::cout << "\nError with number of rows in input file second line\n";
        throw -1; }

    // get/check XLLCORNER
    getline(infile, rec_line);
    try {
        xllcorner = get_dem_param_double(rec_line, "XLLCORNER"); }
    catch (int e) {
        std::cout << "\nError with number of LL x corner in input file third line\n";
        throw -1; }

    // get/check YLLCORNER
    getline(infile, rec_line);
    try {
        yllcorner = get_dem_param_double(rec_line, "YLLCORNER"); }
    catch (int e) {
        std::cout << "\nError with number of LL y corner in input file fourth line\n";
        throw -1; }

    // get/check CELLSIZE
    getline(infile, rec_line);
    try {
        cellsize = get_dem_param_double(rec_line, "CELLSIZE"); }
    catch (int e) {
        std::cout << "\nError with number of cell size value in input file fifth line\n";
        throw -1; }

    // get/check NODATA_VALUE
    getline(infile, rec_line);
    try {
        nodata_value = get_dem_param_double(rec_line, "NODATA_VALUE"); }
    catch (int e) {
        std::cout << "\nError with number of null data value in input file sixth line\n";
        throw -1; }

    // read DEM data values as a vector<double>
    std::vector<double> data_vect;
    while (!infile.eof()) {
        std::string data_line;
        char split_char = ' ';
        std::getline(infile, data_line);
        std::istringstream instr(data_line);
        for (std::string each; std::getline(instr, each, split_char); data_vect.push_back(atof(each.c_str()))); };

    infile.close();

    if (ncols*nrows != data_vect.size()) {
        std::cout << "Error in data count\n";
        throw -1; }

    Point2D pt2d {xllcorner, yllcorner};

    RectangularDomain domain { pt2d, 0.0, cellsize*ncols, cellsize*nrows };

    RectRegularGrid rrgrid { domain, ncols, nrows };

    return DataRRGrid(rrgrid, data_vect, nodata_value);

};

MeshTriangleStrip read_vtk_data_ascii( std::string input_vtk_path ) {

    unsigned int num_points, num_triangle_strips, num_dummy;

    std::ifstream infile;
    infile.open(input_vtk_path.c_str(), std::ios::binary);

    if (infile.fail() ) {
        std::cout << "\nInput file " << input_vtk_path << " not available\n";
        throw -1; }

    // get number of points
    std::string rec_line;
    bool found = false;
    do {
        std::getline(infile, rec_line);
        std::string prefix("POINTS");
        if (rec_line.compare(0, prefix.size(), prefix) == 0) {
            found = true;
            std::string dummy, dummy2;
            std::istringstream instr(rec_line);
            instr >> dummy >> num_points >> dummy2; }
    } while ( ! found );

    std::cout << "num points is " << num_points << "\n";

    std::vector<Point3D> pt3d_vect;
    for (unsigned int n = 0; n < num_points; n++ ) {
        std::getline(infile, rec_line);
        std::istringstream instr(rec_line);
        double x, y, z;
        instr >> x >> y >> z;
        pt3d_vect.push_back( Point3D(x, y, z));
        };

    std::getline(infile, rec_line);

    // TRIANGLE_STRIPS
    std::string prefix("TRIANGLE_STRIPS");

    std::getline(infile, rec_line);
    if (rec_line.compare(0, prefix.size(), prefix) != 0) {
        std::cout << "Error in expected TRIANGLE_STRIPS header line\n";
        throw -1;
    };

    std::istringstream instr(rec_line);
    std::string dummy_str;
    instr >> dummy_str >> num_triangle_strips >> num_dummy;

    std::vector<std::vector<unsigned int> > triangle_strips;
    std::string blank_space = " ";
    for (unsigned int n = 0; n < num_triangle_strips; n++) {

        std::getline(infile, rec_line);

        std::vector<unsigned int> triangle_strip;
        std::size_t found;
        bool line_end = false;
        while ( ! line_end ){
            found = rec_line.find(blank_space);
            if (found != std::string::npos) {
                std::string number = rec_line.substr(0, found);
                triangle_strip.push_back( atoi(number.c_str()) );
                rec_line.erase(0, found+1); }
            else {
                line_end = true;  };   };

        triangle_strips.push_back( triangle_strip );
    };

    infile.close();

    return MeshTriangleStrip(pt3d_vect, triangle_strips);

};


// ciclo sui punti di interesse dell'array del DEM

// definisce doppietta, terzetto o quartetto di punti di interesse

// processamento su unità mesh (triangolo)

// ?fast boundaries (parallelepipedo rettangolo) su DEM triangle

// ?fast boundaries su GAS triangle

// ? do they intersect ?

// check if intersection
// if intersection: set intersection params: local surface attitude


int main() {

    // read DEM data from input file
    std::string input_dem_path = "./test_data/malpi_w4u2n_src.asc";
    DataRRGrid datagrid;
    try {
        datagrid = read_esri_ascii_dem( input_dem_path ); }
    catch (int e) {
        std::cout << "Program will stop\n";
        return -1; }

    // read VTK data from input file
    std::string input_vtk_path = "./test_data/surf3d_sim_01_rot45_04500.vtk";
    MeshTriangleStrip surf3d_mesh;
    try {
        surf3d_mesh = read_vtk_data_ascii( input_vtk_path ); }
    catch (int e) {
        std::cout << "Program will stop\n";
        return -1; }


    std::vector<std::vector<unsigned int> > triangle_strips = surf3d_mesh.trianglestrips();
    unsigned int ndx_strip = 0;
    for(std::vector<std::vector<unsigned int> >::iterator ref_trstr = triangle_strips.begin(); ref_trstr != triangle_strips.end(); ++ref_trstr) {
        std::vector<unsigned int> triangle_strip = *ref_trstr;
        for(std::vector<unsigned int>::iterator ref_ptndx = triangle_strip.begin(); ref_ptndx != triangle_strip.end(); ++ref_ptndx) {
            std::cout << ndx_strip << ": " << *ref_ptndx << "\n";
        }
        ndx_strip++;
    };


    return 0;
};
