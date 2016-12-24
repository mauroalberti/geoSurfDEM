
/*
#
# best_fit_geoplanes.cpp
# by Mauro Alberti - www.malg.eu, alberti.m65@gmail.com
# start 2016-12-10
#
# This program or module is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 2 of the License, or
# version 3 of the License, or (at your option) any later version. It is
# provided for educational purposes and is distributed in the hope that
# it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See
# the GNU General Public License for more details.
#
*/


#include <list>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>  // std::min_element, std::max_element
#include <set>
#include <map>

using namespace std;


bool is_uint_digit (char c) {

    if (c >= '0' && c <='9')
        return true;
    else
        return false; };


bool is_float_digit (char c) {

    if (is_uint_digit(c) or (c == '.') or (c == '-'))
        return true;
    else
        return false; };


bool is_uint(string &str) {

    for (uint i = 0; i < str.size(); i++) {
        if (! is_uint_digit(str[i])) {
            return false; }; };
    return true; };


bool is_float(string &str) {

    for (uint i = 0; i < str.size(); i++) {
        if (! is_float_digit(str[i])) {
            return false; }; };
    return true; };


string read_string(string &line) {

    string readstring;
    istringstream ss(line);
    ss >> readstring;
    return readstring; };


int read_int(string &line) {

    string readstring;
    readstring = read_string(line);
    return std::stoi(readstring);

};

int read_float(string &line) {

    string readstring;
    readstring = read_string(line);
    return std::stof(readstring);

};

void read_params(std::stringstream &param_lines, string &xyz_data_fpth, uint &hdr_ln_num,
                 string &bfp_geoplanes_fpth, string &ptnum_grid_fpth, float &cell_size) {

    string readline;

    getline(param_lines, readline);
    xyz_data_fpth = read_string(readline);

    getline(param_lines, readline);
    readline.erase(readline.end()-1);
    if (is_uint(readline)) {
        hdr_ln_num = read_int(readline);}
    else {
        throw "Number of header lines is not an unsigned integer";};

    getline(param_lines, readline);
    bfp_geoplanes_fpth = read_string(readline);

    getline(param_lines, readline);
    ptnum_grid_fpth = read_string(readline);

    getline(param_lines, readline);
    readline.erase(readline.end()-1);
    if (is_float(readline)) {
        cell_size = read_float(readline);}
    else {
        throw "Cell size is not a float value";};

};



class Point {

private:

    double _x, _y, _z;

public:

    Point( double xval, double yval, double zval) {
       _x = xval;
       _y = yval;
       _z = zval; }

    double x() {
       return _x;}

    double y() {
       return _y;}

    double z() {
       return _z;}
};


int main () {

    //
    // user interaction processes
    //

    // initial screenshot

    std::cout << "\n\nBest Fit Geoplanes\n";
    std::cout << "by M. Alberti - www.malg.eu\n";
    std::cout << "2016-12-17\n\n\n";

    // input of parameter file name and file opening

    std::cout << "Enter input parameter file name: ";
  	std::string param_filename;
    std::cin >> param_filename;

    // start time

    clock_t start_time = clock();

    // input parameter file reading

    std::ifstream param_file(param_filename);
    std::stringstream param_lines;
    param_lines << param_file.rdbuf();

    // reading parameters

    string xyz_data_fpth; // path of input xyz file
    uint hdr_ln_num; // number of header lines in xyz file
    string bfp_geoplanes_fpth; // output geoplanes file path
    string ptnum_grid_fpth;  // output point number grid file path
    float cell_size;  // used cell size

    try {
        read_params(param_lines,
                    xyz_data_fpth,
                    hdr_ln_num,
                    bfp_geoplanes_fpth,
                    ptnum_grid_fpth,
                    cell_size); }
    catch (string error_msg) {
        std::cout << "Parameter input error: " << error_msg;
        return 1;};

    // printout read parameters

    cout << "\nInput parameters\n";
    cout << " - xyz data (input): " << xyz_data_fpth << "\n";
    cout << " - number of header lines in xyz data file: " << hdr_ln_num << "\n";
    cout << " - best fit geoplanes (output): " << bfp_geoplanes_fpth << "\n";
    cout << " - point number grid (output): " << ptnum_grid_fpth << "\n";
    cout << " - cell size: " << cell_size << "\n";

    //
    // internal processings
    //

    cout << "\n\n -- processing started ... please wait\n\n";

    // open xyz input file

    ifstream infile;
    infile.open(xyz_data_fpth.c_str(), ios::binary);
    if (infile.fail())
        {cout << "\n\nUnable to open input file '" << xyz_data_fpth << "'\n";
        return 1;}

    // create output geoplanes file

    ofstream outfile;
    outfile.open(bfp_geoplanes_fpth.c_str(), ios::binary);
    if (outfile.fail())
        {cout << "\n\nUnable to create output file '" << bfp_geoplanes_fpth << "'\n";
         return 1;}

    // create output point number grid file

    ofstream outgridfile;
    outgridfile.open(ptnum_grid_fpth.c_str(), ios::binary);
    if (outgridfile.fail())
        {cout << "\n\nUnable to create output grid file '" <<  ptnum_grid_fpth << "'\n";
         return 1;}

    // read input raw xyz file

    string rec_line;
    list<string> rawdata_list;  // list of raw data strings

    for (uint n = 0; n < hdr_ln_num; n++) {  // read  (and discard) file header
        getline(infile, rec_line);};

    while (! infile.eof()) {  // read raw data
        getline(infile, rec_line);
        if (rec_line.size() > 0) {
            rawdata_list.push_back(rec_line);}}
    infile.close();

    std::cout << rawdata_list.size() << " records read\n";

    // read input xyz data into three double vectors

    int num_recs = rawdata_list.size();  // total number of input points
    vector<double> x(num_recs), y(num_recs), z(num_recs);  // vectors storing the x, y and z coordinates for all points

    char sep;
    int ndx = 0;
    list<string>::iterator string_pos;  // string list iterator
    for (string_pos=rawdata_list.begin(); string_pos!=rawdata_list.end(); string_pos++)
       {istringstream instr(*string_pos);
        instr >> x[ndx] >> sep >> y[ndx] >> sep >> z[ndx];
        ndx++;}

    // calculate spatial ranges of x, y and z

    double x_min = *min_element( x.begin(), x.end() );
    double x_max = *max_element( x.begin(), x.end() );
    double y_min = *min_element( y.begin(), y.end() );
    double y_max = *max_element( y.begin(), y.end() );
    double z_min = *min_element( z.begin(), z.end() );
    double z_max = *max_element( z.begin(), z.end() );
    float x_range = x_max -  x_min;
    float y_range = y_max -  y_min;

    // calculate grid row and column numbers

    int rows = int(y_range/cell_size) + 1;
    int columns = int(x_range/cell_size) + 1;

    // print info on spatial ranges and grid params

    cout << "\nThe spatial range of data is:\n";
    cout << "  x min: " << x_min << " x max: " << x_max << "\n";
    cout << "  y min: " << y_min << " y max: " << y_max << "\n";
    cout << "  z min: " << z_min << " z max: " << z_max << "\n";
    cout << "\nRows number: " << rows << "\nColumn number: " << columns << "\n ";

    // calculate grid indices of input points, together with non-empty grid cells

    vector<int> pt_cellndx_i(num_recs), pt_cellndx_j(num_recs);

    for (int k = 0; k < num_recs; k++) {
        pt_cellndx_i[k] = int((x[k]- x_min)/cell_size);
        pt_cellndx_j[k] = int((y[k]- y_min)/cell_size);}

    // define grid linear indices of input points

    vector<int> pt_cellndx_n(num_recs);
    for (int k = 0; k < num_recs; k++) {
        pt_cellndx_n[k] = (columns * pt_cellndx_i[k]) + pt_cellndx_j[k];}

    // define grid linear indices of cells with at least one point, using a set data structure

    set<int> nonemptycell_ndxs_set;
    for (int k = 0; k < num_recs; k++) {
        nonemptycell_ndxs_set.insert(pt_cellndx_n[k]);}

    // initialize the mapping between each cell grid linear index and associated point indices

    map<int, vector<int> > gridcell_ptndxs_map;
    set<int>::const_iterator pos;
    for(pos = nonemptycell_ndxs_set.begin(); pos != nonemptycell_ndxs_set.end(); ++pos) {  // initialize the mapping between each cell grid linear index and associated point indices
        vector<int> empy_vector;
        gridcell_ptndxs_map[*pos] = empy_vector;}
    for (int k = 0; k < num_recs; k++) {  // insert the point indices into the grid cell linear indices
        gridcell_ptndxs_map[pt_cellndx_n[k]].push_back(k);}

    // iterates on the grid cells, extracting the points for each cell
    // and computing the possible best-fit-plane

    map<int, vector<int> >::iterator map_iter;
    for(map_iter = gridcell_ptndxs_map.begin(); map_iter != gridcell_ptndxs_map.end(); ++map_iter) {
        vector<int> vector_ndx = (*map_iter).second;
        vector<Point> points_in_cell;
        for (uint pi = 0; pi < vector_ndx.size(); pi++) {
            points_in_cell.push_back(Point(x[vector_ndx[pi]], y[vector_ndx[pi]], z[vector_ndx[pi]]));}

        // CALL FORTRAN PASSING FOUND POINTS AND GET ESTIMATED GEOPLANE (DIP-DIR & ANGLE)

        //outfile << id[ndx_median_elevation] << "," << x[ndx_median_elevation] << "," << y[ndx_median_elevation] << "," << z[ndx_median_elevation] << "\n";

    }

    // output of ascii grid storing number of recs per cell

    outgridfile << "ncols " << columns << "\n";
    outgridfile << "nrows " << rows << "\n";
    outgridfile << "xllcorner " << x_min << "\n";
    outgridfile << "yllcorner " << y_min << "\n";
    outgridfile << "cellsize " << cell_size << "\n";
    outgridfile << "nodata_value " << 0 << "\n";

    for (int j=(rows-1);j>=0;j--) {
       for (int i=0;i<columns;i++) {
          int n = (columns*i) + j;
          if (gridcell_ptndxs_map.count(n) == 0) {
             outgridfile << "0 ";}
          else {
             int num_of_recs_in_cell = gridcell_ptndxs_map[n].size();
             outgridfile << num_of_recs_in_cell << " "; } } }

    outgridfile.close();

    // print run time

    float diff_time = ((float)clock() - (float)start_time)/CLOCKS_PER_SEC;  // run time
    printf ("\n\nProcessing completed in %.6lf seconds\n\n", diff_time );

    // appl end

    return 0;
}




