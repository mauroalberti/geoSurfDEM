
/*
#
# best_fit_geoplanes.cpp
# by Mauro Alberti - alberti.m65@gmail.com
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

#include <sys/stat.h>
#include <unistd.h>
#include <ctype.h>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>  // std::min_element, std::max_element
#include <set>
#include <map>

using namespace std;


extern "C" {

    void invert_attitudes(int *num_points, double *pts_array_cform, bool *success, double *dipdir, double *dipang);

}


bool is_emptywhitespace(std::string& s) {

    bool result = true;

    if (s.empty())
        return result;

    std::string::iterator p;
    for (p = s.begin(); p != s.end(); ++p) {
        if (!isspace(*p)) {
            result = false;
            break;
        }
    }

    return result;

}


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


void read_params(std::stringstream &param_lines, string &xyz_data_fpth, int &hdr_ln_num,
                 string &bfp_geoplanes_fpth, string &rep_analysis_fpth, string &ptnum_grid_fpth, float &cell_size) {

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
    rep_analysis_fpth = read_string(readline);

    getline(param_lines, readline);
    ptnum_grid_fpth = read_string(readline);

    getline(param_lines, readline);
    readline.erase(readline.end()-1);
    if (is_float(readline)) {
        cell_size = read_float(readline);}
    else {
        throw "Cell size is not a float value";};

};


class Point2d {

private:

    double _x, _y;

public:

    Point2d( double xval, double yval) {
       _x = xval;
       _y = yval; }

    double x() {
       return _x; }

    double y() {
       return _y; }

};


class Point3d {

private:

    double _x, _y, _z;

public:

    Point3d( double xval, double yval, double zval) {
       _x = xval;
       _y = yval;
       _z = zval; }

    double x() {
       return _x;}

    double y() {
       return _y;}

    double z() {
       return _z;}

    double dist_2d(Point3d another) {

        double dx = _x - another._x;
        double dy = _y - another._y;
        return sqrt(dx*dx + dy*dy); }

    double dist_3d(Point3d another) {

        double dx = _x - another._x;
        double dy = _y - another._y;
        double dz = _z - another._z;
        return sqrt(dx*dx + dy*dy + dz*dz); }

    bool is_coincident_2d(Point3d another, double max_threshold = 1.0e-1) {

        if (dist_2d(another) < max_threshold) {
            return true; }
        else {
            return false; } }

    bool is_coincident_3d(Point3d another, double max_threshold = 1.0e-1) {
        if (dist_3d(another) < max_threshold) {
            return true; }
        else {
            return false; } }

};


std::vector<uint> cellndxn2ij(uint cell_ndx, uint columns, uint rows) {

    std::vector<uint> result;
    uint i = int(cell_ndx / columns);
    uint j = cell_ndx - columns * i;

    result.push_back(i);
    result.push_back(j);

    return result;

}


Point2d cellnxd2centerpoint(uint cell_ndx, double tlc_x, double tlc_y, double cell_size, uint columns, uint rows) {

    std::vector<uint> ij_vect = cellndxn2ij(cell_ndx, columns, rows);
    uint i = ij_vect[0];
    uint j = ij_vect[1];

    double x = tlc_x + cell_size * (j + 0.5);
    double y = tlc_y - cell_size * (i + 0.5);

    return Point2d(x, y);

}


// from: http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
// fourth and fastest function

inline bool exists (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}


int main () {

    //
    // user interaction processes
    //

    // initial screenshot

    std::cout << "\n\nBest Fit Geoplanes\n";
    std::cout << "by M. Alberti\n";
    std::cout << "2017-02-03\n\n";

    // input of parameter file name and file opening

    std::cout << "Enter input parameter file name: ";
  	std::string param_filename;
    std::cin >> param_filename;

    if (! exists(param_filename)) {
        std::cout << "\nError: input parameter file does not exist\n\n";
        return 1;}

    // start time

    clock_t start_time = clock();

    // input parameter file reading

    std::ifstream param_file(param_filename);
    std::stringstream param_lines;
    param_lines << param_file.rdbuf();

    // reading parameters

    string xyz_data_fpth; // path of input xyz file
    int hdr_ln_num; // number of header lines in xyz file
    string bfp_geoplanes_fpth; // output geoplanes file path
    string rep_analysis_fpth; // analysis report file path
    string ptnum_grid_fpth;  // output Point3d number grid file path
    float cell_size;  // used cell size

    try {
        read_params(param_lines,
                    xyz_data_fpth,
                    hdr_ln_num,
                    bfp_geoplanes_fpth,
                    rep_analysis_fpth,
                    ptnum_grid_fpth,
                    cell_size); }
    catch (string error_msg) {
        std::cout << "\nError: parameter input - " << error_msg << "\n\n";
        return 1;};

    if (! exists(xyz_data_fpth)) {
        std::cout << "\nError: input xyz data file does not exist\n\n";
        return 1;}

    // open xyz input file

    ifstream input_xyz_file;
    input_xyz_file.open(xyz_data_fpth.c_str(), ios::binary);
    if (input_xyz_file.fail()) {
        std::cout << "\nError: unable to open input file '" << xyz_data_fpth << "'\n";
        return 1;}

    // create output geoplanes file

    ofstream output_results_file;
    output_results_file.open(bfp_geoplanes_fpth.c_str(), ios::binary);
    if (output_results_file.fail())
        {std::cout << "\nError: unable to create output result file '" << bfp_geoplanes_fpth << "'\n";
         return 1;}

    // create output analysis report

    ofstream output_report_file;
    output_report_file.open(rep_analysis_fpth.c_str(), ios::binary);
    if (output_report_file.fail())
        {std::cout << "\nError: unable to create output report file '" << rep_analysis_fpth << "'\n";
         return 1;}

    // create output Point3d number grid file

    ofstream output_cntgrid_file;
    output_cntgrid_file.open(ptnum_grid_fpth.c_str(), ios::binary);
    if (output_cntgrid_file.fail())
        {std::cout << "\nError: unable to create output grid file '" <<  ptnum_grid_fpth << "'\n";
         return 1;}

    //
    // Calculations
    //

    std::cout << "\nCalculations started ... please wait";

    // printout read parameters

    output_report_file << "BestFitGeoplane analysis report\n";
    output_report_file << "\nInput parameters\n";
    output_report_file << " - xyz data (input): " << xyz_data_fpth << "\n";
    output_report_file << " - number of header lines in xyz data file: " << hdr_ln_num << "\n";
    output_report_file << " - best fit geoplanes (output): " << bfp_geoplanes_fpth << "\n";
    output_report_file << " - analysis report (output): " << rep_analysis_fpth << "\n";
    output_report_file << " - Point3d number grid (output): " << ptnum_grid_fpth << "\n";
    output_report_file << " - cell size: " << cell_size << "\n\n";

    // read input raw xyz file

    string rec_line;
    list<string> rawdata_list;  // list of raw data strings

    for (int n = 0; n < hdr_ln_num; n++) {  // read  (and discard) file header
        getline(input_xyz_file, rec_line);
        output_report_file << "header line: " << rec_line << "\n";
    }

    unsigned long int rd_cnt = 0;
    while (! input_xyz_file.eof()) {  // read raw data
        getline(input_xyz_file, rec_line);
        if (!is_emptywhitespace(rec_line) > 0) {
            rd_cnt++;
            output_report_file << rd_cnt << ": " << rec_line << "\n";
            rawdata_list.push_back(rec_line);
        }
    }
    input_xyz_file.close();

    output_report_file << "\n\nNumber of read records: " << rawdata_list.size() << "\n";

    // read input xyz data into three double vectors

    int num_recs = rawdata_list.size();  // total number of input points
    std::vector<double> x(num_recs), y(num_recs), z(num_recs);  // vectors storing the x, y and z coordinates for all points

    char sep;
    int ndx = 0;
    list<string>::iterator string_pos;  // string list iterator
    for (string_pos=rawdata_list.begin(); string_pos!=rawdata_list.end(); string_pos++) {
        istringstream instr(*string_pos);
        instr >> x[ndx] >> sep >> y[ndx] >> sep >> z[ndx];
        ndx++; }

    if (x.size() != y.size() || x.size() != z.size()) {
        output_report_file << "\nSize of read coordinates is different between x, y or z\n";
        return 1;}
    int num_input_pts = x.size();
    output_report_file << "Number of input points: " << num_input_pts << "\n";

    // calculate spatial ranges of x, y and z

    double x_min = *min_element( x.begin(), x.end() );
    double x_max = *max_element( x.begin(), x.end() );
    double y_min = *min_element( y.begin(), y.end() );
    double y_max = *max_element( y.begin(), y.end() );
    double z_min = *min_element( z.begin(), z.end() );
    double z_max = *max_element( z.begin(), z.end() );
    double x_range = x_max -  x_min;
    double y_range = y_max -  y_min;
    double tlc_x = x_min;
    double tlc_y = y_max;

    // calculate grid row and column numbers

    int rows = int(y_range/cell_size) + 1;
    int columns = int(x_range/cell_size) + 1;

    // print info on spatial ranges and grid params

    output_report_file << "Spatial range of input data :\n";
    output_report_file << "  x min: " << x_min << " x max: " << x_max << "\n";
    output_report_file << "  y min: " << y_min << " y max: " << y_max << "\n";
    output_report_file << "  z min: " << z_min << " z max: " << z_max << "\n";
    output_report_file << "Rows number: " << rows << "\nColumn number: " << columns << "\n";
    output_report_file << "Spatial range of grid :\n";
    output_report_file << "  x min: " << x_min << " x max: " << x_min + (columns * cell_size) << "\n";
    output_report_file << "  y min: " << y_max - (rows * cell_size) << " y max: " << y_max << "\n";

    // creates vector of distinct points

    std::vector<Point3d> distinct_points;
    distinct_points.push_back(Point3d(x[0], y[0], z[0]));
    output_report_file << "Skipped coincident points\n";
    for (int n = 1; n < num_input_pts; n++) {
        Point3d candidate_point = Point3d(x[n], y[n], z[n]);
        bool to_add = true;
        for (uint k = 0; k < distinct_points.size(); k++) {
            Point3d added_pt = distinct_points[k];
            if (candidate_point.is_coincident_2d(added_pt)) {
                to_add = false;
                output_report_file << " " << n + 1 << "\n";
                break; } }
        if (to_add) distinct_points.push_back(candidate_point); }
    int num_distinct_pts = distinct_points.size();

    // write distinct Point3d infos in output result file

    output_report_file << "Number of distinct points: " << num_distinct_pts << "\n";
    output_report_file << "\nNumber of distinct points: " << num_distinct_pts << "\n";

    // calculate grid indices of input points

    std::vector<int> dp_cellndx_i(num_distinct_pts), dp_cellndx_j(num_distinct_pts);
    for (int k = 0; k < num_distinct_pts; k++) {
        dp_cellndx_i[k] = int((y_max - distinct_points[k].y())/cell_size);
        dp_cellndx_j[k] = int((distinct_points[k].x()- x_min)/cell_size);}

    // define grid linear indices of input points

    std::vector<int> dp_cellndx_n(num_distinct_pts);
    for (int k = 0; k < num_distinct_pts; k++) {
        dp_cellndx_n[k] = (columns * dp_cellndx_i[k]) + dp_cellndx_j[k]; }

    // write distinct Point3d infos in output result file

    output_report_file << "\ncnt, x, y, z, i, j, n\n";
    for (uint k = 0; k < distinct_points.size(); k++) {
        Point3d dp = distinct_points[k];
        output_report_file << k << ", " << dp.x() << ", " << dp.y() << ", " << dp.z() << ", " <<
            dp_cellndx_i[k] << ", " << dp_cellndx_j[k] << ", " << dp_cellndx_n[k] << "\n"; }

    // define linear grid indices of cells with at least one Point3d, using a set data structure

    set<int> non_empty_cells;
    for (int k = 0; k < num_distinct_pts; k++) {
        non_empty_cells.insert(dp_cellndx_n[k]); }
    output_report_file << "Number of cells containing at least one point is: " << non_empty_cells.size() << "\n";

    // initialize the mapping between each cell grid linear index and associated Point3d indices

    map<int, std::vector<int> > gridcell_ptndxs_map;
    set<int>::const_iterator pos;

    for(auto const & cell_ndx : non_empty_cells) {  // initialize the mapping between each cell grid linear index and associated Point3d indices
        std::vector<int> empy_vector;
        gridcell_ptndxs_map[cell_ndx] = empy_vector; }

    for (int k = 0; k < num_distinct_pts; k++) {  // insert the Point3d indices into the grid cell linear indices
        gridcell_ptndxs_map[dp_cellndx_n[k]].push_back(k); }

    // write non-empty cell infos in output result file

    output_report_file << "\nNon-empty cell infos\n";
    for (auto const & cell_ndx : non_empty_cells) {
        //int cell_ndx = non_empty_cells[l];
        std::vector<int> dp_indices = gridcell_ptndxs_map[cell_ndx];
        int num_pts = dp_indices.size();
        output_report_file << "grid cell: " << cell_ndx << " - points: " << num_pts << "\n";
        for (auto const & pt_ndx : dp_indices) {
            Point3d dp = distinct_points[pt_ndx];
            output_report_file << pt_ndx << ", " << dp.x() << ", " << dp.y() << ", " << dp.z() << "\n";
        }
    }

    output_report_file.close();

    // iterates on the grid cells, extracting the points for each cell
    // and computing the possible best-fit-plane

    std::cout << "\n\nInverting geological attitudes ...";

    output_results_file << "x, y, pt_num, dip_dir, dip_ang, x_range, y_range, z_range, pseudo-volume\n";
    for (auto const & kv : gridcell_ptndxs_map) {

        int cell_ndx = kv.first;
        std::vector<int> vector_ndxs = kv.second;
        int cell_rec_num = vector_ndxs.size();

        if (cell_rec_num >= 3) {

            double pts_coords[cell_rec_num][3];
            std::vector<double> pts_x, pts_y, pts_z;  // vectors storing the x, y and z coordinates for cell points

            for (int pi = 0; pi < cell_rec_num; pi++) { // fill array and vectors with point coordinates

                double pt_x = distinct_points[vector_ndxs[pi]].x();
                double pt_y = distinct_points[vector_ndxs[pi]].y();
                double pt_z = distinct_points[vector_ndxs[pi]].z();

                pts_coords[pi][0] = pt_x; pts_x.push_back(pt_x);
                pts_coords[pi][1] = pt_y; pts_y.push_back(pt_y);
                pts_coords[pi][2] = pt_z; pts_z.push_back(pt_z); }

            bool success;
            double dipdir, dipang;
            invert_attitudes(&cell_rec_num, &pts_coords[0][0], &success, &dipdir, &dipang);

            if (success) {  // inverted geological attitude

                double pts_x_min = *min_element( pts_x.begin(), pts_x.end() );
                double pts_x_max = *max_element( pts_x.begin(), pts_x.end() );

                double pts_y_min = *min_element( pts_y.begin(), pts_y.end() );
                double pts_y_max = *max_element( pts_y.begin(), pts_y.end() );

                double pts_z_min = *min_element( pts_z.begin(), pts_z.end() );
                double pts_z_max = *max_element( pts_z.begin(), pts_z.end() );

                double range_x = pts_x_max -  pts_x_min;
                double range_y = pts_y_max -  pts_y_min;
                double range_z = pts_z_max -  pts_z_min;

                Point2d cell_center = cellnxd2centerpoint((uint)cell_ndx, tlc_x, tlc_y, cell_size, (uint)columns, (uint)rows);
                output_results_file << cell_center.x() << ", " << cell_center.y() << ", " <<
                                    cell_rec_num << ", " << dipdir << ", " << dipang << ", " <<
                                    range_x << ", " << range_y << ", " << range_z << ", " << range_x*range_y*range_z <<"\n"; } } }

    output_results_file.close();

    // output of ascii grid storing number of recs per cell

    output_cntgrid_file << "ncols " << columns << "\n";
    output_cntgrid_file << "nrows " << rows << "\n";
    output_cntgrid_file << "xllcorner " << x_min << "\n";
    output_cntgrid_file << "yllcorner " << y_max - (rows * cell_size) << "\n";
    output_cntgrid_file << "cellsize " << cell_size << "\n";
    output_cntgrid_file << "nodata_value " << 0 << "\n";

    int num_of_recs_in_cell;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            int n = (columns*i) + j;
            if (gridcell_ptndxs_map.count(n)) {
                num_of_recs_in_cell = gridcell_ptndxs_map[n].size(); }
            else {
                num_of_recs_in_cell = 0; }
            output_cntgrid_file << num_of_recs_in_cell << " "; } }

    output_cntgrid_file.close();

    // print run time

    float diff_time = ((float)clock() - (float)start_time)/CLOCKS_PER_SEC;  // run time
    printf ("\n\nProcessing completed in %.6lf seconds\n\n", diff_time );

    // appl end

    return 0;
}




