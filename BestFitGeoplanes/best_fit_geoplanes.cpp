
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

#include <sys/stat.h>
#include <unistd.h>

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

    double dist_2d(Point another) {

        double dx = _x - another._x;
        double dy = _y - another._y;
        return sqrt(dx*dx + dy*dy); }

    double dist_3d(Point another) {

        double dx = _x - another._x;
        double dy = _y - another._y;
        double dz = _z - another._z;
        return sqrt(dx*dx + dy*dy + dz*dz); }

    bool is_coincident_2d(Point another, double max_threshold = 1.0e-1) {

        if (dist_2d(another) < max_threshold) {
            return true; }
        else {
            return false; } }

    bool is_coincident_3d(Point another, double max_threshold = 1.0e-1) {
        if (dist_3d(another) < max_threshold) {
            return true; }
        else {
            return false; } }

};

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
    std::cout << "by M. Alberti - www.malg.eu\n";
    std::cout << "2016-12-17\n\n\n";

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
        std::cout << "\nParameter input error: " << error_msg << "\n\n";
        return 1;};

    if (! exists(xyz_data_fpth)) {
        std::cout << "\nError: input xyz data file does not exist\n\n";
        return 1;}

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

    cout << "\n-- processing started ... please wait\n\n";

    // open xyz input file

    ifstream input_xyz_file;
    input_xyz_file.open(xyz_data_fpth.c_str(), ios::binary);
    if (input_xyz_file.fail()) {
        cout << "\n\nUnable to open input file '" << xyz_data_fpth << "'\n";
        return 1;}

    // create output geoplanes file

    ofstream output_results_file;
    output_results_file.open(bfp_geoplanes_fpth.c_str(), ios::binary);
    if (output_results_file.fail())
        {cout << "\n\nUnable to create output file '" << bfp_geoplanes_fpth << "'\n";
         return 1;}

    // create output point number grid file

    ofstream output_cntgrid_file;
    output_cntgrid_file.open(ptnum_grid_fpth.c_str(), ios::binary);
    if (output_cntgrid_file.fail())
        {cout << "\n\nUnable to create output grid file '" <<  ptnum_grid_fpth << "'\n";
         return 1;}

    // read input raw xyz file

    string rec_line;
    list<string> rawdata_list;  // list of raw data strings

    for (int n = 0; n < hdr_ln_num; n++) {  // read  (and discard) file header
        getline(input_xyz_file, rec_line); }

    while (! input_xyz_file.eof()) {  // read raw data
        getline(input_xyz_file, rec_line);
        if (rec_line.size() > 0) {
            rawdata_list.push_back(rec_line);}}
    input_xyz_file.close();

    std::cout << "Number of read records: " <<rawdata_list.size() << "\n";

    // read input xyz data into three double vectors

    int num_recs = rawdata_list.size();  // total number of input points
    vector<double> x(num_recs), y(num_recs), z(num_recs);  // vectors storing the x, y and z coordinates for all points

    char sep;
    int ndx = 0;
    list<string>::iterator string_pos;  // string list iterator
    for (string_pos=rawdata_list.begin(); string_pos!=rawdata_list.end(); string_pos++) {
        istringstream instr(*string_pos);
        instr >> x[ndx] >> sep >> y[ndx] >> sep >> z[ndx];
        ndx++; }

    if (x.size() != y.size() || x.size() != z.size()) {
        std::cout << "\nSize of read coordinates is different between x, y or z\n";
        return 1;}
    int num_input_pts = x.size();
    std::cout << "Number of input points: " << num_input_pts << "\n";

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

    cout << "Spatial range of input data :\n";
    cout << "  x min: " << x_min << " x max: " << x_max << "\n";
    cout << "  y min: " << y_min << " y max: " << y_max << "\n";
    cout << "  z min: " << z_min << " z max: " << z_max << "\n";
    cout << "Rows number: " << rows << "\nColumn number: " << columns << "\n";
    cout << "Spatial range of grid :\n";
    cout << "  x min: " << x_min << " x max: " << x_min + (columns * cell_size) << "\n";
    cout << "  y min: " << y_max - (rows * cell_size) << " y max: " << y_max << "\n";

    // creates vector of distinct points

    vector<Point> distinct_points;
    distinct_points.push_back(Point(x[0], y[0], z[0]));
    output_results_file << "Skipped coincident points\n";
    for (int n = 1; n < num_input_pts; n++) {
        Point candidate_point = Point(x[n], y[n], z[n]);
        bool to_add = true;
        for (uint k = 0; k < distinct_points.size(); k++) {
            Point added_pt = distinct_points[k];
            if (candidate_point.is_coincident_2d(added_pt)) {
                to_add = false;
                output_results_file << " - skipped coincident point: " << n + 1 << "\n";
                break; } }
        if (to_add) distinct_points.push_back(candidate_point); }
    int num_distinct_pts = distinct_points.size();

    // write distinct point infos in output result file

    std::cout << "Number of distinct points: " << num_distinct_pts << "\n";
    output_results_file << "\nNumber of distinct points: " << num_distinct_pts << "\n";

    // calculate grid indices of input points

    vector<int> dp_cellndx_i(num_distinct_pts), dp_cellndx_j(num_distinct_pts);
    for (int k = 0; k < num_distinct_pts; k++) {
        dp_cellndx_i[k] = int((y_max - distinct_points[k].y())/cell_size);
        dp_cellndx_j[k] = int((distinct_points[k].x()- x_min)/cell_size);}

    // define grid linear indices of input points

    vector<int> dp_cellndx_n(num_distinct_pts);
    for (int k = 0; k < num_distinct_pts; k++) {
        dp_cellndx_n[k] = (columns * dp_cellndx_i[k]) + dp_cellndx_j[k]; }

    // write distinct point infos in output result file

    output_results_file << "\ncnt, x, y, z, i, j, n\n";
    for (uint k = 0; k < distinct_points.size(); k++) {
        Point dp = distinct_points[k];
        output_results_file << k << ", " << dp.x() << ", " << dp.y() << ", " << dp.z() << ", " <<
            dp_cellndx_i[k] << ", " << dp_cellndx_j[k] << ", " << dp_cellndx_n[k] << "\n"; }

    // define linear grid indices of cells with at least one point, using a set data structure

    set<int> non_empty_cells;
    for (int k = 0; k < num_distinct_pts; k++) {
        non_empty_cells.insert(dp_cellndx_n[k]); }
    std::cout << "Number of cells with at least one point is: " << non_empty_cells.size() << "\n";

    // initialize the mapping between each cell grid linear index and associated point indices

    map<int, vector<int> > gridcell_ptndxs_map;
    set<int>::const_iterator pos;

    for(pos = non_empty_cells.begin(); pos != non_empty_cells.end(); ++pos) {  // initialize the mapping between each cell grid linear index and associated point indices
        vector<int> empy_vector;
        gridcell_ptndxs_map[*pos] = empy_vector; }

    for (int k = 0; k < num_distinct_pts; k++) {  // insert the point indices into the grid cell linear indices
        gridcell_ptndxs_map[dp_cellndx_n[k]].push_back(k); }

    // write non-empty cell infos in output result file

    output_results_file << "\nNon-empty cell infos\n";
    for (auto cell_ndx : non_empty_cells) {
        //int cell_ndx = non_empty_cells[l];
        vector<int> dp_indices = gridcell_ptndxs_map[cell_ndx];
        int num_pts = dp_indices.size();
        output_results_file << "grid cell: " << cell_ndx << " - points: " << num_pts << "\n";
        for (auto pt_ndx : dp_indices) {
            output_results_file << " point: " << pt_ndx << "\n"; } }

    // iterates on the grid cells, extracting the points for each cell
    // and computing the possible best-fit-plane

    map<int, vector<int> >::iterator map_iter;
    output_results_file << "\nInverted geological attitudes\n";
    output_results_file << "cell_ndx, dipdir, dipang\n";
    for(map_iter = gridcell_ptndxs_map.begin(); map_iter != gridcell_ptndxs_map.end(); ++map_iter) {
        int cell_ndx = (*map_iter).first;
        vector<int> vector_ndxs = (*map_iter).second;
        int cell_rec_num = vector_ndxs.size();
        if (cell_rec_num >= 3) {
            double pts_coords[cell_rec_num][3];
            for (int pi = 0; pi < cell_rec_num; pi++) {
                pts_coords[pi][0] = distinct_points[vector_ndxs[pi]].x();
                pts_coords[pi][1] = distinct_points[vector_ndxs[pi]].y();
                pts_coords[pi][2] = distinct_points[vector_ndxs[pi]].z(); }

            bool success;
            double dipdir, dipang;
            invert_attitudes(&cell_rec_num, &pts_coords[0][0], &success, &dipdir, &dipang);

            if (success) {
                output_results_file << cell_ndx << ", " << dipdir << ", " << dipang << "\n"; } } }


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




