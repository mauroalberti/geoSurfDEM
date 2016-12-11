
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

int isfloatpart (char c)
{
    if ((c >= '0' && c <='9') or (c == '.'))
        return 1;
    else
        return 0;
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


int main ()
{

    // initial screenshot

    cout << "\n\nBest Fit Geoplanes\n";
    cout << "by M. Alberti - www.malg.eu\n";
    cout << "2016-12-10\n\n\n";

    // input of parameter file name and file opening

    cout << "Enter input parameter file name: ";
  	string param_filename;
    cin >> param_filename;

    // input parameter file opening

    ifstream param_file;
    param_file.open(param_filename.c_str(), ios::binary);
    if (param_file.fail())
    {
    cout << "\n\nUnable to open parameter file '" << param_filename << "'\n";
    return 1;
    }

    // parameter reading

    string xyzdata_line, bfpgeoplanes_line, ptscntgrd_line, cellsize_line;
    string xyzdata_fpth, bfpgeoplanes_fpth, ptnumgrid_fpth, cellsize_str;

    getline(param_file, xyzdata_line);
    istringstream df(xyzdata_line);
    df >> xyzdata_fpth; // xyz file path

    getline(param_file, bfpgeoplanes_line);
    istringstream fe(bfpgeoplanes_line);
    fe >> bfpgeoplanes_fpth; // output geoplanes file path

    getline(param_file, ptscntgrd_line);
    istringstream ge(ptscntgrd_line);
    ge >> ptnumgrid_fpth; // output point number grid file path

    getline(param_file, cellsize_line);
    istringstream cs(cellsize_line);
    cs >> cellsize_str;
    for (uint i = 0; i < cellsize_str.size(); i++) // tests cell size value being number
        { if (isfloatpart(cellsize_str[i]) == 0)
            {cout << "\n\nCell size input: not correct\n";
             return 1;}};
    float cell_size;
	istringstream instr(cellsize_str);
	instr >> cell_size; // used cell size

    // printout of read parameters

    cout << "\nInput parameters\n";
    cout << " - xyz data (input): " << xyzdata_fpth << "\n";
    cout << " - best fit geoplanes (output): " << bfpgeoplanes_fpth << "\n";
    cout << " - point number grid (output): " << ptnumgrid_fpth << "\n";
    cout << " - cell size: " << cellsize_str << "\n";

    //
    // processing begins
    //

    cout << "\n\n -- processing started ... please wait\n";

    // start time

    clock_t start_time = clock();
    cout << "\nStart time: " << start_time;

    // open xyz input file

    ifstream infile;
    infile.open(xyzdata_fpth.c_str(), ios::binary);
    if (infile.fail())
        {cout << "\n\nUnable to open input file '" << xyzdata_fpth << "'\n";
        return 1;}

    // create output geoplanes file

    ofstream outfile;
    outfile.open(bfpgeoplanes_fpth.c_str(), ios::binary);
    if (outfile.fail())
        {cout << "\n\nUnable to create output file '" << bfpgeoplanes_fpth << "'\n";
         return 1;}

    // create output point number grid file

    ofstream outgridfile;
    outgridfile.open(ptnumgrid_fpth.c_str(), ios::binary);
    if (outgridfile.fail())
        {cout << "\n\nUnable to create output grid file '" <<  ptnumgrid_fpth << "'\n";
         return 1;}

    // vars for reading input raw xyz file

    string rec_line;
    list<string> rawdata_list; //  list of raw data strings

    // reads file header

    getline(infile, rec_line);
    while (!infile.eof())
        {getline(infile, rec_line);
         if (rec_line.size() > 0)
            {rawdata_list.push_back(rec_line);}}
    infile.close();

    // vars for storing xyz input data

    int num_recs = rawdata_list.size();
    list<string>::iterator string_pos; // string list iterator
    vector<int> id(num_recs); vector<double> x(num_recs), y(num_recs), z(num_recs);

    // read input xyz data into three double vectors

    char sep;
    int ndx = 0;
    for (string_pos=rawdata_list.begin(); string_pos!=rawdata_list.end(); string_pos++)
       {istringstream instr(*string_pos);
        instr >> x[ndx] >> sep >> y[ndx] >> sep >> z[ndx];
        ndx++;}

    // define spatial ranges of x, y and z

    double x_min = *min_element( x.begin(), x.end() );
    double x_max = *max_element( x.begin(), x.end() );
    double y_min = *min_element( y.begin(), y.end() );
    double y_max = *max_element( y.begin(), y.end() );
    double z_min = *min_element( z.begin(), z.end() );
    double z_max = *max_element( z.begin(), z.end() );
    float x_range = x_max -  x_min;
    float y_range = y_max -  y_min;

    // define number of rows and columns for grid analysis

    int rows = int(y_range/cell_size) + 1;
    int columns = int(x_range/cell_size) + 1;

    // printout info on spatial ranges and grid params

    cout << "\nThe spatial range of data is:\n";
    cout << "  x min: " << x_min << " x max: " << x_max << "\n";
    cout << "  y min: " << y_min << " y max: " << y_max << "\n";
    cout << "  z min: " << z_min << " z max: " << z_max << "\n";
    cout << "\nRows number: " << rows << "\nColumn number: " << columns << "\n ";

    // define the grid indices for the points, together with non-empty grid cells

   vector<int> pt_cellndx_i(num_recs), pt_cellndx_j(num_recs), pt_cellndx_n(num_recs);
   set<int> nonemptycell_ndxs_set;
    for (int k = 0; k < num_recs; k++)
        {pt_cellndx_i[k] = int((x[k]- x_min)/cell_size);
         pt_cellndx_j[k] = int((y[k]- y_min)/cell_size);
         pt_cellndx_n[k] = (columns * pt_cellndx_i[k]) + pt_cellndx_j[k];
         nonemptycell_ndxs_set.insert(pt_cellndx_n[k]);}

   // initialize the grid cell - point mapping for non empty cells

   map<int, vector<int> > gridcell2pts_map;
   set<int>::const_iterator pos;
   for(pos = nonemptycell_ndxs_set.begin(); pos != nonemptycell_ndxs_set.end(); ++pos)
      {vector<int> empy_vector;
      gridcell2pts_map[*pos] = empy_vector;}

   // insert the point indices into the grid cell references (linear index)

   for (int k = 0; k < num_recs; k++)
      {gridcell2pts_map[pt_cellndx_n[k]].push_back(k);}

    //outfile << "id_rec, x, y, z\n";

    map<int, vector<int> >::iterator map_iter;
    for(map_iter = gridcell2pts_map.begin(); map_iter != gridcell2pts_map.end(); ++map_iter) {
        vector<int> vector_ndx = (*map_iter).second;
        vector<Point> points_in_cell;
        for (uint pi = 0; pi < vector_ndx.size(); pi++) {
            points_in_cell.push_back(Point(x[vector_ndx[pi]], y[vector_ndx[pi]], z[vector_ndx[pi]]));}

        // CALL FORTRAN PASSING FOUND POINTS AND GET ESTIMATED GEOPLANE (DIP-DIR & ANGLE

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
          if (gridcell2pts_map.count(n) == 0) {
             outgridfile << "0 ";}
          else {
             int num_of_recs_in_cell = gridcell2pts_map[n].size();
             outgridfile << num_of_recs_in_cell << " "; } } }

    outgridfile.close();

    /*

    */

    // end time

    clock_t end_time = clock();
    cout << "\nEnd time: " << end_time;
    float diff_time = ((float)end_time - (float)start_time)/CLOCKS_PER_SEC;  // run time

    printf ("\n\nProcessing completed in %.6lf seconds\n\n", diff_time );

    // appl end

    return 0;
}




