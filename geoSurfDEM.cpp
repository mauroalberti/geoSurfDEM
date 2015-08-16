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

    RectangularDomain domain { pt2d, 0.0, cellsize, cellsize, ncols, nrows };

    DataRRGrid datarrgrid = DataRRGrid(domain, data_vect, nodata_value);

    return datarrgrid;

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


std::vector<Triangle3D> extract_triangles_from_mesh(MeshTriangleStrip surf3d_mesh ) {

    std::vector<Point3D> mesh_pts = surf3d_mesh.pts();
    std::vector< std::vector<unsigned int> > triangle_strips = surf3d_mesh.trianglestrips();

    // cycle within triangle strips
    std::vector<Triangle3D> mesh_triangles;
    for(std::vector< std::vector<unsigned int> >::iterator ref_trstr = triangle_strips.begin(); ref_trstr != triangle_strips.end(); ++ref_trstr) {
        std::vector<unsigned int> triangle_strip = *ref_trstr;
        PointTriplet curr_triplet = PointTriplet();

        // cycle within points of a single triangle strip
        for(std::vector<unsigned int>::iterator ref_ptndx = triangle_strip.begin(); ref_ptndx != triangle_strip.end(); ++ref_ptndx) {

            unsigned int curr_pt_ndx = *ref_ptndx;
            Point3D curr_pt = mesh_pts[curr_pt_ndx];
            curr_triplet = curr_triplet.update(curr_pt);
            if (curr_triplet.valid_pts() == 3) {
                Triangle3D geosurf_triagle = Triangle3D(curr_triplet.get(0), curr_triplet.get(1), curr_triplet.get(2));
                mesh_triangles.push_back( geosurf_triagle ); }; }; };

    return mesh_triangles;

};


std::vector<Triangle3D> extract_intersecting_triangles(Space3DPartition dem_vol, std::vector<Triangle3D> mesh_triangles ) {

    std::vector<Triangle3D> intersecting_mesh_triangles;
    for(std::vector<Triangle3D>::iterator ref_ptndx = mesh_triangles.begin(); ref_ptndx != mesh_triangles.end(); ++ref_ptndx) {
        Triangle3D mesh_triangle = *ref_ptndx;
        Space3DPartition mesh_triangle_volume = mesh_triangle.space_volume();
        if (mesh_triangle_volume.intersects( dem_vol )) {
            intersecting_mesh_triangles.push_back( mesh_triangle ); }; };

    return intersecting_mesh_triangles;  };


int vect_ndx(int i, int j, int n_cols) {

    return i * n_cols + j; };


std::vector<int> dem_indices(int ndx, int ncols) {

    int i = ndx / ncols;
    int j = ndx - i*ncols;
    std::vector<int> indices {i, j};

    return indices;

};


std::vector<double> unrotated_coords(int n_rows, int i, int j, double l_size, double m_size) {

    double x_unrot = l_size * ( double(j) + 0.5);
    double y_unrot = m_size * (n_rows - double(i) -0.5 );

    std::vector<double>  unrot_coords {x_unrot, y_unrot};

    return unrot_coords;

};


std::vector<Point3D> create_pts_vector(NumericData grid_data, RectangularDomain grid_geograph) {

    int n_rows = grid_geograph.nrows(), n_cols = grid_geograph.ncols();
    double l_size = grid_geograph.l(), m_size = grid_geograph.m();
    double rotation_angle_degr = grid_geograph.rot_angle();
    Point2D llcorner_pt = grid_geograph.pt();

    std::vector<double> data = grid_data.values();
    double null_data_val = grid_data.nodata_value();

    Matrix2 dem_rot_matrix = angle_to_rot_matrix(rotation_angle_degr);


    std::vector<Point3D> dem_pts_vector;
    int vector_pts_ndx = -1;
    for(std::vector<double>::iterator ref_ptndx = data.begin(); ref_ptndx != data.end(); ++ref_ptndx) {

        vector_pts_ndx++;


        double z = *ref_ptndx;
        bool valid_pt = fabs(z - null_data_val) > 1.0e-12;
        std::vector<int> dem_ndxs = dem_indices(vector_pts_ndx, n_cols);
        int i = dem_ndxs[0], j = dem_ndxs[1];

        std::vector<double> unrot_coords = unrotated_coords(n_rows, i, j, l_size, m_size);
        double x_unrot = unrot_coords[0], y_unrot = unrot_coords[1];
        Point2D unrot_pt2d = Point2D(x_unrot, y_unrot);
        Point2D rot_pt2 = unrot_pt2d.rotateby(dem_rot_matrix);
        Point2D shftd_pt2d = rot_pt2.moveby(llcorner_pt);
        Point3D dem_pts3d = Point3D(shftd_pt2d, z, valid_pt);

        dem_pts_vector.push_back(dem_pts3d);

        };



    return dem_pts_vector;

};


bool check_triangle_points_validity(Point3D pt_a, Point3D pt_b, Point3D pt_c) {

    if (pt_a.is_valid() && pt_b.is_valid() && pt_c.is_valid()) {
        return true; }
    else {
        return false; };

};


std::vector<Triangle3D> create_dem_triangles(std::vector<Point3D> dem_3dpts, int nrows, int ncols ) {

    std::vector<Triangle3D> dem_triangles;

    for (int i = 0; i < nrows-1; ++i) {
        for (int j = 0; j < ncols-1; ++j) {

          int pt_i0_j0_ndx = vect_ndx(i, j, ncols);
          int pt_i0_j1_ndx = vect_ndx(i, j+1, ncols);
          int pt_i1_j0_ndx = vect_ndx(i+1, j, ncols);
          int pt_i1_j1_ndx = vect_ndx(i+1, j+1, ncols);

          Point3D pt_i0_j0 = dem_3dpts[pt_i0_j0_ndx];
          Point3D pt_i0_j1 = dem_3dpts[pt_i0_j1_ndx];
          Point3D pt_i1_j0 = dem_3dpts[pt_i1_j0_ndx];
          Point3D pt_i1_j1 = dem_3dpts[pt_i1_j1_ndx];

          if (check_triangle_points_validity(pt_i0_j0, pt_i0_j1, pt_i1_j0)) {
            dem_triangles.push_back( Triangle3D( pt_i0_j0, pt_i0_j1, pt_i1_j0 )); };

          if (check_triangle_points_validity(pt_i1_j1, pt_i0_j1, pt_i1_j0)) {
            dem_triangles.push_back( Triangle3D( pt_i1_j1, pt_i0_j1, pt_i1_j0 )); }; }; };

    return dem_triangles;

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

    std::cout << "\n ** geoSurfDEM ** \n";
    std::cout << "\nApplication for determining intersections between 3D geosurfaces and DEM topography\n\n";

    // read DEM data from input file
    std::string input_dem_path = "./test_data/malpi_w4u2n_src.asc";
    DataRRGrid datarrgrid = read_esri_ascii_dem( input_dem_path );
    RectangularDomain rect_dom = datarrgrid.rect_domain();

    // solid volume used to check for mesh triangle volume intersection
    Space3DPartition dem_vol = datarrgrid.space_partition();

    // read VTK data from input file
    std::string input_vtk_path = "./test_data/surf3d_sim_01_rot45_04500.vtk";
    MeshTriangleStrip surf3d_mesh;
    try {
        surf3d_mesh = read_vtk_data_ascii( input_vtk_path ); }
    catch (int e) {
        std::cout << "\n ***** Program will stop *****\n";
        return -1; }

    // get triangles (Triangle3D) from mesh
    std::vector<Triangle3D> mesh_triangles = extract_triangles_from_mesh( surf3d_mesh );
    std::cout << "\nnum. total mesh triangles is " << mesh_triangles.size() << "\n";

    // get mesh triangles intersecting with DEM boundaries
    std::vector<Triangle3D> mesh_intersecting_triangles = extract_intersecting_triangles( dem_vol, mesh_triangles );
    std::cout << "\nnum. intersecting mesh triangles is " << mesh_intersecting_triangles.size() << "\n";

    // transform DEM data into a vector of 3D points, valid or invalid
    std::vector<Point3D> dem_3dpts = create_pts_vector(datarrgrid.data(), datarrgrid.rect_domain());
    std::cout << "\nnum. dem 3d pts is " << dem_3dpts.size() << "\n";

    // create vector of valid DEM triangles, for intersecting with the mesh traingles
    std::vector<Triangle3D> dem_triangles = create_dem_triangles( dem_3dpts, datarrgrid.rect_domain().nrows(), datarrgrid.rect_domain().ncols() );
    std::cout << "\nnum. dem 3d triangles is " << dem_triangles.size() << "\n";


    return 0;

};
