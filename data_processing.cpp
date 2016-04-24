#include "data_processing.hpp"


int get_dem_param_int(const std::string& rec_line, const std::string& par_name) {
    // read single key-uint value from row - ESRI ASCII DEM format

    unsigned int par_val;
    std::string str_raw, str;
    std::istringstream instr(rec_line);
    instr >> str_raw >> par_val;
    transform(str_raw.begin(), str_raw.end(), back_inserter(str), ::toupper);
    if (str != par_name)
        throw -1;

    return par_val;
};


double get_dem_param_double(const std::string& rec_line, const std::string& par_name) {
    // read single key-double value from row - ESRI ASCII DEM format

    double par_val;
    std::string str_raw, str;
    std::istringstream instr(rec_line);
    instr >> str_raw >> par_val;
    transform(str_raw.begin(), str_raw.end(), back_inserter(str), ::toupper);
    if (str != par_name)
        throw -1;
    return par_val;
};


DataRRGrid read_esri_ascii_dem( std::string dem_filepath ) {
    // read DEM data provided in ESRI Arc/Info ASCII format

    unsigned int ncols, nrows;
    float xllcorner, yllcorner, cellsize, nodata_value;

    std::ifstream infile;
    infile.open(dem_filepath.c_str(), std::ios::binary);

    if (infile.fail() ) {
        std::cout << "\nError: input file " << dem_filepath << " not available\n";
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

    // check that the number of read values in grid is congruent with the number of rows and cols
    if (ncols*nrows != data_vect.size()) {
        std::cout << "Error in data count\n";
        throw -1; }

    // prepare the data structure to return
    Point2D pt2d {xllcorner, yllcorner};
    RectangularDomain domain { pt2d, 0.0, cellsize, cellsize, ncols, nrows };
    DataRRGrid datarrgrid = DataRRGrid(domain, data_vect, nodata_value);

    return datarrgrid;
};


MeshTriangleStrip read_vtk_data_ascii( std::string output_debmeshfile_path, std::string input_vtk_path ) {

    unsigned int num_points, num_triangle_strips, dummy_num;

    // open input VTK file
    std::ifstream infile;
    infile.open(input_vtk_path.c_str(), std::ios::binary);

    // debmeshfile is a debug file, to be removed
    std::ofstream debmeshfile{output_debmeshfile_path};

    if (infile.fail() ) {
        std::cout << "\nError: input file " << input_vtk_path << " not available\n";
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
        debmeshfile << n << ", " << x << ", " << y << ", " << z << "\n";
        pt3d_vect.push_back( Point3D(x, y, z));
        };

    debmeshfile << "\n\n";
    std::getline(infile, rec_line);

    // read TRIANGLE_STRIPS section
    std::string prefix("TRIANGLE_STRIPS");
    std::getline(infile, rec_line);
    if (rec_line.compare(0, prefix.size(), prefix) != 0) {
        std::cout << "Error in expected TRIANGLE_STRIPS header line\n";
        throw -1;
    };
    std::istringstream instr(rec_line);
    std::string dummy_str;
    instr >> dummy_str >> num_triangle_strips >> dummy_num;

    std::vector<std::vector<unsigned int> > triangle_strips;
    std::string blank_space = " ";
    for (unsigned int n = 0; n < num_triangle_strips; n++) {

        std::getline(infile, rec_line);

        std::vector<unsigned int> triangle_strip;
        std::size_t found;
        bool line_end = false;
        int line_n = -1;
        unsigned int num_line_recs;
        while ( ! line_end ){
            found = rec_line.find(blank_space);
            if (found != std::string::npos) {
                line_n++;
                std::string number = rec_line.substr(0, found);
                if (line_n == 0) {
                    num_line_recs = atoi(number.c_str());}
                else {
                    triangle_strip.push_back( atoi(number.c_str()) );
                    debmeshfile << line_n << " - " << number << "\n";
                };
                rec_line.erase(0, found+1);
                ; }
            else {
                line_end = true;  }; };
        assert (triangle_strip.size() == num_line_recs);
        triangle_strips.push_back( triangle_strip );
    };

    infile.close();
    debmeshfile.close();

    return MeshTriangleStrip(pt3d_vect, triangle_strips);
};


std::vector<Triangle3D> extract_triangles_from_mesh(std::string output_debmeshtrianglesfile_path, MeshTriangleStrip surf3d_mesh ) {

    std::vector<Point3D> mesh_pts = surf3d_mesh.pts();
    std::vector< std::vector<unsigned int> > triangle_strips = surf3d_mesh.trianglestrips();

    // debug file, to be removed
    std::ofstream debmeshtrianglesfile{output_debmeshtrianglesfile_path};

    // cycle within triangle strips
    std::vector<Triangle3D> mesh_triangles;
    unsigned int strip_ndx = 0;
    for(std::vector< std::vector<unsigned int> >::iterator ref_trstr = triangle_strips.begin(); ref_trstr != triangle_strips.end(); ++ref_trstr) {
        strip_ndx++;
        debmeshtrianglesfile << "\nstrip: " << strip_ndx << "\n";
        std::vector<unsigned int> triangle_strip = *ref_trstr;
        PointTriplet curr_triplet = PointTriplet();

        // cycle within points of a single triangle strip
        for(std::vector<unsigned int>::iterator ref_ptndx = triangle_strip.begin(); ref_ptndx != triangle_strip.end(); ++ref_ptndx) {

            unsigned int curr_pt_ndx = *ref_ptndx;
            debmeshtrianglesfile << "ndx: " << curr_pt_ndx << "\n";
            Point3D curr_pt = mesh_pts[curr_pt_ndx];
            curr_triplet = curr_triplet.update(curr_pt);
            if (curr_triplet.valid_pts() == 3) {
                Point3D pt1 = curr_triplet.get(0);
                Point3D pt2 = curr_triplet.get(1);
                Point3D pt3 = curr_triplet.get(2);
                Triangle3D geosurf_triangle = Triangle3D(pt1, pt2, pt3);
                debmeshtrianglesfile << "pt1: " << pt1.x() << "," << pt1.y() << "," << pt1.z() << "\n";
                debmeshtrianglesfile << "pt2: " << pt2.x() << "," << pt2.y() << "," << pt2.z() << "\n";
                debmeshtrianglesfile << "pt3: " << pt3.x() << "," << pt3.y() << "," << pt3.z() << "\n";
                mesh_triangles.push_back( geosurf_triangle ); }; }; };

    return mesh_triangles;
};


std::vector<Triangle3D> extract_intersecting_triangles( std::string output_debmeshinterstrfile_path, Space3DPartition dem_vol, std::vector<Triangle3D> mesh_triangles ) {

    // debug file, to be removed
    std::ofstream debmeshintertrfile{output_debmeshinterstrfile_path};

    std::vector<Triangle3D> intersecting_mesh_triangles;
    unsigned int n_orig_triangle = 0;
    unsigned int n_curr_triangle = 0;
    for(std::vector<Triangle3D>::iterator ref_ptndx = mesh_triangles.begin(); ref_ptndx != mesh_triangles.end(); ++ref_ptndx) {
        n_orig_triangle++;
        Triangle3D mesh_triangle = *ref_ptndx;
        Space3DPartition mesh_triangle_volume = mesh_triangle.space_volume();
        if (mesh_triangle_volume.intersects( dem_vol )) {
            n_curr_triangle++;
            Point3D pt1  = mesh_triangle.pt(0);
            Point3D pt2  = mesh_triangle.pt(1);
            Point3D pt3  = mesh_triangle.pt(2);
            debmeshintertrfile << "ndx triangle orig: " << n_orig_triangle << ", current: " << n_curr_triangle << "\n";
            debmeshintertrfile << "pt1: " << pt1.x() << "," << pt1.y() << "," << pt1.z() << "\n";
            debmeshintertrfile << "pt2: " << pt2.x() << "," << pt2.y() << "," << pt2.z() << "\n";
            debmeshintertrfile << "pt3: " << pt3.x() << "," << pt3.y() << "," << pt3.z() << "\n";

            intersecting_mesh_triangles.push_back( mesh_triangle ); }; };

    return intersecting_mesh_triangles;
};


int vect_ndx(int i, int j, int n_cols) {

    return i * n_cols + j;
};


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
            dem_triangles.push_back( Triangle3D( pt_i1_j0, pt_i0_j1, pt_i1_j1 )); }; }; };

    return dem_triangles;
};


std::tuple<Point3D, bool> intersect_segments(Line3D inters_line, Segment3D dem_segment) {
    bool is_in_segment = false;
    Line3D dem_line = dem_segment.as_line();
    Point3D inters_pt = inters_line.intersect_coplanar(dem_line);
    if (dem_segment.is_point_projection_in_segment(inters_pt)) {
        is_in_segment = true;
    };
    return  std::make_tuple(inters_pt, is_in_segment);
};


std::vector<Point3D> get_inters_pts(Triangle3D mesh_triangle, Triangle3D dem_triangle) {

    std::vector<Point3D> inters_pts;

    CartesianPlane mesh_tr_plane = mesh_triangle.to_cartes_plane();

    CartesianPlane dem_tr_plane = dem_triangle.to_cartes_plane();

    Line3D inters_line = mesh_tr_plane.intersect(dem_tr_plane);

    Point3D pt;
    bool is_in_segment;

    Segment3D dem_segment_a = Segment3D(dem_triangle.pt(0), dem_triangle.pt(1));
    std::tie(pt, is_in_segment) = intersect_segments(inters_line, dem_segment_a);
    if (is_in_segment) {
        inters_pts.push_back(pt); };

    Segment3D dem_segment_b = Segment3D(dem_triangle.pt(1), dem_triangle.pt(2));
    std::tie(pt, is_in_segment) = intersect_segments(inters_line, dem_segment_b);
    if (is_in_segment) {
        inters_pts.push_back(pt); };

    Segment3D dem_segment_c = Segment3D(dem_triangle.pt(2), dem_triangle.pt(0));
    std::tie(pt, is_in_segment) = intersect_segments(inters_line, dem_segment_c);
    if (is_in_segment) {
        inters_pts.push_back(pt); };

    return inters_pts;
};


std::vector<Point3D> intersect_dem_geosurface(std::string output_logfile_path, std::vector<Triangle3D> dem_triangles, std::vector<Triangle3D> mesh_intersecting_triangles) {

    std::ofstream logfile{output_logfile_path, std::ofstream::app};
    //if (!logfile) error("can't open output file", output_logfile_path);

    std::vector<Point3D> intersecting_pts;

    int ndx_curr_dem_triangle = 0;
    for(std::vector<Triangle3D>::iterator dem_ref_ptndx = dem_triangles.begin(); dem_ref_ptndx != dem_triangles.end(); ++dem_ref_ptndx) {
        ++ndx_curr_dem_triangle;
        Triangle3D dem_triangle = *dem_ref_ptndx;

        int ndx_curr_mesh = 0;
        for(std::vector<Triangle3D>::iterator mesh_ref_ptndx = mesh_intersecting_triangles.begin(); mesh_ref_ptndx != mesh_intersecting_triangles.end(); ++mesh_ref_ptndx) {
            ++ndx_curr_mesh;
            Triangle3D mesh_triangle = *mesh_ref_ptndx;

            std::vector<Point3D> inters_pts = get_inters_pts(mesh_triangle, dem_triangle);

            if (inters_pts.size() > 0) {

                int ndx_curr_inters_pt = 0;
                for(std::vector<Point3D>::iterator pt_ref_ptndx = inters_pts.begin(); pt_ref_ptndx != inters_pts.end(); ++pt_ref_ptndx) {

                    ++ndx_curr_inters_pt;

                    Point3D inters_pt = *pt_ref_ptndx;
                    intersecting_pts.push_back(inters_pt); };  };  };

        logfile << "\nmesh: " << ndx_curr_mesh << " - current total intersections: " << intersecting_pts.size() << "\n";

        };

    return intersecting_pts;
};

