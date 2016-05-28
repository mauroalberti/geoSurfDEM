#include "data_processing.hpp"

struct inters_result {
  Point3D inter_pt;
  uint dem_tr_ndx, geosurf_tr_ndx;
  double dist_dem_triangle, dist_geosurf_triangle;
};


struct pt_result {
    Point3D inter_pt;
    double dist_dem_triangle, dist_geosurf_triangle;
};


// modified from: http://www.blackpawn.com/texts/pointinpoly/
bool point_on_same_side(Point3D p1, Point3D p2, Point3D a, Point3D b) {

    Vector3D b_a = Vector3D(a, b);
    Vector3D p1_a = Vector3D(a, p1);
    Vector3D p2_a = Vector3D(a, p2);

    Vector3D cp1 = b_a.vector_prod(p1_a);
    Vector3D cp2 = b_a.vector_prod(p2_a);

    if (cp1.scalar_prod(cp2) >= 0.0) {
        return true; }
    else {
        return false; };
};


// from http://www.blackpawn.com/texts/pointinpoly/
bool point_in_triangle(Point3D p, Point3D a, Point3D b, Point3D c) {

    if (point_on_same_side(p, a, b, c) and
        point_on_same_side(p, b, a, c) and
        point_on_same_side(p, c, a, b) ) {
       return true; }
    else {
       return false; };
};


std::tuple<std::string, std::vector<pt_result> > triangle_pair_inters_pts(Triangle3D dem_triangle, Triangle3D geosurf_triangle) {

    // declares return variables

    std::string msg;
    std::vector<Point3D> inters_pts;
    std::vector<pt_result> inters_pts_dist;

    // check input validity

    if (dem_triangle.area() < 1.0e-10) {
        msg = "degenerate DEM triangle"; }
    else if (geosurf_triangle.area() < 1.0e-10) {
        msg = "degenerate mesh triangle"; }
    else {
         // get cartesian plane from DEM triangle
        CartesianPlane dem_tr_plane = dem_triangle.to_cartes_plane();
        // get cartesian plane from geological surface mesh triangle
        CartesianPlane geosurf_tr_plane = geosurf_triangle.to_cartes_plane();
        // check parallelism/coincidence between the two planes
        bool plane_parallelism = dem_tr_plane.isparallel(geosurf_tr_plane);
        if (plane_parallelism) {
            bool coincident_planes = dem_tr_plane.isequidistant(geosurf_tr_plane);
            if (coincident_planes) {
                msg = "parallel coincident planes";
            }
            else {
                msg = "parallel non-coincident planes";
            };
        }
        else {
            Line3D inters_line = geosurf_tr_plane.intersect(dem_tr_plane);
            Point3D iline_pt = inters_line.orig_pt();
            Vector3D iline_versor = inters_line.versor();
            inters_pts = find_triangle_inters(dem_triangle, inters_line);

            Point3D geosurf_pt_0 = geosurf_triangle.pt(0);
            Point3D geosurf_pt_1 = geosurf_triangle.pt(1);
            Point3D geosurf_pt_2 = geosurf_triangle.pt(2);

            for (uint i = 0; i < inters_pts.size(); i++) {
                Point3D curr_pt = inters_pts[i];

                if (not point_in_triangle(curr_pt, geosurf_pt_0, geosurf_pt_1, geosurf_pt_2)) {
                    continue;
                };

                pt_result curr_result;
                curr_result.inter_pt = curr_pt;
                curr_result.dist_dem_triangle = dem_tr_plane.point_distance(curr_pt);
                curr_result.dist_geosurf_triangle = geosurf_tr_plane.point_distance(curr_pt);
                inters_pts_dist.push_back(curr_result);
            };

            msg = "intersecting planes";
        };
    };

    return  std::make_tuple(msg, inters_pts_dist);

};


std::vector<inters_result> intersect_dem_geosurface(std::vector<Triangle3D> dem_triangles, std::vector<Triangle3D> mesh_intersecting_triangles) {

    std::vector<inters_result> intersections;

    int ndx_curr_dem_triangle = 0;

    for(std::vector<Triangle3D>::iterator dem_ref_ptndx = dem_triangles.begin(); dem_ref_ptndx != dem_triangles.end(); ++dem_ref_ptndx) {

        ++ndx_curr_dem_triangle;
        Triangle3D dem_triangle = *dem_ref_ptndx;

        int ndx_curr_mesh = 0;

        for(std::vector<Triangle3D>::iterator mesh_ref_ptndx = mesh_intersecting_triangles.begin(); mesh_ref_ptndx != mesh_intersecting_triangles.end(); ++mesh_ref_ptndx) {

            ++ndx_curr_mesh;
            Triangle3D mesh_triangle = *mesh_ref_ptndx;

            std::string msg;
            std::vector<pt_result> inters_pts;
            std::tie(msg, inters_pts) = triangle_pair_inters_pts(dem_triangle, mesh_triangle);

            if (msg != "intersecting planes") {
                std::cout << "DEM triangle ndx: " << ndx_curr_dem_triangle << ", mesh triangle ndx: " << ndx_curr_mesh << "; error: " << msg << "\n";
                continue; };

            if (inters_pts.size() > 0) {

                int ndx_curr_inters_pt = 0;
                for(std::vector<pt_result>::iterator result_ref_ptndx = inters_pts.begin(); result_ref_ptndx != inters_pts.end(); ++result_ref_ptndx) {
                    ++ndx_curr_inters_pt;
                    pt_result inters = *result_ref_ptndx;
                    inters_result result;
                    result.inter_pt = inters.inter_pt;
                    result.dist_dem_triangle = inters.dist_dem_triangle;
                    result.dist_geosurf_triangle = inters.dist_geosurf_triangle;
                    result.dem_tr_ndx = ndx_curr_dem_triangle;
                    result.geosurf_tr_ndx = ndx_curr_mesh;
                    intersections.push_back(result);
                };
            };
        };
    };

    return intersections;
};


void print_src_data(std::string outfile_path, std::vector<Triangle3D> src_triangles) {

    std::ofstream outdatafile{outfile_path};
    outdatafile << "nxd,p1_x,p1_y,p1_z,p2_x,p2_y,p2_z,p3_x,p3_y,p3_z\n";
    for (uint i = 0; i < src_triangles.size(); i++) {
        Triangle3D curr_triangle = src_triangles[i];
        Point3D p1, p2, p3;
        p1 = curr_triangle.pt(0);
        p2 = curr_triangle.pt(1);
        p3 = curr_triangle.pt(2);
        outdatafile << i << "," <<  p1.x() << "," << p1.y() << "," << p1.z() << "," <<
                                    p2.x() << "," << p2.y() << "," << p2.z() << "," <<
                                    p3.x() << "," << p3.y() << "," << p3.z() << "\n";
    };
};

int main() {

    // read DEM data from input file

    std::string input_dem_path = "./test_data/malpi_aster_w4u3_rit.asc";

    // input 3D geosurface (VTK format)

    std::string input_vtk_path = "./test_data/geosurf3d_01.vtk";

    // output data file

    std::string output_datafile_path = "./test_data/outdata_10.csv";

    // output source dem data file

    std::string output_srcdem_file_path = "./test_data/src_dem_data_10.csv";
    std::string output_srcgeosurf_file_path = "./test_data/src_geosurface_data_10.csv";

    // program header

    std::cout << "\n*** geoSurfDEM *** \n";
    std::cout << "\nApplication for determining intersections between 3D geosurfaces and DEM topography\n\n";

    // read input DEM data

    DataRRGrid datarrgrid = read_esri_ascii_dem( input_dem_path );
    RectangularDomain rect_dom = datarrgrid.rect_domain();

    // read VTK data from input file

    MeshTriangleStrip surf3d_mesh;
    try {
        surf3d_mesh = read_vtk_data_ascii(input_vtk_path ); }
    catch (int e) {
        std::cout << "\nError: program will stop\n";
        return -1; }

    // get triangles (Triangle3D) from mesh

    std::vector<Triangle3D> mesh_triangles = extract_triangles_from_mesh(surf3d_mesh );
    std::cout << "\nNum. total mesh triangles is " << mesh_triangles.size() << "\n";

    // calculate solid volume used to check for mesh triangle volume intersection

    Space3DPartition dem_vol = datarrgrid.space_partition();

    // get mesh triangles intersecting with DEM boundaries

    std::vector<Triangle3D> mesh_intersecting_triangles = extract_intersecting_triangles(dem_vol, mesh_triangles );
    std::cout << "\nNum. intersecting mesh triangles is " << mesh_intersecting_triangles.size() << "\n";
    print_src_data(output_srcgeosurf_file_path, mesh_intersecting_triangles);

    // transform DEM data into a vector of 3D points, valid or invalid

    std::vector<Point3D> dem_3dpts = create_pts_vector(datarrgrid.data(), datarrgrid.rect_domain());
    std::cout << "\nNum. dem 3d pts is " << dem_3dpts.size() << "\n";

    // create vector of valid DEM triangles, for intersection with triangles mesh

    std::vector<Triangle3D> dem_triangles = create_dem_triangles( dem_3dpts, datarrgrid.rect_domain().nrows(), datarrgrid.rect_domain().ncols() );
    std::cout << "\nNum. dem 3d triangles is " << dem_triangles.size() << "\n";
    print_src_data(output_srcdem_file_path, dem_triangles);

    // get intersection points

    std::vector<inters_result> intersections = intersect_dem_geosurface(dem_triangles, mesh_intersecting_triangles);
    std::cout << "\nNum. intersecting pts is " << intersections.size() << "\n";

    // write intersection points

    if (intersections.size() > 0) {
        std::ofstream outdatafile{output_datafile_path};
        outdatafile << "x,y,z,dem_tr_ndx,dist_dem_triangle,geosurf_tr_ndx,dist_geosurf_triangle\n";
        for(std::vector<inters_result>::iterator intersection_ref = intersections.begin(); intersection_ref != intersections.end(); ++intersection_ref) {
            inters_result intersection = *intersection_ref;
            Point3D inters_pt = intersection.inter_pt;
            outdatafile << inters_pt.x() << "," << inters_pt.y() << "," << inters_pt.z() << "," << intersection.dem_tr_ndx << "," << intersection.dist_dem_triangle << "," << intersection.geosurf_tr_ndx << "," << intersection.dist_geosurf_triangle << "\n";
        };
    };

    return 0;


};







