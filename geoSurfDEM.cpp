#include "data_processing.hpp"


struct inters_result {
  Point3D inter_pt;
  GeologicalPlane geoplane;
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
bool point_in_triangle(Point3D p, Triangle3D tr) {

    Point3D a = tr.pt(0);
    Point3D b = tr.pt(1);
    Point3D c = tr.pt(2);

    if (point_on_same_side(p, a, b, c) and
        point_on_same_side(p, b, a, c) and
        point_on_same_side(p, c, a, b) ) {
       return true; }
    else {
       return false; };
};


std::tuple<std::string, std::vector<Point3D> > triangle_pair_inters_pts(Triangle3D dem_triangle, Triangle3D geosurf_triangle) {

    // declares return variables

    std::string msg;
    std::vector<Point3D> inters_pts;
    std::vector<Point3D> inters_results;

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

            for (uint i = 0; i < inters_pts.size(); i++) {
                Point3D curr_pt = inters_pts[i];
                if (point_in_triangle(curr_pt, geosurf_triangle)) {
                    inters_results.push_back(curr_pt);
                };
            };

            msg = "intersecting planes";
        };
    };

    return  std::make_tuple(msg, inters_results);

};


std::vector<inters_result> intersect_dem_geosurface(std::vector<Triangle3D> dem_triangles, std::vector<geosurf_triangle> mesh_intersecting_triangles) {

    std::vector<inters_result> intersections;

    uint ndx_curr_dem_triangle = 0;

    for(std::vector<Triangle3D>::iterator dem_ref_ptndx = dem_triangles.begin(); dem_ref_ptndx != dem_triangles.end(); ++dem_ref_ptndx) {

        Triangle3D dem_triangle = *dem_ref_ptndx;

        uint ndx_curr_mesh_triangle = 0;

        for(std::vector<geosurf_triangle>::iterator mesh_ref_ptndx = mesh_intersecting_triangles.begin(); mesh_ref_ptndx != mesh_intersecting_triangles.end(); ++mesh_ref_ptndx) {

            geosurf_triangle geosurf_tr = *mesh_ref_ptndx;

            Triangle3D mesh_triangle = geosurf_tr.triangle;

            std::string msg;
            std::vector<Point3D> inters_pts;
            std::tie(msg, inters_pts) = triangle_pair_inters_pts(dem_triangle, mesh_triangle);

            if (msg != "intersecting planes") {
                std::cout << "DEM triangle ndx: " << ndx_curr_dem_triangle << ", mesh triangle ndx: " << ndx_curr_mesh_triangle << "; error: " << msg << "\n";
                continue; };

            if (inters_pts.size() > 0) {

                for(std::vector<Point3D>::iterator result_ref_ptndx = inters_pts.begin(); result_ref_ptndx != inters_pts.end(); ++result_ref_ptndx) {

                    Point3D inters = *result_ref_ptndx;
                    inters_result result;
                    result.inter_pt = inters;
                    result.geoplane = geosurf_tr.geoplane;
                    intersections.push_back(result);
                };
            };

            ndx_curr_mesh_triangle++;
        };

        ndx_curr_dem_triangle++;

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

    std::vector<geosurf_triangle> mesh_triangles = extract_triangles_from_mesh(surf3d_mesh );
    std::cout << "\nNum. total mesh triangles is " << mesh_triangles.size() << "\n";

    // calculate solid volume used to check for mesh triangle volume intersection

    Space3DPartition dem_vol = datarrgrid.space_partition();

    // get mesh triangles intersecting with DEM boundaries

    std::vector<geosurf_triangle> mesh_intersecting_triangles = extract_intersecting_triangles(dem_vol, mesh_triangles );
    std::cout << "\nNum. intersecting mesh triangles is " << mesh_intersecting_triangles.size() << "\n";

    // transform DEM data into a vector of 3D points, valid or invalid

    std::vector<Point3D> dem_3dpts = create_pts_vector(datarrgrid.data(), datarrgrid.rect_domain());
    std::cout << "\nNum. dem 3d pts is " << dem_3dpts.size() << "\n";

    // create vector of valid DEM triangles, for intersection with triangles mesh

    std::vector<Triangle3D> dem_triangles = create_dem_triangles( dem_3dpts, datarrgrid.rect_domain().nrows(), datarrgrid.rect_domain().ncols() );
    std::cout << "\nNum. dem 3d triangles is " << dem_triangles.size() << "\n";

    // get intersection points

    std::vector<inters_result> intersections = intersect_dem_geosurface(dem_triangles, mesh_intersecting_triangles);
    std::cout << "\nNum. intersecting pts is " << intersections.size() << "\n";

    // write intersection points

    if (intersections.size() > 0) {
        std::ofstream outdatafile{output_datafile_path};
        outdatafile << "x,y,z,dipdir,dipang\n";
        for(std::vector<inters_result>::iterator intersection_ref = intersections.begin(); intersection_ref != intersections.end(); ++intersection_ref) {
            inters_result intersection = *intersection_ref;
            Point3D inters_pt = intersection.inter_pt;
            GeologicalPlane geoplane = intersection.geoplane;
            outdatafile << inters_pt.x() << "," << inters_pt.y() << "," << inters_pt.z() << "," << geoplane.dipdir() << "," << geoplane.dipangle() << "\n";
        };
    };

    return 0;


};







