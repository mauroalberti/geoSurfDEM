#include "data_processing.hpp"

struct inters_result {
  Point3D inter_pt;
  uint dem_tr_ndx, geosurf_tr_ndx;
};


std::tuple<std::string, std::vector<Point3D> > triangle_pair_inters_pts(Triangle3D dem_triangle, Triangle3D mesh_triangle) {

    // declares return variables

    std::string msg;
    std::vector<Point3D> inters_pts_dem;

    // check input validity

    if (dem_triangle.area() < 1.0e-10) {
        msg = "degenerate DEM triangle"; }
    else if (mesh_triangle.area() < 1.0e-10) {
        msg = "degenerate mesh triangle"; }
    else {
         // get cartesian plane from DEM triangle
        CartesianPlane dem_tr_plane = dem_triangle.to_cartes_plane();
        // get cartesian plane from geological surface mesh triangle
        CartesianPlane mesh_tr_plane = mesh_triangle.to_cartes_plane();
        // check parallelism/coincidence between the two planes
        bool plane_parallelism = dem_tr_plane.isparallel(mesh_tr_plane);
        if (plane_parallelism) {
            bool coincident_planes = dem_tr_plane.isequidistant(mesh_tr_plane);
            if (coincident_planes) {
                msg = "paralllel coincident planes";
            }
            else {
                msg = "parallel non-coincident planes";
            };
        }
        else {
            Line3D inters_line = mesh_tr_plane.intersect(dem_tr_plane);
            Point3D iline_pt = inters_line.orig_pt();
            Vector3D iline_versor = inters_line.versor();
            inters_pts_dem = find_triangle_inters(dem_triangle, inters_line);
            msg = "intersecting planes";
        };
    };

    return  std::make_tuple(msg, inters_pts_dem);

};

/*
std::string string_intersections()

        Point3D dem_pt1, dem_pt2, dem_pt3, mesh_pt1, mesh_pt2, mesh_pt3;
        dem_pt1 = dem_triangle.pt(0); dem_pt2 = dem_triangle.pt(1); dem_pt3 = dem_triangle.pt(2);
        mesh_pt1 = mesh_triangle.pt(0); mesh_pt2 = mesh_triangle.pt(1); mesh_pt3 = mesh_triangle.pt(2);
        outdebugfile << "dem triangle\n";
        outdebugfile << dem_pt1.x() << ", " << dem_pt1.y() << ", " << dem_pt1.z() << "\n";
        outdebugfile << dem_pt2.x() << ", " << dem_pt2.y() << ", " << dem_pt2.z() << "\n";
        outdebugfile << dem_pt3.x() << ", " << dem_pt3.y() << ", " << dem_pt3.z() << "\n";
        outdebugfile << "mesh triangle\n";
        outdebugfile << mesh_pt1.x() << ", " << mesh_pt1.y() << ", " << mesh_pt1.z() << "\n";
        outdebugfile << mesh_pt2.x() << ", " << mesh_pt2.y() << ", " << mesh_pt2.z() << "\n";
        outdebugfile << mesh_pt3.x() << ", " << mesh_pt3.y() << ", " << mesh_pt3.z() << "\n";
        for (uint i = 0; i < inters_pts_dem.size(); i++) {
            Point3D curr_pt = inters_pts_dem[i];
            outdebugfile << curr_pt.x() << ", " << curr_pt.y() << ", " << curr_pt.z() << "\n"; };
*/


/*
struct inters_result {
  Point3d inter_pt;
  uint dem_tr_ndx, geosurf_tr_ndx;
};
*/


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
            std::vector<Point3D> inters_pts;
            std::tie(msg, inters_pts) = triangle_pair_inters_pts(dem_triangle, mesh_triangle);

            if (msg != "intersecting planes") {
                std::cout << "DEM triangle ndx: " << ndx_curr_dem_triangle << ", mesh triangle ndx: " << ndx_curr_mesh << "; error: " << msg << "\n";
                continue; };

            if (inters_pts.size() > 0) {

                int ndx_curr_inters_pt = 0;
                for(std::vector<Point3D>::iterator pt_ref_ptndx = inters_pts.begin(); pt_ref_ptndx != inters_pts.end(); ++pt_ref_ptndx) {
                    ++ndx_curr_inters_pt;
                    Point3D inters_pt = *pt_ref_ptndx;
                    inters_result result;
                    result.inter_pt = inters_pt;
                    result.dem_tr_ndx = ndx_curr_dem_triangle;
                    result.geosurf_tr_ndx = ndx_curr_mesh;
                    intersections.push_back(result);
                };
            };
        };
    };

    return intersections;
};


int main() {

    // read DEM data from input file

    std::string input_dem_path = "./test_data/malpi_aster_w4u3_rit.asc";

    // input 3D geosurface (VTK format)

    std::string input_vtk_path = "./test_data/geosurf3d_01.vtk";

    // output data file

    std::string output_datafile_path = "./test_data/outdata_08.xyz";

    // output source dem data file

    std::string output_srcdem_file_path = "./test_data/src_dem_data.csv";

    // output source geosurface data file

    std::string output_srcgeosurf_file_path = "./test_data/src_geosurface_data.csv";

    ///////////////////////////

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
        for(std::vector<inters_result>::iterator intersection_ref = intersections.begin(); intersection_ref != intersections.end(); ++intersection_ref) {
            inters_result intersection = *intersection_ref;
            Point3D inters_pt = intersection.inter_pt;
            outdatafile << inters_pt.x() << "," << inters_pt.y() << "," << inters_pt.z() << "," << intersection.dem_tr_ndx << "," << intersection.geosurf_tr_ndx << "\n"; };  };

    return 0;


};







