/*
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <vector>
#include <tuple>
#include <assert.h>

#include <stdlib.h>     // atof
*/


#include "data_processing.hpp"


int main() {

    // read DEM data from input file
    std::string input_dem_path = "./test_data/malpi_aster_w4u3_rit.asc";

    // input 3D geosurface (VTK format)
    std::string input_vtk_path = "./test_data/geosurf3d_01.vtk";

    // debug output mesh file
    std::string output_debmeshfile_path = "./test_data/deb_mesh_06.txt";
    std::string output_debmeshtrianglesfile_path = "./test_data/deb_meshtriangles_06.txt";
    std::string output_debmeshinterstrfile_path = "./test_data/deb_meshinterstr_06.txt";

    // debug output file
    std::string output_logfile_path = "./test_data/log_06.txt";

    // output data file
    std::string output_datafile_path = "./test_data/outdata_06.xyz";

    ///////////////////////////


    std::cout << "\n*** geoSurfDEM *** \n";
    std::cout << "\nApplication for determining intersections between 3D geosurfaces and DEM topography\n\n";

    /*
    Point3D p1 = Point3D(1.0,1.0,1.0);
    Point3D p2 = Point3D(-1.0,1.0,0.0);
    Point3D p3 = Point3D(2.0,0.0,3.0);

    CartesianPlane test = CartesianPlane(p1,p2,p3);
    std::cout << test.a() << ", " << test.b() << ", " << test.c() << ", " << test.d() << "\n";

    assert (test.point_in_plane(p1));
    assert (test.point_in_plane(p2));
    assert (test.point_in_plane(p3));
    */



    // read input DEM data
    DataRRGrid datarrgrid = read_esri_ascii_dem( input_dem_path );
    RectangularDomain rect_dom = datarrgrid.rect_domain();

    // read VTK data from input file
    MeshTriangleStrip surf3d_mesh;
    try {
        surf3d_mesh = read_vtk_data_ascii( output_debmeshfile_path, input_vtk_path ); }
    catch (int e) {
        std::cout << "\nError: program will stop\n";
        return -1; }

    // get triangles (Triangle3D) from mesh
    std::vector<Triangle3D> mesh_triangles = extract_triangles_from_mesh(output_debmeshtrianglesfile_path, surf3d_mesh );
    std::cout << "\nNum. total mesh triangles is " << mesh_triangles.size() << "\n";

    // calculate solid volume used to check for mesh triangle volume intersection
    Space3DPartition dem_vol = datarrgrid.space_partition();

    // get mesh triangles intersecting with DEM boundaries
    std::vector<Triangle3D> mesh_intersecting_triangles = extract_intersecting_triangles( output_debmeshinterstrfile_path, dem_vol, mesh_triangles );
    std::cout << "\nNum. intersecting mesh triangles is " << mesh_intersecting_triangles.size() << "\n";

    // transform DEM data into a vector of 3D points, valid or invalid
    std::vector<Point3D> dem_3dpts = create_pts_vector(datarrgrid.data(), datarrgrid.rect_domain());
    std::cout << "\nNum. dem 3d pts is " << dem_3dpts.size() << "\n";
    assert (dem_3dpts.size()==rect_dom.nrows()*rect_dom.ncols());

    // create vector of valid DEM triangles, for intersection with triangles mesh
    std::vector<Triangle3D> dem_triangles = create_dem_triangles( dem_3dpts, datarrgrid.rect_domain().nrows(), datarrgrid.rect_domain().ncols() );
    std::cout << "\nNum. dem 3d triangles is " << dem_triangles.size() << "\n";
    assert (dem_triangles.size()==2*(rect_dom.nrows()-1)*(rect_dom.ncols()-1));

    // get intersection points
    std::vector<Point3D> intersection_pts = intersect_dem_geosurface(output_logfile_path, dem_triangles, mesh_intersecting_triangles);
    std::cout << "\nNum. intersecting pts is " << intersection_pts.size() << "\n";

    // write intersection points
    if (intersection_pts.size() > 0) {

        std::ofstream outdatafile{output_datafile_path};

        for(std::vector<Point3D>::iterator pt_ref_ptndx = intersection_pts.begin(); pt_ref_ptndx != intersection_pts.end(); ++pt_ref_ptndx) {

            Point3D inters_pt = *pt_ref_ptndx;

            outdatafile << inters_pt.x() << "," << inters_pt.y() << "," << inters_pt.z() << "\n"; };  };

    return 0;
    //

};
