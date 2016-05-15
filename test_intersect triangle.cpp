#include "data_processing.hpp"


std::tuple<std::string, std::vector<Point3D> > triangle_pair_inters_pts(Triangle3D dem_triangle, Triangle3D mesh_triangle) {

    // declares return variables

    std::string msg;
    std::vector<Point3D> inters_pts_dem;

    // get cartesian plane from DEM triangle

    CartesianPlane dem_tr_plane = dem_triangle.to_cartes_plane();

    // get cartesian plane from geological surface mesh triangle

    CartesianPlane mesh_tr_plane = mesh_triangle.to_cartes_plane();

    // check parallelism/coincidence between the two planes

    bool plane_parallelism = dem_tr_plane.isparallel(mesh_tr_plane);
    if (plane_parallelism) {
        bool coincident_planes = dem_tr_plane.isequidistant(mesh_tr_plane);
        if (coincident_planes) {
            msg = "coincident planes";
            return  std::make_tuple(msg, inters_pts_dem); }; };


    Line3D inters_line = mesh_tr_plane.intersect(dem_tr_plane);
    Point3D iline_pt = inters_line.orig_pt();
    Vector3D iline_versor = inters_line.versor();

    inters_pts_dem = find_triangle_inters(dem_triangle, inters_line);

    msg = "intersecting planes";

    return  std::make_tuple(msg, inters_pts_dem);

};


int main() {

    Point3D pt1 = Point3D(0.0, 0.0, 0.0);
    Point3D pt2 = Point3D(1.0, 0.0, 0.0);
    Point3D pt3 = Point3D(0.0, 1.0, 0.0);

    Point3D pt4 = Point3D(-20.0, 0.5, 1.0);
    Point3D pt5 = Point3D(20.0, 0.5, 1.0);
    Point3D pt6 = Point3D(-20.0, 0.5, -20.0);

    //std::cout << "\nDEM triangle\n";

    Triangle3D dem_triangle = Triangle3D(pt1, pt2, pt3);
    //std::cout << " - area is " << dem_triangle.area() << "\n";
    if (dem_triangle.area() < 1.0e-10) {
        std::cout << "Error - degenerate triangle\n";
        exit(1);};

    //std::cout << "\nMesh triangle\n";

    Triangle3D mesh_triangle = Triangle3D(pt4, pt5, pt6);
    //std::cout << " - area is " << mesh_triangle.area() << "\n";
    if (mesh_triangle.area() < 1.0e-10) {
        std::cout << "Error - degenerate triangle\n";
        exit(1);}

    std::vector<Point3D> found_pts;
    std::string msg;

    std::tie(msg, found_pts) = triangle_pair_inters_pts(dem_triangle, mesh_triangle);

    std::cout << msg << "\n";
    std::cout << "num. found pts: " << found_pts.size() << "\n";
    for (uint i = 0; i < found_pts.size(); i++) {
        Point3D curr_pt = found_pts[i];
        std::cout << curr_pt.x() << ", " << curr_pt.y() << ", " << curr_pt.z() << "\n";
    };


};










