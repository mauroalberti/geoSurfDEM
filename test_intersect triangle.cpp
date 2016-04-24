#include "data_processing.hpp"


int main() {

    Point3D pt1 = Point3D(0.0, 0.0, 0.0);
    Point3D pt2 = Point3D(1.0, 0.0, 0.0);
    Point3D pt3 = Point3D(0.0, 1.0, 0.0);

    Point3D pt4 = Point3D(-2.0, 0.5, 1.0);
    Point3D pt5 = Point3D(2.0, 0.5, 1.0);
    Point3D pt6 = Point3D(-2.0, 0.5, -20.0);

    Triangle3D dem_tr = Triangle3D(pt1, pt2, pt3);
    std::cout << "dem triangle area is " << dem_tr.area() << "\n";
    if (dem_tr.area() < 1.0e-10) {
        std::cout << "Error - degenerate dem triangle\n";
        exit(1);};
    CartesianPlane dem_tr_plane = dem_tr.to_cartes_plane();
    std::cout << "dem plane parameters - a: " << dem_tr_plane.a() << " b: " << dem_tr_plane.b() << " c: " << dem_tr_plane.c() << " d: " << dem_tr_plane.d() << "\n";

    Triangle3D mesh_tr = Triangle3D(pt4, pt5, pt6);
    std::cout << "mesh triangle area is " << mesh_tr.area() << "\n";
    if (mesh_tr.area() < 1.0e-10) {
        std::cout << "Error - degenerate mesh triangle\n";
        exit(1);}

    CartesianPlane mesh_tr_plane = mesh_tr.to_cartes_plane();
    std::cout << "mesh plane parameters - a: " << mesh_tr_plane.a() << " b: " << mesh_tr_plane.b() << " c: " << mesh_tr_plane.c() << " d: " << mesh_tr_plane.d() << "\n";

    double plane_angle = dem_tr_plane.angle(mesh_tr_plane);
    std::cout << "plane angle is " << plane_angle << "\n";

    // check parallelism/coincidence between the two planes
    bool plane_parallelism = dem_tr_plane.isparallel(mesh_tr_plane);
    if (plane_parallelism) {
        std::cout << "Warning: parallel planes\n";
        bool coincident_planes = dem_tr_plane.isequidistant(mesh_tr_plane);
        if (coincident_planes) {
            std::cout << "Warning: coincident planes\n";};
        exit(1); };

    Line3D inters_line = mesh_tr_plane.intersect(dem_tr_plane);
    Point3D iline_pt = inters_line.orig_pt();
    Vector3D iline_versor = inters_line.versor();

    std::cout << "point: " << iline_pt.x() << " " << iline_pt.y() << " " << iline_pt.z() << "\n";
    std::cout << "versor: " << iline_versor.x() << " " << iline_versor.y() << " " << iline_versor.z() << "\n";


    //Triangle3D tria2 = Triangle3D(pt4, pt5, pt6);



    //std::vector<Point3D> inter_pts = get_inters_pts(tria2, tria1);

    /*
    std::cout << "intersection point number is " << inter_pts.size() << "\n";
    for (uint i = 0; i < inter_pts.size(); i++) {
        Point3D pt = inter_pts[i];
        std::cout <<  pt.x() << ',' << pt.y() << ',' << pt.z() << '\n';
    };
    */


};










