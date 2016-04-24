#include "data_processing.hpp"


int main() {

    Point3D pt1 = Point3D(0.0, 0.0, 0.0);
    Point3D pt2 = Point3D(1.0, 0.0, 0.0);
    Point3D pt3 = Point3D(0.0, 1.0, 0.0);

    Point3D pt4 = Point3D(-2.0, 0.5, 1.0);
    Point3D pt5 = Point3D(2.0, 0.5, 1.0);
    Point3D pt6 = Point3D(-2.0, 0.5, -20.0);

    Triangle3D tria1 = Triangle3D(pt1, pt2, pt3);
    std::cout << "area 1 is " << tria1.area() << "\n";
    if (tria1.area() < 1.0e-10) {
        std::cout << "Error - degenerate triangle\n";
        exit(1);};

    CartesianPlane plane1 = tria1.to_cartes_plane();
    std::cout << "plane parameters - a: " << plane1.a() << " b: " << plane1.b() << " c: " << plane1.c() << " d: " << plane1.d() << "\n";

    Triangle3D tria2 = Triangle3D(pt4, pt5, pt6);
    std::cout << "area 2 is " << tria2.area() << "\n";
    if (tria2.area() < 1.0e-10) {
        std::cout << "Error - degenerate triangle\n";
        exit(1);}

    CartesianPlane plane2 = tria2.to_cartes_plane();
    std::cout << "plane parameters - a: " << plane2.a() << " b: " << plane2.b() << " c: " << plane2.c() << " d: " << plane2.d() << "\n";

    double plane_angle = plane1.angle(plane2);
    std::cout << "angle is " << plane_angle << "\n";

    // check parallelism/coincidence between the two planes
    bool plane_parallelism = plane1.isparallel(plane2);
    if (plane_parallelism) {
        std::cout << "Warning: parallel planes\n";
        bool coincident_planes = plane1.isequidistant(plane2);
        if (coincident_planes) {
            std::cout << "Warning: coincident planes\n";};
        exit(1); };





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










