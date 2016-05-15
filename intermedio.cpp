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


std::tuple<std::string, std::vector<Point3D> > analize_point_pair(std::vector<Point3D> dem_pt_triplet, std::vector<Point3D> mesh_pt_triplet) {

    std::vector<Point3D> inters_points;
    std::string msg;

    // extraction of DEM triangle points
    Point3D dem_pt_1 = dem_pt_triplet[0];
    Point3D dem_pt_2 = dem_pt_triplet[1];
    Point3D dem_pt_3 = dem_pt_triplet[2];

    // extraction of mesh triangle points
    Point3D mesh_pt_1 = mesh_pt_triplet[0];
    Point3D mesh_pt_2 = mesh_pt_triplet[1];
    Point3D mesh_pt_3 = mesh_pt_triplet[2];


    Triangle3D dem_triangle = Triangle3D(dem_pt_1, dem_pt_2, dem_pt_3);
    if (dem_triangle.area() < 1.0e-10) {
        msg = "degenerate DEM triangle";
        return  std::make_tuple(msg, inters_points); };

    Triangle3D mesh_triangle = Triangle3D(mesh_pt_1, mesh_pt_2, mesh_pt_3);
    if (mesh_triangle.area() < 1.0e-10) {
        msg = "degenerate mesh triangle";
        return  std::make_tuple(msg, inters_points); };

    return triangle_pair_inters_pts(dem_triangle, mesh_triangle);

};










