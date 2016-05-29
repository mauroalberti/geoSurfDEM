#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <vector>
#include <tuple>
#include <assert.h>

#include <stdlib.h>     /* atof */

#include "geographic.hpp"


struct geosurf_triangle {
    Triangle3D triangle;
    GeologicalPlane geoplane;
};


int get_dem_param_int(const std::string&, const std::string&);


double get_dem_param_double(const std::string&, const std::string&);


DataRRGrid read_esri_ascii_dem(std::string);


MeshTriangleStrip read_vtk_data_ascii( std::string);


std::vector<geosurf_triangle> extract_triangles_from_mesh(MeshTriangleStrip);


std::vector<geosurf_triangle> extract_intersecting_triangles(Space3DPartition, std::vector<geosurf_triangle>);


int vect_ndx(int, int, int);


std::vector<int> dem_indices(int, int);


std::vector<double> unrotated_coords(int, int, int, double, double);


std::vector<Point3D> create_pts_vector(NumericData, RectangularDomain);


bool check_triangle_points_validity(Point3D, Point3D, Point3D);


std::vector<Triangle3D> create_dem_triangles(std::vector<Point3D>, int, int);


std::tuple<Point3D, bool> intersect_segments(Line3D, Segment3D);


std::vector<Point3D> get_inters_pts(Triangle3D, Triangle3D);


std::vector<Point3D> intersect_dem_geosurface(std::string, std::vector<Triangle3D>, std::vector<Triangle3D>);


std::vector<Point3D> find_triangle_inters(Triangle3D, Line3D);


