#ifndef CGAL_MESH_H
#define CGAL_MESH_H

#include <vector>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>

namespace PMP = CGAL::Polygon_mesh_processing;

using Point = CGAL::Simple_cartesian<double>::Point_3;
using Vector = CGAL::Simple_cartesian<double>::Vector_3;
using Mesh = CGAL::Surface_mesh<Point>;

namespace cgal
{
    std::vector<double> get_max_point( const Mesh &mesh );
    std::vector<double> get_min_point( const Mesh &mesh );
    bool is_point_in_bounds( const Point& p, const std::vector<double> &maxp, const std::vector<double> &minp, const double threshold=1.0e-10 );
    bool is_vertex_in_bounds( const Mesh::Vertex_index& v, const Mesh &mesh, const std::vector<double> &maxp, const std::vector<double> &minp, const double threshold=1.0e-10 );
    bool is_face_in_bounds( const Mesh::Face_index& f, const Mesh &mesh, const std::vector<double> &maxp, const std::vector<double> &minp, const double threshold=1.0e-10 );
    bool faces_have_identical_coordinates( const Mesh::Face_index& f1, const Mesh &mesh1, const Mesh::Face_index& f2, const Mesh &mesh2, const double threshold=1.0e-10 );
    Mesh extract_contact_surface( const Mesh &mesh1, const Mesh &mesh2, const double threshold=1.0e-10 );

} //namespace cgal


#endif