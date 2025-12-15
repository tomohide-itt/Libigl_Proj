#include <cfloat>
#include "cgal.h"

//メッシュの最大座標を取得する 
std::vector<double> cgal::get_max_point( const Mesh &mesh )
{
    std::vector<double> maxp = { -DBL_MAX, -DBL_MAX, -DBL_MAX };
    for( const auto &v : mesh.vertices() )
    {
        Point p = mesh.point(v);
        if( p.x() > maxp[0] ) maxp[0] = p.x();
        if( p.y() > maxp[1] ) maxp[1] = p.y();
        if( p.z() > maxp[2] ) maxp[2] = p.z();
    }
    return maxp;
}
//メッシュの最小座標を取得する
std::vector<double> cgal::get_min_point( const Mesh &mesh )
{
    std::vector<double> minp = { DBL_MAX, DBL_MAX, DBL_MAX };
    for( const auto &v : mesh.vertices() )
    {
        Point p = mesh.point(v);
        if( p.x() < minp[0] ) minp[0] = p.x();
        if( p.y() < minp[1] ) minp[1] = p.y();
        if( p.z() < minp[2] ) minp[2] = p.z();
    }
    return minp;
}
//ポイントがmaxp, minpで指定した範囲内にあるかを判定
bool cgal::is_point_in_bounds( const Point& p, const std::vector<double> &maxp, const std::vector<double> &minp, const double threshold )
{
    if( (p.x() - minp[0]) < threshold ) return false;
    if( (p.y() - minp[1]) < threshold ) return false;
    if( (p.z() - minp[2]) < threshold ) return false;
    if( (p.x() - maxp[0]) > threshold ) return false;
    if( (p.y() - maxp[1]) > threshold ) return false;
    if( (p.z() - maxp[2]) > threshold ) return false;
    return true;
}
//頂点がmaxp, minpで指定した範囲内にあるかを判定
bool cgal::is_vertex_in_bounds( const Mesh::Vertex_index& v, const Mesh &mesh, const std::vector<double> &maxp, const std::vector<double> &minp, const double threshold )
{
    Point p = mesh.point(v);
    return is_point_in_bounds( p, maxp, minp, threshold );
}
//面がmaxp, minpで指定した範囲内にあるかを判定
bool cgal::is_face_in_bounds( const Mesh::Face_index& f, const Mesh &mesh, const std::vector<double> &maxp, const std::vector<double> &minp, const double threshold )
{
    for( auto v : CGAL::vertices_around_face( mesh.halfedge(f), mesh ) )
    {
      if( !is_vertex_in_bounds(v,mesh,maxp,minp,threshold) ) return false;
    }
    return true;
}
//
bool cgal::faces_have_identical_coordinates( const Mesh::Face_index& f1, const Mesh &mesh1, const Mesh::Face_index& f2, const Mesh &mesh2, const double threshold )
{
    std::vector<Point> pnts1;
    for( const auto& v : CGAL::vertices_around_face( mesh1.halfedge(f1), mesh1 ) )
    {
        pnts1.push_back( mesh1.point(v) );
    }
    std::vector<Point> pnts2;
    for( const auto& v : CGAL::vertices_around_face( mesh2.halfedge(f2), mesh2 ) )
    {
        pnts2.push_back( mesh2.point(v) );
    }
    if( ( CGAL::squared_distance(pnts1[0], pnts2[0]) < threshold && CGAL::squared_distance(pnts1[1], pnts2[1]) < threshold && CGAL::squared_distance(pnts1[2], pnts2[2]) < threshold ) ||
        ( CGAL::squared_distance(pnts1[0], pnts2[0]) < threshold && CGAL::squared_distance(pnts1[1], pnts2[2]) < threshold && CGAL::squared_distance(pnts1[2], pnts2[1]) < threshold ) ||
        ( CGAL::squared_distance(pnts1[0], pnts2[1]) < threshold && CGAL::squared_distance(pnts1[1], pnts2[0]) < threshold && CGAL::squared_distance(pnts1[2], pnts2[2]) < threshold ) ||
        ( CGAL::squared_distance(pnts1[0], pnts2[1]) < threshold && CGAL::squared_distance(pnts1[1], pnts2[2]) < threshold && CGAL::squared_distance(pnts1[2], pnts2[0]) < threshold ) ||
        ( CGAL::squared_distance(pnts1[0], pnts2[2]) < threshold && CGAL::squared_distance(pnts1[1], pnts2[0]) < threshold && CGAL::squared_distance(pnts1[2], pnts2[1]) < threshold ) ||
        ( CGAL::squared_distance(pnts1[0], pnts2[2]) < threshold && CGAL::squared_distance(pnts1[1], pnts2[1]) < threshold && CGAL::squared_distance(pnts1[2], pnts2[0]) < threshold ) )
    {
        return true;
    }
    return false;
}
//2つのメッシュの接面を抽出する
//接面はの面としてmesh1側の面を採用する
Mesh cgal::extract_contact_surface( const Mesh &mesh1, const Mesh &mesh2, const double threshold )
{
  Mesh surface;
  std::map<Mesh::Vertex_index, Mesh::Vertex_index> vmap;
  std::vector<double> maxp2 = get_max_point( mesh2 );
  std::vector<double> minp2 = get_min_point( mesh2 );
  //+++
  //std::cout << "mesh2 maxp: " << maxp2[0] << ", " << maxp2[1] << ", " << maxp2[2] << std::endl;
  //std::cout << "mesh2 minp: " << minp2[0] << ", " << minp2[1] << ", " << minp2[2] << std::endl;
  //---
  for( const auto &f1 : mesh1.faces() )
  {
    //mesh1の面f1がmesh2の範囲内にない場合は，とばす
    bool b = is_face_in_bounds( f1, mesh1, maxp2, minp2, threshold );
    //+++
    if( !b )
    {
        std::cout << "face1: " << f1 << " is_in_bounds_mesh2 = " << std::boolalpha << b << std::endl;
        std::cout << "  mesh2 maxp: " << maxp2[0] << ", " << maxp2[1] << ", " << maxp2[2] << std::endl;
        std::cout << "  mesh2 minp: " << minp2[0] << ", " << minp2[1] << ", " << minp2[2] << std::endl;
        for( const auto &v : CGAL::vertices_around_face( mesh1.halfedge(f1), mesh1 ) )
        {
            std::cout << "  " << v << ": " << mesh1.point(v) << std::endl;
        }
    }
    //---
    if( !is_face_in_bounds( f1, mesh1, maxp2, minp2, threshold ) ) continue;

    std::vector<Point> pnts1;
    std::vector<Mesh::Vertex_index> verts1;
    for( const auto &v : CGAL::vertices_around_face( mesh1.halfedge(f1), mesh1 ) )
    {
      pnts1.push_back(mesh1.point(v));
      verts1.push_back(v);
    }
    //
    for( const auto &f2 : mesh2.faces() )
    {
      //+++
      //std::cout << "face1: " << f1 << " num = " << mesh1.num_faces() << std::endl;
      //std::cout << "  " << verts1[0] << ": " << mesh1.point(verts1[0]) << std::endl;
      //std::cout << "  " << verts1[1] << ": " << mesh1.point(verts1[1]) << std::endl;
      //std::cout << "  " << verts1[2] << ": " << mesh1.point(verts1[2]) << std::endl;
      //std::vector<Point> pnts2;
      //std::vector<Mesh::Vertex_index> verts2;
      //for( const auto &v : CGAL::vertices_around_face( mesh2.halfedge(f2), mesh2 ) )
      //{
      //  pnts2.push_back(mesh2.point(v));
      //  verts2.push_back(v);
      //}
      //std::cout << "face2: " << f2 << " num = " << mesh2.num_faces() << std::endl;
      //std::cout << "  " << verts2[0] << ": " << mesh2.point(verts2[0]) << std::endl;
      //std::cout << "  " << verts2[1] << ": " << mesh2.point(verts2[1]) << std::endl;
      //std::cout << "  " << verts2[2] << ": " << mesh2.point(verts2[2]) << std::endl;
      //---
      if( faces_have_identical_coordinates( f1, mesh1, f2, mesh2, threshold ) )
      {
        //+++
        //std::cout << "same triangle " << std::endl;
        //---
        std::vector<Mesh::Vertex_index> verts(3);
        for( int i=0; i<3; i++ )
        {
            if( vmap.count(verts1[i]) )
            {
                verts[i] = vmap[verts1[i] ];
            }
            else
            {
                verts[i] = surface.add_vertex( pnts1[i] );
                vmap[verts1[i]] = verts[i];
            }
        }
        surface.add_face( verts );
        break;
      }
    }
  }
  return surface;
}